"""
xarray Dataset creation from multiple CTD casts.
"""
import glob
import re
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
from . import config
from .io import get_cast, castaway_cast
from .processing import process_ctd, COMMON_PARAMS, SBE_PARAMS, RBR_PARAMS, CASTAWAY_PARAMS
from .utils import parse_folder_date
from .batch import load_surface_cutoffs


# CF-convention variable metadata: long_name, units, standard_name
VAR_ATTRS: dict[str, dict] = {
    'tv290C':   {'long_name': 'Temperature (ITS-90)',          'units': 'degC',   'standard_name': 'sea_water_temperature'},
    'sal00':    {'long_name': 'Practical Salinity',            'units': 'PSU',    'standard_name': 'sea_water_practical_salinity'},
    'prdM':     {'long_name': 'Pressure',                      'units': 'dbar',   'standard_name': 'sea_water_pressure'},
    'Pressure': {'long_name': 'Pressure',                      'units': 'dbar',   'standard_name': 'sea_water_pressure'},
    'depSM':    {'long_name': 'Depth',                         'units': 'm',      'standard_name': 'depth'},
    'sigma0':   {'long_name': 'Potential Density Anomaly σ₀',  'units': 'kg m-3', 'standard_name': 'sea_water_sigma_theta'},
    'rho0':     {'long_name': 'In Situ Density Anomaly',       'units': 'kg m-3', 'standard_name': 'sea_water_density'},
}

# Processing parameters — derived directly from processing.py constants
PROCESSING_ATTRS: dict = {**COMMON_PARAMS, **SBE_PARAMS, **RBR_PARAMS, **CASTAWAY_PARAMS}


def _match_castaway_repeat(
    base: str,
    date: pd.Timestamp,
    cast_time,
    cast_times: dict,
    already_matched: set,
) -> int:
    """
    Return the repeat index of the reference (RBR or SBE) cast whose timestamp
    is closest to *cast_time* at the same station and date, excluding any
    (base, date, repeat, instrument) keys already claimed by a previous CastAway
    cast in the same day.

    Falls back to the next unused integer repeat when:
      - no reference cast exists at this station/date, or
      - cast_time is None (CastAway header had no timestamp).

    Args:
        base:            Base station name ('S1', 'G1', …).
        date:            Cruise date (pd.Timestamp).
        cast_time:       CastAway cast start time (pd.Timestamp or None).
        cast_times:      The full cast_times dict from create_ctd_dataset().
        already_matched: Set of (base, date, repeat) tuples already assigned to
                         a CastAway cast; updated in-place on successful match.

    Returns: integer repeat index.
    """
    # All reference entries at this (base, date) from non-CastAway instruments
    candidates = [
        (r, t)
        for (b, d, r, inst), t in cast_times.items()
        if b == base and d == date and inst != 'castaway'
        and t is not None
        and (base, date, r) not in already_matched
    ]

    if candidates and cast_time is not None:
        try:
            ct = pd.Timestamp(cast_time)
            closest_repeat, _ = min(
                candidates,
                key=lambda item: abs((pd.Timestamp(item[1]) - ct).total_seconds()),
            )
            already_matched.add((base, date, closest_repeat))
            return closest_repeat
        except Exception:
            pass

    # Fallback: next integer repeat not already used for any instrument
    used_repeats = {r for (b, d, r, _) in cast_times if b == base and d == date}
    used_repeats.update(r for b, d, r in already_matched if b == base and d == date)
    repeat = 0
    while repeat in used_repeats:
        repeat += 1
    already_matched.add((base, date, repeat))
    return repeat


def _parse_station_repeat(raw_station: str) -> tuple[str, int]:
    """
    Split a raw station name into (base_station, repeat_index).

    Primary casts use the plain name ('S1', 'G1', 'B') → repeat 0.
    Repeat casts append a lowercase letter ('S1b', 'S1c', 'G1b') → repeat 1, 2, ...
    """
    m = re.fullmatch(r'([A-Z0-9]+)([a-z]?)', raw_station)
    base   = m.group(1)
    suffix = m.group(2)
    repeat = 0 if not suffix else ord(suffix) - ord('a') + 1
    return base, repeat


def create_ctd_dataset(surface_cutoffs_file=None, castaway_path=None):
    """
    Create xarray Dataset from all CTD files in DATA_PATH and CASTAWAY_PATH.

    Scans the following file types:
      - *.cnv   → SeaBird SBE (processed CTD export)
      - *.xlsx  → RBR Excel export
      - *.hex   → SeaBird SBE raw hex (requires paired XMLCON; skipped if T/C absent)
      - CastAway CSVs from CASTAWAY_PATH/Viola_* subfolders

    Returns: Dataset with dimensions (station, date, repeat, instrument, depth).
      - instrument ∈ ['castaway', 'rbr', 'sbe'] — the source CTD for each profile
      - repeat=0  primary cast, repeat=1 first repeat ('b'), repeat=2 second ('c'), …
      - cast_time coordinate holds the actual timestamp per (station, date, repeat, instrument)
      - Missing combinations and depths outside a cast's range are filled with NaN.

    Use ds.sel(instrument='rbr') vs ds.sel(instrument='castaway') to compare
    profiles from different instruments at the same station and date.

    Args:
        surface_cutoffs_file: Optional path to the CSV produced by
            review_surface_cutoffs().  When provided, each cast's surface
            cutoff is looked up by '{folder}_{station}' key and passed to
            process_ctd() as surface_cutoff_m, overriding the automatic
            stability algorithm for that cast.
        castaway_path: Path to CastAway_profiles folder
            (default: CASTAWAY_PATH from config).
    """
    data_path     = config.DATA_PATH
    castaway_path = Path(castaway_path) if castaway_path else config.CASTAWAY_PATH

    cutoffs: dict = {}
    if surface_cutoffs_file is not None:
        cutoffs = load_surface_cutoffs(surface_cutoffs_file)
        print(f"Loaded {len(cutoffs)} manual surface cutoffs from {surface_cutoffs_file}")

    # key: (base_station, date, repeat_index, instrument)  value: processed DataFrame
    data_dict:  dict[tuple[str, pd.Timestamp, int, str], pd.DataFrame] = {}
    cast_times: dict[tuple[str, pd.Timestamp, int, str], pd.Timestamp] = {}

    # ------------------------------------------------------------------ #
    # 1.  SBE CNV + RBR XLSX — standard dated Data/ subfolders           #
    # ------------------------------------------------------------------ #
    files = (
        glob.glob(str(data_path / '**' / '*.cnv'), recursive=True) +
        glob.glob(str(data_path / '**' / '*.xlsx'), recursive=True)
    )

    for fpath in files:
        folder_date_name = Path(fpath).parent.name
        date = parse_folder_date(folder_date_name)
        if date is None:
            continue

        try:
            casts = get_cast(fpath)
        except Exception as e:
            print(f"Error loading {fpath}: {e}")
            continue

        for cast in casts:
            raw_station = cast.cast_meta.get('station')
            if raw_station is None:
                print(f"Station not found in cast from {fpath}, skipping")
                continue

            instrument = cast.cast_meta.get('instrument_type', 'unknown')
            base, repeat = _parse_station_repeat(raw_station)
            key = (base, date, repeat, instrument)
            label = f"{raw_station} [{instrument}] on {date.date()}"
            print(f"Processing: {label}")

            cutoff_key   = f"{folder_date_name}_{raw_station}"
            cutoff_depth = cutoffs.get(cutoff_key)
            try:
                proc_df = process_ctd(cast, surface_cutoff_m=cutoff_depth)
            except Exception as e:
                print(f"  ERROR processing {label}: {e}")
                continue
            data_dict[key]  = proc_df
            cast_times[key] = cast.cast_meta.get('time')

    # ------------------------------------------------------------------ #
    # 2.  CastAway CSVs from CASTAWAY_PATH/Viola_* subfolders             #
    # ------------------------------------------------------------------ #
    # Per-day tracking: which (base, date, repeat) slots have already been
    # assigned to a CastAway cast, so the greedy matcher doesn't double-assign.
    castaway_matched: set[tuple[str, pd.Timestamp, int]] = set()

    if castaway_path.exists():
        for subfolder in sorted(castaway_path.glob('Viola_*')):
            if not subfolder.is_dir():
                continue
            m = re.match(r'Viola_(\d{1,2})([a-z]{3})(\d{4})$', subfolder.name, re.IGNORECASE)
            if not m:
                continue
            try:
                date = pd.to_datetime(f"{m.group(1)} {m.group(2)} {m.group(3)}", dayfirst=True)
            except Exception:
                continue
            folder_date_name = date.strftime('%Y') + date.strftime('%b') + date.strftime('%d')

            for csv_file in sorted(subfolder.glob('*.csv')):
                try:
                    cast = castaway_cast(str(csv_file))
                except Exception as e:
                    print(f"Error loading {csv_file.name}: {e}")
                    continue

                raw_station = cast.cast_meta.get('station')
                if raw_station is None:
                    continue

                instrument = 'castaway'
                base, _ = _parse_station_repeat(raw_station)

                # Assign the same repeat as the RBR/SBE cast closest in time
                repeat = _match_castaway_repeat(
                    base, date,
                    cast.cast_meta.get('time'),
                    cast_times,
                    castaway_matched,
                )
                key = (base, date, repeat, instrument)

                label = f"{raw_station} [castaway repeat={repeat}] on {date.date()}"
                print(f"Processing: {label}")

                cutoff_key   = f"{folder_date_name}_{raw_station}"
                cutoff_depth = cutoffs.get(cutoff_key)
                try:
                    proc_df = process_ctd(cast, surface_cutoff_m=cutoff_depth)
                except Exception as e:
                    print(f"  ERROR processing {label}: {e}")
                    continue
                data_dict[key]  = proc_df
                cast_times[key] = cast.cast_meta.get('time')

    if not data_dict:
        raise ValueError("No valid CTD data found!")

    # ------------------------------------------------------------------ #
    # 3.  Build 5D arrays (station, date, repeat, instrument, depth)      #
    # ------------------------------------------------------------------ #
    bin_delta  = 0.25
    all_depths = np.concatenate([df.index.values for df in data_dict.values()])
    depth_grid = np.arange(
        np.floor(all_depths.min() * 4) / 4,
        np.ceil(all_depths.max()  * 4) / 4 + bin_delta,
        bin_delta,
    )
    depth_grid = np.round(depth_grid, 4)

    stations    = sorted({s            for s, _, _, _    in data_dict})
    dates       = sorted({d            for _, d, _, _    in data_dict})
    repeats     = sorted({r            for _, _, r, _    in data_dict})
    instruments = sorted({inst         for _, _, _, inst in data_dict})
    n_s, n_d, n_r, n_i, n_z = (
        len(stations), len(dates), len(repeats), len(instruments), len(depth_grid)
    )
    s_idx    = {s: i    for i, s    in enumerate(stations)}
    d_idx    = {d: i    for i, d    in enumerate(dates)}
    r_idx    = {r: i    for i, r    in enumerate(repeats)}
    inst_idx = {inst: i for i, inst in enumerate(instruments)}

    exclude_vars = {'time_local', 'timeS', 'time', 'station', 'instrument', 'cast_date'}
    all_var_names: set[str] = set()
    for proc_df in data_dict.values():
        all_var_names.update(c for c in proc_df.columns if c not in exclude_vars)
    var_names = sorted(all_var_names)

    # 5D arrays — NaN where a cast doesn't exist or doesn't reach a depth
    arrays     = {var: np.full((n_s, n_d, n_r, n_i, n_z), np.nan) for var in var_names}
    time_array = np.full((n_s, n_d, n_r, n_i), np.datetime64('NaT'), dtype='datetime64[ns]')

    for (base, date, repeat, instrument), proc_df in data_dict.items():
        si   = s_idx[base]
        di   = d_idx[date]
        ri   = r_idx[repeat]
        ii   = inst_idx[instrument]
        cast_depths = proc_df.index.values
        in_range    = (depth_grid >= cast_depths.min()) & (depth_grid <= cast_depths.max())
        for var in var_names:
            if var in proc_df.columns:
                arrays[var][si, di, ri, ii, in_range] = np.interp(
                    depth_grid[in_range], cast_depths, proc_df[var].values
                )
        t = cast_times.get((base, date, repeat, instrument))
        if t is not None:
            time_array[si, di, ri, ii] = np.datetime64(pd.Timestamp(t), 'ns')

    # Station lat/lon coords from STATIONS_DF
    station_coords: dict = {}
    if config.STATIONS_DF is not None:
        station_coords = {
            row['name']: (row['lat'], row['lon'])
            for _, row in config.STATIONS_DF.iterrows()
        }

    coords = {
        'station':    stations,
        'date':       dates,
        'repeat':     repeats,
        'instrument': instruments,
        'depth':      depth_grid,
        'cast_time':  (['station', 'date', 'repeat', 'instrument'], time_array),
    }
    if station_coords:
        coords['latitude']  = ('station', [station_coords.get(s, (np.nan,))[0]        for s in stations])
        coords['longitude'] = ('station', [station_coords.get(s, (np.nan, np.nan))[1] for s in stations])

    ds = xr.Dataset(
        {var: (['station', 'date', 'repeat', 'instrument', 'depth'], arr)
         for var, arr in arrays.items()},
        coords=coords,
    )

    # Variable-level metadata
    for var in ds.data_vars:
        if var in VAR_ATTRS:
            ds[var].attrs.update(VAR_ATTRS[var])

    # Dataset-level metadata
    ds.attrs.update(PROCESSING_ATTRS)
    ds.attrs['description'] = (
        'Padilla Bay CTD dataset. '
        'instrument dimension: sbe=SeaBird SBE 19plus, rbr=RBR, castaway=CastAway CTD. '
        'Use ds.sel(instrument="rbr") etc. for cross-instrument comparison.'
    )
    ds.attrs['created'] = pd.Timestamp.now().isoformat()

    return ds
