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
from .io import get_cast
from .processing import process_ctd, COMMON_PARAMS, SBE_PARAMS, RBR_PARAMS
from .utils import parse_folder_date


# CF-convention variable metadata: long_name, units, standard_name
VAR_ATTRS: dict[str, dict] = {
    'tv290C':   {'long_name': 'Temperature (ITS-90)',          'units': 'degC',   'standard_name': 'sea_water_temperature'},
    'sal00':    {'long_name': 'Practical Salinity',            'units': 'PSU',    'standard_name': 'sea_water_practical_salinity'},
    'prdM':     {'long_name': 'Pressure',                      'units': 'dbar',   'standard_name': 'sea_water_pressure'},
    'Pressure': {'long_name': 'Pressure',                      'units': 'dbar',   'standard_name': 'sea_water_pressure'},
    'depSM':    {'long_name': 'Depth',                         'units': 'm',      'standard_name': 'depth'},
    'sigma0':   {'long_name': 'Potential Density Anomaly σ₀',  'units': 'kg m-3', 'standard_name': 'sea_water_sigma_theta'},
    'rho0':     {'long_name': 'In Situ Density Anomaly',       'units': 'kg m-3', 'standard_name': 'sea_water_density'},
    'density':  {'long_name': 'Density Anomaly',               'units': 'kg m-3', 'standard_name': 'sea_water_sigma_theta'},
}

# Processing parameters — derived directly from processing.py constants
PROCESSING_ATTRS: dict = {**COMMON_PARAMS, **SBE_PARAMS, **RBR_PARAMS}


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


def create_ctd_dataset():
    """
    Create xarray Dataset from all CNV (SBE) and XLSX (RBR) files in DATA_PATH.

    Returns: Dataset with dimensions (station, date, repeat, depth).
      - repeat=0  primary cast, repeat=1 first repeat ('b'), repeat=2 second ('c'), …
      - cast_time coordinate holds the actual timestamp for each (station, date, repeat)
      - Missing combinations and depths outside a cast's range are filled with NaN.
    """
    data_path = config.DATA_PATH
    # key: (base_station, date, repeat_index)  value: processed DataFrame
    data_dict: dict[tuple[str, pd.Timestamp, int], pd.DataFrame] = {}
    cast_times: dict[tuple[str, pd.Timestamp, int], pd.Timestamp] = {}

    files = (
        glob.glob(str(data_path / '**' / '*.cnv'), recursive=True) +
        glob.glob(str(data_path / '**' / '*.xlsx'), recursive=True)
    )

    for fpath in files:
        date = parse_folder_date(Path(fpath).parent.name)
        if date is None:
            print(f"Warning: Could not parse date from {fpath}, skipping")
            continue

        try:
            for cast in get_cast(fpath):
                raw_station = cast.cast_meta.get('station')
                if raw_station is None:
                    print(f"Station not found in cast from {fpath}, skipping")
                    continue

                base, repeat = _parse_station_repeat(raw_station)
                key = (base, date, repeat)
                label = raw_station if repeat == 0 else f"{base} repeat {repeat}"
                print(f"Processing: {label} on {date.date()}")

                proc_df = process_ctd(cast)
                data_dict[key]  = proc_df
                cast_times[key] = cast.cast_meta.get('time')

        except Exception as e:
            print(f"Error processing {fpath}: {e}")
            continue

    if not data_dict:
        raise ValueError("No valid CTD data found!")

    # Build common 0.25 m depth grid spanning all casts
    bin_delta  = 0.25
    all_depths = np.concatenate([df.index.values for df in data_dict.values()])
    depth_grid = np.arange(
        np.floor(all_depths.min() * 4) / 4,
        np.ceil(all_depths.max() * 4) / 4 + bin_delta,
        bin_delta,
    )
    depth_grid = np.round(depth_grid, 4)

    stations = sorted({s for s, _, _ in data_dict})
    dates    = sorted({d for _, d, _ in data_dict})
    repeats  = sorted({r for _, _, r in data_dict})
    n_s, n_d, n_r, n_z = len(stations), len(dates), len(repeats), len(depth_grid)
    s_idx = {s: i for i, s in enumerate(stations)}
    d_idx = {d: i for i, d in enumerate(dates)}
    r_idx = {r: i for i, r in enumerate(repeats)}

    exclude_vars = {'time_local', 'timeS', 'time'}
    first_proc   = next(iter(data_dict.values()))
    var_names    = [c for c in first_proc.columns if c not in exclude_vars]

    # 4D arrays — NaN where a cast doesn't exist or doesn't reach a depth
    arrays     = {var: np.full((n_s, n_d, n_r, n_z), np.nan) for var in var_names}
    time_array = np.full((n_s, n_d, n_r), np.datetime64('NaT'), dtype='datetime64[ns]')

    for (base, date, repeat), proc_df in data_dict.items():
        si, di, ri = s_idx[base], d_idx[date], r_idx[repeat]
        cast_depths = proc_df.index.values
        in_range    = (depth_grid >= cast_depths.min()) & (depth_grid <= cast_depths.max())
        for var in var_names:
            if var in proc_df.columns:
                arrays[var][si, di, ri, in_range] = np.interp(
                    depth_grid[in_range], cast_depths, proc_df[var].values
                )
        t = cast_times.get((base, date, repeat))
        if t is not None:
            time_array[si, di, ri] = np.datetime64(pd.Timestamp(t), 'ns')

    # Station lat/lon coords from STATIONS_DF
    station_coords = {}
    if config.STATIONS_DF is not None:
        station_coords = {
            row['name']: (row['lat'], row['lon'])
            for _, row in config.STATIONS_DF.iterrows()
        }

    coords = {
        'station': stations,
        'date':    dates,
        'repeat':  repeats,
        'depth':   depth_grid,
        'cast_time': (['station', 'date', 'repeat'], time_array),
    }
    if station_coords:
        coords['latitude']  = ('station', [station_coords.get(s, (np.nan,))[0]        for s in stations])
        coords['longitude'] = ('station', [station_coords.get(s, (np.nan, np.nan))[1] for s in stations])

    ds = xr.Dataset(
        {var: (['station', 'date', 'repeat', 'depth'], arr) for var, arr in arrays.items()},
        coords=coords,
    )

    # Variable-level metadata
    for var in ds.data_vars:
        if var in VAR_ATTRS:
            ds[var].attrs.update(VAR_ATTRS[var])

    # Dataset-level metadata
    ds.attrs.update(PROCESSING_ATTRS)
    ds.attrs['description'] = 'Padilla Bay CTD dataset'
    ds.attrs['created']     = pd.Timestamp.now().isoformat()

    return ds

