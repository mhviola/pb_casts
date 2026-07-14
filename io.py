"""
Data loading functions for SeaBird, RBR, and CastAway CTD data.
"""
import re
import numpy as np
import pandas as pd
from pathlib import Path
import ctd
from .utils import CastFrame
from . import config


def get_cast(file):
    """
    Load cast(s) from a CNV (SeaBird), XLSX (RBR), or CSV (CastAway) file.

    Always returns a list of CastFrames so callers can iterate uniformly
    regardless of instrument type.

    Args:
        file: Path to CNV, XLSX, or CSV file

    Returns: List of CastFrame objects with cast_meta populated
    """
    if file.endswith('.cnv'):
        cast = sbe_cast(file)
        cast['time'] = cast.cast_meta['time'] + pd.to_timedelta(cast['timeS'], unit='s')
        return [cast]
    elif file.endswith('.xlsx'):
        return rbr_cast(file)
    elif file.endswith('.csv'):
        return [castaway_cast(file)]
    else:
        raise ValueError(f"Unsupported file type: {file}")



def station_from_path(file):
    """
    Infer station name from a file's stem by matching against known station names.

    Handles two naming conventions:
      - Legacy NTS files:  stem starts with the station name (e.g. 'G1NTS' → 'G1')
      - SBE raw files:     stem has a 'Station ' prefix (e.g. 'Station G1' → 'G1')

    Sorts station names longest-first so that multi-character names (e.g. 'G1')
    are matched before single-character ones (e.g. 'B').

    Args:
        file: Path to any CTD file whose stem encodes the station

    Returns: Station name string, or None if no match found.
    """
    stations_sorted = sorted(config.STATIONS_DF['name'], key=len, reverse=True)
    stem = Path(file).stem  # e.g. 'BNTS', 'G1NTS', 'S1', 'Station G1', 'Station B'

    # Strip 'Station ' prefix (case-insensitive) produced by Seasave raw files
    if stem.lower().startswith('station '):
        stem = stem[len('station '):]

    return next((s for s in stations_sorted if stem.lower().startswith(s.lower())), None)


def sbe_cast(cnv_file):
    """
    Load a SeaBird CNV file and return the downcast as a CastFrame.

    Uses python-ctd to parse the CNV, drops instrument-specific channels that
    are not used downstream (transmissometer, PAR, oxygen %, etc.), then splits
    into down/up casts.  The upcast is discarded; a zero-length upcast triggers
    a RuntimeError because it usually means the cast file is truncated or the
    split heuristic failed.

    Args:
        cnv_file: Path to a *NTS.cnv SeaBird file.

    Returns: CastFrame with cast_meta = {station, instrument_type, time}
    """
    cast_df = ctd.from_cnv(cnv_file)
    station = station_from_path(cnv_file)
    start_time = cast_df._metadata.get('time')
    meta = {'station': station, 'instrument_type': 'sbe', 'time': start_time}
    # Keep only the core sensor columns; any new channels added by SBE firmware
    # upgrades or different casts are silently discarded here.
    # depSM is excluded because depth is already the dataset's index/coordinate
    keep_cols = {'timeS', 'tv290C', 'sal00'}
    drop_cols = [c for c in cast_df.columns if c not in keep_cols]
    cast_df.drop(columns=drop_cols, inplace=True)

    down_df, up_df = cast_df.split()

    down_cast = CastFrame(down_df)
    down_cast.cast_meta = meta
    up_cast = CastFrame(up_df)
    up_cast.cast_meta = meta
    if len(up_cast) == 0:
        raise RuntimeError("Cast split failed - check data quality")

    return down_cast


def rbr_cast(excel_file, recasts: dict[str, int] | None = None):
    """
    Split an RBR Excel export into per-station downcast CastFrames.

    The RBR instrument records all casts in a single Excel file.  This function
    reads the 'Data', 'Profile annotation', and 'Metadata' sheets, isolates each
    DOWN cast segment, subtracts atmospheric pressure, and wraps each segment in
    a CastFrame.

    Args:
        excel_file: Path to an RBR Excel export.
        recasts: Dict mapping station names to downcast indices, e.g.
                 {'G1': 0, 'G2': 1, 'S1': 3}.  When None (default), the mapping
                 is looked up automatically from config.CAST_MAP using the parent
                 folder name as the key.  When CAST_MAP has no entry for the
                 folder, all detected DOWN casts are returned as 'cast_0',
                 'cast_1', …

    Returns: List of CastFrames, one per station, with cast_meta populated.
    """

    cols = ['Time', 'Temperature', 'Pressure', 'Depth', 'Salinity']
    rbr_df = pd.read_excel(excel_file, sheet_name='Data', header=1, usecols=cols)
    profs_df = pd.read_excel(excel_file, sheet_name='Profile annotation', header=1, 
                              names=['start_t', 'end_t', 'lab', 'Type'])
    # Get atmospheric pressure from metadata
    metadata_df = pd.read_excel(excel_file, sheet_name='Metadata', header=None)
    indices = np.where(metadata_df == 'Atmospheric pressure')
    atm_pressure = float(metadata_df.iloc[indices[0][0] + 1, indices[1][0]])
    down_prof_df = profs_df[profs_df.Type == 'DOWN'].reset_index(drop=True)
    # rename columns to match sbe format
    rbr_df.columns = ['Time', 'tv290C', 'Pressure', 'depSM', 'sal00']
    # Auto-resolve recasts from CAST_MAP when not explicitly provided.
    # Uses the parent folder name (e.g. '2025Nov18') as the CAST_MAP key.
    # Nested entries (e.g. '2025Oct20') are matched by xlsx filename stem.
    if recasts is None:
        folder_key = Path(excel_file).parent.name
        entry = config.CAST_MAP.get(folder_key)
        if entry is not None:
            if entry and isinstance(next(iter(entry.values())), dict):
                stem = Path(excel_file).stem
                # Sort longest key first so 'allbutG1' matches before 'G1'
                recasts = next(
                    (v for k, v in sorted(entry.items(), key=lambda x: len(x[0]), reverse=True)
                     if stem.endswith(k) or f'_{k}' in stem),
                    None,
                )
            else:
                recasts = entry

    # Normalize recasts into {cast_index: station_name} once before the loop.
    if isinstance(recasts, dict):
        if recasts and isinstance(next(iter(recasts.values())), dict):
            raise ValueError(
                "recasts is a nested dict (e.g. a multi-file CAST_MAP entry). "
                "Pass a specific sub-entry, e.g. CAST_MAP['2025Oct20']['allbutG1']."
            )
        recast_map = {v: k for k, v in recasts.items()}
    elif recasts is not None:
        raise TypeError("recasts must be a dict")
    else:
        recast_map = None

    # Extract downcast segments
    dfs = []
    for i in down_prof_df.index:
        if recast_map is not None and i not in recast_map:
            continue

        df = rbr_df[(rbr_df.Time >= down_prof_df.start_t[i]) & 
                    (rbr_df.Time <= down_prof_df.end_t[i])].copy()
        # Subtract atmospheric pressure so index represents gauge pressure (dbar)
        df['Pressure'] = df['Pressure'] - atm_pressure
        df.index = df.Pressure
        # Keep only core columns; drop Pressure (now the index) and depSM (redundant)
        df = df[['Time', 'tv290C', 'sal00']]

        station = recast_map[i] if recast_map is not None else f'cast_{i}'
        cast = CastFrame(df)
        cast.cast_meta = {'atmospheric_pressure': atm_pressure, 'instrument_type': 'rbr',
                          'time': df.Time.iloc[0], 'station': station}
        dfs.append(cast)

    return dfs


def _parse_castaway_header(filepath):
    """Parse key-value metadata from CastAway CSV header lines (% prefix)."""
    meta = {}
    with open(filepath) as f:
        for line in f:
            if not line.startswith('%'):
                break
            m = re.match(r'% (.+?),(.*)', line.strip())
            if m:
                meta[m.group(1).strip()] = m.group(2).strip()
    return meta


def _nearest_station(lat, lon):
    """Return the name of the nearest station to (lat, lon) using haversine distance."""
    if config.STATIONS_DF is None:
        raise RuntimeError(
            "Station coordinates not loaded. Call pb_casts.set_project_root() first."
        )
    R = 6371.0
    lat1, lon1 = np.radians(lat), np.radians(lon)
    lats2 = np.radians(config.STATIONS_DF['lat'].values)
    lons2 = np.radians(config.STATIONS_DF['lon'].values)
    a = (np.sin((lats2 - lat1) / 2) ** 2
         + np.cos(lat1) * np.cos(lats2) * np.sin((lons2 - lon1) / 2) ** 2)
    distances = 2 * R * np.arcsin(np.sqrt(a))
    return config.STATIONS_DF.iloc[int(np.argmin(distances))]['name']


def castaway_cast(csv_file):
    """
    Load a CastAway CSV file and return the downcast profile as a CastFrame.

    Reads the % comment header to extract GPS coordinates and cast time, finds
    the nearest known station from config.STATIONS_DF using haversine distance,
    then reads the data section and renames columns to the pb_casts standard:
        tv290C  = temperature (°C)
        sal00   = salinity (PSU)
        depSM   = depth (m)

    The DataFrame is indexed by Pressure (Decibar), consistent with sbe_cast()
    and rbr_cast(), so the full process_ctd() pipeline works on CastAway data.

    Args:
        csv_file: Path to a CastAway CSV export.

    Returns: CastFrame indexed by pressure with cast_meta containing:
             station, instrument_type='castaway', time, lat, lon.

    Raises:
        ValueError if the cast header says 'Sample type, Invalid'.
    """
    meta = _parse_castaway_header(csv_file)

    sample_type = meta.get('Sample type', '')
    if 'Invalid' in sample_type:
        raise ValueError(f"Cast is marked Invalid: {csv_file}")

    local_time_str = meta.get('Cast time (local)')
    local_time = pd.to_datetime(local_time_str) if local_time_str else None

    lat_str = meta.get('Start latitude')
    lon_str = meta.get('Start longitude')
    if not lat_str or not lon_str:
        raise ValueError(f"No GPS coordinates found in {csv_file}")
    lat, lon = float(lat_str), float(lon_str)

    station = _nearest_station(lat, lon)

    df = pd.read_csv(csv_file, comment='%')
    col_map = {
        'Pressure (Decibar)':                  'Pressure',
        'Temperature (Celsius)':               'tv290C',
        'Salinity (Practical Salinity Scale)': 'sal00',
    }
    df = df.rename(columns={k: v for k, v in col_map.items() if k in df.columns})
    df = df[['Pressure', 'tv290C', 'sal00']]
    df = df.set_index('Pressure')

    cast = CastFrame(df)
    cast.cast_meta = {
        'station': station,
        'instrument_type': 'castaway',
        'time': local_time,
        'lat': lat,
        'lon': lon,
    }
    return cast


def _castaway_raw_station(base_station, cast_date, cast_time, ds):
    """
    Return the raw station name including any repeat suffix (e.g. 'G1b') by
    matching the CastAway cast time to the closest RBR or SBE cast time stored
    in an xarray Dataset produced by create_ctd_dataset().

    The Dataset stores repeat casts under station=base_station with repeat=1,2,…
    This function finds the repeat whose reference-instrument cast time is closest
    to cast_time and converts that to a station suffix (repeat 1 → 'b', 2 → 'c',
    …).  Falls back to base_station (repeat=0) if no match is found.
    """
    if cast_time is None:
        return base_station

    stations_ds = ds.coords['station'].values if 'station' in ds.coords else []
    if base_station not in stations_ds:
        return base_station

    dates_ds = ds.coords['date'].values if 'date' in ds.coords else []
    date_matches = [d for d in dates_ds if pd.Timestamp(d).date() == cast_date]
    if not date_matches:
        return base_station

    date_val       = date_matches[0]
    cast_time_ts   = pd.Timestamp(cast_time)
    instruments_ds = ds.coords['instrument'].values if 'instrument' in ds.coords else []
    repeats_ds     = ds.coords['repeat'].values     if 'repeat'     in ds.coords else []

    best_repeat = 0
    best_diff   = float('inf')

    for repeat_val in repeats_ds:
        for inst in [i for i in ['rbr', 'sbe'] if i in instruments_ds]:
            try:
                t = ds['cast_time'].sel(
                    station=base_station, date=date_val,
                    repeat=repeat_val, instrument=inst,
                ).values
            except Exception:
                continue
            if pd.isnull(t):
                continue
            diff = abs((pd.Timestamp(t) - cast_time_ts).total_seconds())
            if diff < best_diff:
                best_diff   = diff
                best_repeat = int(repeat_val)

    if best_repeat == 0:
        return base_station
    return base_station + chr(ord('a') + best_repeat - 1)


def _print_no_depth(df: pd.DataFrame) -> None:
    """Print a summary of rows where no CTD depth could be assigned."""
    no_depth = df[df['ctd_source'] == 'no_depth']
    if no_depth.empty:
        return
    print(f"\n{'='*55}")
    print(f"  {len(no_depth)} DOC sample(s) still have no depth:")
    print(f"{'='*55}")
    for _, r in no_depth.iterrows():
        print(f"  {r['date']}  {r['station']}  bottle {int(r['bottle_number'])}")
    print(f"{'='*55}\n")


def build_castaway_doc(castaway_folder, doc_file, invalid_folder=None,
                       bl_doc_path=None, ctd_nc_file=None):
    """
    Combine CastAway profiles with DOC bottle data into a single DataFrame.

    For each valid CastAway cast the function extracts:
      - Surface CTD values (first profile row)  → merged with bottle 1 DOC
      - Bottom  CTD values (last  profile row)  → merged with bottle 5 DOC

    When two casts fall on the same (station, date) — e.g. a full-depth cast
    and a shallower mid-column cast — the deepest cast supplies bottles 1 and 5
    and the shallower cast supplies CTD values for any remaining intermediate
    DOC bottles (e.g. bottle 2).

    Data source priority (highest → lowest):

    **Tier 1 — SBE CNV + BL bottle-firing files** (``bl_doc_path``):
      Exact bottle-firing depth from the SBE scan-match, with SBE sal/temp.
      When available this always wins for that (date, station, bottle_number).

    **Tier 2 — CastAway depth + RBR sal/temp** (``ctd_nc_file``):
      The CastAway rode with the Van Dorn bottles and provides the deployment
      depth.  Salinity and temperature come from the RBR profile (shallowest
      row for bottle 1; nearest-neighbour or deepest row for bottle 5) which
      is more precisely calibrated.  Intermediate bottles keep CastAway
      sal/temp because no clean non-interpolated RBR value exists for a
      mid-column depth.  ``ctd_source = 'castaway'``.

    **Tier 3 — RBR surface/bottom only** (``ctd_nc_file``, no CastAway):
      When no CastAway exists for a (date, station), bottles 1 and 5 are
      filled with the shallowest/deepest valid RBR measurement (depth + sal +
      temp).  Intermediate bottles remain blank.  ``ctd_source = 'rbr'``.

    A ``ctd_source`` column is added to every row:
      - ``'bl_file'``  — Tier 1: SBE bottle-firing file (depth + sal + temp)
      - ``'castaway'`` — Tier 2: depth from CastAway, sal/temp from RBR (or
                         CastAway-only when no RBR is available)
      - ``'rbr'``      — Tier 3: all values from RBR surface/bottom; no CastAway
      - ``'no_depth'`` — no source available for this bottle

    A summary of ``'no_depth'`` rows is printed after processing.

    Args:
        castaway_folder: Path to the root CastAway folder.  The function looks
                         for CSV files in ``Viola_*`` subfolders first; if none
                         are found it searches the folder directly.
        doc_file:        Path to DOC Excel file (sheet 'Processed' with columns:
                         station, date, bottle_number, doc_conc).
        invalid_folder:  Optional path to move Invalid cast files into.  When
                         None (default) invalid casts are silently skipped.
        bl_doc_path:     Optional path to folder of pre-processed BL+DOC CSVs
                         (output of ``process_bl_doc_files()``).  When provided,
                         Tier 1 SBE bottle data is merged and takes priority.
        ctd_nc_file:     Optional path to the NetCDF produced by
                         ``create_ctd_dataset()``.  Enables Tier 2 (CastAway
                         depth + RBR sal/temp) and Tier 3 (RBR surface/bottom
                         for dates with no CastAway).

    Returns: DataFrame with columns:
        doc_conc, date, bottle_number, station, salinity, depth, temperature,
        local_time, ctd_source
    """
    castaway_folder = Path(castaway_folder)

    # Load xarray dataset once up-front (used for both station matching and
    # RBR sal/temp lookup later).
    ds_ctd = None
    if ctd_nc_file is not None:
        import xarray as xr
        ds_ctd = xr.open_dataset(ctd_nc_file)

    csv_files = sorted(castaway_folder.glob('Viola_*/*.csv'))
    if not csv_files:
        csv_files = sorted(castaway_folder.glob('*.csv'))

    rows = []
    for fpath in csv_files:
        header = _parse_castaway_header(str(fpath))
        sample_type = header.get('Sample type', '')
        if 'Invalid' in sample_type:
            if invalid_folder:
                import shutil
                Path(invalid_folder).mkdir(parents=True, exist_ok=True)
                shutil.move(str(fpath), str(Path(invalid_folder) / fpath.name))
            continue

        lat_str = header.get('Start latitude')
        lon_str = header.get('Start longitude')
        if not lat_str or not lon_str:
            continue

        local_time = pd.to_datetime(header.get('Cast time (local)'))
        lat, lon = float(lat_str), float(lon_str)
        base_station = _nearest_station(lat, lon)
        cast_date    = local_time.date()

        # Refine to the correct repeat (e.g. 'G1b') via timestamp matching.
        if ds_ctd is not None:
            station = _castaway_raw_station(base_station, cast_date, local_time, ds_ctd)
        else:
            station = base_station

        df = pd.read_csv(str(fpath), comment='%')
        df = df.rename(columns={
            'Depth (Meter)':                        'depSM',
            'Temperature (Celsius)':                'tv290C',
            'Salinity (Practical Salinity Scale)':  'sal00',
        })

        surface = df.iloc[0]
        bottom  = df.iloc[-1]
        rows.append({
            'file':           fpath.name,
            'local_time':     local_time,
            'date':           cast_date,
            'station':        station,
            'max_depth':      bottom['depSM'],
            'surface_sal':    surface['sal00'],
            'surface_depth':  surface['depSM'],
            'surface_temp':   surface['tv290C'],
            'bottom_sal':     bottom['sal00'],
            'bottom_depth':   bottom['depSM'],
            'bottom_temp':    bottom['tv290C'],
        })

    if not rows:
        return pd.DataFrame()

    cast_summary = pd.DataFrame(rows)

    # Load DOC; average replicates
    doc_df = pd.read_excel(doc_file, sheet_name='Processed')
    doc_df['date'] = pd.to_datetime(doc_df['date']).dt.date
    doc_df = doc_df.groupby(
        ['date', 'station', 'bottle_number'], as_index=False
    )['doc_conc'].mean()

    # Primary cast = deepest on that (station, date); secondary = shallower duplicate(s)
    cast_summary = cast_summary.sort_values('max_depth', ascending=False)
    primary   = cast_summary.drop_duplicates(subset=['date', 'station'], keep='first').copy()
    secondary = cast_summary[
        cast_summary.duplicated(subset=['date', 'station'], keep='first')
    ].copy()

    # ---------------------------------------------------------------------- #
    # Optional: build RBR profile lookup from the CTD NetCDF.               #
    # Key: (date.date(), station str) → dict with full valid depth/sal/temp  #
    # arrays plus convenience surface values.  Uses the first repeat with    #
    # valid data.                                                             #
    # ---------------------------------------------------------------------- #
    rbr_lookup: dict = {}
    if ds_ctd is not None:
        instruments = ds_ctd.coords['instrument'].values if 'instrument' in ds_ctd.coords else []
        stations_ds = ds_ctd.coords['station'].values    if 'station'    in ds_ctd.coords else []
        dates_ds    = ds_ctd.coords['date'].values       if 'date'       in ds_ctd.coords else []
        repeats_ds  = ds_ctd.coords['repeat'].values     if 'repeat'     in ds_ctd.coords else []
        depth_grid  = ds_ctd.coords['depth'].values      if 'depth'      in ds_ctd.coords else np.array([])
        # Support both new CF names (sal/temp) and legacy SBE names (sal00/tv290C)
        sal_var  = 'sal'   if 'sal'   in ds_ctd else ('sal00'  if 'sal00'  in ds_ctd else None)
        temp_var = 'temp'  if 'temp'  in ds_ctd else ('tv290C' if 'tv290C' in ds_ctd else None)
        if 'rbr' in instruments and sal_var and temp_var:
            for station_val in stations_ds:
                for date_val in dates_ds:
                    date_key = pd.Timestamp(date_val).date()
                    for repeat_val in sorted(repeats_ds):
                        prof  = ds_ctd.sel(
                            station=station_val, date=date_val,
                            repeat=repeat_val, instrument='rbr',
                        )
                        n = len(depth_grid)
                        sal  = prof[sal_var].values
                        temp = prof[temp_var].values
                        valid = ~(np.isnan(sal) | np.isnan(temp))
                        if valid.sum() < 1:
                            continue
                        cast_time = None
                        if 'cast_time' in ds_ctd:
                            t = prof['cast_time'].values
                            if not pd.isnull(t):
                                cast_time = pd.Timestamp(t)
                        rbr_lookup[(date_key, str(station_val))] = {
                            'depths':       depth_grid[valid],
                            'sal':          sal[valid],
                            'temp':         temp[valid],
                            'surface_sal':  float(sal[valid][0]),
                            'surface_temp': float(temp[valid][0]),
                            'cast_time':    cast_time,
                        }
                        break  # first repeat with data is enough
        ds_ctd.close()

    # Build CTD lookup keyed by (date, station, bottle_number)
    bottle1 = primary[['date', 'station', 'surface_sal', 'surface_depth',
                        'surface_temp', 'local_time']].rename(
        columns={'surface_sal': 'salinity', 'surface_depth': 'depth',
                 'surface_temp': 'temperature'}
    ).assign(bottle_number=1, ctd_source='castaway')

    bottle5 = primary[['date', 'station', 'bottom_sal', 'bottom_depth',
                        'bottom_temp', 'local_time']].rename(
        columns={'bottom_sal': 'salinity', 'bottom_depth': 'depth',
                 'bottom_temp': 'temperature'}
    ).assign(bottle_number=5, ctd_source='castaway')

    ctd_lookup = pd.concat([bottle1, bottle5], ignore_index=True)

    # Tier 2: override sal/temp with RBR values (depth stays from CastAway).
    # ctd_source remains 'castaway' — the depth still came from the CastAway.
    # Bottle 1: RBR shallowest valid measurement.
    # Bottle 5: if CastAway depth < RBR max depth, nearest-neighbour RBR value
    #           to the CastAway depth; otherwise RBR deepest measurement.
    if rbr_lookup:
        for idx, row in ctd_lookup.iterrows():
            key = (row['date'], str(row['station']))
            if key not in rbr_lookup:
                continue
            rbr = rbr_lookup[key]
            if row['bottle_number'] == 1:
                ctd_lookup.at[idx, 'salinity']    = rbr['surface_sal']
                ctd_lookup.at[idx, 'temperature'] = rbr['surface_temp']
            elif row['bottle_number'] == 5:
                castaway_depth = row['depth']
                rbr_max_depth  = float(rbr['depths'][-1])
                if castaway_depth < rbr_max_depth:
                    nearest = int(np.argmin(np.abs(rbr['depths'] - castaway_depth)))
                    ctd_lookup.at[idx, 'salinity']    = float(rbr['sal'][nearest])
                    ctd_lookup.at[idx, 'temperature'] = float(rbr['temp'][nearest])
                else:
                    ctd_lookup.at[idx, 'salinity']    = float(rbr['sal'][-1])
                    ctd_lookup.at[idx, 'temperature'] = float(rbr['temp'][-1])
            # ctd_source stays 'castaway' — depth is always from the CastAway

    # Secondary casts fill intermediate DOC bottles (typically bottle 2).
    # Sal/temp stay as CastAway; no clean non-interpolated RBR value exists for
    # mid-column depths.
    if not secondary.empty:
        extra_rows = []
        for _, sec in secondary.iterrows():
            mask = (
                (doc_df['date'] == sec['date']) &
                (doc_df['station'] == sec['station']) &
                (~doc_df['bottle_number'].isin([1, 5]))
            )
            for bn in doc_df.loc[mask, 'bottle_number'].values:
                extra_rows.append({
                    'date':        sec['date'],
                    'station':     sec['station'],
                    'bottle_number': bn,
                    'salinity':    sec['bottom_sal'],
                    'depth':       sec['bottom_depth'],
                    'temperature': sec['bottom_temp'],
                    'local_time':  sec['local_time'],
                    'ctd_source':  'castaway',
                })
        if extra_rows:
            ctd_lookup = pd.concat(
                [ctd_lookup, pd.DataFrame(extra_rows)], ignore_index=True
            )

    castaway_result = doc_df.merge(ctd_lookup, on=['date', 'station', 'bottle_number'], how='left')

    # Flag DOC rows with no CastAway cast (depth is NaN after the left-join)
    castaway_result['ctd_source'] = castaway_result['ctd_source'].fillna('no_depth')

    # Tier 3: no CastAway for this (station, date) — fill bottles 1 and 5
    # from the RBR shallowest/deepest measurements (depth + sal + temp).
    # Intermediate bottles have no clean non-interpolated value, so they
    # remain no_depth.
    if rbr_lookup:
        for idx, row in castaway_result[castaway_result['ctd_source'] == 'no_depth'].iterrows():
            m = re.fullmatch(r'([A-Z0-9]+)([a-z]?)', str(row['station']))
            base = m.group(1) if m else str(row['station'])
            key  = (row['date'], base)
            if key not in rbr_lookup:
                continue
            rbr = rbr_lookup[key]
            if row['bottle_number'] == 1:
                castaway_result.at[idx, 'depth']       = float(rbr['depths'][0])
                castaway_result.at[idx, 'salinity']    = float(rbr['surface_sal'])
                castaway_result.at[idx, 'temperature'] = float(rbr['surface_temp'])
                castaway_result.at[idx, 'ctd_source']  = 'rbr'
                if rbr['cast_time'] is not None:
                    castaway_result.at[idx, 'local_time'] = rbr['cast_time']
            elif row['bottle_number'] == 5:
                castaway_result.at[idx, 'depth']       = float(rbr['depths'][-1])
                castaway_result.at[idx, 'salinity']    = float(rbr['sal'][-1])
                castaway_result.at[idx, 'temperature'] = float(rbr['temp'][-1])
                castaway_result.at[idx, 'ctd_source']  = 'rbr'
                if rbr['cast_time'] is not None:
                    castaway_result.at[idx, 'local_time'] = rbr['cast_time']

    # If no BL/SBE data requested, report and return now
    if bl_doc_path is None:
        _print_no_depth(castaway_result)
        return castaway_result

    # Load pre-processed BL+DOC CSVs and normalise column names
    bl_doc_path = Path(bl_doc_path)
    bl_parts = []
    for csv_path in sorted(bl_doc_path.glob('*_DOC.csv')):
        df = pd.read_csv(str(csv_path))
        df['date'] = pd.to_datetime(df['time_local']).dt.date
        df = df.rename(columns={
            'sal00':      'salinity',
            'tv290C':     'temperature',
            'depSM':      'depth',
            'time_local': 'local_time',
        })
        keep = ['bottle_number', 'local_time', 'date', 'station',
                'salinity', 'depth', 'temperature', 'doc_conc']
        df = df[[c for c in keep if c in df.columns]]
        df['ctd_source'] = 'bl_file'
        bl_parts.append(df)

    if not bl_parts:
        return castaway_result

    bl_df_all = pd.concat(bl_parts, ignore_index=True)
    bl_df_all['date'] = pd.to_datetime(bl_df_all['date']).dt.date

    # Merge: stack CastAway/RBR and BL rows, then keep BL row when both cover
    # the same (date, station, bottle_number) — BL depths are most precise.
    combined = pd.concat([castaway_result, bl_df_all], ignore_index=True)
    # Sort so NaN-salinity rows come BEFORE real values so drop_duplicates
    # keep='last' retains the row with actual data.
    combined = combined.sort_values(
        ['date', 'station', 'bottle_number', 'salinity'],
        na_position='first',
    )
    combined = combined.drop_duplicates(
        subset=['date', 'station', 'bottle_number'], keep='last'
    ).sort_values(['date', 'station', 'bottle_number']).reset_index(drop=True)

    _print_no_depth(combined)
    return combined


def load_bottle_file(bl_file, doc_file=None):
    """
    Load bottle (.bl) file, extract CTD values from the paired CNV, and
    optionally merge with DOC data.

    ctd.from_bl() returns only timing information (bottle_number, time_local,
    startscan, endscan). This function additionally opens the paired .cnv file
    and matches each bottle's fire time to the closest CTD scan, extracting
    tv290C (temperature), sal00 (salinity), and depSM (depth) for each bottle.

    Args:
        bl_file:  Path to .bl file; a same-stem .cnv must exist alongside it.
        doc_file: Optional path to DOC Excel file (sheet 'Processed' with
                  columns: station, date, bottle_number, doc_conc).

    Returns: DataFrame with columns:
        bottle_number, time_local, startscan, endscan,
        tv290C, sal00, depSM  (from CNV; NaN if CNV unavailable),
        station, doc_conc     (if doc_file provided)
    """
    bl_df = ctd.from_bl(bl_file)
    # from_bl uses bottle_number as the index; promote it to a regular column
    bl_df = bl_df.reset_index()
    bl_df['bottle_number'] = bl_df['bottle_number'].astype(int)

    # Match each bottle fire time to the nearest CTD scan in the paired CNV.
    # CNV files may be named {stem}.cnv or {stem}NTS.cnv (SeaBird convention).
    bl_path = Path(bl_file)
    cnv_file = bl_path.with_suffix('.cnv')
    if not cnv_file.exists():
        cnv_nts = bl_path.with_name(bl_path.stem + 'NTS.cnv')
        if cnv_nts.exists():
            cnv_file = cnv_nts
    if cnv_file.exists():
        try:
            cnv_df = ctd.from_cnv(str(cnv_file))
            cast_start = cnv_df._metadata.get('time')
            if cast_start is not None and 'timeS' in cnv_df.columns:
                scan_times = (
                    pd.Timestamp(cast_start)
                    + pd.to_timedelta(cnv_df['timeS'].values, unit='s')
                )
                for idx, row in bl_df.iterrows():
                    nearest = int(np.argmin(np.abs(scan_times - pd.Timestamp(row['time_local']))))
                    for col in ('tv290C', 'sal00', 'depSM'):
                        if col in cnv_df.columns:
                            bl_df.at[idx, col] = float(cnv_df.iloc[nearest][col])
        except Exception as e:
            print(f"Warning: Could not extract CTD values from {cnv_file.name}: {e}")

    if doc_file is None:
        bl_df['doc_conc'] = np.nan
        return bl_df

    # Load and filter DOC data
    doc_df = pd.read_excel(doc_file, sheet_name='Processed')
    station_name = station_from_path(bl_file) or Path(bl_file).stem
    bl_date = pd.to_datetime(bl_df.iloc[0]['time_local']).date()

    station_doc = doc_df[doc_df['station'] == station_name].copy()
    if 'date' in station_doc.columns:
        station_doc['date'] = pd.to_datetime(station_doc['date']).dt.date
        station_doc = station_doc[station_doc['date'] == bl_date]

    if not station_doc.empty:
        bl_df = pd.merge(bl_df, station_doc[['bottle_number', 'doc_conc']],
                         on='bottle_number', how='left')
    else:
        bl_df['doc_conc'] = np.nan

    return bl_df

