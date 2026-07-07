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
    # Drop unused sensor columns
    drop_cols = ['CStarAt0', 'CStarTr0', 'par', 'wetStar', 'sbeox0PS', 'v4', 'flag', 'scan', 'c0S/m', 'potemp090C']
    cast_df.drop(columns=[c for c in drop_cols if c in cast_df.columns], inplace=True)

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
        'Pressure (Decibar)':                                   'Pressure',
        'Depth (Meter)':                                        'depSM',
        'Temperature (Celsius)':                                'tv290C',
        'Salinity (Practical Salinity Scale)':                  'sal00',
        'Conductivity (MicroSiemens per Centimeter)':           'conductivity',
        'Specific conductance (MicroSiemens per Centimeter)':   'specific_conductance',
        'Sound velocity (Meters per Second)':                   'sound_velocity',
        'Density (Kilograms per Cubic Meter)':                  'density',
    }
    df = df.rename(columns={k: v for k, v in col_map.items() if k in df.columns})
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


def build_castaway_doc(castaway_folder, doc_file, invalid_folder=None,
                       bl_doc_path=None):
    """
    Combine CastAway profiles with DOC bottle data into a single DataFrame.

    For each valid CastAway cast the function extracts:
      - Surface CTD values (first profile row)  → merged with bottle 1 DOC
      - Bottom  CTD values (last  profile row)  → merged with bottle 5 DOC

    When two casts fall on the same (station, date) — e.g. a full-depth cast
    and a shallower mid-column cast — the deepest cast supplies bottles 1 and 5
    and the shallower cast supplies CTD values for any remaining intermediate
    DOC bottles (e.g. bottle 2).

    When ``bl_doc_path`` is supplied the function also loads the pre-processed
    BL+DOC CSV files produced by ``process_bl_doc_files()`` (columns:
    bottle_number, time_local, sal00, tv290C, depSM, doc_conc, station).
    These SBE-derived rows are combined with the CastAway rows; where both
    sources cover the same (date, station, bottle_number) the BL/SBE row is
    kept because it has a more precise depth from the bottle-firing scan match.

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
                         SBE bottle data is merged in alongside CastAway data.

    Returns: DataFrame with columns:
        doc_conc, date, bottle_number, station, salinity, depth, temperature,
        local_time
    """
    castaway_folder = Path(castaway_folder)

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
        station = _nearest_station(lat, lon)
        cast_date = local_time.date()

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

    # Build CTD lookup keyed by (date, station, bottle_number)
    bottle1 = primary[['date', 'station', 'surface_sal', 'surface_depth',
                        'surface_temp', 'local_time']].rename(
        columns={'surface_sal': 'salinity', 'surface_depth': 'depth',
                 'surface_temp': 'temperature'}
    ).assign(bottle_number=1)

    bottle5 = primary[['date', 'station', 'bottom_sal', 'bottom_depth',
                        'bottom_temp', 'local_time']].rename(
        columns={'bottom_sal': 'salinity', 'bottom_depth': 'depth',
                 'bottom_temp': 'temperature'}
    ).assign(bottle_number=5)

    ctd_lookup = pd.concat([bottle1, bottle5], ignore_index=True)

    # Secondary casts fill intermediate DOC bottles (typically bottle 2)
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
                    'date': sec['date'], 'station': sec['station'],
                    'bottle_number': bn,
                    'salinity':    sec['bottom_sal'],
                    'depth':       sec['bottom_depth'],
                    'temperature': sec['bottom_temp'],
                    'local_time':  sec['local_time'],
                })
        if extra_rows:
            ctd_lookup = pd.concat(
                [ctd_lookup, pd.DataFrame(extra_rows)], ignore_index=True
            )

    castaway_result = doc_df.merge(ctd_lookup, on=['date', 'station', 'bottle_number'], how='left')

    # If no BL/SBE data requested, return CastAway result as-is
    if bl_doc_path is None:
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
        bl_parts.append(df[[c for c in keep if c in df.columns]])

    if not bl_parts:
        return castaway_result

    bl_df_all = pd.concat(bl_parts, ignore_index=True)
    bl_df_all['date'] = pd.to_datetime(bl_df_all['date']).dt.date

    # Merge: stack CastAway and BL rows, then keep BL row when both cover the
    # same (date, station, bottle_number) — BL depths are more precise.
    combined = pd.concat([castaway_result, bl_df_all], ignore_index=True)
    # Sort so NaN-salinity rows (CastAway placeholders) come BEFORE real values;
    # drop_duplicates keep='last' then retains the row with actual data.
    combined = combined.sort_values(
        ['date', 'station', 'bottle_number', 'salinity'],
        na_position='first',
    )
    combined = combined.drop_duplicates(
        subset=['date', 'station', 'bottle_number'], keep='last'
    ).sort_values(['date', 'station', 'bottle_number']).reset_index(drop=True)

    return combined


def load_bottle_file(bl_file, doc_file=None):
    """
    Load bottle (.bl) file and optionally merge with DOC data.
    
    Args:
        bl_file: Path to .bl file
        doc_file: Optional path to DOC excel file (sheet 'Processed' with 
                  columns: station, date, bottle_number, doc_conc)
    
    Returns: DataFrame with bottle times and optionally DOC concentrations
    """
    bl_df = ctd.from_bl(bl_file)
    bl_df['bottle_number'] = bl_df['bottle_number'].astype(int)
    
    if doc_file is None:
        bl_df['doc_conc'] = np.nan
        return bl_df
    
    # Load and filter DOC data
    doc_df = pd.read_excel(doc_file, sheet_name='Processed')
    station_name = station_from_path(bl_file) or Path(bl_file).stem
    bl_date = bl_df.iloc[0]['time_local'].date()
    
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

