"""
Data loading functions for SeaBird and RBR CTD data.
"""
import numpy as np
import pandas as pd
from pathlib import Path
import ctd
from .utils import CastFrame
from . import config


def get_cast(file):
    """
    Load cast(s) from either CNV (SeaBird) or XLSX (RBR) file.

    Always returns a list of CastFrames so callers can iterate uniformly
    regardless of instrument type.

    Args:
        file: Path to CNV or XLSX file

    Returns: List of CastFrame objects with cast_meta populated
    """
    if file.endswith('.cnv'):
        cast = sbe_cast(file)
        cast['time'] = cast.cast_meta['time'] + pd.to_timedelta(cast['timeS'], unit='s')
        return [cast]
    elif file.endswith('.xlsx'):
        return rbr_cast(file)
    else:
        raise ValueError(f"Unsupported file type: {file}")



def station_from_path(file):
    """
    Infer station name from a file's stem by matching against known station names.

    Sorts station names longest-first so that multi-character names (e.g. 'G1')
    are matched before single-character ones (e.g. 'B') when a stem starts with
    both (e.g. 'G1NTS' → 'G1', not 'G').

    Args:
        file: Path to any CTD file whose stem encodes the station (e.g. 'S1NTS.cnv')

    Returns: Station name string, or None if no match found.
    """
    stations_sorted = sorted(config.STATIONS_DF['name'], key=len, reverse=True)
    stem = Path(file).stem  # e.g. 'BNTS', 'G1NTS', 'S1'
    return next((s for s in stations_sorted if stem.startswith(s)), None)


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
    station_name = Path(bl_file).stem  # e.g., 'S1' from 'S1.bl'
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

