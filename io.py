"""
Data loading functions for SeaBird and RBR CTD data.
"""
from traceback import print_tb
import numpy as np
import pandas as pd
from pathlib import Path
import ctd
from .utils import compute_density
from .processing import remove_surface_noise
from .config import CAST_MAP, STATIONS_FILE


def get_cast(file):
    """
    Load a cast from either CNV or XLSX file format.
    
    Args:
        file: Path to CNV (SeaBird) or XLSX (RBR) file
        station_file: CSV file with station names (for RBR processing)
    
    Returns: Cast DataFrame(s)
    """
    stations_df = pd.read_csv(STATIONS_FILE)
    if file.endswith('.cnv'):
        
        return sbe_cast(file, stations_df)
    elif file.endswith('.xlsx'):
        
        return rbr_cast(file, stations_df)
    else:
        raise ValueError(f"Unsupported file type: {file}")





def station_from_path(file, stations_df):
    stations_sorted = sorted(stations_df['name'], key=len, reverse=True)
    print(stations_sorted)
    # e.g. ['G1', 'G2', 'S1', 'S2', 'B']  — multi-char before single-char
    stem = Path(file).stem  # 'BNTS', 'G1NTS', 'BB', 'G1', etc.
    return next((s for s in stations_sorted if stem.startswith(s)), None)

def sbe_cast(cnv_file, stations_df):
    """
    Load and split SeaBird CNV file into downcast.
    
    Returns: (down_df)
    """

    cast_df = ctd.from_cnv(cnv_file)
    station = station_from_path(cnv_file, stations_df)

    cast_df._metadata['station'] = station
    lat = stations_df[stations_df.name == station].lat.values[0]
    lon = stations_df[stations_df.name == station].lon.values[0]
    # Store metadata before processing
    cast_df._metadata['instrument_type'] = 'sbe'
    metadata = cast_df._metadata.copy()
    cast_df = compute_density(cast_df, lat=lat, lon=lon)
    # Drop unused sensor columns
    drop_cols = ['CStarAt0', 'CStarTr0', 'par', 'wetStar', 'sbeox0PS', 'v4', 'flag', 'scan', 'c0S/m', 'potemp090C']
    cast_df.drop(columns=[c for c in drop_cols if c in cast_df.columns], inplace=True)
    
    down_df, up_df = cast_df.split()
    
    # Restore metadata after split
    down_df._metadata = metadata
    up_df._metadata = metadata
    
    if len(up_df) == 0:
        raise RuntimeError("Cast split failed - check data quality")

    return down_df


def rbr_cast(excel_file, stations_df, recasts: dict[str, int] | list[int] = None):
    """
    Split RBR excel export into individual station casts, saving as parquet files.
    
    Args:
        excel_file: RBR excel export path
        stations_df: DataFrame with station names (must have 'name' column)
        recasts: Optional mapping of station names to good cast indices, e.g.
                 {'G1': 0, 'G2': 1, 'S1': 3}.  Can also be a plain list of
                 integer indices.  When provided, only those casts are returned.
    
    Returns: List of cast DataFrames
    """

    cols = ['Time', 'Temperature', 'Pressure', 'Depth', 'Salinity', 'Density anomaly']
    rbr_df = pd.read_excel(excel_file, sheet_name='Data', header=1, usecols=cols)
    profs_df = pd.read_excel(excel_file, sheet_name='Profile annotation', header=1, 
                              names=['start_t', 'end_t', 'lab', 'Type'])
    # Get atmospheric pressure from metadata
    metadata_df = pd.read_excel(excel_file, sheet_name='Metadata', header=None)
    indices = np.where(metadata_df == 'Atmospheric pressure')
    atm_pressure = float(metadata_df.iloc[indices[0][0] + 1, indices[1][0]])
    down_prof_df = profs_df[profs_df.Type == 'DOWN'].reset_index(drop=True)
    # rename columns to match sbe format
    rbr_df.columns = ['Time', 'tv290C', 'Pressure', 'depSM', 'sal00', 'density']
    # Auto-resolve recasts from CAST_MAP when not explicitly provided.
    # Uses the parent folder name (e.g. '2025Nov18') as the CAST_MAP key.
    # Nested entries (e.g. '2025Oct20') are matched by xlsx filename stem.
    if recasts is None:
        folder_key = Path(excel_file).parent.name
        entry = CAST_MAP.get(folder_key)
        if entry is not None:
            if entry and isinstance(next(iter(entry.values())), dict):
                stem = Path(excel_file).stem
                recasts = next(
                    (v for k, v in entry.items() if stem.endswith(k) or f'_{k}' in stem),
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
        df['Pressure'] = df['Pressure'] - atm_pressure
        df.index = df.Pressure
        # Need to remove some data before consistent data

        station = recast_map[i] if recast_map is not None else f'cast_{i}'
        df._metadata = {'atmospheric_pressure': atm_pressure, 'instrument_type': 'rbr',
                        'time': df.Time.iloc[0], 'station': station}
        dfs.append(df)

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


# Backward-compatible alias (formerly update_DOC_file in ctd_step1.py)
update_DOC_file = load_bottle_file

