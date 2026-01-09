"""
Data loading functions for SeaBird and RBR CTD data.
"""
import numpy as np
import pandas as pd
from pathlib import Path
import ctd
from .processing import remove_pump_priming


def get_cast(file, station_file):
    """
    Load a cast from either CNV or XLSX file format.
    
    Args:
        file: Path to CNV (SeaBird) or XLSX (RBR) file
        station_file: CSV file with station names (for RBR processing)
    
    Returns: Cast DataFrame(s)
    """
    if file.endswith('.cnv'):
        return sbe_cast(file)
    elif file.endswith('.xlsx'):
        stations_df = pd.read_csv(station_file)
        return rbr_cast(file, stations_df)
    else:
        raise ValueError(f"Unsupported file type: {file}")


def sbe_cast(cnv_file):
    """
    Load and split SeaBird CNV file into downcast. Removes pump priming artifacts.
    
    Returns: (down_df)
    """
    cast_df = ctd.from_cnv(cnv_file)
    # Store metadata before processing
    cast_df._metadata['instrument_type'] = 'sbe'
    metadata = cast_df._metadata.copy()
    # Drop unused sensor columns
    drop_cols = ['CStarAt0', 'CStarTr0', 'par', 'wetStar', 'sbeox0PS', 'v4', 'flag']
    cast_df.drop(columns=[c for c in drop_cols if c in cast_df.columns], inplace=True)
    
    # Remove pre-pump data using stability method
    cast_df = remove_pump_priming(cast_df, method='stability', 
                                   gradient_threshold=1.0, n_stable=3, min_depth=0.25)
    
    # Preserve metadata before split
    down_df, up_df = cast_df.split()
    
    # Restore metadata after split
    down_df._metadata = metadata
    up_df._metadata = metadata
    
    if len(up_df) == 0:
        raise RuntimeError("Cast split failed - check data quality")

    return down_df


def rbr_cast(excel_file, stations_df, recasts: list[int] = None):
    """
    Split RBR excel export into individual station casts, saving as parquet files.
    
    Args:
        excel_file: RBR excel export path
        stations_df: DataFrame with station names (must have 'name' column)
        recasts: Optional list of cast indices to include (for recast filtering)
    
    Returns: List of cast DataFrames
    """
    cols = ['Time', 'Temperature', 'Pressure', 'Depth', 'Salinity']
    rbr_df = pd.read_excel(excel_file, sheet_name='Data', header=1, usecols=cols)
    profs_df = pd.read_excel(excel_file, sheet_name='Profile annotation', header=1, 
                              names=['start_t', 'end_t', 'lab', 'Type'])
    # Get atmospheric pressure from metadata
    metadata_df = pd.read_excel(excel_file, sheet_name='Metadata', header=None)
    indices = np.where(metadata_df == 'Atmospheric pressure')
    atm_pressure = float(metadata_df.iloc[indices[0][0] + 1, indices[1][0]])
    down_prof_df = profs_df[profs_df.Type == 'DOWN']

    # rename columns to match sbe format
    rbr_df.columns = ['Time', 'tv290C', 'Pressure', 'depSM', 'sal00']

    # Extract downcast segments
    dfs = []
    for i in down_prof_df.index:
        df = rbr_df[(rbr_df.Time >= down_prof_df.start_t[i]) & 
                    (rbr_df.Time <= down_prof_df.end_t[i])].copy()
        df['Pressure'] = df['Pressure'] - atm_pressure
        df.index = df.Pressure
        # Not really 'priming the pump' for rbr but need to remove some data before consistent data
        df = remove_pump_priming(df, method='stability', 
                                   gradient_threshold=1.0, n_stable=3)
        df._metadata = {'atmospheric_pressure': atm_pressure, 'instrument_type': 'rbr', 'time' : df.Time.iloc[0]}
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

