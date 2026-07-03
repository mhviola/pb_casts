"""
xarray Dataset creation from multiple CTD casts.
"""
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
from .config import DATA_PATH
from .io import sbe_cast
from .processing import process_ctd
from .utils import parse_folder_date


def create_ctd_dataset(station_coords_file=None, data_path=None):
    """
    Create xarray Dataset from all CNV files in DATA_PATH.
    
    Args:
        station_coords_file: Path to CSV with station coordinates (columns: station, lat, lon)
        data_path: Path to data directory (default: DATA_PATH from config)
    
    Returns: Dataset with dimensions (station_date) containing ragged depth arrays.
    """
    data_path = Path(data_path) if data_path else DATA_PATH
    
    # Load station coordinates
    coords_file = Path(station_coords_file) if station_coords_file else data_path / "station_coordinates.csv"
    station_coords = {}
    if coords_file.exists():
        station_df = pd.read_csv(coords_file)
        station_coords = {row['station']: (row['lat'], row['lon']) 
                         for _, row in station_df.iterrows()}
    
    cnv_files = list(data_path.glob("**/*.cnv"))
    rbr_files = list(data_path.glob("**/*.xlsx"))
    data_dict = {}
    
    for cnv_file in cnv_files:
        try:
            # Parse station and date
            station = cnv_file.stem.replace('NTS', '')
            date_str = cnv_file.parent.name
            date = parse_folder_date(date_str)
            
            if not station or not date:
                print(f"Warning: Could not parse {cnv_file}, skipping")
                continue
            
            bl_file = cnv_file.with_name(f"{station}.bl")
            if not bl_file.exists():
                print(f"Warning: No BL file for {cnv_file.name}, skipping")
                continue
            
            print(f"Processing: {station} on {date.date()}")
            
            down_df = sbe_cast(cnv_file)
            proc_df = process_ctd(down_df)
            
            data_dict[(station, date)] = {
                'proc': proc_df,
                'station': station,
                'date': date
            }
            
        except Exception as e:
            print(f"Error processing {cnv_file}: {e}")
            continue
    # for rbr_file in rbr_files:
    #     try:
            
            
    if not data_dict:
        raise ValueError("No valid CTD data found!")
    
    # Build xarray Dataset with ragged arrays
    station_date_pairs = list(data_dict.keys())
    station_date_index = pd.MultiIndex.from_tuples(station_date_pairs, names=['station', 'date'])
    
    first_proc = data_dict[station_date_pairs[0]]['proc']
    exclude_vars = ['time_local', 'timeS']
    
    # Collect depth and variable arrays
    depth_arrays = np.empty(len(station_date_pairs), dtype=object)
    data_vars = {var: np.empty(len(station_date_pairs), dtype=object) 
                 for var in first_proc.columns if var not in exclude_vars}
    
    for i, (station, date) in enumerate(station_date_pairs):
        proc = data_dict[(station, date)]['proc']
        depth_arrays[i] = proc.index.values
        for var in data_vars:
            data_vars[var][i] = proc[var].values if var in proc.columns else np.array([])
    
    # Create coordinates
    coords = {
        'station_date': station_date_index,
        'station': ('station_date', [s for s, d in station_date_pairs]),
        'date': ('station_date', [d for s, d in station_date_pairs]),
        'depth': (['station_date'], depth_arrays),
    }
    
    if station_coords:
        coords['latitude'] = ('station_date', [station_coords.get(s, (np.nan,))[0] for s, d in station_date_pairs])
        coords['longitude'] = ('station_date', [station_coords.get(s, (np.nan, np.nan))[1] for s, d in station_date_pairs])
    
    ds = xr.Dataset({k: (['station_date'], v) for k, v in data_vars.items()}, coords=coords)
    ds.attrs['description'] = 'Padilla Bay CTD dataset'
    ds.attrs['created'] = pd.Timestamp.now().isoformat()
    
    return ds

