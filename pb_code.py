"""
Padilla Bay CTD Processing Module

Functions for processing SeaBird and RBR CTD data, bottle files, and DOC data.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
import os
import glob
import ctd
import gsw
import xarray as xr
import pyarrow as pa
import pyarrow.parquet as pq

# Project paths
PROJECT_ROOT = Path("/Users/marisaviola/Library/CloudStorage/OneDrive-WesternWashingtonUniversity/Padilla_Bay_Project")
DATA_PATH = PROJECT_ROOT / "Data"
OUTPUT_PATH = PROJECT_ROOT / "Output"
DOC_FILE = PROJECT_ROOT / "DOC_info" / "DOCdepth_profiles.xlsx"

def detect_precision(series, max_decimals=6):
    """Detect decimal places in a data series. Returns max precision found."""
    sample = series.dropna().head(100)
    if len(sample) == 0:
        return 4
    
    decimals = []
    for val in sample:
        val_str = f"{val:.{max_decimals}f}".rstrip('0')
        if '.' in val_str:
            decimals.append(len(val_str.split('.')[1]))
        else:
            decimals.append(0)
    
    return max(decimals) if decimals else 4

def remove_pump_priming(cast_df, method='stability', depth_threshold=0.5, 
                        gradient_threshold=1.0, n_stable=3, min_depth=0.25):
    """
    Remove pre-pump data before CTD conductivity cell is flushed.
    
    Methods: 'depth' (simple cutoff), 'stability' (n_stable consecutive points 
    with gradient < threshold), 'combined' (both). Always enforces min_depth.
    """
    # Always enforce minimum depth as safety net
    working_df = cast_df[cast_df.index >= min_depth].copy()
    if len(working_df) == 0:
        print(f"Warning: No data below min_depth={min_depth}m")
        return cast_df
    
    if method == 'depth':
        # Simple depth cutoff
        return cast_df[cast_df.index >= depth_threshold]
    
    elif method == 'stability':
        # Find where salinity stabilizes over n_stable consecutive points
        sal_gradient = working_df['sal00'].diff().abs()
        
        # Create rolling check: all of the last n_stable gradients must be below threshold
        is_stable = sal_gradient < gradient_threshold
        # Rolling sum of stable points - when it equals n_stable, we have n_stable consecutive stable readings
        stable_run = is_stable.rolling(window=n_stable, min_periods=n_stable).sum()

        # Find first index where we have n_stable consecutive stable readings
        stable_indices = stable_run[stable_run == n_stable].index
        
        if len(stable_indices) > 0:
            # Start from (n_stable - 1) points before the first fully stable window
            # This is the beginning of the stable region
            first_stable_window_end = stable_indices[0]
            first_stable_window_start_pos = working_df.index.get_loc(first_stable_window_end) - (n_stable - 1)
            first_stable_window_start_pos = max(0, first_stable_window_start_pos)
            first_stable_idx = working_df.index[first_stable_window_start_pos]
            return working_df[working_df.index >= first_stable_idx]
        else:
            # No stable region found - fall back to depth threshold
            print(f"Warning: No stable region found, using depth_threshold={depth_threshold}m")
            return cast_df[cast_df.index >= depth_threshold]
    
    elif method == 'combined':
        # Apply depth threshold first, then check for stability
        depth_filtered = cast_df[cast_df.index >= depth_threshold].copy()
        
        if len(depth_filtered) == 0:
            print(f"Warning: No data below depth_threshold={depth_threshold}m")
            return cast_df
        
        sal_gradient = depth_filtered['sal00'].diff().abs()
        is_stable = sal_gradient < gradient_threshold
        stable_run = is_stable.rolling(window=n_stable, min_periods=n_stable).sum()
        stable_indices = stable_run[stable_run == n_stable].index
        
        if len(stable_indices) > 0:
            first_stable_window_end = stable_indices[0]
            first_stable_window_start_pos = depth_filtered.index.get_loc(first_stable_window_end) - (n_stable - 1)
            first_stable_window_start_pos = max(0, first_stable_window_start_pos)
            first_stable_idx = depth_filtered.index[first_stable_window_start_pos]
            return depth_filtered[depth_filtered.index >= first_stable_idx]
        else:
            # Already depth filtered, return that
            return depth_filtered
    
    else:
        print(f"Warning: Unknown method '{method}'. Using depth threshold only.")
        return cast_df[cast_df.index >= depth_threshold]
def get_cast(file, station_file):
    if file.endswith('.cnv'):
        return sbe_cast(file)
    elif file.endswith('.xlsx'):
        stations_df = pd.read_csv(station_file)
        return rbr_cast(file, stations_df)
    else:
        raise ValueError(f"Unsupported file type: {file}")

def rbr_cast(excel_file, stations_df, recasts: list[int] = None):
    """
    Split RBR excel export into individual station casts, saving as parquet files.
    
    Args:
        excel_file: RBR excel export path
        stations_df: DataFrame with station names (must have 'name' column)
        recasts: Optional list of cast indices to include (for recast filtering)
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
    # # Save to parquet
    # folder = Path(excel_file).parent
    # count = 0
    # for idx, df in enumerate(dfs):
    #     if recasts is None or idx in recasts:
    #         station = stations_df.name.iloc[count]
    #         print(f"{station} ({count})")
    #         df.to_parquet(folder / f"{station}.parquet")
    #         count += 1



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

def process_ctd(df, smooth=True, columns=None):
    """
    Process CTD data: despike, lp_filter, press_check, interpolate, bin, smooth.
    
    Args:
        df: CTD DataFrame with depth/pressure as index
        SHOULD BE IN METADATA:
            instrument_type: 'sbe' (4 Hz) or 'rbr' (8 Hz)
            time: Cast start datetime
        smooth: Apply hanning smoothing (default True)
        columns: Columns to process (default: all numeric except timeS, bottle cols)
    
    Returns: Processed DataFrame with 0.25m depth bins
    """
    # Determine columns to process
    exclude_cols = ['timeS', 'bottle_number', 'doc_conc']
    # saves datetime column if there is one
    date_cols = [col for col in df.columns if np.issubdtype(df[col].dtype, np.datetime64)]
    num_date_cols = len(date_cols)
    date_cols_series = [df[col] for col in date_cols]
    if num_date_cols == 0:
        datetime_df = None
    elif num_date_cols == 1:
        datetime_df = df[date_cols[0]]
    else:
        raise ValueError(f"Expected exactly one datetime column, found {num_date_cols}: {date_cols}")

    if columns is None:
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        columns = [col for col in numeric_cols if col not in exclude_cols]
    
    cols_to_process = [col for col in columns if col in df.columns]
    bottle_cols = ['bottle_number', 'doc_conc']
    
    # Store bottle columns before processing
    bottle_data = {col: df[col].copy() for col in bottle_cols if col in df.columns}
    # Set parameters based on instrument type
    if df._metadata['instrument_type'].lower() == 'rbr':
        # RBR CTD 
        sample_rate = 8.0
        despike_block = 100
        smooth_window = 5
    elif df._metadata['instrument_type'].lower() == 'sbe':
        # SBE 19plus V2 (4 Hz = 0.25s interval)
        sample_rate = 4.0
        despike_block = 75
        smooth_window = 11
    else:
        raise ValueError(f"Invalid instrument type: {df._metadata['instrument_type']}")    
    # Common parameters optimized for Padilla Bay halocline studies
    bin_delta = 0.25
    
    # Process the selected columns (without bindata first to preserve surface data)
    proc_df = (
        df[cols_to_process]
        .despike(n1=2, n2=20, block=despike_block)
        .lp_filter(sample_rate=sample_rate, time_constant=0.15)
        .press_check()
        .interpolate(method="index", limit_direction="both", limit_area="inside")
    )
    
    # Custom binning that preserves near-surface data
    # Create depth grid starting from 0 (or min depth if deeper)
    min_depth = max(0.0, np.floor(proc_df.index.min() * 4) / 4)  # Round down to nearest 0.25m
    max_depth = np.ceil(proc_df.index.max() * 4) / 4  # Round up to nearest 0.25m
    new_index = np.arange(min_depth, max_depth + bin_delta, bin_delta)
    
    # Interpolate to new depth grid
    binned_data = {}
    for col in proc_df.columns:
        binned_data[col] = np.interp(new_index, proc_df.index.values, proc_df[col].values)
    
    proc_df = pd.DataFrame(binned_data, index=new_index)
    proc_df.index.name = df.index.name or 'depth'

    # Apply smoothing if requested
    if smooth:
        proc_df = proc_df.smooth(window_len=smooth_window, window="hanning")
    
    # Round all values to match original data precision (auto-detect from raw data)
    for col in proc_df.columns:
        if col in df.columns:
            # Detect precision from original unprocessed data
            precision = detect_precision(df[col])
            proc_df[col] = np.round(proc_df[col], precision)
    # Add timeS by interpolating from original data (time wasn't processed, just kept separate)
    if 'timeS' in df.columns:
        # Interpolate timeS to new depth grid
        proc_df['timeS'] = np.interp(new_index, df.index.values, df['timeS'].values)
        
        # Round to instrument sampling precision
        time_resolution = 1.0 / sample_rate  # 0.25s for 4Hz, 0.125s for 8Hz
        proc_df['timeS'] = np.round(proc_df['timeS'] / time_resolution) * time_resolution
        proc_df['time_local'] = df._metadata['time'] + pd.to_timedelta(proc_df['timeS'], unit='s')
    elif datetime_df is not None:
        time_resolution = 1.0 / sample_rate
        # Convert datetime to numeric (timestamps), interpolate, then convert back
        datetime_numeric = datetime_df.astype(np.int64) / 1e9  # Convert to seconds since epoch
        interpolated_numeric = np.interp(new_index, df.index.values, datetime_numeric)
        rounded_numeric = np.round(interpolated_numeric / time_resolution) * time_resolution
        proc_df['time_local'] = pd.to_datetime(rounded_numeric, unit='s')
# I had a better way of doing this but lost it somewhere :/ this works

    # Re-add bottle columns using nearest neighbor interpolation
    # This preserves bottle values at their original depths and NaNs elsewhere
    for col in bottle_cols:
        if col in bottle_data:
            # Use nearest neighbor to map bottle data to processed depth grid
            # This preserves values only at bottle depths
            proc_df[col] = bottle_data[col].reindex(proc_df.index, method='nearest')
            # Ensure values are NaN where they weren't originally present
            # Find original non-NaN indices
            orig_non_nan = bottle_data[col].dropna()
            if len(orig_non_nan) > 0:
                # For each processed depth, check if it's close to an original bottle depth
                # If not close enough, set to NaN (threshold: half of bindata delta = 0.125m)
                orig_indices = orig_non_nan.index.values
                proc_indices = proc_df.index.values
                # Vectorized distance calculation
                for i, proc_idx in enumerate(proc_indices):
                    distances = np.abs(orig_indices - proc_idx)
                    min_dist = np.min(distances)
                    # If closest bottle is more than 0.125m away, set to NaN
                    if min_dist > bin_delta / 2:
                        proc_df.iloc[i, proc_df.columns.get_loc(col)] = np.nan
    
        
    return proc_df


def make_parquet(down_df, bl_down_df=None):
    """
    Save CTD data to parquet with metadata (CTD/DOC months included).
    Saves to PROJECT_ROOT/parquets/ with today's date in filename.
    """

    
    parquets_dir = PROJECT_ROOT / "parquets"
    parquets_dir.mkdir(exist_ok=True)
    today = datetime.now().strftime("%Y%m%d")
    
    # Extract months from CTD data
    if 'time_local' in down_df.columns:
        ctd_months = sorted(down_df['time_local'].dt.to_period('M').unique())
        ctd_months_str = [str(m) for m in ctd_months]
    else:
        ctd_months_str = []
    
    # Extract months from DOC data (if bl_down_df provided and has doc_conc with non-null values)
    doc_months_str = []
    if bl_down_df is not None and 'time_local' in bl_down_df.columns:
        if 'doc_conc' in bl_down_df.columns:
            # Only consider rows where doc_conc is not null
            doc_data = bl_down_df[bl_down_df['doc_conc'].notna()]
            if len(doc_data) > 0:
                doc_months = sorted(doc_data['time_local'].dt.to_period('M').unique())
                doc_months_str = [str(m) for m in doc_months]
    
    # Combine dataframes
    if bl_down_df is not None:
        # bl_down_df already contains CTD data at bottle times, so we need to merge carefully
        # Use bl_down_df for rows with bottle data, and down_df for rows without
        # First, identify which rows in down_df don't have corresponding bottle data
        bottle_times = set(bl_down_df['time_local'].values) if 'time_local' in bl_down_df.columns else set()
        down_df_no_bottles = down_df[~down_df['time_local'].isin(bottle_times)] if 'time_local' in down_df.columns else down_df
        
        # Combine: all CTD data without bottles + all bottle data (which includes CTD at those times)
        combined_df = pd.concat([down_df_no_bottles, bl_down_df], ignore_index=True, sort=False)
        # Sort by time_local if available
        if 'time_local' in combined_df.columns:
            combined_df = combined_df.sort_values('time_local').reset_index(drop=True)
    else:
        combined_df = down_df.copy()
    
    # Create filename and save
    filepath = parquets_dir / f"ctd_data_{today}.parquet"
    table = pa.Table.from_pandas(combined_df)
    
    # Create custom metadata
    custom_metadata = {
        b'ctd_months': (','.join(ctd_months_str) if ctd_months_str else 'None').encode('utf-8'),
        b'doc_months': (','.join(doc_months_str) if doc_months_str else 'None').encode('utf-8'),
        b'created_date': today.encode('utf-8'),
        b'description': b'Padilla Bay CTD and bottle data'
    }
    
    # Add metadata to table schema
    existing_metadata = table.schema.metadata or {}
    updated_metadata = {**existing_metadata, **custom_metadata}
    table = table.replace_schema_metadata(updated_metadata)
    
    # Write parquet file
    pq.write_table(table, filepath)
    
    print(f"Saved parquet file: {filepath}")
    print(f"CTD months included: {', '.join(ctd_months_str) if ctd_months_str else 'None'}")
    print(f"DOC months included: {', '.join(doc_months_str) if doc_months_str else 'None'}")
    
    return filepath
    
def plot_ts(down_df, bl_down_df=None):
    """Plot T-S diagram with density contours and salinity vs depth profile."""
    # Calculate axis ranges with padding
    smin = down_df.Salinity.min() * 0.99
    smax = down_df.Salinity.max() * 1.01
    tmin = down_df.Temperature.min() - 0.1 * down_df.Temperature.max()
    tmax = down_df.Temperature.max() * 1.1
    
    xdim = int(round((smax - smin) / 0.1 + 1))
    ydim = int(round(tmax - tmin + 1))
    
    # Create density grid (sigma-t)
    ti = np.linspace(1, ydim - 1, ydim) + tmin
    si = np.linspace(1, xdim - 1, xdim) * 0.1 + smin
    dens = np.zeros((ydim, xdim))
    for j in range(ydim):
        for i in range(xdim):
            dens[j, i] = gsw.rho(si[i], ti[j], 0) - 1000

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 10))
    
    # T-S diagram
    CS = ax1.contour(si, ti, dens, linestyles='dashed', colors='k')
    ax1.clabel(CS, fontsize=12, inline=True)
    ax1.plot(down_df.Salinity, down_df.Temperature, '.')
    ax1.set_xlabel('Salinity (PSU)')
    ax1.set_ylabel('Temperature (°C)')
    ax1.set_title('T-S Diagram')
    
    # Salinity profile
    ax2.plot(down_df.Salinity, down_df.Depth)
    ax2.set_xlabel('Salinity (PSU)')
    ax2.set_ylabel('Depth (m)')
    ax2.set_title('Salinity vs Depth')
    ax2.tick_params(axis='x', rotation=45)
    ax2.invert_yaxis()
    
    # Add bottle annotations
    if bl_down_df is not None:
        ax1.plot(bl_down_df.Salinity, bl_down_df.Temperature, '.')
        ax2.plot(bl_down_df.Salinity, bl_down_df.Depth, '.')
        for _, r in bl_down_df.iterrows():
            for ax, y_col in [(ax1, 'Temperature'), (ax2, 'Depth')]:
                ax.annotate(str(r["bottle_number"]), 
                           xy=(r["Salinity"], r[y_col]),
                           xytext=(5, 5), textcoords="offset points", fontsize=20)
    
    return fig

def from_file_to_plot(cnv_file, bl_file, doc_file=None):
    """Load CNV and BL files, plot T-S diagram. Returns figure."""
    down_df, _ = sbe_cast(cnv_file)
    bl_df = load_bottle_file(bl_file, doc_file)
    return plot_ts(down_df, bl_df)


def plot_bl_files():
    """Process all BL files in Data folder, save T-S plots to Output folder."""
    OUTPUT_PATH.mkdir(exist_ok=True)
    bl_files = list(DATA_PATH.glob("**/*.bl"))
    
    for bl_file in bl_files:
        try:
            cnv_file = bl_file.with_name(bl_file.stem + "NTS.cnv")
            if not cnv_file.exists():
                print(f"Warning: No CNV file for {bl_file.name}")
                continue
            
            output_path = OUTPUT_PATH / f"{bl_file.parent.name}_{bl_file.stem}_plot.png"
            if output_path.exists():
                print(f"Skipping {bl_file.name} - exists")
                continue
            
            print(f"Processing: {bl_file.name}")
            fig = from_file_to_plot(cnv_file, bl_file)
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved: {output_path}")
            
        except Exception as e:
            print(f"Error processing {bl_file}: {e}")
    
    print("Processing complete!")

def process_bl_doc_files():
    """Process all BL files, merge with DOC data, save to CSV."""
    output_dir = OUTPUT_PATH / "bl_doc_files"
    output_dir.mkdir(exist_ok=True)
    bl_files = list(DATA_PATH.glob("**/*.bl"))
    
    for bl_file in bl_files:
        try:
            cnv_file = bl_file.with_name(bl_file.stem + "NTS.cnv")
            if not cnv_file.exists():
                print(f"Warning: No CNV file for {bl_file.name}")
                continue
            
            output_path = output_dir / f"{bl_file.parent.name}_{bl_file.stem}_DOC.csv"
            if output_path.exists():
                print(f"Skipping {bl_file.name} - exists")
                continue
            
            bl_doc = load_bottle_file(bl_file, DOC_FILE)
            bl_doc['station'] = bl_file.stem
            bl_doc.to_csv(output_path, index=False)
            print(f"Saved: {output_path}")
            
        except Exception as e:
            print(f"Error processing {bl_file}: {e}")
    
    print("Processing complete!")

def parse_folder_date(folder_name):
    """Parse date from folder name like '2025Aug20' -> datetime."""
    month_map = {'Jan': '01', 'Feb': '02', 'Mar': '03', 'Apr': '04',
                 'May': '05', 'Jun': '06', 'Jul': '07', 'Aug': '08',
                 'Sep': '09', 'Oct': '10', 'Nov': '11', 'Dec': '12'}
    
    if len(folder_name) >= 9 and folder_name[:4].isdigit():
        year = folder_name[:4]
        month_str = folder_name[4:7]
        day = folder_name[7:9]
        if month_str in month_map:
            return pd.to_datetime(f"{year}-{month_map[month_str]}-{day}")
    return None


def create_ctd_dataset(station_coords_file=None):
    """
    Create xarray Dataset from all CNV files in DATA_PATH.
    
    Returns Dataset with dimensions (station_date) containing ragged depth arrays.
    """
    # Load station coordinates
    coords_file = Path(station_coords_file) if station_coords_file else DATA_PATH / "station_coordinates.csv"
    station_coords = {}
    if coords_file.exists():
        station_df = pd.read_csv(coords_file)
        station_coords = {row['station']: (row['lat'], row['lon']) 
                         for _, row in station_df.iterrows()}
    
    cnv_files = list(DATA_PATH.glob("**/*NTS.cnv"))
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
            
            down_df, start_time = sbe_cast(cnv_file)
            proc_df = process_ctd(down_df, 'sbe', start_time)
            
            data_dict[(station, date)] = {
                'proc': proc_df,
                'station': station,
                'date': date
            }
            
        except Exception as e:
            print(f"Error processing {cnv_file}: {e}")
            continue
    
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
