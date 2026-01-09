"""
Data storage functions for CTD data (parquet format).
"""
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from datetime import datetime
from .config import PROJECT_ROOT


def make_parquet(down_df, bl_down_df=None, output_dir=None):
    """
    Save CTD data to parquet with metadata (CTD/DOC months included).
    Saves to PROJECT_ROOT/parquets/ with today's date in filename.
    
    Args:
        down_df: Processed CTD downcast DataFrame
        bl_down_df: Optional bottle data DataFrame
        output_dir: Optional custom output directory (default: PROJECT_ROOT/parquets)
    
    Returns: Path to saved parquet file
    """
    if output_dir is None:
        parquets_dir = PROJECT_ROOT / "parquets"
    else:
        parquets_dir = output_dir
    
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

