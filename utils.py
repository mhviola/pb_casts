"""
Utility functions for Padilla Bay CTD processing.
"""
import pandas as pd


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

