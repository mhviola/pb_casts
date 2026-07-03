"""
Utility functions for Padilla Bay CTD processing.
"""
import pandas as pd
import gsw


def compute_density(df, lat=48.5, lon=-122.5,
                    sal_col='sal00', temp_col='tv290C',
                    pressure_col=None):
    """
    Compute seawater density from CTD data using TEOS-10 (gsw).

    Converts Practical Salinity → Absolute Salinity and in-situ temperature
    → Conservative Temperature before computing density.  Both the in situ
    density (rho, kg/m³) and the potential density anomaly referenced to the
    surface (sigma0, kg/m³ − 1000) are returned so the result is directly
    comparable to the 'Density anomaly' column that RBR instruments export.

    Parameters
    ----------
    df : pandas.DataFrame
        Cast DataFrame.  Pressure is taken from the index (as returned by
        python-ctd / pb_casts.sbe_cast) unless *pressure_col* is given.
    lat : float, optional
        Station latitude in decimal degrees North.  Default 48.5 (Padilla Bay).
    lon : float, optional
        Station longitude in decimal degrees East (negative = West).
        Default -122.5 (Padilla Bay).
    sal_col : str, optional
        Column name for Practical Salinity [PSU].  Default 'sal00'.
    temp_col : str, optional
        Column name for in-situ temperature ITS-90 [°C].  Default 'tv290C'.
    pressure_col : str or None, optional
        Column name for pressure [dbar].  When None (default) the DataFrame
        index is used, which is the convention for python-ctd DataFrames.

    Returns
    -------
    pandas.DataFrame
        A copy of *df* with two new columns added:

        ``rho``
            In situ density [kg/m³].
        ``sigma0``
            Potential density anomaly referenced to 0 dbar [kg/m³ − 1000].
            Matches the scale of the RBR 'Density anomaly' / 'density' column.

    Examples
    --------
    >>> down = pb_casts.sbe_cast('station_S1.cnv')
    >>> down = pb_casts.compute_density(down)
    >>> down[['sal00', 'tv290C', 'rho', 'sigma0']].head()
    """
    df = df.copy()

    p  = df[pressure_col] if pressure_col else df.index.to_series()
    SP = df[sal_col]
    t  = df[temp_col]

    SA = gsw.SA_from_SP(SP, p, lon, lat)
    CT = gsw.CT_from_t(SA, t, p)
    df['rho'] = gsw.rho(SA, CT, p)
    df['sigma0'] = gsw.sigma0(SA, CT)

    return df


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

