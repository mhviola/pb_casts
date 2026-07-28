"""
Utility functions for Padilla Bay CTD processing.
"""
import re
import numpy as np
import pandas as pd
import gsw
from . import config


class CastFrame(pd.DataFrame):
    """
    DataFrame subclass for CTD cast data that preserves instrument metadata
    through pandas operations (slicing, copy, resample, etc.).

    Metadata is stored in the ``cast_meta`` attribute as a plain dict and is
    automatically propagated by pandas via the ``_metadata`` mechanism.

    Keys used by pb_casts:
        station          – station name string (e.g. 'S1')
        instrument_type  – 'sbe' or 'rbr'
        time             – cast start timestamp (pandas Timestamp)
        atmospheric_pressure – float, dbar (RBR casts only)
    """
    _metadata = ['cast_meta']

    @property
    def _constructor(self):
        return CastFrame

    @property
    def cast_meta(self):
        return self.__dict__.get('_cast_meta', {})

    @cast_meta.setter
    def cast_meta(self, value):
        self.__dict__['_cast_meta'] = value if value is not None else {}

    def __finalize__(self, other, method=None, **kwargs):
        super().__finalize__(other, method=method, **kwargs)
        if hasattr(other, 'cast_meta'):
            self.cast_meta = other.cast_meta
        return self


def compute_density(df,
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

        ``rho0``
            In situ density anomaly [kg/m³ − 1000]  (i.e. gsw.rho − 1000).
        ``sigma0``
            Potential density anomaly referenced to 0 dbar [kg/m³ − 1000].
            Matches the scale of the RBR 'Density anomaly' / 'density' column.

    Examples
    --------
    >>> down = pb_casts.sbe_cast('station_S1.cnv')
    >>> down = pb_casts.compute_density(down)
    >>> down[['sal00', 'tv290C', 'rho', 'sigma0']].head()
    """
    meta = df.cast_meta if isinstance(df, CastFrame) else (
           df._metadata if isinstance(df._metadata, dict) else {})
    station = meta.get('station')
    if station is None:
        raise ValueError(
            "DataFrame has no 'station' metadata. "
            "Load data via pb_casts.sbe_cast() or pb_casts.rbr_cast()."
        )
    if config.STATIONS_DF is None:
        raise RuntimeError(
            "Station coordinates not loaded. Call pb_casts.set_project_root('/path/to/project') first."
        )
    # Strip repeat suffix (e.g. 'S1b', 'G1c' → 'S1', 'G1') before coordinate lookup
    base_station = re.sub(r'[a-z]+$', '', station)
    matches = config.STATIONS_DF[config.STATIONS_DF.name == base_station]
    if matches.empty:
        raise ValueError(
            f"Station '{station}' (base: '{base_station}') not found in station_coordinates.csv. "
            f"Available: {config.STATIONS_DF.name.tolist()}"
        )
    lat = matches.lat.values[0]
    lon = matches.lon.values[0]
    df = df.copy()
    p  = df[pressure_col] if pressure_col else df.index.to_series()
    SP = df[sal_col]
    t  = df[temp_col]

    SA = gsw.SA_from_SP(SP, p, lon, lat)
    CT = gsw.CT_from_t(SA, t, p)
    df['rho0'] = gsw.rho(SA, CT, p)-1000
    df['sigma0'] = gsw.sigma0(SA, CT)

    return df


def detect_precision(series, max_decimals=6):
    """
    Infer the number of decimal places used in a numeric Series.

    Inspects up to 100 non-NaN values and returns the maximum number of
    significant decimal digits found.  Used by process_ctd() to round
    processed output to match the original sensor resolution.

    Args:
        series:       Numeric pandas Series.
        max_decimals: Upper bound on decimal places to check (default 6).

    Returns: int — decimal precision (0 if all integers, 4 if series is empty).
    """
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
    """
    Parse a date from a folder name in the format 'YYYYMon##' (e.g. '2025Aug20').

    Returns a pandas Timestamp on success, or None if the name does not match
    the expected pattern.  Used by batch functions to attach a cast date to
    processed output.
    """
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


def godin_filter(data, dt_hours=1):
    """
    Apply Godin filter: three moving averages (24h, 24h, 25h).
    Uses centered rolling means so edge values are NaN rather than
    artificially low from zero-padding.

    Parameters:
    - data: array-like or pandas Series of data to filter
    - dt_hours: time step in hours (default: 1 hour)
    """
    N24 = int(24 / dt_hours)
    N25 = int(25 / dt_hours)

    s = pd.Series(data.values if hasattr(data, 'values') else np.asarray(data))
    filtered = s.rolling(N24, center=True, min_periods=N24).mean()
    filtered = filtered.rolling(N24, center=True, min_periods=N24).mean()
    filtered = filtered.rolling(N25, center=True, min_periods=N25).mean()

    return filtered.values
