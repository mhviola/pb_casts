"""
CTD data processing functions for Padilla Bay.
"""
import numpy as np
import pandas as pd
from .utils import detect_precision, compute_density, CastFrame

# ── Processing constants ──────────────────────────────────────────────────────
# Edit these to change behaviour; dataset attrs are derived from them directly.

COMMON_PARAMS = {
    'depth_bin_size_m':                 0.25,
    'surface_noise_method':             'stability',
    'surface_noise_gradient_threshold': 0.02,
    'surface_noise_n_stable':           10,
    'surface_noise_min_depth_m':        0.5,
    'despike_n1':                       3,
    'despike_n2':                       20,
    'lp_filter_time_constant_s':        0.15,
    'smoothing_window':                 'hanning',
    'density_standard':                 'TEOS-10 (gsw)',
}

SBE_PARAMS = {
    'sbe_sample_rate_hz':       4.0,
    'sbe_despike_block':        75,
    'sbe_smooth_window_pts':    11,
}

RBR_PARAMS = {
    'rbr_sample_rate_hz':       8.0,
    'rbr_despike_block':        100,
    'rbr_smooth_window_pts':    5,
}
# ─────────────────────────────────────────────────────────────────────────────


def remove_surface_noise(cast_df, method='stability', depth_threshold=0.5,
                        gradient_threshold=0.02, n_stable=5, min_depth=0.25,
                        label=''):
    """
    Remove surface noise from CTD data.
    
    Methods: 'depth' (simple cutoff), 'stability' (n_stable consecutive points 
    with gradient < threshold), 'combined' (both). Always enforces min_depth.
    """
    # Always enforce minimum depth as safety net
    working_df = cast_df[cast_df.index >= min_depth].copy()
    if len(working_df) == 0:
        print(f"Warning: No data below min_depth={min_depth}m")
        return cast_df
    
    if method == 'depth':
        return cast_df[cast_df.index >= depth_threshold]
    
    elif method == 'stability':
        # Find where density stabilizes over n_stable consecutive points.
        # Check sign as well as magnitude: a density reversal (negative gradient)
        # should never count as "stable" even if its magnitude is small.
        den_gradient = working_df['sigma0'].diff()
        
        # Stable = small positive change: density non-decreasing AND not changing
        # too fast. Lower bound is 0 (not -threshold) so any inversion fails,
        # even a gradual one that would have passed the original abs() < threshold check.
        is_stable = (den_gradient >= 0) & (den_gradient <= gradient_threshold)

        stable_run = is_stable.rolling(window=n_stable, min_periods=n_stable).sum()

        # Second condition: the NET density change over the window must be positive.
        # This catches noisy surface zones where N consecutive tiny positive steps are
        # found by chance even though the overall density trend is flat or reversing.
        net_increase = working_df['sigma0'].diff(n_stable) > 0

        stable_indices = stable_run[(stable_run == n_stable) & net_increase].index
        
        if len(stable_indices) > 0:
            # Walk back to the start of the first fully stable window
            first_stable_window_end = stable_indices[0]
            loc = working_df.index.get_loc(first_stable_window_end)
            # get_loc returns int, slice, or bool array when index has duplicates
            if isinstance(loc, slice):
                end_pos = loc.start
            elif isinstance(loc, np.ndarray):
                end_pos = int(loc.nonzero()[0][0])
            else:
                end_pos = int(loc)
            first_stable_window_start_pos = max(0, end_pos - (n_stable - 1))
            first_stable_idx = working_df.index[first_stable_window_start_pos]
            return working_df[working_df.index >= first_stable_idx]
        else:
            # No stable region found - fall back to depth threshold
            prefix = f"[{label}] " if label else ""
            print(f"Warning: {prefix}No stable region found, using depth_threshold={depth_threshold}m")
            return cast_df[cast_df.index >= depth_threshold]
    
    elif method == 'combined':
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
            loc = depth_filtered.index.get_loc(first_stable_window_end)
            if isinstance(loc, slice):
                end_pos = loc.start
            elif isinstance(loc, np.ndarray):
                end_pos = int(loc.nonzero()[0][0])
            else:
                end_pos = int(loc)
            first_stable_window_start_pos = max(0, end_pos - (n_stable - 1))
            first_stable_idx = depth_filtered.index[first_stable_window_start_pos]
            return depth_filtered[depth_filtered.index >= first_stable_idx]
        else:
            # Already depth filtered, return that
            return depth_filtered
    
    else:
        print(f"Warning: Unknown method '{method}'. Using depth threshold only.")
        return cast_df[cast_df.index >= depth_threshold]


def process_ctd(df, smooth=True, columns=None, surface_cutoff_m=None):
    """
    Full CTD processing pipeline: surface-noise removal → despike → low-pass
    filter → pressure check → interpolation → 0.25 m depth binning → smoothing.

    Processing parameters (sample rate, despike block, smooth window) are
    selected automatically from SBE_PARAMS / RBR_PARAMS based on
    df.cast_meta['instrument_type'].  The cast must have been loaded via
    sbe_cast() or rbr_cast() for the metadata to be present.

    Args:
        df:               CastFrame with pressure/depth as index and cast_meta populated.
        smooth:           Apply Hanning-window smoothing after binning (default True).
        columns:          Columns to run through the QC pipeline.  Defaults to all
                          numeric columns except timeS, bottle_number, and doc_conc.
        surface_cutoff_m: If given, skip the automatic stability algorithm and simply
                          remove all data shallower than this depth (metres).  Use
                          cutoffs from review_surface_cutoffs() for manual overrides.

    Returns: Processed CastFrame with 0.25 m depth bins and cast_meta preserved.
    """
    df = compute_density(df)

    exclude_cols = ['timeS', 'bottle_number', 'doc_conc']

    # Stash any datetime column so it can be re-interpolated after binning
    # (datetime columns cannot go through the numeric QC chain).
    date_cols = [col for col in df.columns if np.issubdtype(df[col].dtype, np.datetime64)]
    num_date_cols = len(date_cols)
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
    bottle_data = {col: df[col].copy() for col in bottle_cols if col in df.columns}

    meta = df.cast_meta if isinstance(df, CastFrame) else (
           df._metadata if isinstance(df._metadata, dict) else {})
    instrument_type = meta.get('instrument_type')
    if instrument_type is None:
        raise ValueError(
            "DataFrame has no 'instrument_type' metadata. "
            "Load data via pb_casts.sbe_cast() or pb_casts.rbr_cast()."
        )
    if instrument_type.lower() == 'rbr':
        sample_rate   = RBR_PARAMS['rbr_sample_rate_hz']
        despike_block = RBR_PARAMS['rbr_despike_block']
        smooth_window = RBR_PARAMS['rbr_smooth_window_pts']
    elif instrument_type.lower() == 'sbe':
        sample_rate   = SBE_PARAMS['sbe_sample_rate_hz']
        despike_block = SBE_PARAMS['sbe_despike_block']
        smooth_window = SBE_PARAMS['sbe_smooth_window_pts']
    else:
        raise ValueError(f"Invalid instrument type: {instrument_type}")

    bin_delta = COMMON_PARAMS['depth_bin_size_m']
    if surface_cutoff_m is not None:
        df_clean = df[df.index >= surface_cutoff_m].copy()
    else:
        cast_label = f"{meta.get('station', '?')} {meta.get('time', '')}"
        df_clean = remove_surface_noise(
            df,
            method=COMMON_PARAMS['surface_noise_method'],
            gradient_threshold=COMMON_PARAMS['surface_noise_gradient_threshold'],
            n_stable=COMMON_PARAMS['surface_noise_n_stable'],
            min_depth=COMMON_PARAMS['surface_noise_min_depth_m'],
            label=cast_label,
        )
    if len(df_clean) == 0:
        raise ValueError(
            f"Cast at {meta.get('station','?')} on {meta.get('time','?')} is empty after "
            "surface-noise removal — check the raw data or relax remove_surface_noise parameters."
        )
    # Cap block size to the actual data length — despike raises an error if
    # len(df_clean) < block. Short casts still go through the full pipeline; a
    # smaller block just means tighter local statistics for spike detection.
    if despike_block > len(df_clean):
        print(f"Warning: short cast at {meta.get('station','?')} on {meta.get('time','?')} "
              f"({len(df_clean)} rows); capping despike block from {despike_block} to {len(df_clean)}.")
        despike_block = len(df_clean    )

    proc_df = (
        df_clean[cols_to_process]
        .despike(n1=COMMON_PARAMS['despike_n1'], n2=COMMON_PARAMS['despike_n2'], block=despike_block)
        .lp_filter(sample_rate=sample_rate, time_constant=COMMON_PARAMS['lp_filter_time_constant_s'])
        .press_check()
        .interpolate(method="index", limit_direction="both", limit_area="inside")
    )
    
    # Custom binning: build depth grid then interpolate to preserve near-surface data
    min_depth = max(0.0, np.floor(proc_df.index.min() * 4) / 4)
    max_depth = np.ceil(proc_df.index.max() * 4) / 4
    new_index = np.arange(min_depth, max_depth + bin_delta, bin_delta)

    binned_data = {}
    for col in proc_df.columns:
        binned_data[col] = np.interp(new_index, proc_df.index.values, proc_df[col].values)
    
    proc_df = pd.DataFrame(binned_data, index=new_index)
    proc_df.index.name = df.index.name or 'depth'

    if smooth:
        proc_df = proc_df.smooth(window_len=smooth_window, window="hanning")
    
    for col in proc_df.columns:
        if col in df.columns:
            precision = detect_precision(df[col])
            proc_df[col] = np.round(proc_df[col], precision)

    if 'timeS' in df.columns:
        proc_df['timeS'] = np.interp(new_index, df.index.values, df['timeS'].values)
        time_resolution = 1.0 / sample_rate  # 0.25s for 4Hz, 0.125s for 8Hz
        proc_df['timeS'] = np.round(proc_df['timeS'] / time_resolution) * time_resolution
        proc_df['time_local'] = meta.get('time') + pd.to_timedelta(proc_df['timeS'], unit='s')
    elif datetime_df is not None:
        time_resolution = 1.0 / sample_rate
        # Convert datetime to numeric (timestamps), interpolate, then convert back
        datetime_numeric = datetime_df.astype(np.int64) / 1e9  # Convert to seconds since epoch
        interpolated_numeric = np.interp(new_index, df.index.values, datetime_numeric)
        rounded_numeric = np.round(interpolated_numeric / time_resolution) * time_resolution
        proc_df['time_local'] = pd.to_datetime(rounded_numeric, unit='s')

    # Re-add bottle columns: snap each value to the nearest depth bin,
    # then NaN-out any bin that is more than half a bin-width from a real bottle depth.
    for col in bottle_cols:
        if col in bottle_data:
            proc_df[col] = bottle_data[col].reindex(proc_df.index, method='nearest')
            orig_non_nan = bottle_data[col].dropna()
            if len(orig_non_nan) > 0:
                orig_indices = orig_non_nan.index.values
                proc_indices = proc_df.index.values
                for i, proc_idx in enumerate(proc_indices):
                    if np.min(np.abs(orig_indices - proc_idx)) > bin_delta / 2:
                        proc_df.iloc[i, proc_df.columns.get_loc(col)] = np.nan

    result = CastFrame(proc_df)
    result.cast_meta = meta
    return result

