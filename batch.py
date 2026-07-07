"""
Batch processing functions for multiple CTD files.
"""
import re
from pathlib import Path
import numpy as np
import pandas as pd
from . import config
from .io import sbe_cast, rbr_cast, castaway_cast, load_bottle_file, station_from_path
from .processing import process_ctd
from .utils import parse_folder_date


def load_surface_cutoffs(cutoffs_file):
    """
    Load surface cutoffs from a JSON or CSV file.

    JSON format (produced by the interactive notebook cell):
        { "2025Oct20_S1": 0.75, "2025Nov18_G1": 0.5, ... }

    CSV format (produced by review_surface_cutoffs()):
        folder, station, cutoff_depth_m
        2025Oct20, S1, 0.75

    Returns a dict keyed by '{folder}_{station}' → cutoff depth (float, metres).
    Entries with missing or non-numeric values are silently skipped.
    """
    import json
    cutoffs_file = Path(cutoffs_file)
    cutoffs = {}

    if cutoffs_file.suffix.lower() == '.json':
        with open(cutoffs_file) as f:
            raw = json.load(f)
        for key, val in raw.items():
            try:
                cutoffs[key] = float(val)
            except (ValueError, TypeError):
                pass
    else:
        df = pd.read_csv(cutoffs_file)
        for _, row in df.iterrows():
            key = f"{row['folder']}_{row['station']}"
            try:
                val = float(row['cutoff_depth_m'])
                if not np.isnan(val):
                    cutoffs[key] = val
            except (ValueError, TypeError):
                pass

    return cutoffs


def review_surface_cutoffs(data_path=None, output_path=None, review_depth_m=5.0):
    """
    Generate a surface-review figure and a CSV template for manual cutoff entry.

    For every cast found in the dated folders this function:
      1. Loads the raw (unprocessed) cast.
      2. Plots the top *review_depth_m* of salinity and temperature side-by-side.
      3. Saves a multi-page PDF (one page per cast) to
         ``output_path/surface_review.pdf``.
      4. Saves a CSV template to ``output_path/surface_cutoffs.csv`` with columns
         ``folder``, ``station``, ``cutoff_depth_m`` (blank — fill these in).

    After calling this function:
      - Open ``surface_review.pdf`` and decide on a cutoff depth for each cast.
      - Fill in ``surface_cutoffs.csv`` with those depths.
      - Pass the CSV path to ``batch_process_all(surface_cutoffs_file=...)`` to
        reprocess with your manual cutoffs.

    Args:
        data_path:      Path to Data folder (default: DATA_PATH from config).
        output_path:    Where to save output files (default: OUTPUT_PATH from config).
        review_depth_m: How many metres from the surface to plot (default 5).
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    data_path   = Path(data_path)   if data_path   else config.DATA_PATH
    output_path = Path(output_path) if output_path else config.OUTPUT_PATH

    output_path.mkdir(parents=True, exist_ok=True)
    pdf_path = output_path / "surface_review.pdf"
    csv_path = output_path / "surface_cutoffs.csv"

    folder_pattern = re.compile(r'^\d{4}[A-Z][a-z]{2}\d{2}$')
    dated_folders  = sorted([
        f for f in data_path.iterdir()
        if f.is_dir() and folder_pattern.match(f.name)
    ])

    from matplotlib.ticker import MultipleLocator, ScalarFormatter

    def _percentile_xlim(series, margin=0.05):
        """50th–95th percentile limits to focus on real water values."""
        lo, hi = series.quantile(0.50), series.quantile(0.95)
        pad = max((hi - lo) * margin, 0.05)
        return lo - pad, hi + pad

    def _fmt_xaxis(ax, data_range):
        if data_range <= 1.0:
            step = 0.1
        elif data_range <= 3.0:
            step = 0.2
        elif data_range <= 6.0:
            step = 0.5
        else:
            step = 1.0
        ax.xaxis.set_major_locator(MultipleLocator(step))
        fmt = ScalarFormatter(useOffset=False)
        fmt.set_scientific(False)
        ax.xaxis.set_major_formatter(fmt)
        ax.tick_params(axis='x', labelrotation=45)

    # ── Pass 1: load all casts, compute global per-variable x limits ──────────
    SAL_COLS  = ['sal00', 'salinity', 'Salinity']
    TEMP_COLS = ['tv290C', 'temperature', 'Temperature']

    all_casts = []   # list of (folder_name, station, surface_df)
    rows      = []
    global_sal_lims  = []   # (lo, hi) per cast
    global_temp_lims = []

    for folder in dated_folders:
        cnv_files  = sorted([
            f for f in folder.glob('*NTS.cnv')
            if 'do not use' not in f.name.lower() and 'redone' not in f.stem.lower()
        ])
        xlsx_files = sorted(folder.glob('*.xlsx'))

        cast_list = []
        for cnv_file in cnv_files:
            station = cnv_file.stem.replace('NTS', '')
            try:
                cast_list.append((station, sbe_cast(str(cnv_file))))
            except Exception as e:
                print(f"[{folder.name}] Could not load {cnv_file.name}: {e}")

        for xlsx_file in xlsx_files:
            try:
                casts = rbr_cast(str(xlsx_file))
                for cast_df in casts:
                    station = cast_df.cast_meta.get('station', 'unknown')
                    cast_list.append((station, cast_df))
            except Exception as e:
                print(f"[{folder.name}] Could not load {xlsx_file.name}: {e}")

        for station, raw_df in cast_list:
            rows.append({'folder': folder.name, 'station': station, 'cutoff_depth_m': ''})
            surface = raw_df[raw_df.index <= review_depth_m]
            all_casts.append((folder.name, station, surface))

            sal_col  = next((c for c in SAL_COLS  if c in surface.columns), None)
            temp_col = next((c for c in TEMP_COLS if c in surface.columns), None)
            if sal_col:
                global_sal_lims.append(_percentile_xlim(surface[sal_col]))
            if temp_col:
                global_temp_lims.append(_percentile_xlim(surface[temp_col]))

    # Global range WIDTH: widest cast across all dates sets the common scale.
    # Each cast is then centered on its own data but spans this same width,
    # so 1 PSU (or 1°C) is always the same physical size on every page.
    g_sal_range  = max(hi - lo for lo, hi in global_sal_lims)  if global_sal_lims  else 1.0
    g_temp_range = max(hi - lo for lo, hi in global_temp_lims) if global_temp_lims else 1.0

    def _centred_xlim(series, common_range):
        mid = (series.quantile(0.50) + series.quantile(0.95)) / 2
        return mid - common_range / 2, mid + common_range / 2

    # Fixed vertical scale: 1.2 inches per metre of depth so 0.25 m gridlines
    # are always the same physical size regardless of how deep each cast goes.
    INCHES_PER_M = 1.2
    LABEL_PAD    = 1.8  # extra inches for title, x-labels, tick labels

    # ── Pass 2: plot every cast using the global x limits ─────────────────────
    with PdfPages(pdf_path) as pdf:
        for folder_name, station, surface in all_casts:
            sal_col  = next((c for c in SAL_COLS  if c in surface.columns), None)
            temp_col = next((c for c in TEMP_COLS if c in surface.columns), None)

            n_panels = sum(c is not None for c in [sal_col, temp_col])
            if n_panels == 0:
                print(f"[{folder_name}] {station}: no salinity or temperature column, skipping plot.")
                continue

            # Y range: surface (0) to the deepest data point, capped at review_depth_m
            depth_max  = min(float(surface.index.max()), review_depth_m)
            fig_height = depth_max * INCHES_PER_M + LABEL_PAD

            fig, axes = plt.subplots(1, n_panels, figsize=(4 * n_panels, fig_height), sharey=True)
            if n_panels == 1:
                axes = [axes]

            panel = 0
            if sal_col:
                axes[panel].plot(surface[sal_col], surface.index, color='steelblue', lw=1, marker='o', markersize=2)
                axes[panel].set_xlabel('Salinity (PSU)')
                axes[panel].set_xlim(*_centred_xlim(surface[sal_col], g_sal_range))
                axes[panel].set_ylabel('Depth (m)')
                axes[panel].grid(True, alpha=0.3)
                _fmt_xaxis(axes[panel], g_sal_range)
                panel += 1
            if temp_col:
                axes[panel].plot(surface[temp_col], surface.index, color='tomato', lw=1, marker='o', markersize=2)
                axes[panel].set_xlabel('Temperature (°C)')
                axes[panel].set_xlim(*_centred_xlim(surface[temp_col], g_temp_range))
                axes[panel].grid(True, alpha=0.3)
                _fmt_xaxis(axes[panel], g_temp_range)

            axes[0].set_ylim(depth_max, 0)
            axes[0].yaxis.set_major_locator(MultipleLocator(0.25))

            fig.suptitle(f"{folder_name}  |  {station}", fontsize=11, fontweight='bold')
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    pd.DataFrame(rows).to_csv(csv_path, index=False)
    print(f"Review PDF saved to:      {pdf_path}")
    print(f"Cutoffs CSV template at:  {csv_path}")
    print(f"Fill in 'cutoff_depth_m' for each cast, then run batch_process_all(surface_cutoffs_file='{csv_path}').")


def plot_cutoff_check(surface_cutoffs_file, data_path=None, output_path=None, review_depth_m=3.0):
    """
    Generate a QC PDF showing raw vs processed data for every cast.

    For each cast the plot shows:
      - Raw salinity and temperature (thin grey line, full surface data)
      - Processed / smoothed data on top (coloured line)
      - A dashed horizontal line at the surface cutoff depth

    Casts with a manual cutoff in *surface_cutoffs_file* use that depth;
    others fall back to the automatic stability algorithm.

    Args:
        surface_cutoffs_file: Path to JSON or CSV produced by the interactive
                              notebook cell or review_surface_cutoffs().
        data_path:            Path to Data folder (default: DATA_PATH from config).
        output_path:          Where to save the PDF (default: OUTPUT_PATH from config).
        review_depth_m:       Depth range to plot (default 5 m).
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.ticker import MultipleLocator, ScalarFormatter

    data_path   = Path(data_path)   if data_path   else config.DATA_PATH
    output_path = Path(output_path) if output_path else config.OUTPUT_PATH
    output_path.mkdir(parents=True, exist_ok=True)

    cutoffs = load_surface_cutoffs(surface_cutoffs_file)
    pdf_path = output_path / "cutoff_check.pdf"

    SAL_COLS  = ['sal00', 'salinity', 'Salinity']
    TEMP_COLS = ['tv290C', 'temperature', 'Temperature']

    INCHES_PER_M = 1.2
    LABEL_PAD    = 1.8

    folder_pattern = re.compile(r'^\d{4}[A-Z][a-z]{2}\d{2}$')
    dated_folders  = sorted([
        f for f in data_path.iterdir()
        if f.is_dir() and folder_pattern.match(f.name)
    ])

    def _fmt_xaxis(ax, data_range):
        if data_range <= 1.0:
            step = 0.1
        elif data_range <= 3.0:
            step = 0.2
        elif data_range <= 6.0:
            step = 0.5
        else:
            step = 1.0
        ax.xaxis.set_major_locator(MultipleLocator(step))
        fmt = ScalarFormatter(useOffset=False)
        fmt.set_scientific(False)
        ax.xaxis.set_major_formatter(fmt)
        ax.tick_params(axis='x', labelrotation=45)

    # ── Pass 1: load all casts, collect global x ranges ───────────────────────
    all_casts        = []   # (folder_name, station, raw_df, proc_df, cutoff_depth)
    global_sal_lims  = []
    global_temp_lims = []

    for folder in dated_folders:
        cnv_files  = sorted([
            f for f in folder.glob('*NTS.cnv')
            if 'do not use' not in f.name.lower() and 'redone' not in f.stem.lower()
        ])
        xlsx_files = sorted(folder.glob('*.xlsx'))

        cast_list = []
        for cnv_file in cnv_files:
            station = cnv_file.stem.replace('NTS', '')
            try:
                cast_list.append((station, sbe_cast(str(cnv_file))))
            except Exception as e:
                print(f"[{folder.name}] Could not load {cnv_file.name}: {e}")

        for xlsx_file in xlsx_files:
            try:
                for cast_df in rbr_cast(str(xlsx_file)):
                    station = cast_df.cast_meta.get('station', 'unknown')
                    cast_list.append((station, cast_df))
            except Exception as e:
                print(f"[{folder.name}] Could not load {xlsx_file.name}: {e}")

        for station, raw_df in cast_list:
            cutoff_key   = f"{folder.name}_{station}"
            cutoff_depth = cutoffs.get(cutoff_key)
            try:
                proc_df = process_ctd(raw_df, surface_cutoff_m=cutoff_depth)
            except Exception as e:
                print(f"[{folder.name}] {station}: processing failed — {e}")
                continue

            surface = raw_df[raw_df.index <= review_depth_m]
            proc_surface = proc_df[proc_df.index <= review_depth_m]
            all_casts.append((folder.name, station, surface, proc_surface, cutoff_depth))

            sal_col  = next((c for c in SAL_COLS  if c in surface.columns), None)
            temp_col = next((c for c in TEMP_COLS if c in surface.columns), None)
            # Global range from PROCESSED data only — avoids extreme surface
            # noise (near-zero salinity in air) blowing out the x scale.
            if sal_col and sal_col in proc_surface.columns:
                s = proc_surface[sal_col].dropna()
                if len(s):
                    lo, hi = s.min(), s.max()
                    pad = max((hi - lo) * 0.05, 0.05)
                    global_sal_lims.append((lo - pad, hi + pad))
            if temp_col and temp_col in proc_surface.columns:
                s = proc_surface[temp_col].dropna()
                if len(s):
                    lo, hi = s.min(), s.max()
                    pad = max((hi - lo) * 0.05, 0.05)
                    global_temp_lims.append((lo - pad, hi + pad))

    g_sal_range  = max(hi - lo for lo, hi in global_sal_lims)  if global_sal_lims  else 1.0
    g_temp_range = max(hi - lo for lo, hi in global_temp_lims) if global_temp_lims else 1.0

    def _centred_xlim(raw_series, proc_series, common_range):
        # Centre on the processed data midpoint (clean water-column values),
        # then extend the left edge if raw surface values fall outside the window.
        if len(proc_series.dropna()):
            mid = (proc_series.quantile(0.10) + proc_series.quantile(0.90)) / 2
        else:
            mid = raw_series.quantile(0.50)
        left  = mid - common_range / 2
        right = mid + common_range / 2
        # Pull left edge to include raw surface dip (robust low percentile)
        raw_lo = raw_series.quantile(0.05)
        if raw_lo < left:
            left = raw_lo - max((right - raw_lo) * 0.02, 0.05)
        return left, right

    # ── Pass 2: plot ───────────────────────────────────────────────────────────
    with PdfPages(pdf_path) as pdf:
        for folder_name, station, surface, proc_surface, cutoff_depth in all_casts:
            sal_col  = next((c for c in SAL_COLS  if c in surface.columns), None)
            temp_col = next((c for c in TEMP_COLS if c in surface.columns), None)

            n_panels = sum(c is not None for c in [sal_col, temp_col])
            if n_panels == 0:
                continue

            depth_max  = min(float(surface.index.max()), review_depth_m)
            fig_height = depth_max * INCHES_PER_M + LABEL_PAD
            fig, axes  = plt.subplots(1, n_panels, figsize=(4 * n_panels, fig_height), sharey=True)
            if n_panels == 1:
                axes = [axes]

            cutoff_label = f"cutoff = {cutoff_depth} m" if cutoff_depth is not None else "cutoff = auto"

            panel = 0
            for col, color, xlabel, g_range in [
                (sal_col,  'steelblue', 'Salinity (PSU)',    g_sal_range),
                (temp_col, 'tomato',    'Temperature (°C)',  g_temp_range),
            ]:
                if col is None:
                    continue
                ax = axes[panel]

                # Raw data — thin grey
                ax.plot(surface[col], surface.index,
                        color='grey', lw=0.8, alpha=0.5, zorder=1, label='raw')

                # Processed data — coloured, on top
                if col in proc_surface.columns:
                    ax.plot(proc_surface[col], proc_surface.index,
                            color=color, lw=1.5, zorder=2, label='processed')

                # Cutoff line
                if cutoff_depth is not None:
                    ax.axhline(cutoff_depth, color='black', lw=1.2,
                               linestyle='--', zorder=3,
                               label=cutoff_label if panel == 0 else '_')

                proc_col = proc_surface[col] if col in proc_surface.columns else pd.Series(dtype=float)
                ax.set_xlim(*_centred_xlim(surface[col], proc_col, g_range))
                ax.set_xlabel(xlabel)
                ax.grid(True, alpha=0.3)
                _fmt_xaxis(ax, g_range)
                if panel == 0:
                    ax.set_ylabel('Depth (m)')
                    ax.legend(fontsize=7, loc='lower left')
                panel += 1

            axes[0].set_ylim(depth_max, 0)
            axes[0].yaxis.set_major_locator(MultipleLocator(0.25))

            fig.suptitle(f"{folder_name}  |  {station}  |  {cutoff_label}",
                         fontsize=11, fontweight='bold')
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    print(f"Cutoff check PDF saved to: {pdf_path}")


def process_bl_doc_files(data_path=None, output_path=None, doc_file=None):
    """
    Process all BL files, merge with DOC data, save to CSV.
    
    Args:
        data_path: Path to data directory (default: DATA_PATH from config)
        output_path: Path to output directory (default: OUTPUT_PATH from config)
        doc_file: Path to DOC excel file (default: DOC_FILE from config)
    """
    data_path = Path(data_path) if data_path else config.DATA_PATH
    output_path = Path(output_path) if output_path else config.OUTPUT_PATH
    doc_file = Path(doc_file) if doc_file else config.DOC_FILE
    
    output_dir = output_path / "bl_doc_files"
    output_dir.mkdir(exist_ok=True)
    bl_files = list(data_path.glob("**/*.bl"))
    
    for bl_file in bl_files:
        try:
            cnv_file = bl_file.with_name(bl_file.stem + ".cnv")
            if not cnv_file.exists():
                print(f"Warning: No CNV file for {bl_file.name}")
                continue
            
            output_file = output_dir / f"{bl_file.parent.name}_{bl_file.stem}_DOC.csv"
            if output_file.exists():
                print(f"Skipping {bl_file.name} - exists")
                continue
            
            bl_doc = load_bottle_file(bl_file, doc_file)
            bl_doc['station'] = station_from_path(str(bl_file)) or bl_file.stem
            bl_doc.to_csv(output_file, index=False)
            print(f"Saved: {output_file}")
            
        except Exception as e:
            print(f"Error processing {bl_file}: {e}")
    
    print("Processing complete!")


def _find_castaway_files_for_date(date, castaway_path):
    """
    Return all valid CastAway CSV files in castaway_path that match a given date.

    CastAway files live in subfolders named 'Viola_DDMONYYYY' (e.g. 'Viola_18nov2025').
    Each subfolder maps to a single cruise date.  Only 'Cast' (not 'Invalid') files
    are returned.

    Args:
        date:          pd.Timestamp (or datetime) for the cruise date.
        castaway_path: Path to the root CastAway_profiles/ folder.

    Returns: sorted list of Path objects for matching CSV files.
    """
    from .io import _parse_castaway_header
    castaway_path = Path(castaway_path)
    if not castaway_path.exists():
        return []

    matched = []
    for subfolder in castaway_path.glob('Viola_*'):
        if not subfolder.is_dir():
            continue
        # Parse date from folder name: 'Viola_18nov2025' → 2025-11-18
        m = re.match(r'Viola_(\d{1,2})([a-z]{3})(\d{4})$', subfolder.name, re.IGNORECASE)
        if not m:
            continue
        try:
            folder_date = pd.to_datetime(f"{m.group(1)} {m.group(2)} {m.group(3)}", dayfirst=True)
        except Exception:
            continue
        if date is None or folder_date.date() != pd.Timestamp(date).date():
            continue
        for csv_file in sorted(subfolder.glob('*.csv')):
            try:
                header = _parse_castaway_header(str(csv_file))
                if 'Invalid' not in header.get('Sample type', ''):
                    matched.append(csv_file)
            except Exception:
                pass

    return matched


def plot_multi_instrument_pdf(nc_file, output_pdf=None):
    """
    Generate a PDF comparing profiles from multiple instruments at the same cast.

    Scans every (station, date, repeat) combination in an xarray Dataset
    loaded from *nc_file* and creates one page per combination where two or
    more instruments have valid data.  Each page shows three panels:
    salinity, density anomaly (rho0), and temperature side-by-side with
    depth (inverted) on the y-axis and one line per instrument.

    Args:
        nc_file:    Path to the NetCDF file produced by create_ctd_dataset().
        output_pdf: Path for the output PDF.  Defaults to the same directory
                    as *nc_file* with the stem suffixed ``_multi_inst.pdf``.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import xarray as xr

    nc_file = Path(nc_file)
    if output_pdf is None:
        output_pdf = nc_file.parent / (nc_file.stem + '_multi_inst.pdf')
    output_pdf = Path(output_pdf)

    ds = xr.open_dataset(nc_file)

    INST_COLORS  = {'castaway': '#1f77b4', 'rbr': '#ff7f0e', 'sbe': '#2ca02c'}
    INST_MARKERS = {'castaway': 'o', 'rbr': None, 'sbe': 's'}
    DEFAULT_COLOR  = '#9467bd'
    DEFAULT_MARKER = None

    panels = [
        ('sal00',  'Salinity (PSU)'),
        ('tv290C', 'Temperature (°C)'),
        ('rho0',   'Density Anomaly (kg m⁻³)'),
    ]

    pages = 0
    with PdfPages(output_pdf) as pdf:
        for station in ds.station.values:
            for date in ds.date.values:
                for repeat in ds.repeat.values:
                    sel = ds.sel(station=station, date=date, repeat=repeat)

                    # Identify instruments with at least one non-NaN value
                    active_insts = []
                    for inst in ds.instrument.values:
                        prof = sel.sel(instrument=inst)
                        # Check if any panel variable has data
                        has_data = any(
                            var in ds and not np.all(np.isnan(prof[var].values))
                            for var, _ in panels
                        )
                        if has_data:
                            active_insts.append(inst)

                    if len(active_insts) < 2:
                        continue

                    depth = ds.depth.values

                    fig, axes = plt.subplots(1, 3, figsize=(10, 10), sharey=True)

                    for var, xlabel in zip([p[0] for p in panels], [p[1] for p in panels]):
                        ax = axes[[p[0] for p in panels].index(var)]
                        for inst in active_insts:
                            prof   = sel.sel(instrument=inst)
                            if var not in ds:
                                continue
                            values = prof[var].values
                            mask   = ~np.isnan(values)
                            if not mask.any():
                                continue
                            color  = INST_COLORS.get(inst, DEFAULT_COLOR)
                            marker = INST_MARKERS.get(inst, DEFAULT_MARKER)
                            ax.plot(
                                values[mask], depth[mask],
                                color=color, marker=marker,
                                markersize=2, linewidth=1.2,
                                label=inst,
                            )
                        ax.set_xlabel(xlabel)
                        ax.grid(True, alpha=0.3)

                    axes[0].invert_yaxis()
                    axes[0].set_ylabel('Depth (m)')
                    axes[0].legend(fontsize=8)

                    date_str   = pd.Timestamp(date).strftime('%Y-%m-%d')
                    repeat_str = '' if repeat == 0 else f'  repeat {repeat}'
                    fig.suptitle(
                        f'{station}  |  {date_str}{repeat_str}',
                        fontsize=12, fontweight='bold',
                    )
                    fig.tight_layout()
                    pdf.savefig(fig)
                    plt.close(fig)
                    pages += 1

    print(f"Saved {pages} page(s) to {output_pdf}")
    return output_pdf


def batch_process_all(data_path=None, output_path=None, castaway_path=None,
                      overwrite=False, surface_cutoffs_file=None):
    """
    Process all CTD casts in yearMON## folders (e.g. 2025Aug20).
    Handles SBE CNV, SBE raw HEX, RBR XLSX, and CastAway CSV data.
    Saves one parquet file per cast to output_path/{folder}/{station}_{instrument}.parquet.

    Each output parquet includes an 'instrument' column ('sbe', 'rbr', or 'castaway')
    so multi-instrument casts on the same date can be distinguished.

    Args:
        data_path:            Path to Data folder (default: DATA_PATH from config)
        output_path:          Where to save output parquet files (default: OUTPUT_PATH/casts)
        castaway_path:        Path to CastAway_profiles folder (default: CASTAWAY_PATH from config)
        overwrite:            Re-process and overwrite existing parquet files (default: False)
        surface_cutoffs_file: Optional path to CSV produced by review_surface_cutoffs().
                              When provided, manual cutoff depths override the automatic
                              stability algorithm for each matched cast.
    """
    data_path     = Path(data_path)     if data_path     else config.DATA_PATH
    output_path   = Path(output_path)   if output_path   else config.OUTPUT_PATH / "casts"
    castaway_path = Path(castaway_path) if castaway_path else config.CASTAWAY_PATH

    cutoffs: dict = {}
    if surface_cutoffs_file is not None:
        cutoffs = load_surface_cutoffs(surface_cutoffs_file)
        print(f"Loaded {len(cutoffs)} manual surface cutoffs from {surface_cutoffs_file}")

    folder_pattern = re.compile(r'^\d{4}[A-Z][a-z]{2}\d{2}$')
    dated_folders = sorted([
        f for f in data_path.iterdir()
        if f.is_dir() and folder_pattern.match(f.name)
    ])

    if not dated_folders:
        print("No dated folders found.")
        return

    print(f"Found {len(dated_folders)} folder(s): {[f.name for f in dated_folders]}\n")
    total_processed = 0
    total_skipped = 0
    total_errors = 0

    for folder in dated_folders:
        date = parse_folder_date(folder.name)

        # SBE: any *.cnv in the folder or its subdirectories
        cnv_files = sorted([
            f for f in folder.rglob('*.cnv')
            if 'do not use' not in f.name.lower() and 'redone' not in f.stem.lower()
        ])
        # RBR: .xlsx files at the top level of the dated folder
        xlsx_files = sorted(folder.glob('*.xlsx'))
        # CastAway: CSVs from the matching Viola_* subfolder
        castaway_files = _find_castaway_files_for_date(date, castaway_path)

        if not cnv_files and not xlsx_files and not castaway_files:
            print(f"[{folder.name}] No CTD files found, skipping.\n")
            continue

        out_dir = output_path / folder.name
        out_dir.mkdir(parents=True, exist_ok=True)

        # ref_cast_times: base_station → [(station_with_suffix, timestamp), …]
        # Populated by SBE/RBR loops so CastAway can time-match to the right repeat.
        ref_cast_times: dict[str, list[tuple[str, pd.Timestamp]]] = {}

        def _register_ref(cast_df):
            """Store (station_tag, time) in ref_cast_times keyed by base station."""
            tag = cast_df.cast_meta.get('station')
            t   = cast_df.cast_meta.get('time')
            if tag and t is not None:
                base = re.match(r'([A-Z0-9]+)', tag).group(1)
                ref_cast_times.setdefault(base, []).append((tag, pd.Timestamp(t)))

        def _resolve_castaway_station_tag(cast_df, matched_tags):
            """
            Return the station tag (with repeat suffix if needed) that best matches
            the CastAway cast time to an already-loaded reference cast.

            Uses greedy closest-time assignment: each reference tag can only be
            claimed once.  Falls back to the base station name when no reference
            exists or the cast has no timestamp.
            """
            base      = cast_df.cast_meta.get('station', '')
            cast_time = cast_df.cast_meta.get('time')
            candidates = [
                (tag, t)
                for tag, t in ref_cast_times.get(base, [])
                if tag not in matched_tags
            ]
            if candidates and cast_time is not None:
                try:
                    ct = pd.Timestamp(cast_time)
                    closest_tag, _ = min(
                        candidates,
                        key=lambda x: abs((pd.Timestamp(x[1]) - ct).total_seconds()),
                    )
                    matched_tags.add(closest_tag)
                    return closest_tag
                except Exception:
                    pass
            return base  # fallback: no reference or no timestamp

        # Track how many times each station_tag has been used for non-reference
        # instruments so auto-suffixing still works when there's no reference.
        station_counts: dict[str, int] = {}

        def _save_cast(cast_df, instrument_tag, label, station_tag_override=None):
            """Process one cast and write a parquet; returns True on success."""
            if station_tag_override is not None:
                station_tag = station_tag_override
            else:
                station = cast_df.cast_meta.get('station')
                if station is None:
                    print(f"[{folder.name}] Skipping {label} — could not identify station.")
                    return False
                count = station_counts.get(station, 0)
                suffix = '' if count == 0 else chr(ord('b') + count - 1)
                station_tag = f"{station}{suffix}"
                station_counts[station] = count + 1

            out_file = out_dir / f"{station_tag}_{instrument_tag}.parquet"
            if out_file.exists() and not overwrite:
                print(f"[{folder.name}] {station_tag} ({instrument_tag}): already exists, skipping.")
                return None  # skipped
            try:
                print(f"[{folder.name}] Processing {instrument_tag.upper()} {station_tag}...", end=' ')
                cutoff_depth = cutoffs.get(f"{folder.name}_{station_tag}")
                proc_df = process_ctd(cast_df, surface_cutoff_m=cutoff_depth)
                proc_df['station'] = station_tag
                proc_df['instrument'] = instrument_tag
                if date is not None:
                    proc_df['cast_date'] = date.date()
                proc_df.to_parquet(out_file)
                print(f"saved ({len(proc_df)} depth bins) -> {out_file.name}")
                return True
            except Exception as e:
                print(f"ERROR: {e}")
                return False

        # --- SBE CNV processing ---
        for cnv_file in cnv_files:
            try:
                cast_df = sbe_cast(str(cnv_file))
            except Exception as e:
                print(f"[{folder.name}] ERROR loading {cnv_file.name}: {e}")
                total_errors += 1
                continue
            _register_ref(cast_df)
            result = _save_cast(cast_df, 'sbe', cnv_file.name)
            if result is True:
                total_processed += 1
            elif result is False:
                total_errors += 1
            else:
                total_skipped += 1

        # --- RBR processing ---
        for xlsx_file in xlsx_files:
            try:
                casts = rbr_cast(str(xlsx_file))
            except Exception as e:
                print(f"[{folder.name}] ERROR loading {xlsx_file.name}: {e}")
                total_errors += 1
                continue
            for i, cast_df in enumerate(casts):
                _register_ref(cast_df)
                result = _save_cast(cast_df, 'rbr', f"{xlsx_file.name}[{i}]")
                if result is True:
                    total_processed += 1
                elif result is False:
                    total_errors += 1
                else:
                    total_skipped += 1

        # --- CastAway processing ---
        # Each CastAway cast is matched to the closest-time RBR/SBE cast at
        # the same base station.  The matched reference tag (e.g. 'S1b') is
        # used as the CastAway station tag so repeat indices align across
        # instruments in the final dataset.
        castaway_matched_tags: set[str] = set()
        for csv_file in castaway_files:
            try:
                cast_df = castaway_cast(str(csv_file))
            except Exception as e:
                print(f"[{folder.name}] ERROR loading {csv_file.name}: {e}")
                total_errors += 1
                continue
            station_tag = _resolve_castaway_station_tag(cast_df, castaway_matched_tags)
            result = _save_cast(cast_df, 'castaway', csv_file.name,
                                station_tag_override=station_tag)
            if result is True:
                total_processed += 1
            elif result is False:
                total_errors += 1
            else:
                total_skipped += 1

        print()

    print(f"Done: {total_processed} processed, {total_skipped} skipped, {total_errors} errors.")


def fill_manual_depths(csv_file, ctd_nc_file, output_file=None):
    """
    For rows in a build_castaway_doc() CSV that were manually given a depth
    after being flagged as 'no_depth', look up the nearest-neighbour RBR
    salinity and temperature at that depth from the xarray dataset and fill
    them in.

    Only rows where ``ctd_source == 'no_depth'`` AND ``depth`` is not NaN are
    touched.  All other rows are left unchanged.  ``ctd_source`` is updated to
    ``'rbr_manual_depth'`` for every row that is successfully filled.

    Rows that still have no depth after manual editing remain as ``'no_depth'``
    and are printed as a reminder.

    Args:
        csv_file:    Path to the CSV produced by build_castaway_doc() after
                     manual depth entry.
        ctd_nc_file: Path to the NetCDF produced by create_ctd_dataset().
        output_file: Where to save the result.  Defaults to overwriting
                     csv_file in-place.

    Returns: updated DataFrame.
    """
    import xarray as xr

    csv_file = Path(csv_file)
    df = pd.read_csv(csv_file)
    df['date'] = pd.to_datetime(df['date']).dt.date

    needs_fill  = (df['ctd_source'] == 'no_depth') & df['depth'].notna()
    still_empty = (df['ctd_source'] == 'no_depth') & df['depth'].isna()

    if still_empty.any():
        print(f"\n{'='*55}")
        print(f"  {still_empty.sum()} row(s) still have no depth — skipped:")
        print(f"{'='*55}")
        for _, r in df[still_empty].iterrows():
            print(f"  {r['date']}  {r['station']}  bottle {int(r['bottle_number'])}")
        print(f"{'='*55}\n")

    n = needs_fill.sum()
    if n == 0:
        print("No manually-filled depths found — nothing to update.")
        return df

    print(f"Filling {n} row(s) with manually-provided depth from RBR profiles...")

    ds = xr.open_dataset(ctd_nc_file)
    depth_grid     = ds.coords['depth'].values      if 'depth'      in ds.coords else np.array([])
    stations_ds    = ds.coords['station'].values    if 'station'    in ds.coords else []
    dates_ds       = ds.coords['date'].values       if 'date'       in ds.coords else []
    repeats_ds     = ds.coords['repeat'].values     if 'repeat'     in ds.coords else []
    instruments_ds = ds.coords['instrument'].values if 'instrument' in ds.coords else []

    updated = 0
    for idx, row in df[needs_fill].iterrows():
        # Parse base station and repeat index from raw station name
        # e.g. 'G1b' → base='G1', repeat=1;  'G1' → base='G1', repeat=0
        m = re.fullmatch(r'([A-Z0-9]+)([a-z]?)', str(row['station']))
        base_station = m.group(1) if m else str(row['station'])
        suffix       = m.group(2) if m else ''
        repeat_idx   = 0 if not suffix else ord(suffix) - ord('a') + 1

        if base_station not in stations_ds:
            print(f"  SKIP {row['station']} {row['date']} bottle {int(row['bottle_number'])}: "
                  f"station '{base_station}' not in dataset")
            continue

        date_matches = [d for d in dates_ds if pd.Timestamp(d).date() == row['date']]
        if not date_matches:
            print(f"  SKIP {row['station']} {row['date']} bottle {int(row['bottle_number'])}: "
                  f"date not in dataset")
            continue
        date_val = date_matches[0]

        # Try RBR first, then SBE — own repeat first, then others.
        sal_vals = temp_vals = depths_valid = None
        inst_used = None
        repeat_order = [repeat_idx] + [r for r in sorted(repeats_ds) if r != repeat_idx]
        for inst_try in [i for i in ['rbr', 'sbe'] if i in instruments_ds]:
            for rv in repeat_order:
                if rv not in repeats_ds:
                    continue
                try:
                    prof = ds.sel(station=base_station, date=date_val,
                                  repeat=rv, instrument=inst_try)
                except Exception:
                    continue
                sal  = prof['sal00'].values  if 'sal00'  in ds else np.full(len(depth_grid), np.nan)
                temp = prof['tv290C'].values if 'tv290C' in ds else np.full(len(depth_grid), np.nan)
                valid = ~(np.isnan(sal) | np.isnan(temp))
                if valid.sum() < 1:
                    continue
                sal_vals     = sal[valid]
                temp_vals    = temp[valid]
                depths_valid = depth_grid[valid]
                inst_used    = inst_try
                break
            if sal_vals is not None:
                break

        if sal_vals is None:
            print(f"  SKIP {row['station']} {row['date']} bottle {int(row['bottle_number'])}: "
                  f"no RBR or SBE data in dataset")
            continue

        manual_depth = float(row['depth'])
        if manual_depth == 0:
            # Surface placeholder: take all three values from the shallowest
            # RBR/SBE measurement rather than using depth=0 literally.
            i = 0
            fill_depth = float(depths_valid[i])
        else:
            i = int(np.argmin(np.abs(depths_valid - manual_depth)))
            fill_depth = manual_depth

        df.at[idx, 'depth']       = fill_depth
        df.at[idx, 'salinity']    = float(sal_vals[i])
        df.at[idx, 'temperature'] = float(temp_vals[i])
        df.at[idx, 'ctd_source']  = f'{inst_used}_manual_depth'
        updated += 1
        print(f"  Filled {row['station']}  {row['date']}  bottle {int(row['bottle_number'])}: "
              f"depth={fill_depth:.2f} m  "
              f"sal={sal_vals[i]:.4f}  temp={temp_vals[i]:.4f}  [{inst_used}]")

    ds.close()
    print(f"\n{updated} of {n} row(s) updated.")

    out_path = Path(output_file) if output_file else csv_file
    df.to_csv(out_path, index=False)
    print(f"Saved to {out_path}")
    return df

