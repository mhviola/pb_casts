"""
Batch processing functions for multiple CTD files.
"""
import re
from pathlib import Path
import numpy as np
import pandas as pd
from . import config
from .io import sbe_cast, rbr_cast, load_bottle_file, station_from_path
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
            print(cnv_file)
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


def batch_process_all(data_path=None, output_path=None, overwrite=False, surface_cutoffs_file=None):
    """
    Process all CTD casts in yearMON## folders (e.g. 2025Aug20).
    Handles SBE (CNV) and RBR (XLSX) data automatically.
    Saves one parquet file per cast to output_path/{folder}/{station}.parquet.

    Args:
        data_path:            Path to Data folder (default: DATA_PATH from config)
        output_path:          Where to save output parquet files (default: OUTPUT_PATH/casts)
        overwrite:            Re-process and overwrite existing parquet files (default: False)
        surface_cutoffs_file: Optional path to CSV produced by review_surface_cutoffs().
                              When provided, manual cutoff depths override the automatic
                              stability algorithm for each matched cast.
    """
    data_path = Path(data_path) if data_path else config.DATA_PATH
    output_path = Path(output_path) if output_path else config.OUTPUT_PATH / "casts"

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

        # SBE: any *.cnv in the folder or its subdirectories (e.g. "Station B.cnv"
        # from raw Seasave output or legacy "*NTS.cnv" files).
        # Files marked "DO NOT USE" or "redone" are excluded.
        cnv_files = sorted([
            f for f in folder.rglob('*.cnv')
            if 'do not use' not in f.name.lower() and 'redone' not in f.stem.lower()
        ])
        # RBR: .xlsx files
        xlsx_files = sorted(folder.glob('*.xlsx'))

        if not cnv_files and not xlsx_files:
            print(f"[{folder.name}] No CTD files found, skipping.\n")
            continue

        out_dir = output_path / folder.name
        out_dir.mkdir(parents=True, exist_ok=True)

        # --- SBE processing ---
        for cnv_file in cnv_files:
            station = station_from_path(str(cnv_file))
            if station is None:
                print(f"[{folder.name}] Skipping {cnv_file.name} — could not identify station.")
                total_skipped += 1
                continue
            out_file = out_dir / f"{station}.parquet"

            if out_file.exists() and not overwrite:
                print(f"[{folder.name}] {station}: already exists, skipping.")
                total_skipped += 1
                continue

            try:
                print(f"[{folder.name}] Processing SBE {station}...", end=' ')
                down_df = sbe_cast(str(cnv_file))
                cutoff_depth = cutoffs.get(f"{folder.name}_{station}")
                proc_df = process_ctd(down_df, surface_cutoff_m=cutoff_depth)
                proc_df['station'] = station
                if date is not None:
                    proc_df['cast_date'] = date.date()
                proc_df.to_parquet(out_file)
                print(f"saved ({len(proc_df)} depth bins) -> {out_file.name}")
                total_processed += 1
            except Exception as e:
                print(f"ERROR: {e}")
                total_errors += 1

        # --- RBR processing ---
        if xlsx_files:
            for xlsx_file in xlsx_files:
                try:
                    casts = rbr_cast(str(xlsx_file))
                except Exception as e:
                    print(f"[{folder.name}] ERROR loading {xlsx_file.name}: {e}")
                    total_errors += 1
                    continue

                for i, cast_df in enumerate(casts):
                    station = cast_df.cast_meta.get('station', f'cast_{i}')
                    out_file = out_dir / f"{station}.parquet"

                    if out_file.exists() and not overwrite:
                        print(f"[{folder.name}] {station}: already exists, skipping.")
                        total_skipped += 1
                        continue

                    try:
                        print(f"[{folder.name}] Processing RBR {station}...", end=' ')
                        cutoff_depth = cutoffs.get(f"{folder.name}_{station}")
                        proc_df = process_ctd(cast_df, surface_cutoff_m=cutoff_depth)
                        proc_df['station'] = station
                        if date is not None:
                            proc_df['cast_date'] = date.date()
                        proc_df.to_parquet(out_file)
                        print(f"saved ({len(proc_df)} depth bins) -> {out_file.name}")
                        total_processed += 1
                    except Exception as e:
                        print(f"ERROR: {e}")
                        total_errors += 1

        print()

    print(f"Done: {total_processed} processed, {total_skipped} skipped, {total_errors} errors.")

