
"""
Batch processing functions for multiple CTD files.
"""
import re
import pandas as pd
from pathlib import Path
from .config import DATA_PATH, OUTPUT_PATH, DOC_FILE
from .io import sbe_cast, rbr_cast, load_bottle_file
from .processing import process_ctd
from .utils import parse_folder_date


def process_bl_doc_files(data_path=None, output_path=None, doc_file=None):
    """
    Process all BL files, merge with DOC data, save to CSV.
    
    Args:
        data_path: Path to data directory (default: DATA_PATH from config)
        output_path: Path to output directory (default: OUTPUT_PATH from config)
        doc_file: Path to DOC excel file (default: DOC_FILE from config)
    """
    data_path = Path(data_path) if data_path else DATA_PATH
    output_path = Path(output_path) if output_path else OUTPUT_PATH
    doc_file = Path(doc_file) if doc_file else DOC_FILE
    
    output_dir = output_path / "bl_doc_files"
    output_dir.mkdir(exist_ok=True)
    bl_files = list(data_path.glob("**/*.bl"))
    
    for bl_file in bl_files:
        try:
            cnv_file = bl_file.with_name(bl_file.stem + "NTS.cnv")
            if not cnv_file.exists():
                print(f"Warning: No CNV file for {bl_file.name}")
                continue
            
            output_file = output_dir / f"{bl_file.parent.name}_{bl_file.stem}_DOC.csv"
            if output_file.exists():
                print(f"Skipping {bl_file.name} - exists")
                continue
            
            bl_doc = load_bottle_file(bl_file, doc_file)
            bl_doc['station'] = bl_file.stem
            bl_doc.to_csv(output_file, index=False)
            print(f"Saved: {output_file}")
            
        except Exception as e:
            print(f"Error processing {bl_file}: {e}")
    
    print("Processing complete!")


def batch_process_all(data_path=None, output_path=None, station_file=None, overwrite=False):
    """
    Process all CTD casts in yearMON## folders (e.g. 2025Aug20).
    Handles SBE (CNV) and RBR (XLSX) data automatically.
    Saves one parquet file per cast to output_path/{folder}/{station}.parquet.

    Args:
        data_path: Path to Data folder (default: DATA_PATH from config)
        output_path: Where to save output parquet files (default: OUTPUT_PATH/casts)
        station_file: Path to station coordinates CSV (default: data_path/station_coordinates.csv)
        overwrite: Re-process and overwrite existing parquet files (default: False)
    """
    data_path = Path(data_path) if data_path else DATA_PATH
    output_path = Path(output_path) if output_path else OUTPUT_PATH / "casts"
    station_file = Path(station_file) if station_file else data_path / "station_coordinates.csv"

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

        # SBE: *NTS.cnv files (skip anything marked "DO NOT USE" or "redone")
        cnv_files = sorted([
            f for f in folder.glob('*NTS.cnv')
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
            station = cnv_file.stem.replace('NTS', '')
            out_file = out_dir / f"{station}.parquet"

            if out_file.exists() and not overwrite:
                print(f"[{folder.name}] {station}: already exists, skipping.")
                total_skipped += 1
                continue

            try:
                print(f"[{folder.name}] Processing SBE {station}...", end=' ')
                down_df = sbe_cast(str(cnv_file))
                proc_df = process_ctd(down_df)
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
            if station_file.exists():
                stations_df = pd.read_csv(station_file)
            else:
                print(f"[{folder.name}] Warning: station_coordinates.csv not found.")
                stations_df = pd.DataFrame({'name': []})

            for xlsx_file in xlsx_files:
                try:
                    casts = rbr_cast(str(xlsx_file), stations_df)
                except Exception as e:
                    print(f"[{folder.name}] ERROR loading {xlsx_file.name}: {e}")
                    total_errors += 1
                    continue

                for i, cast_df in enumerate(casts):
                    station = cast_df._metadata.get('station', f'cast_{i}')
                    out_file = out_dir / f"{station}.parquet"

                    if out_file.exists() and not overwrite:
                        print(f"[{folder.name}] {station}: already exists, skipping.")
                        total_skipped += 1
                        continue

                    try:
                        print(f"[{folder.name}] Processing RBR {station}...", end=' ')
                        proc_df = process_ctd(cast_df)
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

