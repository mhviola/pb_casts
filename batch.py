
"""
Batch processing functions for multiple CTD files.
"""
from pathlib import Path
from .config import DATA_PATH, OUTPUT_PATH, DOC_FILE
from .io import load_bottle_file


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

