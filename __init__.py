"""
Padilla Bay CTD Processing Package

A Python package for processing SeaBird (CNV) and RBR (Excel) CTD data from
Padilla Bay field campaigns, including bottle files and DOC measurements.

Typical workflow
----------------
1. Call set_project_root() once to point the package at your data directory
   and load station coordinates (required for density computation).

2. Load a cast with sbe_cast() or rbr_cast() (or get_cast() for either type).

3. Process it with process_ctd() → 0.25 m binned CastFrame.

4. Optionally merge bottle data with load_bottle_file() and visualise with
   plot_ts(), or export to parquet with make_parquet().

5. Use batch_process_all() to run the full pipeline on every dated folder,
   or create_ctd_dataset() to assemble an xarray Dataset across all casts.

Quick start
-----------
>>> import pb_casts
>>> pb_casts.set_project_root('/path/to/PadillaBay')
>>>
>>> cast = pb_casts.sbe_cast('Data/2025Aug20/S1NTS.cnv')
>>> processed = pb_casts.process_ctd(cast)
>>>
>>> bottles = pb_casts.load_bottle_file(
...     'Data/2025Aug20/S1.bl',
...     'DOC_info/DOCdepth_profiles.xlsx'
... )
>>> fig = pb_casts.plot_ts(processed, bottles)
>>> fig.savefig('Output/S1_ts.png')
"""

__version__ = '0.1.0'
__author__ = 'Marisa Viola'

# Import key functions for convenient access
from .io import sbe_cast, rbr_cast, get_cast, load_bottle_file, station_from_path
from .processing import process_ctd, remove_surface_noise
from .plotting import plot_ts, from_file_to_plot, plot_bl_files
from .storage import make_parquet
from .batch import process_bl_doc_files, batch_process_all, review_surface_cutoffs, load_surface_cutoffs, plot_cutoff_check
from .dataset import create_ctd_dataset
from .utils import detect_precision, parse_folder_date, compute_density, CastFrame
from .config import set_project_root, PROJECT_ROOT, DATA_PATH, OUTPUT_PATH, DOC_FILE

__all__ = [
    # IO functions
    'sbe_cast',
    'rbr_cast',
    'get_cast',
    'load_bottle_file',
    'station_from_path',
    # Processing functions
    'process_ctd',
    'remove_surface_noise',
    
    # Plotting functions
    'plot_ts',
    'from_file_to_plot',
    'plot_bl_files',
    
    # Storage functions
    'make_parquet',
    
    # Batch processing
    'process_bl_doc_files',
    'batch_process_all',
    'review_surface_cutoffs',
    'load_surface_cutoffs',
    'plot_cutoff_check',
    
    # Dataset creation
    'create_ctd_dataset',
    
    # Utilities
    'detect_precision',
    'parse_folder_date',
    'compute_density',
    'CastFrame',
    
    # Configuration
    'set_project_root',
    'PROJECT_ROOT',
    'DATA_PATH',
    'OUTPUT_PATH',
    'DOC_FILE',
]

