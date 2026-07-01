"""
Padilla Bay CTD Processing Package

A Python package for processing SeaBird and RBR CTD data, bottle files, 
and DOC data from Padilla Bay field campaigns.

Example usage:
    >>> import pb_casts
    >>> 
    >>> # Load a cast
    >>> cast = pb_casts.sbe_cast('data/station_S1.cnv')
    >>> 
    >>> # Process CTD data
    >>> processed = pb_casts.process_ctd(cast)
    >>> 
    >>> # Load bottle data with DOC
    >>> bottles = pb_casts.load_bottle_file('data/station_S1.bl', 'DOC_info/DOCdepth_profiles.xlsx')
    >>> 
    >>> # Make a T-S plot
    >>> fig = pb_casts.plot_ts(processed, bottles)
    >>> fig.savefig('output/ts_diagram.png')
"""

__version__ = '0.1.0'
__author__ = 'Marisa Viola'

# Import key functions for convenient access
from .io import sbe_cast, rbr_cast, get_cast, load_bottle_file
from .processing import process_ctd, remove_pump_priming
from .plotting import plot_ts, from_file_to_plot, plot_bl_files
from .storage import make_parquet
from .batch import process_bl_doc_files, batch_process_all
from .dataset import create_ctd_dataset
from .utils import detect_precision, parse_folder_date
from .config import set_project_root, PROJECT_ROOT, DATA_PATH, OUTPUT_PATH, DOC_FILE

__all__ = [
    # IO functions
    'sbe_cast',
    'rbr_cast',
    'get_cast',
    'load_bottle_file',
    
    # Processing functions
    'process_ctd',
    'remove_pump_priming',
    
    # Plotting functions
    'plot_ts',
    'from_file_to_plot',
    'plot_bl_files',
    
    # Storage functions
    'make_parquet',
    
    # Batch processing
    'process_bl_doc_files',
    'batch_process_all',
    
    # Dataset creation
    'create_ctd_dataset',
    
    # Utilities
    'detect_precision',
    'parse_folder_date',
    
    # Configuration
    'set_project_root',
    'PROJECT_ROOT',
    'DATA_PATH',
    'OUTPUT_PATH',
    'DOC_FILE',
]

