"""
Configuration and path management for Padilla Bay project.
"""
from pathlib import Path
import os

# Default project paths - can be overridden by user
_default_root = Path(__file__).parent.parent
PROJECT_ROOT = Path(os.environ.get('PB_PROJECT_ROOT', _default_root))
DATA_PATH = PROJECT_ROOT / "Data"
OUTPUT_PATH = PROJECT_ROOT / "Output"
DOC_FILE = PROJECT_ROOT / "DOC_info" / "DOCdepth_profiles.xlsx"


def set_project_root(path):
    """Set custom project root path."""
    global PROJECT_ROOT, DATA_PATH, OUTPUT_PATH, DOC_FILE
    PROJECT_ROOT = Path(path)
    DATA_PATH = PROJECT_ROOT / "Data"
    OUTPUT_PATH = PROJECT_ROOT / "Output"
    DOC_FILE = PROJECT_ROOT / "DOC_info" / "DOCdepth_profiles.xlsx"

