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
STATIONS_FILE = PROJECT_ROOT / "Data" / "stations.csv"

CAST_MAP = {
    '2025Oct20': {
        'G1':       {'G1': 0},
        'allbutG1': {'G2': 0, 'S2': 2, 'S1': 3, 'B': 5},
    },
    '2025Nov18': {'G1': 0, 'G2': 1, 'S2': 3, 'S1': 5, 'B': 6},
    '2025Dec10': {'G1': 0, 'G2': 1, 'S2': 3, 'S1': 4, 'B': 5},
    '2026Jan12': {'G1': 0, 'G2': 1, 'S2': 2, 'S1': 3, 'B': 4},
    '2026Jan26': {'G1': 0, 'G2': 1, 'B': 2, 'S2': 3, 'S1': 4, 'S1b': 5, 'S1c': 6, 'S1d': 7, 'S1e': 8},
    '2026Feb17': {'G1': 0, 'G2': 5, 'S2': 1, 'S1': 2, 'B': 3},
    '2026Mar16': {'G1': 0, 'G2': 1, 'S2': 3, 'S1': 4, 'B': 5},
    '2026Apr16': {'G1': 0, 'G2': 2, 'S2': 3, 'S1': 4, 'B': 5},
    '2026May15': {'G1': 0, 'G2': 1},
    '2026Jun02': {'G1': 0, 'G2': 2, 'S2': 3, 'S1': 5, 'B': 7},
    '2026Jun09': {'G1': 0, 'G2': 1, 'S2': 2, 'S1': 3, 'B': 4},
    '2026Jun17': {'G1': 1, 'G2': 2, 'B': 3, 'S2': 4, 'S1': 5, 'S1b': 6, 'S1c': 7, 'S1d': 8, 'Bb': 9, 'G1b': 10},
}


def set_project_root(path):
    """Set custom project root path."""
    global PROJECT_ROOT, DATA_PATH, OUTPUT_PATH, DOC_FILE, STATIONS_FILE
    PROJECT_ROOT = Path(path)
    DATA_PATH = PROJECT_ROOT / "Data"
    OUTPUT_PATH = PROJECT_ROOT / "Output"
    DOC_FILE = PROJECT_ROOT / "DOC_info" / "DOCdepth_profiles.xlsx"
    STATIONS_FILE = PROJECT_ROOT / "Data" / "stations.csv"

