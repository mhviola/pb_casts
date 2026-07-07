"""
Configuration and path management for Padilla Bay project.

Path constants default to the directory one level above this package so they
work out-of-the-box when the repo is cloned alongside the Data/Output folders.
Call set_project_root() (or export PB_PROJECT_ROOT) to point elsewhere.
"""
from pathlib import Path
import os
import pandas as pd

# ---------------------------------------------------------------------------
# Project paths — override via set_project_root() or PB_PROJECT_ROOT env var
# ---------------------------------------------------------------------------
_default_root = Path(__file__).parent.parent
PROJECT_ROOT   = Path(os.environ.get('PB_PROJECT_ROOT', _default_root))
DATA_PATH      = PROJECT_ROOT / "Data"
OUTPUT_PATH    = PROJECT_ROOT / "Output"
DOC_FILE       = PROJECT_ROOT / "DOC_info" / "DOCdepth_profiles.xlsx"
CASTAWAY_PATH  = PROJECT_ROOT / "Data" / "CastAway_profiles"

# STATIONS_FILE is only used internally by set_project_root to load STATIONS_DF.
# STATIONS_DF is what the rest of the package actually reads.
STATIONS_FILE = PROJECT_ROOT / "Data" / "station_coordinates.csv"
STATIONS_DF   = None  # populated by set_project_root(); required before computing density

# ---------------------------------------------------------------------------
# RBR cast index map — tells rbr_cast() which DOWN-cast index corresponds to
# which station for each field day.
#
# Format:  folder_name → {station_name: cast_index, ...}
#
# When a single day used two separate xlsx files (e.g. 2025Oct20 where G1 was
# recorded separately from the rest), the value is a nested dict keyed by a
# substring that appears in the xlsx filename stem:
#   folder_name → {filename_key: {station: cast_index, ...}, ...}
# ---------------------------------------------------------------------------
CAST_MAP = {
    '2025Oct20': {
        'G1':       {'G1': 0},
        'allbutG1': {'G2': 0, 'S2': 2, 'S1': 3, 'B': 5},
    },
    '2025Nov18': {'G1': 0, 'G2': 1, 'S2': 3, 'S1': 5, 'B': 6},
    '2025Dec10': {'G1': 0, 'G2': 1, 'S2': 3, 'S1': 4, 'B': 5},
    '2026Jan12': {'G1': 0, 'G2': 1, 'S2': 2, 'S1': 3, 'B': 4},
    '2026Jan26': {'G1': 0, 'G2': 1, 'B': 3, 'S2': 4, 'S1': 5, 'S1b': 6, 'S1c': 7, 'S1d': 8, 'S1e': 9},
    '2026Feb17': {'G1': 0, 'G2': 5, 'S2': 1, 'S1': 2, 'B': 3},
    '2026Mar16': {'G1': 0, 'G2': 1, 'S2': 3, 'S1': 4, 'B': 5},
    '2026Apr16': {'G1': 0, 'G2': 2, 'S2': 3, 'S1': 5, 'B': 6},
    '2026May15': {'G1': 0, 'G2': 1},
    '2026Jun02': {'G1': 0, 'G2': 2, 'S2': 3, 'S1': 5, 'B': 7},
    '2026Jun09': {'G1': 0, 'G2': 1, 'S2': 2, 'S1': 3, 'B': 4},
    '2026Jun17': {'G1': 1, 'G2': 2, 'B': 3, 'S2': 4, 'S1': 5, 'S1b': 6, 'S1c': 7, 'S1d': 8, 'Bb': 9, 'G1b': 10},
}


def set_project_root(path):
    """
    Override all project paths and load station coordinates.

    Must be called before any function that computes density (compute_density,
    process_ctd) because those functions look up station lat/lon from STATIONS_DF.

    Args:
        path: New project root directory (string or Path).  Expected layout:
              <path>/Data/               — raw CTD and bottle files
              <path>/Output/             — processed output
              <path>/DOC_info/           — DOC Excel file
              <path>/Data/CastAway_profiles/  — CastAway CSV files (Viola_* subfolders)
              <path>/Data/station_coordinates.csv  — station lat/lon table
    """
    global PROJECT_ROOT, DATA_PATH, OUTPUT_PATH, DOC_FILE, CASTAWAY_PATH, STATIONS_FILE, STATIONS_DF
    PROJECT_ROOT  = Path(path)
    DATA_PATH     = PROJECT_ROOT / "Data"
    OUTPUT_PATH   = PROJECT_ROOT / "Output"
    DOC_FILE      = PROJECT_ROOT / "DOC_info" / "DOCdepth_profiles.xlsx"
    CASTAWAY_PATH = PROJECT_ROOT / "Data" / "CastAway_profiles"
    STATIONS_FILE = PROJECT_ROOT / "Data" / "station_coordinates.csv"
    STATIONS_DF   = pd.read_csv(STATIONS_FILE)

