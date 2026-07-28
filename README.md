# pb_casts — Padilla Bay CTD Processing

A Python package for processing SeaBird and RBR CTD data from Padilla Bay field campaigns, including bottle data and DOC (Dissolved Organic Carbon) measurements.

## Features

- **Multi-instrument support**: Process both SeaBird (CNV) and RBR (Excel) CTD data
- **Automated quality control**: Remove surface artifacts, despike, filter
- **Bottle data integration**: Merge CTD profiles with discrete bottle samples
- **DOC analysis**: Integrate DOC measurements with salinity profiles
- **Batch processing**: Process every cast in every dated folder in one call
- **Visualization**: T-S diagrams with density contours and depth profiles
- **Data export**: Parquet files with embedded month metadata

## Installation

### Local development

```bash
git clone https://github.com/yourusername/pb_casts.git
cd pb_casts
pip install -e .
```

`python-ctd` is not on PyPI — install it via conda or directly from GitHub before running `pip install`:

```bash
conda install -c conda-forge ctd
# or
pip install git+https://github.com/pyoceans/python-ctd.git
```

### Requirements

- Python 3.10+
- numpy, pandas, matplotlib, gsw, xarray, pyarrow, openpyxl
- python-ctd (conda-forge or GitHub)

## Setup

Call `set_project_root()` once at the top of any script.  This sets all path
constants and loads `station_coordinates.csv` — required before density is
computed.

```python
import pb_casts

pb_casts.set_project_root('/path/to/PadillaBay')

# Paths are now available:
print(pb_casts.PROJECT_ROOT)   # /path/to/PadillaBay
print(pb_casts.DATA_PATH)      # /path/to/PadillaBay/Data
print(pb_casts.OUTPUT_PATH)    # /path/to/PadillaBay/Output

# Or set via environment variable (before import):
# export PB_PROJECT_ROOT=/path/to/PadillaBay
```

Expected directory layout:

```
PadillaBay/
├── Data/
│   ├── station_coordinates.csv
│   ├── 2025Aug20/
│   │   ├── S1NTS.cnv
│   │   ├── S1.bl
│   │   └── ...
│   └── 2026Jun09/
│       └── jun09_2026.xlsx
├── DOC_info/
│   └── DOCdepth_profiles.xlsx
└── Output/
```

## Quick Start

```python
import pb_casts

pb_casts.set_project_root('/path/to/PadillaBay')

# Load and process a SeaBird cast
cast = pb_casts.sbe_cast('Data/2025Aug20/S1NTS.cnv')
processed = pb_casts.process_ctd(cast)

# Load bottle data with DOC
bottles = pb_casts.load_bottle_file(
    'Data/2025Aug20/S1.bl',
    'DOC_info/DOCdepth_profiles.xlsx'
)

# T-S diagram
fig = pb_casts.plot_ts(processed, bottles)
fig.savefig('Output/S1_ts.png')

# Save to parquet
pb_casts.make_parquet(processed, bottles)
```

## Usage Examples

### Load an RBR Excel export

The RBR records all stations in one file.  Cast indices are mapped to station
names in `config.CAST_MAP`; if your folder is already in the map, `rbr_cast()`
resolves the mapping automatically.

```python
import pb_casts

pb_casts.set_project_root('/path/to/PadillaBay')

# Returns a list of CastFrames, one per station
casts = pb_casts.rbr_cast('Data/2025Nov18/nov18_2025.xlsx')

for cast in casts:
    station = cast.cast_meta['station']
    processed = pb_casts.process_ctd(cast)
    print(f"{station}: {len(processed)} depth bins, "
          f"max depth {processed.index.max():.1f} m")
```

If the folder is not yet in `CAST_MAP`, pass a `recasts` dict explicitly:

```python
casts = pb_casts.rbr_cast(
    'Data/2026Jun09/jun09_2026.xlsx',
    recasts={'G1': 0, 'G2': 1, 'S2': 2, 'S1': 3, 'B': 4}
)
```

### Auto-detect instrument type

`get_cast()` returns a list of CastFrames regardless of whether the file is a
CNV or Excel, so you can loop over mixed datasets uniformly:

```python
for fpath in ['Data/2025Aug20/S1NTS.cnv', 'Data/2026Jun09/jun09_2026.xlsx']:
    for cast in pb_casts.get_cast(fpath):
        processed = pb_casts.process_ctd(cast)
        print(cast.cast_meta['station'], len(processed))
```

### Batch process all casts

Scans every `YYYYMon##` folder under `DATA_PATH`, processes SBE and RBR files,
and writes one `.parquet` per cast under `OUTPUT_PATH/casts/`:

```python
import pb_casts

pb_casts.set_project_root('/path/to/PadillaBay')
pb_casts.batch_process_all()

# Re-run and overwrite existing files:
pb_casts.batch_process_all(overwrite=True)
```

### Batch process bottle + DOC files

```python
import pb_casts

pb_casts.set_project_root('/path/to/PadillaBay')
pb_casts.process_bl_doc_files()
# Saves CSV files to OUTPUT_PATH/bl_doc_files/
```

### Create a multi-cast xarray Dataset

Assembles all casts into a 4D Dataset with dimensions
`(station, date, repeat, depth)`:

```python
import pb_casts

pb_casts.set_project_root('/path/to/PadillaBay')
ds = pb_casts.create_ctd_dataset()

print(ds)
# Dimensions: station, date, repeat, depth
# Data vars:  tv290C, sal00, sigma0, rho0, ...

ds.to_netcdf('Output/padilla_bay_ctd.nc')
```

## Package Structure

```
pb_casts/
├── __init__.py      # Public API
├── config.py        # Paths, CAST_MAP, set_project_root()
├── io.py            # Data loading: sbe_cast, rbr_cast, load_bottle_file
├── processing.py    # QC pipeline: process_ctd, remove_surface_noise
├── plotting.py      # Visualization: plot_ts, plot_bl_files
├── storage.py       # Export: make_parquet
├── batch.py         # Batch runners: batch_process_all, process_bl_doc_files
├── dataset.py       # xarray builder: create_ctd_dataset
└── utils.py         # CastFrame, compute_density, helpers
```

## Data Processing Pipeline

`process_ctd()` runs the following steps in order:

1. **Density computation** — TEOS-10 `sigma0` and `rho0` via gsw
2. **Surface noise removal** — find first stable density gradient window
3. **Despiking** — block-based MAD spike removal (python-ctd)
4. **Low-pass filter** — Butterworth filter, 0.15 s time constant
5. **Pressure check** — enforce monotonic depth
6. **Interpolation** — fill small internal gaps
7. **Depth binning** — resample to 0.25 m grid
8. **Hanning smoothing** — optional (default on)
9. **Precision rounding** — round to original sensor resolution

All tuning parameters live in `COMMON_PARAMS`, `SBE_PARAMS`, `RBR_PARAMS` at
the top of `processing.py`.

## API Reference

### Data Loading

| Function | Description |
|----------|-------------|
| `sbe_cast(cnv_file)` | Load SeaBird CNV file → downcast CastFrame |
| `rbr_cast(excel_file, recasts=None)` | Load RBR Excel export → list of CastFrames |
| `get_cast(file)` | Auto-detect CNV/XLSX and load → list of CastFrames |
| `load_bottle_file(bl_file, doc_file=None)` | Load bottle data, optionally merge DOC |
| `station_from_path(file)` | Infer station name from filename |

### Processing

| Function | Description |
|----------|-------------|
| `process_ctd(df, smooth=True, columns=None)` | Full QC pipeline → binned CastFrame |
| `remove_surface_noise(cast_df, method='stability', ...)` | Trim unstable surface data |

### Plotting

| Function | Description |
|----------|-------------|
| `plot_ts(down_df, bl_down_df=None)` | T-S diagram + salinity profile |
| `from_file_to_plot(cnv_file, bl_file, doc_file=None)` | Load files and plot |
| `plot_bl_files(data_path=None, output_path=None)` | Batch plot all BL files |

### Storage

| Function | Description |
|----------|-------------|
| `make_parquet(down_df, bl_down_df=None, output_dir=None)` | Save to parquet with metadata |

### Batch Processing

| Function | Description |
|----------|-------------|
| `batch_process_all(data_path=None, output_path=None, overwrite=False)` | Process all dated folders |
| `process_bl_doc_files(data_path=None, output_path=None, doc_file=None)` | Batch BL + DOC → CSV |

### Dataset Creation

| Function | Description |
|----------|-------------|
| `create_ctd_dataset()` | Build xarray Dataset from all casts in DATA_PATH |

### Utilities

| Function | Description |
|----------|-------------|
| `compute_density(df, ...)` | Add `rho0` and `sigma0` columns via TEOS-10 |
| `detect_precision(series)` | Infer decimal precision of a Series |
| `parse_folder_date(folder_name)` | Parse `2025Aug20` → pandas Timestamp |
| `godin_filter(data, dt_hours)` | Apply Godin (24h–24h–25h) low-pass filter to a time series |
| `CastFrame` | DataFrame subclass that preserves `cast_meta` through pandas ops |

### Configuration

| Symbol | Description |
|--------|-------------|
| `set_project_root(path)` | Set project root and load station coordinates |
| `PROJECT_ROOT` | Current project root Path |
| `DATA_PATH` | `PROJECT_ROOT / "Data"` |
| `OUTPUT_PATH` | `PROJECT_ROOT / "Output"` |
| `DOC_FILE` | Path to DOC Excel file |

## Project Background

This package supports CTD data collection and analysis for research in Padilla Bay, Washington, with a focus on halocline dynamics, freshwater–seawater mixing, and dissolved organic carbon distributions in an estuarine environment.

## Acknowledgments

Built on [python-ctd](https://github.com/pyoceans/python-ctd) by Filipe Fernandes and contributors, and the [Gibbs Seawater (gsw)](https://teos-10.github.io/GSW-Python/) TEOS-10 toolbox.

## Contact

Marisa Viola — Western Washington University
