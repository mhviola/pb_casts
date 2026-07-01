# pb_casts - Padilla Bay CTD Processing

A Python package for processing SeaBird and RBR CTD data from Padilla Bay field campaigns, including bottle data and DOC (Dissolved Organic Carbon) measurements.

## Features

- **Multi-instrument support**: Process both SeaBird (CNV) and RBR (Excel) CTD data
- **Automated quality control**: Remove pump priming artifacts, despike, filter
- **Bottle data integration**: Merge CTD profiles with discrete bottle samples
- **DOC analysis**: Integrate DOC measurements with salinity profiles
- **Batch processing**: Process multiple casts automatically
- **Visualization**: Create T-S diagrams and depth profiles
- **Data export**: Save processed data in Parquet format with metadata

## Installation

### From Git Repository

```bash
pip install git+https://github.com/yourusername/pb_casts.git
```

### Local Development Installation

```bash
git clone https://github.com/yourusername/pb_casts.git
cd pb_casts
pip install -e .
```

### Requirements

- Python 3.8+
- numpy
- pandas
- matplotlib
- ctd (python-ctd)
- gsw (gibbs seawater)
- xarray
- pyarrow
- openpyxl

## Quick Start

```python
import pb_casts

# Load and process a SeaBird cast
cast = pb_casts.sbe_cast('Data/2025Aug20/S1NTS.cnv')
processed = pb_casts.process_ctd(cast)

# Load bottle data with DOC measurements
bottles = pb_casts.load_bottle_file(
    'Data/2025Aug20/S1.bl',
    'DOC_info/DOCdepth_profiles.xlsx'
)

# Create T-S diagram
fig = pb_casts.plot_ts(processed, bottles)
fig.savefig('output/ts_diagram.png')

# Save to parquet
pb_casts.make_parquet(processed, bottles)
```

## Package Structure

```
pb_casts/
├── __init__.py      # Main package interface
├── config.py        # Configuration and paths
├── io.py            # Data loading (CNV, RBR, bottle files)
├── processing.py    # CTD processing pipeline
├── plotting.py      # Visualization functions
├── storage.py       # Data export (parquet)
├── batch.py         # Batch processing utilities
├── dataset.py       # xarray Dataset creation
└── utils.py         # Utility functions
```

## Usage Examples

### Process RBR Excel export

```python
import pb_casts
import pandas as pd

# Load station list
stations = pd.read_csv('Data/station_coordinates.csv')

# Process RBR file
casts = pb_casts.rbr_cast('Data/2025Nov18/nov_18_2025.xlsx', stations)

# Process each cast
for i, cast in enumerate(casts):
    processed = pb_casts.process_ctd(cast)
    print(f"Cast {i}: {len(processed)} depth bins")
```

### Batch process all bottle files

```python
import pb_casts

# Process all BL files in Data folder
pb_casts.plot_bl_files()

# Create CSV files with DOC data
pb_casts.process_bl_doc_files()
```

### Create multi-cast xarray Dataset

```python
import pb_casts

# Create dataset from all CNV files
ds = pb_casts.create_ctd_dataset(
    station_coords_file='Data/station_coordinates.csv'
)

# Save to NetCDF
ds.to_netcdf('output/padilla_bay_ctd.nc')
```

### Custom project paths

```python
import pb_casts

# Set custom project root
pb_casts.set_project_root('/path/to/your/project')

# Or use environment variable
import os
os.environ['PB_PROJECT_ROOT'] = '/path/to/your/project'
import pb_casts
```

## Data Processing Pipeline

The CTD processing pipeline includes:

1. **Pump priming removal**: Automatically detect and remove unstable surface data
2. **Despiking**: Remove sensor spikes using configurable block averaging
3. **Low-pass filtering**: Apply Butterworth filter with 0.15s time constant
4. **Pressure checks**: Validate monotonic depth progression
5. **Interpolation**: Fill small gaps in data
6. **Depth binning**: Resample to 0.25m depth bins
7. **Smoothing**: Apply Hanning window (optional, default=True)
8. **Precision matching**: Round to original sensor precision

## Configuration

Default paths are automatically detected based on package location. Override using:

```python
import pb_casts

# Set custom paths
pb_casts.set_project_root('/custom/path')

# Access current paths
print(pb_casts.PROJECT_ROOT)
print(pb_casts.DATA_PATH)
print(pb_casts.OUTPUT_PATH)
```

## API Reference

### Data Loading

- `sbe_cast(cnv_file)` - Load SeaBird CNV file
- `rbr_cast(excel_file, stations_df)` - Load RBR Excel export
- `get_cast(file, station_file)` - Auto-detect and load cast
- `load_bottle_file(bl_file, doc_file)` - Load bottle data with DOC

### Processing

- `process_ctd(df, smooth=True, columns=None)` - Process CTD data
- `remove_pump_priming(cast_df, method='stability')` - Remove surface artifacts

### Plotting

- `plot_ts(down_df, bl_down_df=None)` - Create T-S diagram
- `from_file_to_plot(cnv_file, bl_file, doc_file)` - Load and plot
- `plot_bl_files()` - Batch process all BL files

### Storage

- `make_parquet(down_df, bl_down_df=None)` - Save to parquet with metadata

### Batch Processing

- `process_bl_doc_files()` - Batch process bottle files with DOC

### Dataset Creation

- `create_ctd_dataset(station_coords_file)` - Create xarray Dataset from multiple casts

### Utilities

- `detect_precision(series)` - Auto-detect decimal precision
- `parse_folder_date(folder_name)` - Parse date from folder name

## Project Background

This package supports CTD data collection and analysis for research in Padilla Bay, Washington. The focus is on studying halocline dynamics, freshwater-seawater mixing, and dissolved organic carbon distributions in an estuarine environment.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

For academic and research use. Please cite appropriately if used in publications.

## Contact

Marisa Viola - Western Washington University

## Acknowledgments

This package builds on the excellent [python-ctd](https://github.com/pyoceans/python-ctd) library by Filipe Fernandes and contributors.

