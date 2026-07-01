"""
Plotting functions for CTD data visualization.
"""
import numpy as np
import matplotlib.pyplot as plt
import gsw
from pathlib import Path
from .config import DATA_PATH, OUTPUT_PATH
from .io import sbe_cast, load_bottle_file


def plot_ts(down_df, bl_down_df=None):
    """Plot T-S diagram with density contours and salinity vs depth profile."""
    # Calculate axis ranges with padding
    smin = down_df.Salinity.min() * 0.99
    smax = down_df.Salinity.max() * 1.01
    tmin = down_df.Temperature.min() - 0.1 * down_df.Temperature.max()
    tmax = down_df.Temperature.max() * 1.1
    
    xdim = int(round((smax - smin) / 0.1 + 1))
    ydim = int(round(tmax - tmin + 1))
    
    # Create density grid (sigma-t)
    ti = np.linspace(1, ydim - 1, ydim) + tmin
    si = np.linspace(1, xdim - 1, xdim) * 0.1 + smin
    dens = np.zeros((ydim, xdim))
    for j in range(ydim):
        for i in range(xdim):
            dens[j, i] = gsw.rho(si[i], ti[j], 0) - 1000

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 10))
    
    # T-S diagram
    CS = ax1.contour(si, ti, dens, linestyles='dashed', colors='k')
    ax1.clabel(CS, fontsize=12, inline=True)
    ax1.plot(down_df.Salinity, down_df.Temperature, '.')
    ax1.set_xlabel('Salinity (PSU)')
    ax1.set_ylabel('Temperature (°C)')
    ax1.set_title('T-S Diagram')
    
    # Salinity profile
    ax2.plot(down_df.Salinity, down_df.Depth)
    ax2.set_xlabel('Salinity (PSU)')
    ax2.set_ylabel('Depth (m)')
    ax2.set_title('Salinity vs Depth')
    ax2.tick_params(axis='x', rotation=45)
    ax2.invert_yaxis()
    
    # Add bottle annotations
    if bl_down_df is not None:
        ax1.plot(bl_down_df.Salinity, bl_down_df.Temperature, '.')
        ax2.plot(bl_down_df.Salinity, bl_down_df.Depth, '.')
        for _, r in bl_down_df.iterrows():
            for ax, y_col in [(ax1, 'Temperature'), (ax2, 'Depth')]:
                ax.annotate(str(r["bottle_number"]), 
                           xy=(r["Salinity"], r[y_col]),
                           xytext=(5, 5), textcoords="offset points", fontsize=20)
    
    return fig


def from_file_to_plot(cnv_file, bl_file, doc_file=None):
    """Load CNV and BL files, plot T-S diagram. Returns figure."""
    down_df = sbe_cast(cnv_file)
    bl_df = load_bottle_file(bl_file, doc_file)
    return plot_ts(down_df, bl_df)


def plot_bl_files(data_path=None, output_path=None):
    """Process all BL files in Data folder, save T-S plots to Output folder."""
    data_path = Path(data_path) if data_path else DATA_PATH
    output_path = Path(output_path) if output_path else OUTPUT_PATH
    output_path.mkdir(exist_ok=True)
    
    bl_files = list(data_path.glob("**/*.bl"))
    
    for bl_file in bl_files:
        try:
            cnv_file = bl_file.with_name(bl_file.stem + "NTS.cnv")
            if not cnv_file.exists():
                print(f"Warning: No CNV file for {bl_file.name}")
                continue
            
            output_file = output_path / f"{bl_file.parent.name}_{bl_file.stem}_plot.png"
            if output_file.exists():
                print(f"Skipping {bl_file.name} - exists")
                continue
            
            print(f"Processing: {bl_file.name}")
            fig = from_file_to_plot(cnv_file, bl_file)
            fig.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved: {output_file}")
            
        except Exception as e:
            print(f"Error processing {bl_file}: {e}")
    
    print("Processing complete!")

