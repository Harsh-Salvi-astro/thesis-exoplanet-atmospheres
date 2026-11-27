"""
GJ 1132 b: opacity contributions plot from the TauREx retrieval.

This script is the version used for the MSc thesis.

What it does
1. Reads the model stored by TauREx in the retrieval output file
   (the "Spectra" group of solution0): the total best-fit model and the
   contribution of each molecule, cloud and scattering term.
2. Smooths each spectrum with a 400-point moving average. The same smoothing
   is used in the retrieval fit plot, so that the two figures match.
3. Crops the wavelength range to the first and last observed data point and
   plots the contributions in ppm on a log wavelength axis. No data points
   are drawn.
4. Saves gj1132b_opacity_contributions.png.

How to run (from the main repository folder)
   python scripts/gj1132b/04_opacity_plot.py

Inputs
   data/gj1132b/gj1132b_retrieval.h5
   (TauREx output; made with the command in taurex_params/run_commands.txt.
   Output files are not included in this repository because of their size.)
   data/gj1132b/gj1132b_spectrum.dat
   (made by 01_combine_data.py)
Output
   data/gj1132b/gj1132b_opacity_contributions.png
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.ticker as ticker

# ============================================================
# 1. PATHS
# ============================================================
h5_file = "data/gj1132b/gj1132b_retrieval.h5"
obs_file = "data/gj1132b/gj1132b_spectrum.dat"

# ============================================================
# 2. DEFINE BOUNDARIES (Strictly cutting at first/last data point)
# ============================================================
obs = pd.read_csv(obs_file, sep='\s+', names=['wl','depth','err','id'])
obs = obs.sort_values('wl')
start_wl = obs['wl'].min()
end_wl = obs['wl'].max()

# ============================================================
# 3. EXTRACT AND SMOOTH (Identical Logic to Fitting Plot)
# ============================================================
contrib_dict = {}
window_size = 400  # Must match the value in your fitting script

with h5py.File(h5_file, 'r') as f:
    sol_spectra = f['Output/Solutions/solution0/Spectra']
    native_wl = sol_spectra['native_wlgrid'][()]
    
    idx = np.argsort(native_wl)
    wl_sorted = native_wl[idx]
    
    # 1. Total Model (The Black Line)
    total_raw = sol_spectra['native_spectrum'][()]
    total_smooth = np.convolve(total_raw[idx], np.ones(window_size)/window_size, mode='same')
    
    # Apply Crop Mask to Total
    mask = (wl_sorted >= start_wl) & (wl_sorted <= end_wl)
    contrib_dict['Total Best-Fit'] = (wl_sorted[mask], total_smooth[mask])

    # 2. Individual Species
    base_path = 'Output/Solutions/solution0/Spectra/Contributions'
    
    def get_clean_contribs(name, obj):
        if 'native_spectrum' in name and isinstance(obj, h5py.Dataset):
            parts = name.split('/')
            if len(parts) > 2:
                # Label format matching your approved style
                clean_name = f"{parts[1]} {parts[0]}"
                
                raw_data = obj[()][idx]
                # Apply identical convolution
                smooth_data = np.convolve(raw_data, np.ones(window_size)/window_size, mode='same')
                
                # Apply identical crop
                m = (wl_sorted >= start_wl) & (wl_sorted <= end_wl)
                contrib_dict[clean_name] = (wl_sorted[m], smooth_data[m])

    f[base_path].visititems(get_clean_contribs)

# ============================================================
# 4. PLOTTING (Pure Opacity - NO POINTS)
# ============================================================
plt.figure(figsize=(15, 8))

# Define colors (Tab10 provides clear distinction)
colors = plt.cm.tab10(np.linspace(0, 1, len(contrib_dict)))

# Plot Total Model line in thick Black for the baseline
tx, ty = contrib_dict.pop('Total Best-Fit')
plt.plot(tx, ty * 1e6, color='black', lw=3, label='Total Best-Fit (Model)', zorder=10)

# Plot each individual species
for i, (name, (x, y)) in enumerate(sorted(contrib_dict.items())):
    plt.plot(x, y * 1e6, lw=2, label=name, color=colors[i], alpha=0.8)

# ============================================================
# 5. STYLING (Strictly matching your fitting plot axis)
# ============================================================
plt.xscale('log')

# CROP THE AXIS
plt.xlim(start_point if 'start_point' in locals() else start_wl, 
         end_point if 'end_point' in locals() else end_wl)

# Y zoom
plt.ylim(2150, 3100)

# Axis labels (size 16)
plt.xlabel("Wavelength ($\mu$m)", fontsize=16)
plt.ylabel("Transit Depth (ppm)", fontsize=16)

# Axis tick numbers (size 14)
plt.tick_params(axis='both', which='major', labelsize=14)

# Legend box (size 14)
plt.legend(loc='upper center',
           bbox_to_anchor=(0.5, -0.15),
           ncol=3,
           frameon=True,
           shadow=True,
           fontsize=14)

plt.grid(True, which='both', linestyle='--', alpha=0.2)
plt.tight_layout()


# ============================================================
# 6. SAVE
# ============================================================
plt.savefig("data/gj1132b/gj1132b_opacity_contributions.png", dpi=300, bbox_inches='tight')
plt.show()
