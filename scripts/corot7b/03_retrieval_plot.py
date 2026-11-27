"""
CoRoT-7 b: TauREx best-fit spectrum compared with the synthetic observation
(Fig. 4.2).

This script is the version used for the MSc thesis.

What it does
1. Reads the binned model spectrum from the TauREx retrieval output file.
2. Reads the synthetic noisy spectrum made by 01_synthetic_spectrum_prt.py and
   averages it in groups of 5 points to make the plot easier to read (the
   error bar of a group is the mean error divided by sqrt(5)). The retrieval
   itself used the unbinned spectrum.
3. Plots the binned synthetic observation with the model and saves
   corot7b_retrieval_fit.pdf.

How to run (from the main repository folder)
   python scripts/corot7b/03_retrieval_plot.py

Inputs
   data/corot7b/corot7b_retrieval.h5
   (TauREx output; made with the command in taurex_params/run_commands.txt.
   Output files are not included in this repository because of their size.)
   data/corot7b/corot7b_synthetic_spectrum.dat
Output
   data/corot7b/corot7b_retrieval_fit.pdf
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.ticker as ticker

# ============================================================
# 1. PATHS
# ============================================================
h5_file = "data/corot7b/corot7b_retrieval.h5"
obs_file = "data/corot7b/corot7b_synthetic_spectrum.dat"

# ============================================================
# 2. EXTRACT TAUREX BINNED MODEL
# ============================================================
with h5py.File(h5_file, 'r') as f:
    sol = f['Output/Solutions/solution0/Spectra']
    
    model_wl = sol['binned_wlgrid'][()]
    model_depth = sol['binned_spectrum'][()]

# ============================================================
# 3. LOAD SYNTHETIC OBSERVATION
# ============================================================
obs = pd.read_csv(obs_file, sep='\s+', names=['wl','depth','err'])

wl = obs['wl'].values
depth = obs['depth'].values
err = obs['err'].values




# ============================================================
# BIN OBSERVATION FOR CLEAN VISUALIZATION
# ============================================================

bin_size = 5

n_bins = len(wl) // bin_size

wl_binned = wl[:n_bins*bin_size].reshape(n_bins, bin_size).mean(axis=1)
depth_binned = depth[:n_bins*bin_size].reshape(n_bins, bin_size).mean(axis=1)

# noise reduces when binning
err_binned = err[:n_bins*bin_size].reshape(n_bins, bin_size).mean(axis=1) / np.sqrt(bin_size)





# ============================================================
# 5. PLOT
# ============================================================
plt.figure(figsize=(18,7))

# binned observation
plt.errorbar(wl_binned,
             depth_binned*1e6,
             yerr=err_binned*1e6,
             fmt='o',
             color='#1f77b4',
             label='Synthetic Observation (binned)',
             alpha=0.7,
             markersize=6,
             elinewidth=1,
             capsize=2)

# TauREx best-fit model
plt.plot(model_wl,
         model_depth*1e6,
         color='black',
         lw=2.5,
         label='TauREx Best-Fit Model',
         zorder=10)

# ============================================================
# 6. AXIS SETTINGS
# ============================================================
# AXIS SETTINGS
plt.xscale('log')

plt.xlabel("Wavelength ($\mu$m)", fontsize=16)
plt.ylabel("Transit Depth (ppm)", fontsize=16)

plt.tick_params(axis='both', which='major', labelsize=14)

plt.xlim(np.min(wl_binned), np.max(wl_binned))

data_min, data_max = np.min(depth_binned)*1e6, np.max(depth_binned)*1e6
plt.ylim(data_min - 15, data_max + 15)

plt.legend(fontsize=14, frameon=True)

plt.grid(alpha=0.3)
# ============================================================
# 7. SAVE FIGURE
# ============================================================

plt.legend(frameon=True,fontsize=14)
plt.grid(alpha=0.3)

plt.tight_layout()

# SAVE FIGURE
plt.savefig(
    "data/corot7b/corot7b_retrieval_fit.pdf",
    bbox_inches='tight'
)
plt.show()
