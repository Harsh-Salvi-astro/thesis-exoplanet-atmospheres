"""
GJ 1132 b: TauREx best-fit spectrum compared with the observed data,
and the reduced chi-square quoted in the thesis (Fig. 4.6).

This script is the version used for the MSc thesis.

What it does
1. Reads the observed spectrum and uses its first and last wavelength
   to crop the model.
2. Reads the model stored by TauREx in the retrieval output file
   (the "Spectra" group of solution0) and smooths it with a 400-point
   moving average. The smoothing is a plotting choice. The smoothed model is
   also the one used in the chi-square calculation below.
3. Interpolates the smoothed model onto the observed wavelengths and computes
   chi-square. The reduced chi-square is chi-square divided by (number of data
   points - 1). It prints "Reduced Chi-square" and shows it in the plot title.
4. Plots the observed data with the fitted model and saves
   gj1132b_retrieval_fit.png.

How to run (from the main repository folder)
   python scripts/gj1132b/03_retrieval_plot.py

Inputs
   data/gj1132b/gj1132b_retrieval.h5
   (TauREx output; made with the command in taurex_params/run_commands.txt.
   Output files are not included in this repository because of their size.)
   data/gj1132b/gj1132b_spectrum.dat
   (made by 01_combine_data.py)
Output
   data/gj1132b/gj1132b_retrieval_fit.png
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.ticker as ticker

# ============================================================
# 1. PATHS
# ============================================================
h5_file = (
    "data/gj1132b/"
    "gj1132b_retrieval.h5"
)

obs_file = (
    "data/gj1132b/"
    "gj1132b_spectrum.dat"
)

# ============================================================
# 2. LOAD OBSERVATION (To define the crop range)
# ============================================================
obs = pd.read_csv(obs_file, sep='\s+', names=['wl','depth','err','id'])
obs = obs.sort_values('wl')
start_wl = obs['wl'].min()
end_wl = obs['wl'].max()

# ============================================================
# 3. EXTRACT & CROP MODEL (To start/end exactly on data)
# ============================================================
with h5py.File(h5_file, 'r') as f:
    sol = f['Output/Solutions/solution0/Spectra']
    native_wl = sol['native_wlgrid'][()]
    native_depth = sol['native_spectrum'][()]

    # Sort
    idx = np.argsort(native_wl)
    native_wl, native_depth = native_wl[idx], native_depth[idx]

    # SMOOTHING (To remove high-res wiggles and match 'Binned' look)
    window = 400 # Increased for a very smooth look
    smooth_depth = np.convolve(native_depth, np.ones(window)/window, mode='same')

    # CROP THE MODEL: Keep only points between the first and last data point
    mask = (native_wl >= start_wl) & (native_wl <= end_wl)
    final_wl = native_wl[mask]
    final_depth = smooth_depth[mask]


# ============================================================
# 4. CALCULATE REDUCED CHI-SQUARE
# ============================================================

# Interpolate smoothed model onto observation wavelengths
model_interp = np.interp(obs['wl'], final_wl, final_depth)

# Chi-square
chi2 = np.sum(((obs['depth'] - model_interp) / obs['err'])**2)

# Degrees of freedom
dof = len(obs['depth']) - 1

# Reduced chi-square
chi2_red = chi2 / dof

print("Reduced Chi-square =", chi2_red)


# ============================================================
# 4. PLOTTING
# ============================================================
plt.figure(figsize=(14, 7))

# 1. Observed Data (Blue points, thin grey error bars)
plt.errorbar(obs['wl'], obs['depth']*1e6, yerr=obs['err']*1e6, 
             fmt='o', color='#1f77b4', label='Observed Data', 
             alpha=0.7, markersize=6, elinewidth=1, capsize=2)

# 2. Fitted Model (Black Line - NO WIGGLES, NO OVERFLOW)
# This line now starts and ends exactly where your points are
plt.plot(final_wl, final_depth*1e6, color='black', lw=2.5, 
         label='Fitted model', zorder=10)

# ============================================================
# 5. STYLING & AXIS
# ============================================================
plt.xscale('log')


# Tight crop on the axis to match the points
plt.xlim(start_wl - 0.05, end_wl + 0.1)

# Transit Depth in PPM
plt.ylim(2150, 3100) 

plt.xlabel("Wavelength ($\mu$m)", fontsize=16, )
plt.ylabel("Transit Depth (ppm)", fontsize=16,)

plt.tick_params(axis='both', which='major', labelsize=14)

# Legend box (size 14)
plt.legend(loc='upper center',
           bbox_to_anchor=(0.5, -0.15),
           ncol=3,
           frameon=True,
           shadow=True,
           fontsize=14)

plt.legend(frameon=True, shadow=True, loc='upper right', fontsize=14)
plt.grid(True, which='both', linestyle='--', alpha=0.3)
plt.title(f"GJ 1132 b: Reduced $\\chi^2$ = {chi2_red:.2f}", fontsize=16)
plt.tight_layout()

# ============================================================
# 6. SAVE
# ============================================================
plt.savefig("data/gj1132b/"
    "gj1132b_retrieval_fit.png", dpi=300)
plt.show()
