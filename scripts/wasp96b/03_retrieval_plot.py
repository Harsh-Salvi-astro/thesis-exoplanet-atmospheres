"""
WASP-96 b: TauREx best-fit spectrum compared with the observed data, and the
reduced chi-square quoted in the thesis (Fig. 4.10).

This script is the version used for the MSc thesis.

What it does
1. Reads the binned model spectrum from the TauREx retrieval output file.
2. Reads the combined observed spectrum made by 01_combine_data.py.
3. Interpolates the model onto the observed wavelengths and computes
   chi-square. The reduced chi-square is chi-square divided by (number of data
   points - 1). It prints "Reduced Chi-square" and shows it in the plot title.
4. Plots the observed data with the model and saves wasp96b_retrieval_fit.png.

How to run (from the main repository folder)
   python scripts/wasp96b/03_retrieval_plot.py

Inputs
   data/wasp96b/wasp96b_retrieval.h5
   (TauREx output; made with the command in taurex_params/run_commands.txt.
   Output files are not included in this repository because of their size.)
   data/wasp96b/wasp96b_spectrum.dat   (made by 01_combine_data.py)
Output
   data/wasp96b/wasp96b_retrieval_fit.png
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.ticker as ticker

# ============================================================
# 1. PATHS
# ============================================================
h5_file = "data/wasp96b/wasp96b_retrieval.h5"
obs_file = "data/wasp96b/wasp96b_spectrum.dat"

# ============================================================
# 2. EXTRACT BINNED MODEL FROM TAUREX
# ============================================================
with h5py.File(h5_file, 'r') as f:
    sol = f['Output/Solutions/solution0/Spectra']
    model_wl = sol['binned_wlgrid'][()]
    model_depth = sol['binned_spectrum'][()]

# ============================================================
# 3. LOAD OBSERVATIONAL DATA
# ============================================================
obs = pd.read_csv(obs_file, sep='\s+', names=['wl','depth','err'])

wl = obs['wl'].values
depth = obs['depth'].values
err = obs['err'].values

# ============================================================
# 4. CALCULATE REDUCED CHI-SQUARE
# ============================================================

# Interpolate TauREx model to observation wavelengths
model_interp = np.interp(wl, model_wl, model_depth)

# Chi-square
chi2 = np.sum(((depth - model_interp) / err)**2)

# Degrees of freedom
dof = len(depth) - 1

# Reduced chi-square
chi2_red = chi2 / dof

print("Reduced Chi-square =", chi2_red)

# ============================================================
# 5. PLOTTING
# ============================================================
plt.figure(figsize=(14,7))

# Observed data
plt.errorbar(obs['wl'], obs['depth']*1e6,
             yerr=obs['err']*1e6,
             fmt='o',
             color='#1f77b4',
             label='Observed Data',
             alpha=0.5,
             markersize=4,
             elinewidth=1,
             capsize=2)

# TauREx binned model
plt.plot(model_wl,
         model_depth*1e6,
         color='black',
         lw=2,
         label='TauREx Binned Prediction',
         zorder=10)



# ============================================================
# 6. AXIS SETTINGS
# ============================================================
plt.xscale('log')

# Clean log tick labels
plt.gca().xaxis.set_major_formatter(ticker.ScalarFormatter())

# Axis limits
data_min, data_max = np.min(obs['depth'])*1e6, np.max(obs['depth'])*1e6
plt.ylim(data_min - 1000, data_max + 1500)
plt.xlim(0.35, 5.0)

# Axis labels
plt.xlabel("Wavelength ($\mu$m)", fontsize=16)
plt.ylabel("Transit Depth (ppm)", fontsize=16)

# Axis numbers size
plt.tick_params(axis='both', which='major', labelsize=14)

# ============================================================
# 7. TITLE WITH REDUCED CHI-SQUARE
# ============================================================
plt.title("WASP-96 b: Reduced $\\chi^2$ = {:.2f}".format(chi2_red), fontsize=16)

# Legend size
plt.legend(frameon=True, shadow=True, fontsize=14)

plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("data/wasp96b/wasp96b_retrieval_fit.png", dpi=300)
plt.show()
