"""
WASP-96 b: contribution of each opacity source to the TauREx best-fit
spectrum, plotted with the total best-fit model (Fig. 4.11).

This script is the version used for the MSc thesis.

What it does
1. Reads the binned model spectrum (total best fit) and the binned spectrum of
   each contribution (absorption by each molecule or atom, collision-induced
   absorption, Rayleigh scattering, clouds) from the TauREx output file.
2. Plots the total as a thick black line and each contribution as a coloured
   line, and saves wasp96b_opacity_contributions.png.

How to run (from the main repository folder)
   python scripts/wasp96b/04_opacity_plot.py

Input
   data/wasp96b/wasp96b_retrieval.h5
   (TauREx output; made with the command in taurex_params/run_commands.txt.
   Output files are not included in this repository because of their size.)
Output
   data/wasp96b/wasp96b_opacity_contributions.png
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# 1. PATH
h5_file = "data/wasp96b/wasp96b_retrieval.h5"

# 2. EXTRACT ONLY THE FINAL SOLUTION CONTRIBUTIONS
contrib_dict = {}
with h5py.File(h5_file, 'r') as f:

    sol_spectra = f['Output/Solutions/solution0/Spectra']

    # wavelength grid
    wl = sol_spectra['binned_wlgrid'][()]

    # -------------------------------
    # ADD TOTAL BEST-FIT MODEL
    # -------------------------------
    total_model = sol_spectra['binned_spectrum'][()]

    # We only look inside the Final Solutions folder
    base_path = 'Output/Solutions/solution0/Spectra/Contributions'
    
    def get_clean_contribs(name, obj):
        if 'binned_spectrum' in name and isinstance(obj, h5py.Dataset):
            parts = name.split('/')
            if len(parts) > 2:
                clean_name = f"{parts[1]} {parts[0]}"
                contrib_dict[clean_name] = obj[()]

    f[base_path].visititems(get_clean_contribs)

# 3. GENERATE THE CLEAR PLOT
plt.figure(figsize=(14, 8))

# ------------------------------------------------
# PLOT TOTAL BEST FIT (BLACK LINE)
# ------------------------------------------------
plt.plot(wl, total_model * 1e6,
         color='black',
         lw=3,
         label='Total Best-Fit Model',
         zorder=10)

# Define a nice color map so lines are easy to distinguish
colors = plt.cm.tab10(np.linspace(0, 1, len(contrib_dict)))

# Individual species (unchanged)
for i, (name, spectrum) in enumerate(contrib_dict.items()):
    plt.plot(wl, spectrum * 1e6, lw=2, label=name, color=colors[i])

# 4. STYLING (UNCHANGED)
plt.xscale('log')

plt.xlim(0.35, 5.0)

plt.ylim(13800, 15300)

plt.xlabel("Wavelength ($\mu$m)", fontsize=16, )
plt.ylabel("Transit Depth (ppm)", fontsize=16, )


plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3,
           frameon=True, shadow=True, fontsize=14)

plt.grid(True, which='both', linestyle='--', alpha=0.3)
plt.tight_layout()
plt.tick_params(axis='both', which='major', labelsize=14)

plt.savefig("data/wasp96b/wasp96b_opacity_contributions.png", dpi=300)
plt.show()
