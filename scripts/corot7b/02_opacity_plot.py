"""
CoRoT-7 b: contribution of each opacity source to the TauREx best-fit
spectrum of the synthetic observation, plotted with the total best-fit model
(Fig. 4.3).

This script is the version used for the MSc thesis.

What it does
1. Reads the binned model spectrum (total best fit) and the binned spectrum of
   each contribution from the TauREx retrieval output file.
2. Plots the total as a thick black line and each contribution as a coloured
   line, and saves corot7b_opacity_contributions.pdf.

How to run (from the main repository folder)
   python scripts/corot7b/02_opacity_plot.py

Input
   data/corot7b/corot7b_retrieval.h5
   (TauREx output; made with the command in taurex_params/run_commands.txt.
   Output files are not included in this repository because of their size.)
Output
   data/corot7b/corot7b_opacity_contributions.pdf
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ============================================================
# 1. PATH
# ============================================================
h5_file = "data/corot7b/corot7b_retrieval.h5"

# ============================================================
# 2. EXTRACT DATA
# ============================================================
contrib_dict = {}

with h5py.File(h5_file, 'r') as f:
    # Get the binned wavelength grid
    sol_spectra = f['Output/Solutions/solution0/Spectra']
    wl = sol_spectra['binned_wlgrid'][()] 
    
    # NEW: Extract the TOTAL Best-Fit line (it is located outside the Contributions folder)
    total_model = sol_spectra['binned_spectrum'][()]

    # Extract individual species
    base_path = 'Output/Solutions/solution0/Spectra/Contributions'
    
    def get_clean_contribs(name, obj):
        if 'binned_spectrum' in name and isinstance(obj, h5py.Dataset):
            parts = name.split('/')
            if len(parts) > 2:
                # Example: Absorption/SiO/binned_spectrum → SiO Absorption
                clean_name = f"{parts[1]} {parts[0]}"
                if "Observed" not in clean_name:
                    contrib_dict[clean_name] = obj[()]

    f[base_path].visititems(get_clean_contribs)

# ============================================================
# 3. GENERATE THE PLOT
# ============================================================
plt.figure(figsize=(14, 8))

# Define colors for species
colors = plt.cm.tab10(np.linspace(0, 1, len(contrib_dict)))

# 1. Plot the TOTAL line in thick Black (Matches GJ 1132 b style)
plt.plot(wl, total_model * 1e6, color='black', lw=3, label='Total Best-Fit', zorder=10, alpha=0.9)

# 2. Plot individual species
for i, (name, spectrum) in enumerate(sorted(contrib_dict.items())):
    plt.plot(wl, spectrum * 1e6, lw=1.5, label=name, color=colors[i], alpha=0.8)

# ============================================================
# 4. STYLING (Labels set to normal weight as requested)
# ============================================================
plt.xscale('log')


# Set ticks for the range used in CoRoT-7 b

plt.xlim(np.min(wl), np.max(wl))

# Set Y-axis to match CoRoT-7 b signal (roughly 260 to 275 ppm)
plt.ylim(261, 275)

# Normal font weight for labels
plt.xlabel("Wavelength ($\mu$m)", fontsize=16)
plt.ylabel("Transit Depth (ppm)", fontsize=16)

plt.grid(True, which='both', linestyle='--', alpha=0.3)

plt.legend(
    loc='upper center',
    bbox_to_anchor=(0.5, -0.15),
    ncol=3,
    frameon=True,
    fontsize=14
)

plt.tight_layout()
plt.tick_params(axis='both', which='major', labelsize=14)
# ============================================================
# 5. SAVE FIGURE
# ============================================================
plt.savefig(
    "data/corot7b/corot7b_opacity_contributions.pdf",
    bbox_inches='tight'
)
print("Min wavelength:", np.min(wl))
print("Max wavelength:", np.max(wl))
plt.show()
