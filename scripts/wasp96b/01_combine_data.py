"""
WASP-96 b: combine the published transmission spectra into one file for the
TauREx retrieval.

This script is the version used for the MSc thesis.

What it does
1. Reads three published tables (CSV files from the NASA Exoplanet Archive):
   VLT/FORS2 (Nikolov et al. 2018), JWST/NIRISS (Radica et al. 2023) and
   Spitzer/IRAC (Nikolov et al. 2022). Depths and errors are given in percent,
   so they are divided by 100.
2. Shifts the VLT and Spitzer points vertically so that their median transit
   depth matches the JWST median (JWST is the reference). This is done so that
   an offset between instruments is not read as a chemical signal.
3. Removes the VLT points that lie inside the JWST wavelength range.
4. Saves wasp96b_spectrum.dat with three columns: wavelength (micron),
   transit depth, error. No header.
5. Plots the combined spectrum.

How to run (from the main repository folder)
   python scripts/wasp96b/01_combine_data.py

Inputs (the published tables, not included in this repository because they
belong to the original authors; see data/README.md for the references)
   data/wasp96b/table_WASP-96-b-Nikolov-et-al.-2018 0.3-0.79.csv
   data/wasp96b/table_WASP-96-b-Radica-et-al.-2023._0.6-2.82.csv
   data/wasp96b/table_WASP-96-b-Nikolov-et-al.-2022 3.6-4.5.csv
Output
   data/wasp96b/wasp96b_spectrum.dat
"""

import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# 1. LOAD AND SCALE DATA (1.4% -> 0.014)
# ============================================================
# ID 0: VLT/FORS2 (Nikolov 2018) - Blue end/Sodium
vlt = pd.read_csv("data/wasp96b/table_WASP-96-b-Nikolov-et-al.-2018 0.3-0.79.csv")[['CENTRALWAVELNG', 'PL_TRANDEP', 'PL_TRANDEPERR1']]
vlt['id'] = 0

# ID 2: JWST/NIRISS (Radica 2023) - The Master Anchor
jwst = pd.read_csv("data/wasp96b/table_WASP-96-b-Radica-et-al.-2023._0.6-2.82.csv")[['CENTRALWAVELNG', 'PL_TRANDEP', 'PL_TRANDEPERR1']]
jwst['id'] = 2

# ID 3: Spitzer IRAC (Nikolov 2022) - Red end anchors
spitzer = pd.read_csv("data/wasp96b/table_WASP-96-b-Nikolov-et-al.-2022 3.6-4.5.csv")[['CENTRALWAVELNG', 'PL_TRANDEP', 'PL_TRANDEPERR1']]
spitzer['id'] = 3

# Combine and convert to Fractional Transit Depth
df = pd.concat([vlt, jwst, spitzer])
df.columns = ['wl', 'depth', 'err', 'id']
df['depth'] /= 100
df['err'] /= 100

# ============================================================
# 2. CALCULATE AND APPLY OFFSETS (Anchor = JWST ID 2)
# ============================================================
medians = df.groupby('id')['depth'].median()
jwst_anchor = medians[2]

offsets = {
    0: jwst_anchor - medians[0],
    2: 0.0,
    3: jwst_anchor - medians[3]
}

# Apply vertical shifts to align with JWST baseline
df_corrected = df.copy()
for inst in offsets:
    df_corrected.loc[df_corrected['id'] == inst, 'depth'] += offsets[inst]

# ============================================================
# 3. REMOVE ONLY TRUE OVERLAP (Dynamic Boundary)
# ============================================================

# Find JWST minimum wavelength
jwst_min = df_corrected[df_corrected['id'] == 2]['wl'].min()

# Keep VLT only below the first JWST point
vlt_clean = df_corrected[
    (df_corrected['id'] == 0) &
    (df_corrected['wl'] < jwst_min)
]

# Keep all JWST
jwst_clean = df_corrected[df_corrected['id'] == 2]

# Keep all Spitzer
spitzer_clean = df_corrected[df_corrected['id'] == 3]

# Merge final spectrum
master_final = pd.concat([
    vlt_clean,
    jwst_clean,
    spitzer_clean
]).sort_values('wl')

# ============================================================
# 4. SAVE FOR TAUREX (3 Columns: wl, depth, err)
# ============================================================
master_final[['wl', 'depth', 'err']].to_csv(
    "data/wasp96b/wasp96b_spectrum.dat", 
    sep=' ', 
    index=False, 
    header=False
)

print("SUCCESS: wasp96b_spectrum.dat created.")
print(f"Offsets Applied: VLT={offsets[0]*1e6:+.1f}ppm, Spitzer={offsets[3]*1e6:+.1f}ppm")

# ============================================================
# 5. VISUALIZATION (Professional MSc Plot)
# ============================================================
plt.figure(figsize=(14, 7))

# Define colors and labels
colors = {0: '#1f77b4', 2: '#ff7f0e', 3: '#d62728'}
labels = {0: 'VLT/FORS2 (Nikolov 2018)', 2: 'JWST/NIRISS (Radica 2023)', 3: 'Spitzer (Nikolov 2022)'}

for inst_id in [0, 2, 3]:
    subset = master_final[master_final['id'] == inst_id]
    plt.errorbar(subset['wl'], subset['depth']*1e6, yerr=subset['err']*1e6, 
                 fmt='o', color=colors[inst_id], label=labels[inst_id], 
                 alpha=0.2, markersize=4, elinewidth=1)



plt.xscale('log')
plt.xlabel("Wavelength ($\mu$m)", fontsize=12)
plt.ylabel("Transit Depth (ppm)", fontsize=12)
plt.title("WASP-96 b: Combined Full-Picture Spectrum (Offset Corrected)", fontsize=14, pad=20)
plt.legend(loc='upper right', frameon=True, shadow=True)
plt.grid(True, which='both', linestyle='--', alpha=0.4)
plt.tight_layout()

# Show the plot
plt.show()
