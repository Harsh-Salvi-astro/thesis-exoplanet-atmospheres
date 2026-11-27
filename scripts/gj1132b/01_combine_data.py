"""
GJ 1132 b: combine the published transmission spectra into one file for the
forward model and the TauREx retrieval.

This script is the version used for the MSc thesis. It joins three notebook
cells, in the order I ran them.

What it does
1. Reads six sets of published data (transit depths are converted to a
   fraction, not percent):
     id 0  optical, Diamond-Lowe et al. (2018)
     id 1  HST/WFC3, Libby-Roberts et al. (2022) (table in MRT format; the first
           32 lines are the description of the columns and are skipped)
     id 2  JWST, Bennett et al. (2025)
     id 3  Spitzer 4.5 micron point, Dittmann et al. (2017)
     id 4  K-band 2.15 micron point, Southworth et al. (2017)
     id 5  TESS point, Libby-Roberts et al. (2022) Table 1
   The Spitzer, K-band and TESS points were typed in by hand from the papers
   (planet-to-star radius ratio, squared, with its error propagated).
   Saves gj1132b_spectrum_raw.dat (wavelength in micron, depth, error, id).
2. Computes the median depth of each data set and the shift needed to bring
   each one to the JWST median (JWST is the reference).
3. Applies the shift to every set except JWST and the K-band point, and
   saves gj1132b_spectrum.dat.

How to run (from the main repository folder)
   python scripts/gj1132b/01_combine_data.py

Inputs (the published tables, not included in this repository because they
belong to the original authors; see data/README.md for the references)
   data/gj1132b/table_GJ-1132-b-Bennett-et-al.-2025.csv
   data/gj1132b/ajac75det2_Libby-Roberts et al 2022_mrt.txt
   data/gj1132b/table_GJ-1132-b-Diamond-Lowe-et-al.-2018.csv
Outputs
   data/gj1132b/gj1132b_spectrum_raw.dat
   data/gj1132b/gj1132b_spectrum.dat
"""

import pandas as pd
import numpy as np

# ============================================================
# 1. PROCESS THE 6 DATASETS (CONVERT ALL TO FRACTION)
# ============================================================

# --- A. Bennett 2025 (JWST) - ID: 2 ---
df_ben = pd.read_csv("data/gj1132b/table_GJ-1132-b-Bennett-et-al.-2025.csv")
bennett = pd.DataFrame({
    'wl': df_ben['CENTRALWAVELNG'],
    'depth': df_ben['PL_TRANDEP'] / 100,  # % to Fraction
    'err': df_ben['PL_TRANDEPERR1'] / 100,
    'id': 2
})

# --- B. Libby-Roberts 2022 (HST) - ID: 1 ---
df_lr = pd.read_csv("data/gj1132b/ajac75det2_Libby-Roberts et al 2022_mrt.txt", sep='\s+', skiprows=32,
                    names=['WS','WE','V1','eV1','V2','eV2','V3','eV3','V4','eV4','V5','eV5','RpRs','eRpRs'])
hst = pd.DataFrame({
    'wl': (df_lr['WS'] + df_lr['WE']) / 2,
    'depth': df_lr['RpRs']**2,
    'err': 2 * df_lr['RpRs'] * df_lr['eRpRs'],
    'id': 1
})

# --- C. Diamond-Lowe 2018 (Optical) - ID: 0 ---
df_dl = pd.read_csv("data/gj1132b/table_GJ-1132-b-Diamond-Lowe-et-al.-2018.csv")
optical = pd.DataFrame({
    'wl': df_dl['CENTRALWAVELNG'],
    'depth': df_dl['PL_TRANDEP'] / 100,
    'err': df_dl['PL_TRANDEPERR1'] / 100,
    'id': 0
})

# --- D. Dittmann 2017 (Spitzer Anchor) - ID: 3 ---
# Rp/Rs = 0.0492 +/- 0.0008 at 4.5um
spitzer = pd.DataFrame({'wl': [4.5], 'depth': [0.0492**2], 'err': [2*0.0492*0.0008], 'id': [3]})

# --- E. Southworth 2017 (The K-band Bridge) - ID: 4 ---
# Rp/Rs = 0.0526 +/- 0.0039 at 2.15um
bridge = pd.DataFrame({'wl': [2.15], 'depth': [0.0526**2], 'err': [2*0.0526*0.0039], 'id': [4]})

# --- F. TESS Broadband Point (Libby-Roberts 2022 Table 1) - ID: 5 ---
# Rp/Rs = 0.04913 +/- 0.00130 at 0.786um
tess = pd.DataFrame({
    'wl': [0.786],
    'depth': [0.04913**2],
    'err': [2*0.04913*0.00130],
    'id': [5]
})

# ============================================================
# 2. COMBINE AND SAVE (NO OFFSETS APPLIED YET)
# ============================================================
master_obs = pd.concat([optical, hst, bennett, spitzer, bridge, tess]).sort_values('wl')

# Save as .dat (Space separated, no header) - The Official TauREx/pRT format
master_obs.to_csv("data/gj1132b/gj1132b_spectrum_raw.dat", sep=' ', index=False, header=False)

print("SUCCESS: gj1132b_spectrum_raw.dat created with 6 datasets (including TESS)(0.47 to 5.2 um)")


# ============================================================
# 3. OFFSETS RELATIVE TO JWST
# ============================================================

# Load raw master file (now includes ID 5)
df = pd.read_csv(
    "data/gj1132b/gj1132b_spectrum_raw.dat",
    sep='\s+',
    names=['wl','depth','err','id']
)

# Convert to ppm for clarity
df['depth_ppm'] = df['depth'] * 1e6

# 1. Compute median per instrument
medians = df.groupby('id')['depth_ppm'].median()

# 2. DEFINE THE ANCHOR FIRST (Move this line up!)
jwst_median = medians[2]

# 3. Now calculate offsets
offsets = {}
for inst in medians.index:
    offsets[inst] = jwst_median - medians[inst]

print("Instrument Offsets (ppm relative to JWST):")
print(offsets)

print("Median depths (ppm):")
print(medians)


# ============================================================
# 4. APPLY OFFSETS AND SAVE
# ============================================================

# Apply offsets (convert ppm back to fractional depth)

df_corrected = df.copy()

for inst in offsets:
    if inst != 2 and inst != 4:   # keep JWST fixed, ignore K-band
        shift_fraction = offsets[inst] / 1e6
        df_corrected.loc[df_corrected['id'] == inst, 'depth'] += shift_fraction

# Save corrected file
df_corrected[['wl','depth','err','id']].to_csv(
    "data/gj1132b/gj1132b_spectrum.dat",
    sep=' ',
    index=False,
    header=False
)

print("Offset-corrected file created.")
