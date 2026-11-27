# Data

The data files are not stored in this repository, because most of them are
tables published by other authors, and the TauREx output files are very large.
The one exception is the noisy CoRoT-7 b spectrum (`corot7b/corot7b_synthetic_spectrum.dat`),
which I made myself and which cannot be regenerated exactly (the noise is random).

To run the scripts, make the sub-folders below and put the files in them.
Scripts are run from the main repository folder, and paths in the scripts start
with `data/`.

## CoRoT-7 b (`data/corot7b/`)
CoRoT-7 b has no measured transmission spectrum. The spectrum used is a model
with added noise (see `scripts/corot7b/01_synthetic_spectrum_prt.py`).
- `corot7b_synthetic_spectrum.dat`: included in this repository.
- `corot7b_retrieval.h5`: TauREx output (not included).

## GJ 1132 b (`data/gj1132b/`)
Published data used by `scripts/gj1132b/01_combine_data.py`:
- `table_GJ-1132-b-Bennett-et-al.-2025.csv`: Bennett et al. (2025), JWST
- `ajac75det2_Libby-Roberts et al 2022_mrt.txt`: Libby-Roberts et al. (2022), HST/WFC3
- `table_GJ-1132-b-Diamond-Lowe-et-al.-2018.csv`: Diamond-Lowe et al. (2018), ground-based optical
- Three single points typed into the script from the papers: Dittmann et al.
  (2017) Spitzer 4.5 micron, Southworth et al. (2017) K-band 2.15 micron, and
  the TESS point from Table 1 of Libby-Roberts et al. (2022)

The script makes `gj1132b_spectrum_raw.dat` (before the offsets) and
`gj1132b_spectrum.dat` (offset-corrected, used in the fits).
TauREx output: `gj1132b_retrieval.h5` (not included).

## WASP-96 b (`data/wasp96b/`)
Published tables used by `scripts/wasp96b/01_combine_data.py`
(CSV files from the NASA Exoplanet Archive):
- `table_WASP-96-b-Nikolov-et-al.-2018 0.3-0.79.csv`: Nikolov et al. (2018), VLT/FORS2
- `table_WASP-96-b-Radica-et-al.-2023._0.6-2.82.csv`: Radica et al. (2023), JWST/NIRISS
- `table_WASP-96-b-Nikolov-et-al.-2022 3.6-4.5.csv`: Nikolov et al. (2022), Spitzer/IRAC

The script makes `wasp96b_spectrum.dat` (wavelength in micron, transit depth
as a fraction, error). TauREx output: `wasp96b_retrieval.h5` (not included).

## TauREx cross-section files (`taurex_input/`)
The TauREx parameter files read the molecular cross-section (xsec) and
collision-induced absorption (CIA) files from `taurex_input/<planet>/xsec` and
`taurex_input/<planet>/cia`. These are large files from public databases (molecular cross-sections from
ExoMol, CIA from HITRAN, as in Section 3.4 of the thesis) and are not included here.

## References
Please cite the original papers if you use these data. The full references are
in the reference list of the thesis.
