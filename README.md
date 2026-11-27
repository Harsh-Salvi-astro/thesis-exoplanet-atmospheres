# Atmospheric characterisation of exoplanets: MSc thesis analysis code

Code from my MSc thesis, "Exploration of Properties of Exoplanetary
Atmospheres" (MSc Astrophysics, St Xavier's College, Mumbai; supervisor
Dr Katherine Rawlins). The thesis studies three planets with transmission
spectroscopy, using forward models (petitRADTRANS) and Bayesian retrievals
(TauREx 3 with MultiNest):

| Planet | Role in the thesis |
|---|---|
| CoRoT-7 b | Synthetic benchmark. There is no measured spectrum, so a model spectrum with added noise (20 ppm) is used to test whether the retrieval recovers a known atmosphere. |
| GJ 1132 b | Archival HST, JWST, ground-based, Spitzer and TESS data combined and fitted. |
| WASP-96 b | Archival VLT, JWST and Spitzer data combined and fitted. |

The code was written during the thesis. It is a student project, and some
scripts keep the style and comments of the notebooks they came from.
The release tagged `thesis-submission` is the code as submitted with the thesis.

## Results quoted in the thesis (for reference)
- **GJ 1132 b:** retrieved mean molecular weight about 28 amu; reduced chi-square
  1.21 for the best forward model and 1.02 for the retrieval fit.
- **WASP-96 b:** retrieved mean molecular weight 2.32 amu; H2O, Na and a cloud
  deck found; reduced chi-square about 1.75 for the best forward model and 1.86
  for the retrieval fit. The full parameter table is Table 4.3 of the thesis.
- **CoRoT-7 b:** a synthetic benchmark, not a measurement. Script 01 prints a
  signal-to-noise of about 0.46 for the added 20 ppm noise, so this is a hard test.

## What is in the repository
```
scripts/
  corot7b/   01 synthetic spectrum, 02 opacity plot, 03 retrieval-fit plot
  gj1132b/   01 combine data, 02 forward model, 03 retrieval-fit plot, 04 opacity plot
  wasp96b/   01 combine data, 02 forward model, 03 retrieval-fit plot, 04 opacity plot
taurex_params/   TauREx parameter files and the commands used (run_commands.txt)
taurex_input/    TauREx cross-section files (not included, see data/README.md)
plots/thesis_figures/   the figures as printed in the thesis
data/            see data/README.md (data files are not included)
```
Scripts are run from the main repository folder, for example
`python scripts/wasp96b/03_retrieval_plot.py`. The header of each script says
what it reads and what it writes.

## Software
- petitRADTRANS 3.2.0 for the forward models (Python 3.10, `requirements_petitradtrans.txt`)
- TauREx 3.3.2 with MultiNest for the retrievals (Python 3.11, `requirements_taurex.txt`)
- numpy, scipy, matplotlib, pandas, h5py

I used two separate environments, one for each. Line lists as run in
petitRADTRANS: GJ 1132 b used H2O (HITEMP), CO2 (UCL-4000), CH4 (HITEMP), CO
(HITEMP) and N2-N2 collision-induced absorption. CoRoT-7 b used SiO (SiOUVenIR),
Na (NewAllard), K (Allard), Fe (Kurucz), Mg (Kurucz), CO (HITEMP) and CO2
(UCL-4000). The cross-section files used by TauREx are set in the parameter
files in `taurex_params/`.

## Things to know when reading the results
- **Re-running will not give identical numbers.** MultiNest is stochastic. The
  CoRoT-7 b noise is random and unseeded, which is why the noisy spectrum used
  in the thesis is included in the repository.
- **WASP-96 b.** The numbers and figures in the thesis (Table 4.3, Figs 4.10-4.12)
  come from a retrieval run with 1000 live points (February 2026), the setting
  in `taurex_params/wasp96b_params.txt`. A new run will not reproduce them
  exactly.
- **Reduced chi-square.** In the retrieval-fit scripts, chi-square is divided by
  (number of points - 1). In the forward-model scripts it is divided by the
  number of points, because those models were adjusted by hand and not fitted.
  The usual definition divides by (points - number of fitted parameters), so the
  values quoted should be read as indicative.
- **GJ 1132 b retrieval figure.** The model curve in the fit and opacity plots is
  smoothed with a 400-point moving average to make the plot readable. The
  reduced chi-square of 1.02 in the title of the fit figure is computed from
  that smoothed curve.
- **GJ 1132 b abundances.** CO2 was a free parameter in the retrieval, but no CO2 absorption appears in the best-fit opacity contributions.
- **CoRoT-7 b.** The synthetic spectrum includes Fe and Mg, but TauREx
  cross-section files for them were not available on the ExoMol website, where the cross-sections for the other species came from (see Section 3.4 of the thesis), so they are left out of the retrieval chemistry (SiO, Na, K, CO and CO2 only). The retrieval model is therefore not
  identical to the model that made the data.
- **Plots.** The plot scripts are the notebook cells used for the thesis figures.
  The figures as printed in the thesis are in `plots/thesis_figures/`.

## Updates after submission
(None yet. Any later change will be added as a new commit and described here.
The tag `thesis-submission` always points to the submitted code.)

## Licence
MIT, see `LICENSE`. If you use the code, please cite the thesis or this
repository.
