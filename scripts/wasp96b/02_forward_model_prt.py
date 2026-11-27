"""
WASP-96 b: forward model with petitRADTRANS, compared with the observed spectrum.

This script is the version used for the MSc thesis (Fig. 4.9).

What it does
1. Builds an isothermal atmosphere at T = 1350 K on a pressure grid of
   1e-6 to 1e2 bar (130 layers). Planet and star values follow McGruder et al.
   (2023): 13.79 R_earth, 149 M_earth, stellar radius 1.055 R_sun.
2. Chemistry: hydrogen/helium atmosphere (H2:He = 0.83:0.17 of the remaining
   gas) with H2O = 10^-4.05, Na = 10^-4.0, K = 10^-7, CO2 = 10^-6.5,
   CO = 10^-4.5, CH4 = 1e-10, H2S = 1e-7, HCN = 1e-9 (volume mixing ratios).
   The abundances and the cloud were adjusted by hand over several trials to
   lower the chi-square (the model is labelled "Tuned" in the plot). This is a
   manual trial, not a formal fit.
3. Adds a grey cloud (kappa = 5e-2) at pressures of 0.02 bar and higher.
   (The comment next to that line in the code mentions 1e4; that is an old
   comment from an earlier trial. The value used is 5e-2.)
4. Computes the transmission spectrum with petitRADTRANS 3.2.0.
5. Compares the model with the observed spectrum and prints
   "Reduced Chi-Square". This is chi-square divided by the number of data
   points (not by points minus fitted parameters), because no formal fit
   was done here.

How to run (from the main repository folder)
   python scripts/wasp96b/02_forward_model_prt.py

Inputs
   data/wasp96b/wasp96b_spectrum.dat
   (made by 01_combine_data.py; the published tables are not in this repository)
   petitRADTRANS opacity files (see the main README).
Output (folder results/wasp96b_forward/)
   forward_model_vs_data.png
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from petitRADTRANS.radtrans import Radtrans
import matplotlib.ticker as ticker

# ============================================================
# 1. FIXED PARAMETERS (McGruder et al. 2023)
# ============================================================
R_earth = 6.371e8
M_earth = 5.972e27
R_sun   = 6.957e10
G       = 6.67430e-8

planet_radius  = 13.79 * R_earth   
planet_mass    = 149.0 * M_earth
stellar_radius = 1.055 * R_sun
teq_val        = 1350.0 

gravity = G * planet_mass / planet_radius**2
pressures = np.logspace(-6, 2, 130)
temperature = np.full_like(pressures, teq_val)


from pathlib import Path
from petitRADTRANS.plotlib import plot_opacity_contributions, plot_radtrans_opacities

OUT = Path("results/wasp96b_forward")
OUT.mkdir(parents=True, exist_ok=True)

# ============================================================
# 2. CHEMISTRY (Sub-solar H2O for better fit - Radica 2023)
# ============================================================
vmr = {
    'H2O': 10**-4.05,  
    'Na':  10**-4.0,   
    'K':   10**-7.0,
    'CO2': 10**-6.5,
    'CO':  10**-4.5,
    'CH4': 10**-10.0, 
    'H2S': 10**-7.0,  
    'HCN': 10**-9.0   
}

vmr_rem = 1.0 - sum(vmr.values())
vmr['H2'], vmr['He'] = vmr_rem * 0.83, vmr_rem * 0.17

weights = {
    'H2O': 18.01, 
    'Na':  22.99, 
    'K':   39.10, 
    'CO2': 44.01, 
    'CO':  28.01, 
    'CH4': 16.04, 
    'H2S': 34.08, 
    'HCN': 27.02, 
    'H2':  2.016, 
    'He':  4.002
}
MMW_val = sum(vmr[sp] * weights[sp] for sp in vmr)
MMW = np.full_like(pressures, MMW_val)
mass_fractions = {sp: np.full_like(pressures, vmr[sp] * weights[sp] / MMW_val) for sp in vmr}

# ============================================================
# 3. INITIALIZE & CLOUD DECK
# ============================================================
atmosphere = Radtrans(
    pressures=pressures,
    line_species=['H2O', 'Na', 'K', 'CO2', 'CO', 'CH4', 'H2S', 'HCN'], 
    rayleigh_species=['H2', 'He'],
    gas_continuum_contributors=['H2--H2', 'H2--He'],
    wavelength_boundaries=[0.3, 5.0]
)

P_cloud = 0.02   # 50 mbar
kappa_cloud = np.zeros_like(pressures)
# HIGH VALUE (1e4) creates the "Floor" needed to lower chi-square
kappa_cloud[pressures >= P_cloud] = 5e-2 

# CALCULATE
wavelength_cm, transit_radius_cm, _ = atmosphere.calculate_transit_radii(
    temperatures=temperature, mass_fractions=mass_fractions,
    mean_molar_masses=MMW, reference_gravity=gravity,
    planet_radius=planet_radius, 
    reference_pressure=0.01, 
    gray_opacity=kappa_cloud
)

wl_um = wavelength_cm * 1e4
model_depth = (transit_radius_cm / stellar_radius)**2

# ============================================================
# 4. LOAD DATA & PLOT
# ============================================================
data = pd.read_csv("data/wasp96b/wasp96b_spectrum.dat", sep='\s+', names=['wl','depth','err'])
data['id'] = data['wl'].apply(lambda x: 0 if x < 0.599 else (2 if x < 3.2 else 3))

interp_model = interp1d(wl_um, model_depth, bounds_error=False, fill_value="extrapolate")
model_on_obs = interp_model(data['wl'].values)
chi2_red = np.sum((data['depth'] - model_on_obs)**2 / data['err']**2) / len(data)

plt.figure(figsize=(15, 7))

plt.plot(wl_um, model_depth*1e6,
         color='black',zorder=10,lw=2,
         label='Tuned pRT Model', alpha=0.9)

colors = {0: '#1f77b4', 2: '#ff7f0e', 3: '#d62728'}
labels = {0: 'VLT', 2: 'JWST', 3: 'Spitzer'}

for i in [0, 2, 3]:
    sub = data[data['id'] == i]
    plt.errorbar(sub['wl'], sub['depth']*1e6,
                 yerr=sub['err']*1e6,
                 fmt='o',
                 color=colors[i],
                 label=labels[i],
                 alpha=0.4,
                 markersize=6,
                 elinewidth=1,
                 capsize=2)

# Log wavelength axis
plt.xscale('log')




# Axis label sizes
plt.xlabel("Wavelength ($\\mu$m)", fontsize=16)
plt.ylabel("Transit Depth (ppm)", fontsize=16)

# Axis number sizes
plt.tick_params(axis='both', which='major', labelsize=14)

# Fixed Y-limit
all_depths = np.concatenate([model_depth, data['depth']])
plt.ylim(all_depths.min()*1e6 - 300,
         all_depths.max()*1e6 + 500)

plt.title(
    f"WASP-96 b: Reduced $\\chi^2$ = {chi2_red:.2f}",
    fontsize=16
)

# Legend size
plt.legend(fontsize=14)

plt.grid(alpha=0.25)

plt.tight_layout()
plt.savefig(OUT/"forward_model_vs_data.png", dpi=300)
plt.show()
plt.close()

print(f"Mean Molecular Weight: {MMW_val:.3f} amu")
print(f"Surface Gravity: {gravity/100:.2f} m/s^2")
print(f"Reduced Chi-Square: {chi2_red:.4f}")


