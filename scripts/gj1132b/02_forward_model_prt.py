"""
GJ 1132 b: forward model with petitRADTRANS, compared with the observed spectrum.

This script is the version used for the MSc thesis (Fig. 4.5).

What it does
1. Builds an isothermal atmosphere at T = 583.8 K on a pressure grid of
   1e-6 to 1e2 bar (120 layers). Planet and star values follow Xue et al. (2024).
2. Chemistry: nitrogen-dominated atmosphere (N2 is the remaining gas) with
   H2O = 7e-5 and CO2 = CH4 = CO = 1e-8 (volume mixing ratios).
   Several sets of abundances were tried by hand and the set with the lowest
   chi-square was kept. This is a manual trial, not a formal fit.
3. Adds a grey cloud (kappa = 1e-6) at pressures of 0.01 bar and higher.
4. Computes the transmission spectrum with petitRADTRANS 3.2.0.
5. Compares the model with the observed spectrum and prints
   "Reduced Chi-square". This is chi-square divided by the number of data
   points (not by points minus fitted parameters), because no formal fit
   was done here.
6. Saves the figures listed below.

How to run (from the main repository folder)
   python scripts/gj1132b/02_forward_model_prt.py

Inputs
   data/gj1132b/gj1132b_spectrum.dat
   (made by 01_combine_data.py; the published tables are not in this repository)
   petitRADTRANS opacity files (see the main README).
Outputs (folder results/gj1132b_forward/)
   forward_model_vs_data.png, TP_profile.png, abundances_vs_pressure.png,
   mmw_vs_pressure.png, opacity_contributions_transmission.png, line_opacities.png
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from petitRADTRANS.radtrans import Radtrans
import matplotlib.ticker as ticker
from pathlib import Path
from scipy.ndimage import gaussian_filter1d
from petitRADTRANS.plotlib import plot_opacity_contributions, plot_radtrans_opacities


# ============================================================
# 1. CONSTANTS & PARAMETERS (Xue et al. 2024)
# ============================================================
R_earth = 6.371e8
M_earth = 5.972e27
R_sun   = 6.957e10
G       = 6.67430e-8

planet_radius  = 1.192 * R_earth
planet_mass    = 1.84 * M_earth
stellar_radius = 0.2211 * R_sun
teq_val        = 583.8

gravity = G * planet_mass / planet_radius**2
pressures = np.logspace(-6, 2, 120)
temperature = np.full_like(pressures, teq_val)

OUT = Path("results/gj1132b_forward")
OUT.mkdir(parents=True, exist_ok=True)


# ============================================================
# 2. CHEMISTRY
# ============================================================
vmr = {
    'H2O': 7e-5,
    'CO2': 1e-8,
    'CH4': 1e-8,
    'CO':  1e-8
}

vmr['N2'] = 1.0 - sum(vmr.values())

weights = {
    'H2O': 18.015,
    'CO2': 44.01,
    'CH4': 16.04,
    'CO':  28.01,
    'N2':  28.014
}

MMW_val = sum(vmr[sp] * weights[sp] for sp in vmr)
MMW = np.full_like(pressures, MMW_val)

mass_fractions = {}
for sp in vmr:
    mass_fractions[sp] = np.full_like(pressures, vmr[sp] * weights[sp] / MMW_val)

print(f"Mean Molecular Weight (Rocky): {MMW_val:.3f}")
print(f"Surface Gravity: {gravity/100:.2f} m/s^2")


# ============================================================
# 3. CLOUD DECK
# ============================================================
P_cloud = 1e-2
kappa_cloud_value = 1e-6
kappa_cloud = np.zeros_like(pressures)
kappa_cloud[pressures >= P_cloud] = kappa_cloud_value


# ============================================================
# 4. INITIALIZE pRT
# ============================================================
atmosphere = Radtrans(
    pressures=pressures,
    line_species=['H2O', 'CO2', 'CH4', 'CO'],
    rayleigh_species=['N2'],
    gas_continuum_contributors=['N2--N2'],
    wavelength_boundaries=[0.6, 6.0]
)

wavelength_cm, transit_radius_cm, _ = atmosphere.calculate_transit_radii(
    temperatures=temperature,
    mass_fractions=mass_fractions,
    mean_molar_masses=MMW,   
    reference_gravity=gravity,
    planet_radius=planet_radius,
    reference_pressure=1.0,
    gray_opacity=kappa_cloud,
    return_contribution=True
)

wl_um = wavelength_cm * 1e4
model_depth = (transit_radius_cm / stellar_radius)**2


# ============================================================
# 5. LOAD DATA
# ============================================================
data = pd.read_csv(
    "data/gj1132b/gj1132b_spectrum.dat",
    sep='\s+',
    names=['wl','depth','err','id']
)

obs_wl = data['wl'].values
obs_depth = data['depth'].values
obs_err = data['err'].values


# ============================================================
# 6. CHI-SQUARE
# ============================================================
interp_model = interp1d(wl_um, model_depth, bounds_error=False, fill_value="extrapolate")
model_on_obs = interp_model(obs_wl)

chi2 = np.sum((obs_depth - model_on_obs)**2 / obs_err**2)
chi2_red = chi2 / len(obs_depth)

print(f"Reduced Chi-square: {chi2_red:.4f}")


# ============================================================
# 7. TRANSMISSION SPECTRUM PLOT
# ============================================================
plt.figure(figsize=(15,7))

plt.plot(
    wl_um,
    gaussian_filter1d(model_depth,3)*1e6,
    color='black',
    lw=1.5,
    label='pRT Model'
)

colors = {0:'blue',1:'cyan',2:'orange',3:'red',4:'green',5:'magenta'}

labels = {
0:'Optical (Diamond-Lowe 2018)',
1:'HST (Libby-Roberts 2022)',
2:'JWST (Bennett 2025)',
3:'Spitzer (Dittmann 2017)',
4:'K-band (Southworth 2017)',
5:'TESS (Libby-Roberts 2022)'
}

for inst_id in sorted(data['id'].unique()):
    subset=data[data['id']==inst_id]
    plt.errorbar(
        subset['wl'],
        subset['depth']*1e6,
        yerr=subset['err']*1e6,
        fmt='o',
        color=colors[inst_id],
        label=labels[inst_id],
        alpha=0.7,
        markersize=6,
        elinewidth=1,
        capsize=2
    )

# Molecular Bands
ax = plt.gca()

ax.axvspan(1.1,1.7,color='cyan',alpha=0.2,label='H$_2$O')
ax.axvspan(2.5,3.2,color='cyan',alpha=0.2)
ax.axvspan(3.2,3.5,color='purple',alpha=0.2,label='CH$_4$')
ax.axvspan(4.1,4.5,color='red',alpha=0.2,label='CO$_2$')
ax.axvspan(4.5,4.8,color='brown',alpha=0.2,label='CO')

# Log scale with 10^0 style
plt.xscale('log')


# Axis limits
plt.ylim(obs_depth.min()*1e6-200,obs_depth.max()*1e6+200)

# Axis labels (size 16)
plt.xlabel("Wavelength ($\mu$m)", fontsize=16)
plt.ylabel("Transit Depth (ppm)", fontsize=16)

# Axis numbers size
plt.tick_params(axis='both', which='major', labelsize=14)

# Title
plt.title(f"GJ 1132 b: Reduced $\\chi^2$ = {chi2_red:.2f}", fontsize=16)

# Legend size
plt.legend(ncol=2, fontsize=14)

plt.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(OUT/"forward_model_vs_data.png",dpi=300)

plt.show()
plt.close()


# ============================================================
# T–P PROFILE
# ============================================================
plt.figure(figsize=(5,6))

plt.semilogy(temperature,pressures,color='black',linewidth=2)

plt.gca().invert_yaxis()

plt.xlabel("Temperature [K]")
plt.ylabel("Pressure [bar]")

plt.title("Temperature–Pressure Profile")

plt.grid(alpha=0.3)

plt.tight_layout()

plt.savefig(OUT/"TP_profile.png",dpi=300)

plt.close()


# ============================================================
# ABUNDANCES
# ============================================================
plt.figure(figsize=(6,7))

for sp in vmr:
    if sp!='N2':
        plt.plot(mass_fractions[sp],pressures,label=sp,linewidth=2)

plt.xscale("log")
plt.yscale("log")

plt.gca().invert_yaxis()

plt.xlabel("Mass Fraction")
plt.ylabel("Pressure [bar]")

plt.title("Species Mass Fractions vs Pressure")

plt.legend()

plt.grid(alpha=0.3)

plt.tight_layout()

plt.savefig(OUT/"abundances_vs_pressure.png",dpi=300)

plt.close()


# ============================================================
# MMW PROFILE
# ============================================================
plt.figure(figsize=(5,6))
plt.semilogx(MMW, pressures, color='purple', linewidth=2)
plt.gca().invert_yaxis()
plt.xlabel("Mean Molecular Weight [amu]", fontsize=12)
plt.ylabel("Pressure [bar]", fontsize=12)
plt.title("Mean Molecular Weight vs Pressure", fontsize=13)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(OUT/"mmw_vs_pressure.png", dpi=300)
plt.close()


# ============================================================
# Opacity Contributions (Transmission)
# ============================================================
plot_opacity_contributions(
    atmosphere,
    mode="transmission",
    temperatures=temperature,
    mass_fractions=mass_fractions,
    mean_molar_masses=MMW,
    reference_gravity=gravity,
    planet_radius=planet_radius,
    reference_pressure=1.0
)

plt.gcf().set_size_inches(10,5)
plt.tight_layout()
plt.savefig(OUT/"opacity_contributions_transmission.png", dpi=300)
plt.close()


# ============================================================
# Line Opacities
# ============================================================
fig, ax = plt.subplots(figsize=(8,5))

plot_radtrans_opacities(
    atmosphere,
    atmosphere.line_species,
    temperature=teq_val,
    pressure_bar=0.1
)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Wavelength [μm]")
ax.set_ylabel("Cross-section [cm²]")
ax.set_title("Line Opacities at T_eq, P=0.1 bar")
plt.tight_layout()
plt.savefig(OUT/"line_opacities.png", dpi=300)
plt.close()
