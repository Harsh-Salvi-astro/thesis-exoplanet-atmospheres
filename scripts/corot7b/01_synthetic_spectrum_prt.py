"""
CoRoT-7 b: synthetic transmission spectrum with petitRADTRANS (benchmark case).

This script is the version used for the MSc thesis (Fig. 4.1).
CoRoT-7 b has no observed transmission spectrum. The spectrum made here is
a model with added random noise. It is a test of the method, not a measurement.

What it does
1. Builds an isothermal atmosphere at T = 1756 K on a pressure grid of
   1e-6 to 1e2 bar (130 layers). Planet 1.528 R_earth and 6.056 M_earth,
   star 0.864 R_sun.
2. Chemistry (volume mixing ratios): SiO 0.30, Na 0.08, K 0.005, Fe 0.02,
   Mg 0.02, CO 0.05, CO2 0.01, and O2 as the remaining gas.
3. Computes the transmission spectrum (0.3-15 micron), smooths it to R = 1500,
   adds Gaussian noise of 20 ppm, and smooths the noisy spectrum to R = 500.
4. Saves plots and the noisy spectrum as corot7b_synthetic_spectrum.dat
   (wavelength in micron, transit depth as a fraction, error 20e-6).
5. Prints scale height, signal size and signal/noise (about 0.46).

The noise is random and no seed is set, so running the script again gives a
different noisy spectrum. The file used in the thesis retrievals is included
in this repository as data/corot7b/corot7b_synthetic_spectrum.dat.
This script writes to results/corot7b_forward/ so that it does not overwrite it.

How to run (from the main repository folder)
   python scripts/corot7b/01_synthetic_spectrum_prt.py

Inputs
   petitRADTRANS opacity files (see the main README).
Outputs (folder results/corot7b_forward/)
   transmission_spectrum.png, TP_profile.png, abundances_vs_pressure.png,
   mmw_profile.png, opacity_contributions_transmission.png, line_opacities.png,
   corot7b_synthetic_spectrum.dat
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from petitRADTRANS.radtrans import Radtrans
from petitRADTRANS.plotlib import plot_opacity_contributions, plot_radtrans_opacities
from scipy.ndimage import gaussian_filter1d
import matplotlib.ticker as ticker

# ============================================================
# OUTPUT DIRECTORY
# ============================================================
OUT = Path("results/corot7b_forward")
OUT.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. PHYSICAL PARAMETERS (Calculated Dynamically)
# ============================================================
R_earth = 6.371e8
M_earth = 5.972e27
R_sun   = 6.96e10
G       = 6.67430e-8

planet_radius  = 1.528 * R_earth
planet_mass    = 6.056 * M_earth
stellar_radius = 0.864 * R_sun
teq_val        = 1756.0 

gravity = (G * planet_mass) / (planet_radius**2) 

pressures = np.logspace(-6, 2, 130)
temperature = np.full_like(pressures, teq_val)

# ============================================================
# 2. CHEMISTRY (Lava World Mix)
# ============================================================
vmr = {'SiO': 0.30, 'Na': 0.08, 'K': 0.005, 'Fe': 0.02, 'Mg': 0.02, 'CO': 0.05, 'CO2': 0.01}
vmr['O2'] = 1.0 - sum(vmr.values())
weights = {'SiO': 44.08, 'Na': 22.99, 'K': 39.10, 'Fe': 55.84, 'Mg': 24.30, 'CO': 28.01, 'CO2': 44.01, 'O2': 31.99}
MMW_val = sum(vmr[sp] * weights[sp] for sp in vmr)
MMW = np.full_like(pressures, MMW_val)
mass_fractions = {sp: np.full_like(pressures, vmr[sp] * weights[sp] / MMW_val) for sp in vmr}

# ============================================================
# 3. pRT CALCULATION
# ============================================================
atmosphere = Radtrans(
    pressures=pressures,
    line_species=['SiO','Na','K','Fe','Mg','CO','CO2'],
    rayleigh_species=['O2'],
    gas_continuum_contributors=['O2--O2'],
    wavelength_boundaries=[0.3, 15.0]
)

wl_cm, transit_radius_cm, _ = atmosphere.calculate_transit_radii(
    temperatures=temperature, mass_fractions=mass_fractions,
    mean_molar_masses=MMW, reference_gravity=gravity,
    planet_radius=planet_radius, reference_pressure=1.0
)

wl_um = wl_cm * 1e4
depth_clean_ppm = (transit_radius_cm / stellar_radius)**2 * 1e6

# ============================================================
# 4. CONVOLUTION & NOISE
# ============================================================
def convolve_log(wl, spec, R=1500):
    loglam = np.log(wl)
    sigma = (1.0 / R) / np.mean(np.diff(loglam))
    return gaussian_filter1d(spec, sigma=sigma)

depth_conv_ppm = convolve_log(wl_um, depth_clean_ppm)

noise_ppm = 20.0
raw_noisy_ppm = depth_conv_ppm + np.random.normal(0, noise_ppm, size=len(depth_conv_ppm))
depth_noisy_conv_ppm = convolve_log(wl_um, raw_noisy_ppm, R=500) 

# ============================================================
# 5. MAIN TRANSMISSION PLOT (SAVED ONLY)
# ============================================================
fig, ax = plt.subplots(figsize=(12, 6))

ax.axvspan(0.58, 0.60, color='orange', alpha=0.2, label='Na')
ax.axvspan(0.76, 0.78, color='green', alpha=0.2, label='K')
ax.axvspan(1.8, 2.8,   color='red', alpha=0.2, label='Fe/Mg')
ax.axvspan(4.15, 4.45, color='magenta', alpha=0.2, label='CO2')
ax.axvspan(4.6, 5.0,   color='brown', alpha=0.2, label='CO')
ax.axvspan(8.0, 11.5,  color='blue', alpha=0.2, label='SiO')

ax.plot(wl_um, depth_conv_ppm, color='black', lw=2, label='Clean Model (R=1500)')
ax.plot(wl_um, depth_noisy_conv_ppm, color='#ff7f0e', lw=1, alpha=0.8, label='Simulated Trace (20 ppm noise)')

ax.set_xscale('log')
ax.set_xlim(0.3, 15.0)

mean_d = np.mean(depth_conv_ppm)
ax.set_ylim(mean_d - 60, mean_d + 60)

# Axis labels (size 16)
ax.set_xlabel("Wavelength ($\mu$m)", fontsize=16)
ax.set_ylabel("Transit Depth (ppm)", fontsize=16)

# Axis tick numbers (size 14)
ax.tick_params(axis='both', which='major', labelsize=14)

# Legend box (size 14)
ax.legend(ncol=3, fontsize=14)

ax.grid(alpha=0.2)

plt.tight_layout()
plt.savefig(OUT/"transmission_spectrum.png", dpi=300)
plt.close()

# ============================================================
# 6. T-P PROFILE
# ============================================================
plt.figure(figsize=(5,6))
plt.semilogy(temperature, pressures)
plt.gca().invert_yaxis()
plt.xlabel("Temperature [K]")
plt.ylabel("Pressure [bar]")
plt.title("Temperature–Pressure Profile")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(OUT/"TP_profile.png", dpi=300)
plt.close()

# ============================================================
# 7. ABUNDANCES VS PRESSURE
# ============================================================
plt.figure(figsize=(6,7))
for sp in vmr:
    plt.plot(mass_fractions[sp], pressures, label=sp)

plt.xscale("log")
plt.yscale("log")
plt.gca().invert_yaxis()
plt.xlabel("Mass Fraction")
plt.ylabel("Pressure [bar]")
plt.legend(fontsize=8)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(OUT/"abundances_vs_pressure.png", dpi=300)
plt.close()

# ============================================================
# 8. MMW PROFILE
# ============================================================
plt.figure(figsize=(5,6))
plt.semilogx(MMW, pressures)
plt.gca().invert_yaxis()
plt.xlabel("Mean Molecular Weight (amu)")
plt.ylabel("Pressure [bar]")
plt.title("MMW Profile")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(OUT/"mmw_profile.png", dpi=300)
plt.close()

# ============================================================
# 9. OPACITY CONTRIBUTIONS (TRANSMISSION)
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
# 10. LINE OPACITIES
# ============================================================
fig, ax = plt.subplots(figsize=(8,5))

plot_radtrans_opacities(
    atmosphere,
    atmosphere.line_species,
    temperature=2000.0,
    pressure_bar=0.1
)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Wavelength [μm]")
ax.set_ylabel("Cross-section [cm²]")
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(OUT/"line_opacities.png", dpi=300)
plt.close()

# ============================================================
# 11. SAVE .dat FILE (NO HEADER, 3 COLUMNS)
# ============================================================

# Convert ppm back to fraction for TauREx
transit_depth_fraction = depth_noisy_conv_ppm * 1e-6

# Constant error (same noise level everywhere)
error_fraction = np.full_like(transit_depth_fraction, noise_ppm * 1e-6)

# Stack columns
output_array = np.column_stack([
    wl_um,
    transit_depth_fraction,
    error_fraction
])

# Save WITHOUT header
np.savetxt(
    OUT / "corot7b_synthetic_spectrum.dat",
    output_array,
    fmt="%.8e"
)

print("Transmission .dat file saved successfully.")
print("All plots saved successfully to:", OUT)

# Calculate Scale Height H = (kB * T) / (mu * m_u * g)
k_b = 1.3807e-16 # erg/K
m_u = 1.6605e-24 # g
H_cm = (k_b * teq_val) / (MMW_val * m_u * gravity)
H_km = H_cm / 1e5

signal_size = np.max(depth_conv_ppm) - np.min(depth_conv_ppm)

print("-" * 30)
print("DIAGNOSTICS SUMMARY:")
print("-" * 30)
print(f"Calculated Gravity: {gravity/100:.2f} m/s^2")
print(f"Mean Molecular Weight: {MMW_val:.2f} amu")
print(f"Scale Height: {H_km:.2f} km")
print(f"Atmosphere Signal (Max-Min): {signal_size:.2f} ppm")
print(f"Noise Level: {noise_ppm:.2f} ppm")
print(f"SNR (Signal/Noise): {signal_size/noise_ppm:.2f}")
print("-" * 30)