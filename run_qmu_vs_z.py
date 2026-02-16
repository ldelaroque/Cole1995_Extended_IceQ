# ======================================================
# main.py – Attenuation model sensitivity with depth
# ======================================================

# --- Imports librairies ---
import numpy as np
import matplotlib.pyplot as plt

# --- Internal imports ---
from src.materials import saline_ice
from src.attenuation_models import Qinv_Cammarano, Qinv_Tobie
from src.cole_sas_maxwell import combined_mechanisms

# ======================================================
# =========== PARAMETER RANGES TO EXPLORE ==============
# ======================================================

# Conversion factor. In Cole (1995), data are given in eV. So need a conversion
# for D and GBS activation energies
conv = 96485.33212  # eV → J/mol 

param_sets = {
    # Dislocation (D)
    "alpha_d"  : np.linspace(0.50, 0.60, 30),
    "E_d"      : np.linspace(0.55 * conv, 0.65 * conv, 30),
    "deltaD_d" : np.linspace(1.7e-10, 2.23e-8, 30),
    "B0"       : np.linspace(1.205e-11, 1.205e-9, 30),
    "D_u_d"    : np.linspace(5.0e-11, 5.0e-10, 30),

    # Grain boundary sliding (GBS)
    "alpha_gb" : np.linspace(0.50, 0.60, 30),
    "E_gb"     : np.linspace(1.26 * conv, 1.38 * conv, 30),
    "deltaD_gb": np.linspace(3.0e-12, 9.0e-12, 30),
    "tau_gb0"  : np.linspace(8.0e-28, 2.7e-27, 30),
    "D_u_gb"   : np.linspace(1.5e-10, 5.0e-10, 30),

    # Proton reorientation (PR)
    "c_0"       : np.linspace(5e-8, 1e-7, 30),
    "tau_pr_0"  : np.linspace(1.7e-16, 6.9e-16, 30),
    "Em"        : np.linspace(25e3, 35e3, 30), 
    "Ep"        : np.linspace(24e3, 48e3, 30), 
}

# Nominal PR values
PR_nominal = dict(c_0=5e-8, tau_pr_0=6.9e-16, Em=25e3, Ef=25e3, Ep=25e3)

# --- 1) Reading model ---

data = np.loadtxt("data/TitanProfile_MgSO4_10WtPct_Zb88km_pp.txt", skiprows=1)
P_profile = data[:, 0]*1e3      # in Pa, pressure
T_profile = data[:, 1]          # in K, temperature
r_profile = data[:, 2]          # in km, radius
params = saline_ice()

f = 1           # in Hz, frequency

R_surface = 2.57473e3  
depth = R_surface - r_profile

# --- 2) Calcul Q⁻¹(T) ---
curves = []

for param_name, values in param_sets.items():
    print(f"→ Sensitivity test: {param_name}")

    for val in values:
        params = saline_ice()
        params[param_name] = val

        qmu = []

        for i in range(len(T_profile)):

        # --- Mear-surface pressure ---
            Qmu, D1t, D2t = combined_mechanisms(
                P=P_profile[i], 
                T=T_profile[i], 
                f=f, params=params, 
                mechanisms=("D", "GBS", "PR")
            )

            Q_total = np.abs(1.0 / Qmu) 
            qmu.append(Q_total)
        qmu = np.array(qmu)
        curves.append(qmu)

curves = np.array(curves)

curve_mean = np.mean(curves, axis=0)
curve_min  = np.min(curves, axis=0)
curve_max  = np.max(curves, axis=0)

# Autres modèles
Qinv_cam   = Qinv_Cammarano(T_profile, f)
Qinv_tobie = Qinv_Tobie(T_profile, f)


# --- Figure ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 5), sharey=True)

# --- (1) Q⁻¹ ---
ax1.axhspan(0, 29.5, color='#e6ffff', alpha=0.2, zorder=0)
ax1.axhspan(29.5, 90, color='lightblue', alpha=0.4, zorder=0)
ax1.fill_betweenx(depth, curve_min, curve_max, color='#fff0e6', alpha=0.4, label=r"$Q_\mu^{-1} \in (\min,\max)$", zorder=1)

for i, c in enumerate(curves):
    ax1.plot(c, depth, color="gray", alpha=0.2, label=r"$Q_\mu^{-1}$(z,T,P)" if i == 0 else None)

ax1.plot(curve_mean, depth, color="r", linewidth=2, label=r'$\langle Q_\mu^{-1}\rangle(z,T,P)$', zorder=2)

ax1.plot(Qinv_cam, depth, '--', label=f'Cammarano et al. 2006', color='#00cc00', zorder=2)
ax1.plot(Qinv_tobie, depth, '--', label='Tobie et al. 2025', color='#00ccff', zorder=2)
ax1.set_ylabel("Depth below surface [km]", fontsize=14)
ax1.set_xlim([0, 0.1])
ax1.set_ylim([0, R_surface - 2.48670e3])
# ax1.invert_yaxis()
ax1.legend(frameon=False, fontsize=9)
ax1.secondary_xaxis('top').set_xlabel(r"Shear attenuation $Q_\mu^{-1}$", fontsize=14)

# --- (2) Q ---
ax2.axhspan(0, 29.5, color='#e6ffff', alpha=0.2, zorder=0)
ax2.axhspan(29.5, 90, color='lightblue', alpha=0.4, zorder=0)
ax2.fill_betweenx(depth, 1/curve_min, 1/curve_max, color='#fff0e6', alpha=0.4, label=r"$Q_\mu \in (\min,\max)$", zorder=1)

for i, c in enumerate(curves):
    ax2.semilogx(1/c, depth, color="gray", alpha=0.2, label=r"$Q_\mu$(z,T,P)" if i == 0 else None)

ax2.semilogx(1/curve_mean, depth, color="r", linewidth=2, label=r'$\langle Q_\mu\rangle(z,T,P)$', zorder=2)

ax2.semilogx(1/Qinv_cam, depth, '--', color='#00cc00', zorder=2)
ax2.semilogx(1/Qinv_tobie, depth, '--', color='#00ccff', zorder=2)
ax2.set_xlim([1, 1e4])
ax2.set_ylim([0, R_surface - 2.48670e3])
# ax1.invert_yaxis()
ax2.legend(frameon=False, loc='upper left')
ax2.secondary_xaxis('top').set_xlabel(r"Shear quality factor $Q_\mu$", fontsize=14)

plt.gca().invert_yaxis()
ax1.xaxis.set_visible(False)
ax2.xaxis.set_visible(False)
plt.tight_layout()
# plt.savefig("./qmu_vs_depth.pdf", format="PDF")
plt.show()
