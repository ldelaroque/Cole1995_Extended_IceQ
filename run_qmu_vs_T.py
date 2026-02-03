# ======================================================
# main.py – Attenuation model sensitivity study
# ======================================================

# --- Imports librairies ---
import numpy as np
import matplotlib.pyplot as plt

# --- Internal imports ---
from src.constants import R
from src.materials import saline_ice
from src.grids import temperature_grid
from src.attenuation_models import Qinv_Cammarano, Qinv_Tobie
from src.cole_sas_maxwell import combined_mechanisms
from src.rheology import get_mu   


# ======================================================
# ================ CONFIGURATION =======================
# ======================================================

# Temperature and frequency grids
T = temperature_grid(93.15, 268.15, 500)
f = 1.0                    # Hz
omega = 2 * np.pi * f

# Pressure names are defined such that:
# 1. Pressure at near-surface (0.1 MPa)
# 2. Pressure representative of ice-ocean boundary (~100 MPa)
P_low = 0.1   # MPa
P_high = 100  # MPa

# Conversion factor
conv = 96485.33212  # eV → J/mol 

# ======================================================
# =========== PARAMETER RANGES TO EXPLORE ==============
# ======================================================

param_sets = {
    # Dislocation (D)
    "alpha_d"  : np.linspace(0.50, 0.60, 30),
    "Q_d"      : np.linspace(0.55 * conv, 0.65 * conv, 30),
    "deltaD_d" : np.linspace(1.7e-10, 2.23e-8, 30),
    "B0"       : np.linspace(1.205e-11, 1.205e-9, 30),
    "D_u_d"    : np.linspace(5.0e-11, 5.0e-10, 30),

    # Grain boundary sliding (GBS)
    "alpha_gb" : np.linspace(0.50, 0.60, 30),
    "Q_gb"     : np.linspace(1.26 * conv, 1.38 * conv, 30),
    "deltaD_gb": np.linspace(3.0e-12, 9.0e-12, 30),
    "tau_gb0"  : np.linspace(8.0e-28, 2.7e-27, 30),
    "D_u_gb"   : np.linspace(1.5e-10, 5.0e-10, 30),

    # Proton reorientation (PR)
    "c_0"       : np.linspace(5e-8, 1e-7, 30),
    "tau_pr_0"  : np.linspace(1.7e-16, 6.9e-16, 30),
    "Em"        : np.linspace(30e3, 35e3, 30),
}

# Nominal PR values
PR_nominal = dict(c_0=5e-8, tau_pr_0=6e-16, Em=30e3)


# ======================================================
# ======= Cammarano & Tobie reference models ===========
# ======================================================

# --- Cammarano et al. (2006) ---
Qinv_Cammarano = Qinv_Cammarano(T, f)

# --- Tobie et al. (2025) / Cole simplified ---
Qinv_Tobie     = Qinv_Tobie(T, f)

# ======================================================
# ================ MAIN COMPUTATION ====================
# ======================================================

curves_LP = []
curves_HP = []

for param_name, values in param_sets.items():
    print(f"→ Sensitivity test: {param_name}")

    for val in values:
        params = saline_ice()
        params[param_name] = val

        # --- Mear-surface pressure ---
        Qmu, D1t, D2t = combined_mechanisms(P=0.1*1e6, T=T, f=f, params=params, mechanisms=("D", "GBS", "PR"))

        PR_vals = PR_nominal.copy()
        if param_name in PR_vals:
            PR_vals[param_name] = val

        Q_total = np.abs(1.0 / Qmu) 
        curves_LP.append(Q_total)

        # --- Ice-ocean boundary pressure ---
        Qmu_hp, D1t_hp, D2t_hp = combined_mechanisms(P=100*1e6, T=T, f=f, params=params, mechanisms=("D", "GBS", "PR"))
        Q_total_hp = np.abs(1.0 / Qmu_hp) 
        curves_HP.append(Q_total_hp)


curves_LP = np.array(curves_LP)
curves_HP = np.array(curves_HP)

mean_curve_lp = np.vstack([curves_LP]).mean(axis=0)
mean_curve_hp = np.vstack([curves_HP]).mean(axis=0)


# ================
# === PLOTTING ===
# ================

fig, ax = plt.subplots(figsize=(6, 4))

for i, c in enumerate(curves_LP):
    ax.plot(T, c, color="gray", alpha=0.2, label=r"$Q_\mu^{-1}$(P=0.1 MPa)" if i == 0 else None)

for i, c in enumerate(curves_HP):
    ax.plot(T, c, color="lightgray", alpha=0.2, linestyle="--", label=r"$Q_\mu^{-1}$(P=100 MPa)" if i == 0 else None)

ax.plot(T, mean_curve_lp, color="red", lw=2, label=r"$\langle Q_\mu^{-1} \rangle$ (P=0.1 MPa)")
ax.plot(T, mean_curve_hp, color="red", lw=2, ls="--", label=r"$\langle Q_\mu^{-1} \rangle$ (P=100 MPa)")

ax.plot(T, Qinv_Cammarano, "--", label="Cammarano et al. 2006", color="#00cc00")
ax.plot(T, Qinv_Tobie, "--", label="Tobie et al. 2025", color="#00ccff")

ax.set_xlabel("Temperature [K]", fontsize=14)
ax.set_ylabel(r"Shear attenuation $Q_\mu^{-1}$", fontsize=14)
ax.set_xlim([93.15, 268.15])
ax.set_ylim([0, 0.1])
ax.invert_xaxis()
ax.grid(True, linestyle=":", color='lightgray')
ax.legend(frameon=False)

plt.tight_layout()
plt.show()
