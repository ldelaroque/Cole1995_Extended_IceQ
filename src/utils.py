# -- Import librairies --
import numpy as np
import matplotlib.pyplot as plt
from math import pi

eta   = 1e22         # in Pa s, shear viscosity
V     = 10e-6        # in m3/mol, activation volume
P     = 100e6        # in Pa, pressure at the ice Ih-ocean interface
k_B   = 8.61e-5      # in eV/K, Boltzmann constant
R     = 8.314        # in J/(mol K), gas constant
# Constants and parameters (from the paper / Table 1)
t0 = 6.9e-16         # s (Tatibouet et al. 1981)
Em = 30e3            # J/mol, Bjerrum detects'energies of migratio)
Ef = 25e3            # J/mol, Bjerrum detects'energies of formation
Ep = 25e3            # J/mol (for tanδ_max empirical)
c0 = 1e-7            # reference concentration of extrinsic defects (we use 1e-7 as in Table 1; you can try 5e-8 as alternative)


mu_data = np.loadtxt("./data/TitanProfile_MgSO4_10WtPct_Zb88km_pp.txt", skiprows=1)
T_mu    = mu_data[:, 1]  
mu_vals = mu_data[:, 8] * 1e9

def Tm_p_Hirschmann2000(P_MPa):
    """
    Compute melting temperature (K) from pressure (MPa)
    using Hirschmann (2000) parameterization.
    """
    P_GPa = 1e-3 * P_MPa
    coeffs = [-5.140, 132.899, 1120.661]  # coefficients for polynomial in GPa
    # np.polyval takes highest degree first (like MatLab)
    Tm_oC = np.polyval(coeffs, P_GPa)
    return 273.15 + Tm_oC  # convert to Kelvin

def Tm_ice_linear(P_MPa):
    # Local linear approximation (valid for low pressures typical of glaciers)
    # Feistel & Wagner (2005) "A New Equation of State for H2O Ice Ih"
    # Clausius-Clapeyron slope
    clapeyron_slope = -0.074  # K / MPa (approx.), see eqn 21 of the paper
    return 273.15 + clapeyron_slope * P_MPa


def cole_model_compliances(T_C, f, params):
    """
    Calculate the compliances D1 and D2 for a given temperature (in °C) 
    and frequency f (Hz), according to the Cole (1995) model.
    """
    T = T_C + 273.15        # in K, temperature
    omega = 2 * np.pi * f   # in rad/s, angular frequency

    # µ depends now on temperature
    kE = get_mu(T)

    # Temperature-dependant dislocation drag
    B_T_d = params["B0"] * np.exp(params["Q_d"] / (R * T))   # partial second units
    # Relaxation time for dislocation sliding
    tau_d = B_T_d / params["K"]

    # Relaxation time for grain-boundary sliding
    tau_gb = params["tau_gb_0"] * np.exp(params["Q_gb"] / (R * T))

    # FIRST MECHANISM - DISLOCATION
    # Dislocation-based relaxation mechanism
    s_d = np.log(tau_d * omega)
    # Storage compliance
    D1_d = params["δD_d"] * (1 - (2/np.pi) * np.arctan(np.exp(params["α_d"] * s_d))) #+ 1/kE
    # Loss compliance
    D2_d = params["α_d"] * params["δD_d"] / (np.exp(params["α_d"] * s_d) + np.exp(-params["α_d"] * s_d)) #+ 1/(eta*omega)

    # SECOND MECHANISM - GBS
    # Grain boundary-based relaxation mechanism
    s_gb = np.log(tau_gb * omega)
    # Storage compliance (in Pa-1)
    D1_gb =  params["δD_gb"] * (1 - (2/np.pi) * np.arctan(np.exp(params["α_gb"] * s_gb))) #+ 1/kE
    # Loss compliance (in Pa-1)
    D2_gb = params["α_gb"] * params["δD_gb"] / (np.exp(params["α_gb"] * s_gb) + np.exp(-params["α_gb"] * s_gb)) #+ 1/(eta*omega)
    D1_t = D1_d + D1_gb + 1/kE
    D2_t = D2_d + D2_gb - 1/(eta*omega)

    return D1_t, D2_t, D1_d, D1_gb, D2_d, D2_gb

def cole_model_compliances_HP(T_C, f, params):
    """
    Calculate the compliances D1 and D2 for a given temperature (in °C) 
    and frequency f (Hz), according to the Cole (1995) model.
    """
    T = T_C + 273.15        # in K, temperature
    omega = 2 * np.pi * f   # in rad/s, angular frequency

    # µ(T) - temperature dpdt
    kE = get_mu(T)

    # Temperature-dependant dislocation drag
    B_T_d = params["B0"] * np.exp((params["Q_d"]+P*V) / (R * T))   # partial second units
    # Relaxation time for dislocation sliding
    tau_d = B_T_d / params["K"]

    # Relaxation time for GBS
    tau_gb = params["tau_gb_0"] * np.exp((params["Q_gb"]+P*V)/ (R * T))

    # FIRST MECHANISM - DISLOCATION
    # Dislocation-based relaxation mechanism
    s_d = np.log(tau_d * omega)
    # Storage compliance
    D1_d =  params["δD_d"] * (1 - (2/np.pi) * np.arctan(np.exp(params["α_d"] * s_d))) #+ 1/kE
    # Loss compliance
    D2_d = params["α_d"] * params["δD_d"] / (np.exp(params["α_d"] * s_d) + np.exp(-params["α_d"] * s_d)) #+ 1/(eta*omega)

    # SECOND MECHANISM - GBS
    # Grain boundary-based relaxation mechanism
    s_gb = np.log(tau_gb * omega)
    # Storage compliance (in Pa-1)
    D1_gb =  params["δD_gb"] * (1 - (2/np.pi) * np.arctan(np.exp(params["α_gb"] * s_gb))) #+ 1/kE
    # Loss compliance (in Pa-1)
    D2_gb = params["α_gb"] * params["δD_gb"] / (np.exp(params["α_gb"] * s_gb) + np.exp(-params["α_gb"] * s_gb)) #+ 1/(eta*omega)
    D1_t = D1_d + D1_gb + 1/kE
    D2_t = D2_d + D2_gb - 1/(eta*omega)

    return D1_t, D2_t, D1_d, D1_gb, D2_d, D2_gb

# ---------------------
# ICE Ih MELTING CURVE
# ---------------------
def Tm_ice_linear(P_MPa):
    """Linear approximation of the melting curve (Feistel & Wagner, 2005)."""
    clapeyron_slope = -0.074  # K / MPa
    return 273.15 + clapeyron_slope * P_MPa

# ------------------------------------
# CAMMARANO'S MODEL (adapted for ice)
# ------------------------------------
def Qinv_cammarano(T, Ba, gamma_c, ga, R, P_MPa, omega):
    """Quality factor according to Cammarano (2006)."""
    Tm_K = Tm_ice_linear(P_MPa)
    Hp = ga * Tm_K
    Qs = Ba * np.exp(gamma_c * Hp / T)   
    Qinv = 1.0 / (Qs * omega**gamma_c)
    return Qinv

# Cole's (1995) model
def mu_cole(omega, mu, eta, alpha_g, tau_d, deltaD):
    i = 1j
    s = np.log(tau_d * omega)
    term_real = 1.0/mu + deltaD * (1.0 - 2.0/pi * np.arctan(np.exp(alpha_g * s)))
    term_imag = deltaD * (alpha_g / (np.exp(alpha_g * s) + np.exp(-alpha_g * s))) + 1.0/(omega * eta)
    D = term_real - i * term_imag
    return 1.0 / D

# Cole's model
def mu_cole_T(omega, mu, eta, alpha_g, deltaD, T):
    i = 1j
    B = 1.205e-9 * np.exp(0.54 / (8.61e-5 * T))   # partial second units
    tau_d = B / 0.07
    s = np.log(tau_d * omega)
    term_real = 1.0/mu + deltaD * (1.0 - 2.0/pi * np.arctan(np.exp(alpha_g * s)))
    term_imag = deltaD * (alpha_g / (np.exp(alpha_g * s) + np.exp(-alpha_g * s))) + 1.0/(omega * eta)
    D = term_real - i * term_imag
    return 1.0 / D

def get_mu(T):
    """
    Linear interpolation of the shear modulus µ(T)
    for a given temperature T (in K).
    If T > max(T_mu), µ is fixed to 3.5 GPa.
    """
    mu_interp = np.interp(T, T_mu, mu_vals, left=mu_vals[0], right=mu_vals[-1])

    mu_interp = np.where(mu_interp <= 0, 3.5e9, mu_interp)

    return mu_interp

# -------------------------------
# PROTON REORIENTATION MECHANISM
# -------------------------------

def tP_of_T(T, c0):
    """Equation (58): relaxation time for proton reorientation (Vassouille 1974, Tatibouret 1981)"""
    return t0 * np.exp(Em/(R*T)) / (c0 + np.exp(-Ef/(R*T)))

def tan_delta_debye(omega, tP, T):
    """Debye-like peak using empirical tanδ_max from eq. (59): tanδ_max = Ep/(R T)"""
    tan_delta_max = Ep / (R * T)
    x = omega * tP
    return 2.0 * tan_delta_max * (x / (1.0 + x**2))

def tP_of_T_HP(T, c0):
    """Equation (58): relaxation time for proton reorientation (Vassouille 1974, Tatibouret 1981)"""
    return t0 * np.exp((Em+P*V)/(R*T)) / (c0 + np.exp(-Ef/(R*T)))

def tan_delta_debye_HP(omega, tP, T):
    """Debye-like peak using empirical tanδ_max from eq. (59): tanδ_max = Ep/(R T)"""
    tan_delta_max = (Ep+P*V) / (R * T)
    x = omega * tP
    return 2.0 * tan_delta_max * (x / (1.0 + x**2))

def Q_inv_PR(T, f, c0, tau0, Em_val):
    """
    Proton reorientation attenuation Q^-1_PR(T,f) using Debye peak.
    Same equations you used earlier.
    """
    R = 8.314
    omega = 2 * np.pi * f
    
    # Relaxation time tP(T)
    tP = tau0 * np.exp(Em_val/(R*T)) / (c0 + np.exp(-Ef/(R*T)))

    # tanδ_max empirical from Ep/(RT)
    tan_delta_max = Ep / (R * T)
    
    x = omega * tP
    tan_delta = 2.0 * tan_delta_max * (x / (1.0 + x**2))
    
    return tan_delta/1e3
    
def Q_inv_PR_hp(T, f, c0, tau0, Em):
    """
    Proton reorientation attenuation Q^-1_PR(T,f) using Debye peak.
    Same equations you used earlier.
    """
    R = 8.314
    omega = 2 * np.pi * f
    
    # Relaxation time tP(T)
    tP = tau0 * np.exp((Em+P*6e-6)/(R*T)) / (c0 + np.exp(-(Ef+P*V)/(R*T)))

    # tanδ_max empirical from Ep/(RT)
    tan_delta_max = (Ep+P*6e-6) / (R * T)
    
    x = omega * tP
    tan_delta = 2.0 * tan_delta_max * (x / (1.0 + x**2))
    
    return tan_delta/1e3
