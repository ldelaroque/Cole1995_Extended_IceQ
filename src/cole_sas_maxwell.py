# cole_sas_maxwell.py

# -- Import librairies --
import numpy as np
import matplotlib.pyplot as plt
from .rheology import get_mu
from .constants import R
from .proton_reorientation import compute_Qmu_PR

# -----------------
# COLE-SAS-MAXWELL
# -----------------

def compute_qmu(P, T, f, params):
    """
    Calculate the compliances D1 and D2 for a given temperature (in °C) 
    and frequency f (Hz), according to the Cole (1995) model.
    """
    omega = 2*np.pi*f   # in rad/s, angular frequency

    # µ(T) - temperature dpdt
    kE = get_mu(T)

    # Temperature-dependant dislocation drag
    B_T_d = params["B0"] * np.exp((params["E_d"]+P*params["V"]) / (R * T))   # partial second units
    # Relaxation time for dislocation sliding
    tau_d = B_T_d / params["K"]

    # Relaxation time for GBS
    tau_gb = params["tau_gb_0"] * np.exp((params["E_gb"]+P*params["V"])/ (R * T))

    # FIRST MECHANISM - DISLOCATION
    # Dislocation-based relaxation mechanism
    s_d = np.log(tau_d * omega)
    # Storage compliance
    D1_d = params["dD_d"] * (1 - (2/np.pi) * np.arctan(np.exp(params["a_d"] * s_d))) 
    # Loss compliance
    D2_d = params["a_d"] * params["dD_d"] / (np.exp(params["a_d"] * s_d) + np.exp(-params["a_d"] * s_d)) 

    # SECOND MECHANISM - GBS
    # Grain boundary-based relaxation mechanism
    s_gb = np.log(tau_gb * omega)
    # Storage compliance (in Pa-1)
    D1_gb =  params["dD_gb"] * (1 - (2/np.pi) * np.arctan(np.exp(params["a_gb"] * s_gb))) 
    # Loss compliance (in Pa-1)
    D2_gb = params["a_gb"] * params["dD_gb"] / (np.exp(params["a_gb"] * s_gb) + np.exp(-params["a_gb"] * s_gb)) 

    # Adding the Maxwell correction
    D1_t = D1_d + D1_gb + 1/kE
    D2_t = D2_d + D2_gb - 1/(params["eta"]*omega)

    return D1_t, D2_t, D1_d, D1_gb, D2_d, D2_gb


def combined_mechanisms(P, T, f, params, mechanisms=("D", "GBS", "PR")):
    """
    Combine anelastic relaxation mechanisms into total shear attenuation.

    Parameters
    ----------
    P : float
        Pressure [Pa]
    T_C : float or ndarray
        Temperature [°C]
    f : float or ndarray
        Frequency [Hz]
    params : dict
        Rheological parameters
    mechanisms : tuple
        Active mechanisms ("D", "GBS", "PR")

    Returns
    -------
    Qmu : ndarray
        Shear quality factor
    D1_tot, D2_tot : ndarray
        Total storage and loss compliances
    """

    # Elastic compliance
    omega = 2 * np.pi * f
    mu_E = get_mu(T)

    # Elastic compliance
    D1_el = 1.0 / mu_E
    D1_tot = D1_el.copy()
    D2_tot = np.zeros_like(D1_el)

    # Cole-SAS-Maxwell
    if "D" in mechanisms or "GBS" in mechanisms:
        D1_t, D2_t, D1_d, D1_gb, D2_d, D2_gb = compute_qmu(P, T, f, params)

        if "D" in mechanisms:
            D1_tot += D1_d
            D2_tot += D2_d

        if "GBS" in mechanisms:
            D1_tot += D1_gb
            D2_tot += D2_gb

    # Proton reorientation (additive in Q-1)
    if "PR" in mechanisms:
        Qinv_pr = compute_Qmu_PR(P, T, f, params)
        D2_tot += Qinv_pr * D1_el

    Qmu = D1_tot / D2_tot
    return Qmu, D1_tot, D2_tot
