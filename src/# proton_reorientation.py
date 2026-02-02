# proton_reorientation.py

# -- Import librairies --
import numpy as np
from .attenuation_models import *
from .constants import R


def tP_of_T(P, T, params):
    """Equation (58): relaxation time for proton reorientation (Vassouille 1974, Tatibouret 1981)"""
    return params["tau0"] * np.exp((params["Em"]+P*params["V"])/(R*T)) / (params["c0"] + np.exp(-(params["Ef"]+P*params["V"])/(R*T)))


def tan_delta_debye(T, P, f, params, tP):
    """Debye-like peak using empirical tanδ_max from eq. (59): tanδ_max = Ep/(R T)"""
    omega = 2*np.pi*f

    tan_delta_max = (params["Ep"]+P*params["V"]) / (R*T)
    x = omega * tP

    return 2.0 * tan_delta_max * (x / (1.0 + x**2))


def compute_Qmu_PR(P, T, f, params):
    """
    Proton reorientation attenuation Q^-1_PR(T,f) using Debye peak.
    Same equations you used earlier.
    """
    omega = 2*np.pi*f
    
    # relaxation time tP(T)
    tP = params["tau0"] * np.exp((params["Em"]+P*params["V"])/(R*T)) / (params["c0"] + np.exp(-(params["Ef"]+P*params["V"])/(R*T)))

    # tanδ_max empirical from Ep/(RT)
    tan_delta_max = (params["Ep"]+P*params["V"]) / (R * T)
    x = omega * tP
    tan_delta = 2.0 * tan_delta_max * (x / (1.0 + x**2))
    
    return tan_delta/1e3
