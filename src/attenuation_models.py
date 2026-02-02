# attenuation_models.py

# -- Import librairies --
import numpy as np
from .thermodynamics import Tm_ice_linear


# ------------------------------------
# CAMMARANO'S MODEL (adapted for ice)
# ------------------------------------
def Qinv_Cammarano(T, f):
    """Intrinsic attenuation model adapted from Cammarano et al. (2006).

    Parameters
    ----------
    T : float or ndarray
        Temperature [K]
    f : float or ndarray
        Frequency [Hz]
    P_MPa : float
        Pressure [MPa]
    Ba : float
        Pre-exponential factor
    gamma : float
        Frequency exponent
    g : float
        Scaling coefficient for homologous temperature

    Returns
    -------
    Qinv : ndarray
        Inverse shear quality factor (= shear attenuation)
        """
    omega = 2*np.pi*f
    
    Ba      = 0.56
    gamma_c = 0.2
    ga      = 22
    P_MPa   = 2.0  

    Tm_K  = Tm_ice_linear(P_MPa)

    Hp = ga * Tm_K
    Qs = Ba * np.exp(gamma_c * Hp / T)   

    return 1.0 / (Qs * omega**gamma_c)

# ----------------------------------------------
# TOBIE'S MODEL (simple Cole - dislocation only)
# ----------------------------------------------

def mu_cole(omega, mu_E, eta, alpha_d, tau_d, deltaD):
    """
    Complex shear modulus following Cole (1995).

    Parameters
    ----------
    omega : float or ndarray
        Angular frequency [rad/s]
    mu_E : float or ndarray
        Elastic shear modulus [Pa]
    eta : float
        Maxwell viscosity [Pa s]
    alpha : float
        Cole distribution parameter
    tau : float or ndarray
        Relaxation time [s]
    deltaD : float
        Relaxation strength [Pa^-1]

    Returns
    -------
    mu_star : complex ndarray
        Complex shear modulus
    """
    s = np.log(tau_d * omega)

    D1 = 1.0/mu_E + deltaD * (1.0 - 2.0/np.pi * np.arctan(np.exp(alpha_d * s)))                  # Real part
    D2 = deltaD * (alpha_d / (np.exp(alpha_d * s) + np.exp(-alpha_d * s))) + 1.0/(omega * eta)   # Imaginary part
    
    return 1.0 / (D1 - 1j * D2)


def Qinv_Tobie(T, f):
    """
    Intrinsic attenuation model following Tobie et al. (Cole-type rheology).

    Parameters
    ----------
    T : float or ndarray
        Temperature [K]
    f : float or ndarray
        Frequency [Hz]
    (other params : already defined as in the article).

    Returns
    -------
    Qinv : ndarray
        Inverse shear quality factor (= shear attenuation)
    """
    omega = 2*np.pi*f
    B     = 1.205e-9 * np.exp(0.54 / (8.61e-5 * T))   # partial second units
    tau_d = B / 0.07

    # Parameters used from Tobie et al. 2025
    mu_star = mu_cole(
        omega   = omega,
        mu_E    = 3.3e9,
        eta     = 1e22,
        alpha_d = 0.53,
        tau_d   = tau_d,
        deltaD  = 1.4e-9,
    )

    return np.abs(np.imag(mu_star)) / np.real(mu_star)

