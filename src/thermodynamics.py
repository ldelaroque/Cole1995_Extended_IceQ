# thermodynamics.py

# -- Import librairies --
import numpy as np

def Tm_p_Hirschmann2000(P_MPa):
    """
    Compute melting temperature (K) from pressure (MPa)
    using Hirschmann (2000) parameterization.
    """
    P_GPa  = 1e-3 * P_MPa
    coeffs = [-5.140, 132.899, 1120.661]  # coefficients for polynomial in GPa
    # np.polyval takes highest degree first (like MatLab)
    Tm_oC  = np.polyval(coeffs, P_GPa)
    return 273.15 + Tm_oC  # convert to Kelvin


# -- Ice Ih melting curve --
def Tm_ice_linear(P_MPa):
    """
    Local linear approximation of the melting curve (valid for low pressures typical of glaciers)
    from Feistel & Wagner, 2005 ("A New Equation of State for H2O Ice Ih").
    """
    # Clausius-Clapeyron slope
    clapeyron_slope = -0.074  # K / MPa (approx.), see eqn 21 of the paper
    return 273.15 + clapeyron_slope * P_MPa
