# rheology.py

# -- Import librairies --
import numpy as np

mu_data = np.loadtxt("./data/TitanProfile_MgSO4_10WtPct_Zb88km_pp.txt", skiprows=1)
T_mu    = mu_data[:, 1]  
mu_vals = mu_data[:, 8] * 1e9

def get_mu(T):
    """
    Linear interpolation of the shear modulus µ(T)
    for a given temperature T (in K).
    If T > max(T_mu), µ is fixed to 3.5 GPa.
    """
    mu_interp = np.interp(T, T_mu, mu_vals, left=mu_vals[0], right=mu_vals[-1])

    mu_interp = np.where(mu_interp <= 0, 3.5e9, mu_interp)

    return mu_interp
