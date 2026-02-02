# materials.py

import numpy as np

conv = 96485.33212

# According to parameters in Cole (1995) - sea water ice with 2-3% of NaCl
def saline_ice():
    return dict(
        eta = 1e22,                   # in Pa·s, shear viscosity
        V   = 10e-6,                  # in m3/mol, activation volume
        # -- Dislocation parameters --
        B0    = 1.205e-9,             # in Pa·s, dislocation drag pre-exponential
        Q_d   = 0.54*conv,            # in eV, activation energy for dislocation (NB. valid from -10ºC to -50ºC, see Fig. A1)
        K     = 0.07,                 # in Pa, restoring stress
        dD_d  = 1.4e-9,               # in Pa⁻¹, max ampli
        a_d   = 0.54,                 # Debye-type peak width at LF
        D_u_d = 1.2e-10,
        b     = 4.52e-10,             # in m, Burger's vector. NB. depends on the crystalline structure
        # -- Grain boundary parameters --
        tau_gb_0 = 8e-28,             # in s, relaxation time
        Q_gb     = 1.32*conv,         # in eV, activation energy for grain-boundary. Provient orginellement de l'article de Tatibouret (1987) 
                                      # sur de la glace pure, car pas de valeur fiable pour la glace saline selon Cole (1995)
        dD_gb    = 3e-11,             # in Pa⁻¹, max amplil
        a_gb     = 0.6,               # Debye-type peak width at HF
        D_u_gb   = 1.5e-10,
        # -- PR mechanism parameters --
        Em       = 35e3,              # J/mol, Bjerrum detects'energies of migration
        Ef       = 25e3,              # J/mol, Bjerrum detects'energies of formation
        Ep       = 25e3,              # in J/mol (for tanδ_max empirical)
        tau_pr_0 = 6e-16,             # in s [Tatibouet et al. 1981]
        c0       = 1e-7               # reference concentration of extrinsic defects
    )
