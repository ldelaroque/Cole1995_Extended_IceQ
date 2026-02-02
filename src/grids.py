# grids.py

# -- Import librairies --
import numpy as np


def temperature_grid(Tmin=94.0, Tmax=266.0, nT=300, spacing="linear"):
    """
    Generate a temperature grid for ice attenuation calculations.

    Parameters
    ----------
    Tmin : float
        Minimum temperature [K]
    Tmax : float
        Maximum temperature [K]
    nT : int
        Number of temperature points
    spacing : str
        'linear' or 'log'

    Returns
    -------
    T : ndarray
        Temperature array [K]
    """
    if spacing == "linear":
        return np.linspace(Tmin, Tmax, nT)

    elif spacing == "log":
        return np.logspace(np.log10(Tmin), np.log10(Tmax), nT)

    else:
        raise ValueError("spacing must be 'linear' or 'log'")


def frequency_grid(fmin=1e-3, fmax=100.0, nf=200, spacing="log"):
    """
    Generate a frequency grid for attenuation calculations.

    Parameters
    ----------
    fmin : float
        Minimum frequency [Hz]
    fmax : float
        Maximum frequency [Hz]
    nf : int
        Number of frequency points
    spacing : str
        'log' or 'linear'

    Returns
    -------
    f : ndarray
        Frequency array [Hz]
    """
    if spacing == "log":
        return np.logspace(np.log10(fmin), np.log10(fmax), nf)

    elif spacing == "linear":
        return np.linspace(fmin, fmax, nf)

    else:
        raise ValueError("spacing must be 'log' or 'linear'")
