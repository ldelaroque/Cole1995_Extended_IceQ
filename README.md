# Cole-SAS-Maxwell Attenuation Model

## Overview

This reprository contains a reference implementation of the **Cole-SAS-Maxwell anelastic attenuation model**, developed to computed frequency-, temperature-, and pressure-dependent shear attenuation $Q_\mu^{-1}(\omega; T,P)$ in **water ice (ice Ih)**.

The model integrated the dominant anelastic relaxation mechanisms relevant to planetary cryospheres, including:
- **Proton reorientation** (H-bond),
- **Dislocation-based anelasticity**,
- **Viscoelastic background response (Maxwell)**.
It is designed for applications spanning planetary seismology / wave propagation in icy shells and investigating attenuation-controlled processes at relevant conditions.

## Context

Intrinsic attenution models decribe the frequency-dependent viscoelastic response of materials and are widely used across disciplines such as seismology, geodynamics, and planetary science.
This implementation targets conditions relevant to **icy satellites** and **planetary cryospheres**, bridging seismic (Hz) to tidal (mHz) frequency regimes.

This model is consistent with laboratory- and field-derived attenuation measurements and reproduces published quality factors over a wide range of thermodynamic conditions.

## Features

Compute the (seismic) shear quality factor $Q_\mu^{-1}(\omega; T,P)$ for ice Ih.

Validity:
- Frequency range: 0.001-100 Hz
- Temperature: 94-260 K
- Pressure range: 0.1 (near-surface)-100 MPa (ice-ocean)

Suitable for:
- synthetic seismogram modeling,
- depth-dependent attenuation profiles,
- tidal dissipation studies in icy shells

## Usage

Run the model located in the main directory with:

```ruby
python run_qmu_vs_T.py
```

## References
 If you use this code, please cite:
 - Delaroque et al. 2026 - original developement of the Cole-SAS-Maxwell model

Other usefule references:
 - [Cole et al. (1995)](https://doi.org/10.1080/01418619508239592) - background attenuation model
 - [Castillo et al. (2011)](https://doi.org/10.1029/2010JE003664) – proton orientation formulation included in the Cole-SAS-Maxwell model
 - [Cammarano et al. (2006)](https://doi.org/10.1029/2006JE002710); [Tobie et al. (2025)](https://doi.org/10.1007/s11214-025-01136-y) - comparative attenuation models

## Limitations
The model is currently calibrated for **pure water ice (ice Ih)**. Some parameters associated with high-pressure ice polymorphs (e.g.., III, V, VI) cannot be idrectly extrapolated to thermodynamic conditions comparable to those of ice Ih. 

Grain size, impurities, and partial melt effects are not explicitly included. 

Pressure dependence is parametrerized using activation volumes.

## Contact
Developed by **Lorraine Delaroque** \
*Institut de Physique du Globe de Paris* \
E-mail: delaroque .at. ipgp.fr
