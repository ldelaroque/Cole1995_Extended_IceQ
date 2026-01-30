# Cole-SAS-Maxwell Attenuation Model

## Overview

This reprository contains a reference implementation of the **Cole-SAS-Maxwell anelastic attenuation mode**, developed to computed frequency-, temperature-, and pressure-dependent shear attenuation $Q_\mu^{-1}(\omega; T,P)$ in **water ice (ice Ih)**.

The model integrated the dominant anelastic relaxation mechanisms relevant to planetary cryospheres, including:
- **Proton reorientation**,
- **Dislocation-based anelasticity**,
- **Viscoelastic background response (Maxwell)**.
It is designed for applications spanning planetary seismology / wave propagation in icy shells and investigating attenuation-controlled processes at relevant conditions.

## Context

Intrinsic attenution models decribe the frequency-dependent viscoelastic response of materials and are widely used across disciplines such as seismology, geodynamics, and planetary science.
This implementation targets conditions relevant to **icy satellites** and **planetary cryospheres**, bridging seismic (Hz) to tidal (mHz) frequency regimes.

This model is consistent with laboratory- and field-derived attenuation measurements and reproduces published quality factors over a wide range of thermodynamic conditions.
