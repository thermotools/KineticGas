---
layout: default
version: 2.0.0
title: Known Bugs
permalink: /v2.0.0/bugs.html
---

# Known bugs in v2.0.0

This page contains an overview of known bugs in KineticGas 2.0.0, along with information regarding patch status and 
possible workarounds.

## ImportError cyclic import
In some cases, for some installations, importing `pykingas` causes an `ImportError` to be raised. This appears to 
occur most commonly when `pykingas` has been installed on the system level with Anaconda. 

### Workaround
* The bug has never been observed when installing to a properly set up virtual environment.

If you are running a system level interpreter, uninstall `pykingas`, set up a virtual environment, and reinstall `pykingas`
to that virtual environment.

### Status
* Newer versions have a more robust import system

## Ideal gas `thermal_diffusion_ratio`
Models initialised with `is_idealgas=True` will produce incorrect thermal diffusion ratios ($k_{T,i}$). This bug 
does not effect other properties.

### Status
* Fixed in later versions

## Multicomponent (> 2) thermal properties
For multicomponent (more than two species) mixtures, properties relying on the the coefficients $\ell_i^{(0)}$ and $\ell_i^{(1)}$ 
(see: Eq. (6), (13), (15) and (24) in [Revised Enskog Theory for Mie Fluids](https://pubs.aip.org/aip/jcp/article/158/22/224101/2895227/Revised-Enskog-theory-for-Mie-fluids-Prediction-of))
are incorrect. This bug effects the methods
* `thermal_conductivity` / `thermal_coductivity_tp`
* `thermal_diffusion_coeff` / `thermal_diffusion_coeff_tp`
* `thermal_diffusion_factor` / `thermal_diffusion_factor_tp`
* `thermal_diffusion_ratio` / `thermal_diffusion_ratio_tp`
* `soret_coefficient` / `soret_coefficient_tp`

### Status
* Fixed in later versions