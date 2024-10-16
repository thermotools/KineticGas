---
layout: default
version: 
title: Installing KineticGas
permalink: /vcurrent/units.html
---

All IO in the KineticGas package is in SI units, on a molar basis.

To support working in reduced units, `pykingas` uses the `Units` struct, found in `pykingas/units.py`. Models may implement the method `get_reducing_units`, which returns an instance of this struct, which holds conversion factors for quickly computing non-dimensional quantities from the SI-quantities returned by `pykingas` models.

The `Units` struct has the attributes

* `T` : Temperature (K)
* `E` : Energy (J)
* `L` : Length (m)
* `m` : Mass (kg)
* `V` : Volume (m$^3$)
* `t` : Time (s)
* `F` : Force (N)
* `speed` : Speed (m / s)
* `rho` : Density (mol / m$^3$)
* `D` : Diffusion coefficient (m$^2$ / s)
* `p` : Pressure (Pa)
* `visc` : Viscosity (Pa s)
* `tcond` : Thermal conductivity (W / K m)
  
Different models may use different scaling parameters to deduce the reducing factors. The `MieKinGas` model uses the Mie potential parameters, such that e.g. the dimensionless viscosity may be computed as

```python
from pykingas.MieKinGas import MieKinGas
mie = MieKinGas('LJF') # RET-Mie for Lennard-Jones fluid
unt = mie.get_reducing_units()
T, p = 300, 1e5 # Temperature and pressure
visc = mie.viscosity_tp(T, p, [1], N=2) # Viscosity in SI-units (Pa s)

visc_red = visc / unt.visc # Viscosity in Lennard-Jones units

mie2 = MieKinGas('CO2,H2') # RET-Mie for CO2 / H2 mixture
unt1 = mie2.get_reducing_units(0) # Use potential parameters of first component for reducing units
unt2 = mie2.get_reducing_units(1) # Use potential parameters of second component
```

An overview of what reducing units different models use is found below
* `HardSphere` : No reducing units (has no energy unit)
* `MieKinGas` : Mie potential parameters
* `ModTangToennies` : Potential root (`sigma`) and well depth (`epsilon`)