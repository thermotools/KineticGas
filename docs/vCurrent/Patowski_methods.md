---
layout: default
version: 
title: Methods in the Patowski class
permalink: /vcurrent/patowski_methods.html
---

<!--- 
Generated at: 2025-11-14T17:17:20.220042
This is an auto-generated file, generated using the script at KineticGas/pyUtils/markdown_from_docstrings.py
The file is created by parsing the docstrings of the methods in the 
Patowski class. For instructions on how to use the parser routines, see the
file KineticGas/pyUtils/markdown_from_docstrings.py--->

The `Patowski` class, found in `pykingas/multiparam.py`, inherrits from the py_KineticGas class, and  is the interface to the 
Patowski Model. This class implements utility methods to access mixing parameters etc.

## Table of contents
  * [Constructor](#constructor)
    * [\_\_init\_\_](#__init__self-comps-quantum_activetrue)
  * [Utility methods](#utility-methods)
    * [JKWB_phase_shift](#jkwb_phase_shiftself-i-j-l-e)
    * [JKWB_upper_E_limit](#jkwb_upper_e_limitself-i0-jnone)
    * [cross_section](#cross_sectionself-i-j-l-e-reducedfalse)
    * [de_broglie_wavelength](#de_broglie_wavelengthself-i-t)
    * [get_de_boer](#get_de_boerself-inone-jnone)
    * [get_quantum_active](#get_quantum_activeself)
    * [get_r_min](#get_r_minself-i-j)
    * [get_reducing_units](#get_reducing_unitsself-i0-jnone)
    * [omega](#omegaself-i-j-n-s-t)
    * [phase_shift](#phase_shiftself-i-j-l-e)
    * [potential](#potentialself-r)
    * [potential_dn](#potential_dnself-r-n)
    * [potential_r](#potential_rself-r)
    * [potential_rr](#potential_rrself-r)
    * [quantum_omega](#quantum_omegaself-i-j-n-s-t)
    * [set_de_boer_mass](#set_de_boer_massself-i-de_boer)
    * [set_quantum_active](#set_quantum_activeself-active)
    * [vdw_alpha](#vdw_alphaself)
    * [wave_function](#wave_functionself-i-j-l-e-r_end-dr01)

## Constructor

Methods to initialise Patowski model.

### Table of contents
  * [Constructor](#constructor)
    * [\_\_init\_\_](#__init__self-comps-quantum_activetrue)


### `__init__(self, comps, quantum_active=True)`
Initialize Patowski potential (used for H2)
 

## Utility methods

Set- and get methods for interaction parameters, mixing parameters ...

### Table of contents
  * [Utility methods](#utility-methods)
    * [JKWB_phase_shift](#jkwb_phase_shiftself-i-j-l-e)
    * [JKWB_upper_E_limit](#jkwb_upper_e_limitself-i0-jnone)
    * [cross_section](#cross_sectionself-i-j-l-e-reducedfalse)
    * [de_broglie_wavelength](#de_broglie_wavelengthself-i-t)
    * [get_de_boer](#get_de_boerself-inone-jnone)
    * [get_quantum_active](#get_quantum_activeself)
    * [get_r_min](#get_r_minself-i-j)
    * [get_reducing_units](#get_reducing_unitsself-i0-jnone)
    * [omega](#omegaself-i-j-n-s-t)
    * [phase_shift](#phase_shiftself-i-j-l-e)
    * [potential](#potentialself-r)
    * [potential_dn](#potential_dnself-r-n)
    * [potential_r](#potential_rself-r)
    * [potential_rr](#potential_rrself-r)
    * [quantum_omega](#quantum_omegaself-i-j-n-s-t)
    * [set_de_boer_mass](#set_de_boer_massself-i-de_boer)
    * [set_quantum_active](#set_quantum_activeself-active)
    * [vdw_alpha](#vdw_alphaself)
    * [wave_function](#wave_functionself-i-j-l-e-r_end-dr01)


### `JKWB_phase_shift(self, i, j, l, E)`
Compute the phase shift for a collision with angular momentum quantum number `l` and energy `E`, using the JKWB approximation
Args:
i, j (int): Species indices
l (int): Angular momentum quantum number
E (float): Total energy (J)
Returns:
float: The relative phase shift $(- \pi / 2, \pi / 2)$
 

### `JKWB_upper_E_limit(self, i=0, j=None)`
Get the upper energy limit for when the JKWB approximation is automatically applied.
 

### `cross_section(self, i, j, l, E, reduced=False)`
Calculate the collision cross section. If `reduced=True`, return the cross section divided by the hard-sphere cross section.
 

### `de_broglie_wavelength(self, i, T)`
Get the de Broglie wavelength of species `i` at temperature `T`.
 

### `get_de_boer(self, i=None, j=None)`
Get the de Boer parameter
 

### `get_quantum_active(self)`
Get the current quantum_active state.
 

### `get_r_min(self, i, j)`
Compute the position of the potential minimum.
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **i, j (int):** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Species indices

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **float :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  r_min (m) 

### `get_reducing_units(self, i=0, j=None)`
See `py_KineticGas`.
 

### `omega(self, i, j, n, s, T)`
Calculate the collision integral $\Omega^{(n, s)}$ as defined in The Limits of the Feynman-Hibbs corrections ... paper (see cite page).

This method uses quantum mechanical or classical calculation based on whether `self.get_quantum_active()` is `True`
 

### `phase_shift(self, i, j, l, E)`
Compute the phase shift for a collision with angular momentum quantum number `l` and energy `E`
Args:
i, j (int): Species indices
l (int): Angular momentum quantum number
E (float): Total energy (J)
Returns:
float: The relative phase shift $(- \pi / 2, \pi / 2)$
 

### `potential(self, r)`
Evaluate potential
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **r (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Distance (m) 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **float :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Potential (J) 

### `potential_dn(self, r, n)`
Calculate the `n`'th derivative of the potential wrt. distance.
Args:
r (float) : Distance (m)

Returns:
float : Potential n'th derivative wrt. distance (J / m^n)
 

### `potential_r(self, r)`
Evaluate potential derivative wrt. distance
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **r (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Distance (m) 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **float :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Potential derivative wrt. distance (N) 

### `potential_rr(self, r)`
Evaluate potential second derivative wrt. distance
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **r (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Distance (m) 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **float :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Potential second derivative wrt. distance (N / m) 

### `quantum_omega(self, i, j, n, s, T)`
Calculate the quantal collision integral $\Omega^{(n, s)}$ as defined in The Limits of the Feynman-Hibbs corrections ... paper (see cite page).
 

### `set_de_boer_mass(self, i, de_boer)`
Set the particle mass to get the specified de Boer parameter
 

### `set_quantum_active(self, active)`
Activate/deactivate quantum mechanical calculation of things.
 

### `vdw_alpha(self)`
Get the dimensionless Van der Waals alpha-parameter
 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **float :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  alpha (-) 

### `wave_function(self, i, j, l, E, r_end, dr=0.1)`
Solve the Schrödinger equation for the two-particle wave function at energy E, out to the distance `r_end`
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **i, j (int):** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Species indices

&nbsp;&nbsp;&nbsp;&nbsp; **l (int):** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Angular momentum quantum number

&nbsp;&nbsp;&nbsp;&nbsp; **E (float):** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Total energy (J)

&nbsp;&nbsp;&nbsp;&nbsp; **r_end (float):** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Maximum particle separation (m)

&nbsp;&nbsp;&nbsp;&nbsp; **dr (float):** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Step size

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **list[float] :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  The two-particle non-normalized wave function out to `r_end`. 

