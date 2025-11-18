---
layout: default
version: 
title: Methods in the Quantum class
permalink: /vcurrent/quantum_methods.html
---

<!--- 
Generated at: 2025-11-14T17:17:20.220931
This is an auto-generated file, generated using the script at KineticGas/pyUtils/markdown_from_docstrings.py
The file is created by parsing the docstrings of the methods in the 
Quantum class. For instructions on how to use the parser routines, see the
file KineticGas/pyUtils/markdown_from_docstrings.py--->

The Quantum class is an abstract base class that contains all the interfaces to quantum mechanical calculations for Spherical potentials.Any spherical potential model can be made to inherrit from the Quantum class instead of inherriting directly from `py_KineticGas` in order to gain access toquantum mechanical calculation of collision integrals and other stuff. **NOTE**: That the `Quantum` class has the attribute `quantum_active` which is used toturn quantal calculations on/off. So even if a class inherits from `Quantum`, the classical calculations can be accessed by turning off the quantal calculations.## Table of contents
  * [Constructor](#constructor)
    * [\_\_init\_\_](#__init__self-comps-is_idealgastrue)
  * [Utility methods](#utility-methods)
    * [JKWB_phase_shift](#jkwb_phase_shiftself-i-j-l-e)
    * [JKWB_upper_E_limit](#jkwb_upper_e_limitself-i0-jnone)
    * [cross_section](#cross_sectionself-i-j-l-e-reducedfalse)
    * [de_broglie_wavelength](#de_broglie_wavelengthself-i-t)
    * [get_de_boer](#get_de_boerself-inone-jnone)
    * [get_quantum_active](#get_quantum_activeself)
    * [get_reducing_units](#get_reducing_unitsself-i0-jnone)
    * [omega](#omegaself-i-j-n-s-t)
    * [phase_shift](#phase_shiftself-i-j-l-e)
    * [potential](#potentialself-i-j-r)
    * [potential_r](#potential_rself-i-j-r)
    * [potential_rr](#potential_rrself-i-j-r)
    * [quantum_omega](#quantum_omegaself-i-j-n-s-t)
    * [set_de_boer_mass](#set_de_boer_massself-i-de_boer)
    * [set_quantum_active](#set_quantum_activeself-active)
    * [wave_function](#wave_functionself-i-j-l-e-r_end-dr01)

## Constructor

Methods to initialise Quantum model.

### Table of contents
  * [Constructor](#constructor)
    * [\_\_init\_\_](#__init__self-comps-is_idealgastrue)


### `__init__(self, comps, is_idealgas=True)`
Interface to quantum mechanical stuff.
 

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
    * [get_reducing_units](#get_reducing_unitsself-i0-jnone)
    * [omega](#omegaself-i-j-n-s-t)
    * [phase_shift](#phase_shiftself-i-j-l-e)
    * [potential](#potentialself-i-j-r)
    * [potential_r](#potential_rself-i-j-r)
    * [potential_rr](#potential_rrself-i-j-r)
    * [quantum_omega](#quantum_omegaself-i-j-n-s-t)
    * [set_de_boer_mass](#set_de_boer_massself-i-de_boer)
    * [set_quantum_active](#set_quantum_activeself-active)
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
 

### `potential(self, i, j, r)`
Potential
 

### `potential_r(self, i, j, r)`
Potential derivative
 

### `potential_rr(self, i, j, r)`
Potential second derivative
 

### `quantum_omega(self, i, j, n, s, T)`
Calculate the quantal collision integral $\Omega^{(n, s)}$ as defined in The Limits of the Feynman-Hibbs corrections ... paper (see cite page).
 

### `set_de_boer_mass(self, i, de_boer)`
Set the particle mass to get the specified de Boer parameter
 

### `set_quantum_active(self, active)`
Activate/deactivate quantum mechanical calculation of things.
 

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

