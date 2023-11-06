---
layout: default
version: 
title: Methods in the MieKinGas class
permalink: /vcurrent/MieKinGas_methods.html
---

<!--- 
Generated at: 2023-11-06T12:02:00.296867
This is an auto-generated file, generated using the script at KineticGas/pyUtils/markdown_from_docstrings.py
The file is created by parsing the docstrings of the methods in the 
MieKinGas class. For instructions on how to use the parser routines, see the
file KineticGas/pyUtils/markdown_from_docstrings.py--->

The `MieKinGas` class, found in `pykingas/MieKinGas.py`, inherrits from the MieType class, and  is the interface to the 
RET-Mie Model. This class implements utility methods to access mixing parameters etc.

## Table of contents
  * [Constructor](#constructor)
    * [\_\_init\_\_](#__init__self-comps-mole_weightsnone-sigmanone-eps_div_knone-lanone-lrnone-lij0-kij0-n4-is_idealgasfalse-use_eosnone-parameter_refdefault)

## Constructor

Methods to initialise RET-Mie model.

### Table of contents
  * [Constructor](#constructor)
    * [\_\_init\_\_](#__init__self-comps-mole_weightsnone-sigmanone-eps_div_knone-lanone-lrnone-lij0-kij0-n4-is_idealgasfalse-use_eosnone-parameter_refdefault)


### `__init__(self, comps, mole_weights=None, sigma=None, eps_div_k=None, la=None, lr=None, lij=0, kij=0, N=4, is_idealgas=False, use_eos=None, parameter_ref='default')`
If parameters are explicitly supplied through optional arguments, these will be used instead of those in the database.
To supply specific parameters for only some components, give `None` for the components that should use the database
value
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **comps (str) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Comma-separated list of components

&nbsp;&nbsp;&nbsp;&nbsp; **mole_weights (1D array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar weights [g/mol]

&nbsp;&nbsp;&nbsp;&nbsp; **sigma (1D array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  sigma-parameters [m]

&nbsp;&nbsp;&nbsp;&nbsp; **eps_div_k (1D array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  epsilon parameter / Boltzmann constant [-]

&nbsp;&nbsp;&nbsp;&nbsp; **la, lr (1D array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  attractive and repulsive exponent of the pure components [-]

&nbsp;&nbsp;&nbsp;&nbsp; **lij (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Mixing parameter for sigma (lij > 0 => smaller sigma_12, lij < 0 => larger sigma_12)

&nbsp;&nbsp;&nbsp;&nbsp; **kij (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Mixing parameter for epsilon (kij > 0 => favours mixing, kij < 0 => favours separation)

&nbsp;&nbsp;&nbsp;&nbsp; **use_eos :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  (thermopack eos object, optional) EoS to use (initialized), defaults to `saftvrmie` 

