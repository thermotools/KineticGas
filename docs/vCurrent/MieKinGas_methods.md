---
layout: default
version: 
title: Methods in the MieKinGas class
permalink: /vcurrent/MieKinGas_methods.html
---

<!--- 
Generated at: 2024-04-18T17:40:59.602391
This is an auto-generated file, generated using the script at KineticGas/pyUtils/markdown_from_docstrings.py
The file is created by parsing the docstrings of the methods in the 
MieKinGas class. For instructions on how to use the parser routines, see the
file KineticGas/pyUtils/markdown_from_docstrings.py--->

The `MieKinGas` class, found in `pykingas/MieKinGas.py`, inherrits from the MieType class, and  is the interface to the 
RET-Mie Model. This class implements utility methods to access mixing parameters etc.

## Table of contents
  * [Constructor](#constructor)
    * [\_\_init\_\_](#__init__self-comps-mole_weightsnone-sigmanone-eps_div_knone-lanone-lrnone-lij0-kij0-n4-is_idealgasfalse-use_eosnone-parameter_refdefault-use_default_eos_paramfalse)
  * [Utility methods](#utility-methods)
    * [set_eps_div_k](#set_eps_div_kself-eps_div_k-update_eostrue)
    * [set_la](#set_laself-la-update_eostrue)
    * [set_lr](#set_lrself-lr-update_eostrue)
    * [set_sigma](#set_sigmaself-sigma-update_eostrue)
  * [Internal methods](#internal-methods)
    * [\_\_update_cpp_kingas_param\_\_](#__update_cpp_kingas_param__self)

## Constructor

Methods to initialise RET-Mie model.

### Table of contents
  * [Constructor](#constructor)
    * [\_\_init\_\_](#__init__self-comps-mole_weightsnone-sigmanone-eps_div_knone-lanone-lrnone-lij0-kij0-n4-is_idealgasfalse-use_eosnone-parameter_refdefault-use_default_eos_paramfalse)


### `__init__(self, comps, mole_weights=None, sigma=None, eps_div_k=None, la=None, lr=None, lij=0, kij=0, N=4, is_idealgas=False, use_eos=None, parameter_ref='default', use_default_eos_param=False)`
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

&nbsp;&nbsp;&nbsp;&nbsp; **use_eos (thermopack eos object, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  EoS to use (initialized), defaults to `saftvrmie`

&nbsp;&nbsp;&nbsp;&nbsp; **use_default_eos_param (bool) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  If `False` (default), ensure that the EoS and RET-model use the same parameters(if applicable). If `False`, do not forward specified parameters to the EoS.  

## Utility methods

Set- and get methods for interaction parameters, mixing parameters ...

### Table of contents
  * [Utility methods](#utility-methods)
    * [set_eps_div_k](#set_eps_div_kself-eps_div_k-update_eostrue)
    * [set_la](#set_laself-la-update_eostrue)
    * [set_lr](#set_lrself-lr-update_eostrue)
    * [set_sigma](#set_sigmaself-sigma-update_eostrue)


### `set_eps_div_k(self, eps_div_k, update_eos=True)`
See MieType
 

### `set_la(self, la, update_eos=True)`
See MieType
 

### `set_lr(self, lr, update_eos=True)`
See MieType
 

### `set_sigma(self, sigma, update_eos=True)`
See MieType
 

## Internal methods

Internal methods are not intended for use by end-users.

### Table of contents
  * [Internal methods](#internal-methods)
    * [\_\_update_cpp_kingas_param\_\_](#__update_cpp_kingas_param__self)


### `__update_cpp_kingas_param__(self)`
See MieType
 

