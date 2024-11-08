---
layout: default
version: 
title: Methods in the ModTangToennies class
permalink: /vcurrent/modtangtoennies_methods.html
---

<!--- 
Generated at: 2024-10-16T17:34:53.825965
This is an auto-generated file, generated using the script at KineticGas/pyUtils/markdown_from_docstrings.py
The file is created by parsing the docstrings of the methods in the 
ModTangToennies class. For instructions on how to use the parser routines, see the
file KineticGas/pyUtils/markdown_from_docstrings.py--->

The `ModTangToennies` class, found in `pykingas/multiparam.py`, inherrits from the py_KineticGas class, and  is the interface to the 
Modified Tang-Toennies Model. This class implements utility methods to access mixing parameters etc.

## Table of contents
  * [Constructor](#constructor)
    * [\_\_init\_\_](#__init__self-comps-parameter_refdefault)
  * [Utility methods](#utility-methods)
    * [potential](#potentialself-r)
    * [potential_r](#potential_rself-r)
    * [potential_rr](#potential_rrself-r)
    * [second_virial](#second_virialself-t)
    * [vdw_alpha](#vdw_alphaself)

## Constructor

Methods to initialise Modified Tang-Toennies model.

### Table of contents
  * [Constructor](#constructor)
    * [\_\_init\_\_](#__init__self-comps-parameter_refdefault)


### `__init__(self, comps, parameter_ref='default')`
Initialize modified Tang-Toennies potential
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **comps (str) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Single component identifier

&nbsp;&nbsp;&nbsp;&nbsp; **parameter_ref (str, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Identifier for parameter set to use 

## Utility methods

Set- and get methods for interaction parameters, mixing parameters ...

### Table of contents
  * [Utility methods](#utility-methods)
    * [potential](#potentialself-r)
    * [potential_r](#potential_rself-r)
    * [potential_rr](#potential_rrself-r)
    * [second_virial](#second_virialself-t)
    * [vdw_alpha](#vdw_alphaself)


### `potential(self, r)`
Evaluate potential
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **r (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Distance (m) 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **float :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Potential (J) 

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

### `second_virial(self, T)`
Compute second virial coefficient
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature (K) 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **float :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Second virial coefficient 

### `vdw_alpha(self)`
Get the dimensionless Van der Waals alpha-parameter
 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **float :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  alpha (-) 

