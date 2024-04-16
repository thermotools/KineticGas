---
layout: default
version: 
title: Methods in the HardSphere class
permalink: /vcurrent/HardSphere_methods.html
---

<!--- 
Generated at: 2024-04-16T15:51:32.962670
This is an auto-generated file, generated using the script at KineticGas/pyUtils/markdown_from_docstrings.py
The file is created by parsing the docstrings of the methods in the 
HardSphere class. For instructions on how to use the parser routines, see the
file KineticGas/pyUtils/markdown_from_docstrings.py--->

The `HardSphere` class, found in `pykingas/HardSphere.py`, inherrits from the py_KineticGas class, and  is the interface to the 
RET-HS Model. This class implements utility methods to access mixing parameters etc.

## Table of contents
  * [Constructor](#constructor)
    * [\_\_init\_\_](#__init__self-comps-mole_weightsnone-sigmanone-n4-is_idealgasfalse-parameter_refdefault)

## Constructor

Methods to initialise RET-HS model.

### Table of contents
  * [Constructor](#constructor)
    * [\_\_init\_\_](#__init__self-comps-mole_weightsnone-sigmanone-n4-is_idealgasfalse-parameter_refdefault)


### `__init__(self, comps, mole_weights=None, sigma=None, N=4, is_idealgas=False, parameter_ref='default')`
If parameters are explicitly supplied through optional arguments, these will be used instead of those in the database.
To supply specific parameters for only some components, give `None` for the components that should use the database
value
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **comps (str) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Comma-separated list of components

&nbsp;&nbsp;&nbsp;&nbsp; **mole_weights (1D array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar weights [g/mol]

&nbsp;&nbsp;&nbsp;&nbsp; **sigma (1D array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  hard-sphere diameters [m]

&nbsp;&nbsp;&nbsp;&nbsp; **parameter_ref (str) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Id for parameter set to use 

