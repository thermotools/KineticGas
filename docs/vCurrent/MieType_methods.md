---
layout: default
version: 
title: Methods in the MieType class
permalink: /vcurrent/MieType_methods.html
---

<!--- 
Generated at: 2024-02-20T16:59:43.575853
This is an auto-generated file, generated using the script at KineticGas/pyUtils/markdown_from_docstrings.py
The file is created by parsing the docstrings of the methods in the 
MieType class. For instructions on how to use the parser routines, see the
file KineticGas/pyUtils/markdown_from_docstrings.py--->

The `MieType` class, found in `pykingas/MieType.py`, inherrits from the py_KineticGas class, and  is the interface to the 
Mie-Type Model. This class implements utility methods to access mixing parameters etc.

## Table of contents
  * [Constructor](#constructor)
    * [\_\_init\_\_](#__init__self-comps-potential-mole_weightsnone-sigmanone-eps_div_knone-lanone-lrnone-lij0-kij0-n4-is_idealgasfalse-parameter_refdefault)
  * [Utility methods](#utility-methods)
    * [get_BH_diameters](#get_bh_diametersself-t)
    * [get_epsilon_matrix](#get_epsilon_matrixself-eps_div_k-kij)
    * [get_lambda_matrix](#get_lambda_matrixself-lambdas-lij)
    * [get_sigma_matrix](#get_sigma_matrixself-sigma)
    * [potential](#potentialself-i-j-r)
    * [potential_r](#potential_rself-i-j-r)
    * [potential_rr](#potential_rrself-i-j-r)

## Constructor

Methods to initialise Mie-Type model.

### Table of contents
  * [Constructor](#constructor)
    * [\_\_init\_\_](#__init__self-comps-potential-mole_weightsnone-sigmanone-eps_div_knone-lanone-lrnone-lij0-kij0-n4-is_idealgasfalse-parameter_refdefault)


### `__init__(self, comps, potential, mole_weights=None, sigma=None, eps_div_k=None, la=None, lr=None, lij=0, kij=0, N=4, is_idealgas=False, parameter_ref='default')`
If optional parameters are supplied, these are used instead of the parameters found in the database. To supply specific parameters for only some components, give `None` for the components that should use the database
value
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **comps (str) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Comma-separated list of components

&nbsp;&nbsp;&nbsp;&nbsp; **mole_weights (optional, 1D array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar masses [g/mol]

&nbsp;&nbsp;&nbsp;&nbsp; **sigma (optional, 1D array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  hard-sphere diameters [m]

&nbsp;&nbsp;&nbsp;&nbsp; **eps_div_k (optional, 1D array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  epsilon parameter / Boltzmann constant [-]

&nbsp;&nbsp;&nbsp;&nbsp; **la, lr (optional, 1D array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  attractive and repulsive exponent of the pure components [-]

&nbsp;&nbsp;&nbsp;&nbsp; **lij (optional, float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Mixing parameter for sigma (lij > 0 => smaller sigma_12, lij < 0 => larger sigma_12)

&nbsp;&nbsp;&nbsp;&nbsp; **kij (optional, float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Mixing parameter for epsilon (kij > 0 => favours mixing, kij < 0 => favours separation) 

## Utility methods

Set- and get methods for interaction parameters, mixing parameters ...

### Table of contents
  * [Utility methods](#utility-methods)
    * [get_BH_diameters](#get_bh_diametersself-t)
    * [get_epsilon_matrix](#get_epsilon_matrixself-eps_div_k-kij)
    * [get_lambda_matrix](#get_lambda_matrixself-lambdas-lij)
    * [get_sigma_matrix](#get_sigma_matrixself-sigma)
    * [potential](#potentialself-i-j-r)
    * [potential_r](#potential_rself-i-j-r)
    * [potential_rr](#potential_rrself-i-j-r)


### `get_BH_diameters(self, T)`
Compute Barker-Henderson diameters

Args:
T (float) : Temperature [K]

Returns:
2d array : Barker-Henderson Diameters, indexed by component pair [m]
 

### `get_epsilon_matrix(self, eps_div_k, kij)`
Compute matrix of well-depths, given well depth of each component
Warning: Use of mixing parameters is not thouroughly tested.
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **eps_div_k (1d array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Well depth parameter of each component

&nbsp;&nbsp;&nbsp;&nbsp; **kij (2d array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Not in use, internal parameter `self.kij` is used for mixing. 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **2d array :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Well depth for each interaction pair. 

### `get_lambda_matrix(self, lambdas, lij)`
Compute pair-interaction $\lambda_r$ parameters, apply mixing parameter.
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **lambdas (1d array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Repulsive exponents for each pure-component interaction potential

&nbsp;&nbsp;&nbsp;&nbsp; **lij (1d array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Mixing parameters 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **2d array :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Repulsive exponent for each pair-interaction. 

### `get_sigma_matrix(self, sigma)`
Compute interaction parameter $sigma$ for each particle pair, applying mixing parameters given by `self.lij`.
Warning: Use of mixing parameters is not thouroughly tested.
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **sigma (1D array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  sigma-parameters [m] 

&nbsp;&nbsp;&nbsp;&nbsp; **Retunrs:** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; 

&nbsp;&nbsp;&nbsp;&nbsp; **2d array :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  N x N matrix of sigma parameters, where sigma_ij = 0.5 * (sigma_i + sigma_j), if self.lij = 0.

&nbsp;&nbsp;&nbsp;&nbsp; **Warning:** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Use of mixing parameters is not thouroughly tested.  

### `potential(self, i, j, r)`
Evaluate the interaction potential between types i and j at distance r

Args:
i, j (int) : Component indices
r (float) : Distance [m]
Returns:
float : Interaction potential [J]
 

### `potential_r(self, i, j, r)`
Evaluate the derivative of the interaction potential between types i and j at distance r

Args:
i, j (int) : Component indices
r (float) : Distance [m]
Returns:
float : First derivative of interaction potential [N]
 

### `potential_rr(self, i, j, r)`
Evaluate the second derivative of the interaction potential between types i and j at distance r

Args:
i, j (int) : Component indices
r (float) : Distance [m]
Returns:
float : Second derivative of interaction potential [N / m]
 

