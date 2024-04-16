---
layout: default
version: 
title: Methods in the py_KineticGas class
permalink: /vcurrent/py_KineticGas_methods.html
---

<!--- 
Generated at: 2024-04-16T15:54:30.853573
This is an auto-generated file, generated using the script at KineticGas/pyUtils/markdown_from_docstrings.py
The file is created by parsing the docstrings of the methods in the 
py_KineticGas class. For instructions on how to use the parser routines, see the
file KineticGas/pyUtils/markdown_from_docstrings.py--->

The `py_KineticGas` class, found in `pykingas/py_KineticGas.py`, is the core of the KineticGas Python interface. All Revised Enskog Theory models inherit from `py_KineticGas`. This is the class that contains the interface to all practical calculations that can be done from the python-side of KineticGas. Derived classes only implement specific functions for parameter handling etc.

## Table of contents
  * [The constructor](#the-constructor)
    * [\_\_init\_\_](#__init__self-comps-mole_weightsnone-n3-is_idealgasfalse)
  * [TV-property interfaces](#tv-property-interfaces)
    * [bulk_viscosity](#bulk_viscosityself-t-vm-x-nnone)
    * [conductivity_matrix](#conductivity_matrixself-t-vm-x-n2-formulationt-psi-frame_of_referencecom-use_thermal_conductivitynone)
    * [interdiffusion](#interdiffusionself-t-vm-x-nnone-use_independenttrue-dependent_idxnone-frame_of_referencecon-use_binarytrue-solvent_idxnone)
    * [interdiffusion_general](#interdiffusion_generalself-t-vm-x-nnone)
    * [resistivity_matrix](#resistivity_matrixself-t-vm-x-n2-formulationt-psi-frame_of_referencecom-use_thermal_conductivitynone)
    * [soret_coefficient](#soret_coefficientself-t-vm-x-nnone)
    * [thermal_conductivity](#thermal_conductivityself-t-vm-x-nnone)
    * [thermal_diffusion_coeff](#thermal_diffusion_coeffself-t-vm-x-nnone-use_independentfalse-dependent_idxnone-frame_of_referencecon-solvent_idxnone)
    * [thermal_diffusion_factor](#thermal_diffusion_factorself-t-vm-x-nnone)
    * [thermal_diffusion_ratio](#thermal_diffusion_ratioself-t-vm-x-nnone)
    * [viscosity](#viscosityself-t-vm-x-nnone)
  * [Tp-property interfaces](#tp-property-interfaces)
    * [interdiffusion_tp](#interdiffusion_tpself-t-p-x-nnone-use_independenttrue-dependent_idxnone-frame_of_referencecon-use_binarytrue-solvent_idxnone)
    * [thermal_coductivity_tp](#thermal_coductivity_tpself-t-p-x-nnone)
    * [thermal_diffusion_coeff_tp](#thermal_diffusion_coeff_tpself-t-p-x-nnone-use_independentfalse-dependent_idxnone-frame_of_referencecon-solvent_idxnone)
    * [thermal_diffusion_factor_tp](#thermal_diffusion_factor_tpself-t-p-x-nnone)
    * [viscosity_tp](#viscosity_tpself-t-p-x-nnone)
  * [Frame of Reference transformations](#frame-of-reference-transformations)
    * [get_com_2_con_matr](#get_com_2_con_matrself-x)
    * [get_com_2_cov_matr](#get_com_2_cov_matrself-t-vm-x)
    * [get_com_2_for_matr](#get_com_2_for_matrself-t-vm-x-for-**kwargs)
    * [get_com_2_solv_matr](#get_com_2_solv_matrself-x-solvent_idx)
    * [get_solv_2_solv_matr](#get_solv_2_solv_matrself-x-prev_solv_idx-new_solv_idx)
    * [get_zarate_W_matr](#get_zarate_w_matrself-x-dependent_idx)
    * [get_zarate_X_matr](#get_zarate_x_matrself-x-dependent_idx)
  * [Interfaces to C++ methods](#interfaces-to-c++-methods)
    * [get_collision_diameters](#get_collision_diametersself-particle_density-t-x)
    * [get_conductivity_matrix](#get_conductivity_matrixself-particle_density-t-mole_fracs-nnone)
    * [get_conductivity_vector](#get_conductivity_vectorself-particle_density-t-mole_fracs-n)
    * [get_diffusion_vector](#get_diffusion_vectorself-particle_density-t-mole_fracs-nnone)
    * [get_rdf](#get_rdfself-particle_density-t-x)
  * [Utility methods](#utility-methods)
    * [check_valid_composition](#check_valid_compositionself-x)
    * [compute_cond_vector](#compute_cond_vectorself-particle_density-t-mole_fracs-nnone)
    * [compute_diffusion_coeff_vector](#compute_diffusion_coeff_vectorself-particle_density-t-mole_fracs-nnone)
    * [compute_dth_vector](#compute_dth_vectorself-particle_density-t-mole_fracs-nnone)
    * [compute_visc_vector](#compute_visc_vectorself-t-particle_density-mole_fracs-nnone)
    * [get_Eij](#get_eijself-vm-t-x)
    * [get_P_factors](#get_p_factorsself-vm-t-x)
    * [reshape_diffusion_coeff_vector](#reshape_diffusion_coeff_vectorself-d)

## The constructor

The constructor

### Table of contents
  * [The constructor](#the-constructor)
    * [\_\_init\_\_](#__init__self-comps-mole_weightsnone-n3-is_idealgasfalse)


### `__init__(self, comps, mole_weights=None, N=3, is_idealgas=False)`
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **comps (str):** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Comma-separated list of components, following ThermoPack-convention

&nbsp;&nbsp;&nbsp;&nbsp; **mole_weights (1d array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Mole weights [g/mol]. Will be used instead of database values if provided

&nbsp;&nbsp;&nbsp;&nbsp; **is_idealgas (bool) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  If true, radial distribution function is unity, if false use radial distribution function of modelIn addition, several density-dependent factors are set to zero, to ensure consistency with the ideal gas law / Gibbs-Duhem for an ideal gas.  

## TV-property interfaces

Computing properties as a function of temperature and volume.

### Table of contents
  * [TV-property interfaces](#tv-property-interfaces)
    * [bulk_viscosity](#bulk_viscosityself-t-vm-x-nnone)
    * [conductivity_matrix](#conductivity_matrixself-t-vm-x-n2-formulationt-psi-frame_of_referencecom-use_thermal_conductivitynone)
    * [interdiffusion](#interdiffusionself-t-vm-x-nnone-use_independenttrue-dependent_idxnone-frame_of_referencecon-use_binarytrue-solvent_idxnone)
    * [interdiffusion_general](#interdiffusion_generalself-t-vm-x-nnone)
    * [resistivity_matrix](#resistivity_matrixself-t-vm-x-n2-formulationt-psi-frame_of_referencecom-use_thermal_conductivitynone)
    * [soret_coefficient](#soret_coefficientself-t-vm-x-nnone)
    * [thermal_conductivity](#thermal_conductivityself-t-vm-x-nnone)
    * [thermal_diffusion_coeff](#thermal_diffusion_coeffself-t-vm-x-nnone-use_independentfalse-dependent_idxnone-frame_of_referencecon-solvent_idxnone)
    * [thermal_diffusion_factor](#thermal_diffusion_factorself-t-vm-x-nnone)
    * [thermal_diffusion_ratio](#thermal_diffusion_ratioself-t-vm-x-nnone)
    * [viscosity](#viscosityself-t-vm-x-nnone)


### `bulk_viscosity(self, T, Vm, x, N=None)`
Not implemented

Raises:
NotImplementedError
 

### `conductivity_matrix(self, T, Vm, x, N=2, formulation='T-psi', frame_of_reference='CoM', use_thermal_conductivity=None)`
Compute the conductivity matrix $L$, for use in NET calculations. The Flux/Force formulation used in the NET
model is selected using the `formulation` kwarg. Currently implemented formulations are:
----------------------------------------------------------------------------------------------

`'T-psi'`:

----------------------------------------------------

$$ J_q = L_{qq} \nabla (1 / T) - \sum_{i=1}^{N_c-1}(1 / T) L_{qi} \nabla_T \Psi_i $$

$$ J_i = L_{iq} \nabla (1 / T) - \sum_{j=1}^{N_c-1}(1 / T) L_{ij} \nabla_T \Psi_j $$

Where

$$ \Psi_i = \mu_i - \mu_{N_c} $$

and $N_c$ denotes the number of components. The last component is used as the dependent component.
The fluxes in this formulation are on a *mass basis*, and the formulation is implemented in the centre of
mass frame of reference. The formulation is only implemented for ideal gases.
----------------------------------------------------------------------------------------------
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature [K]

&nbsp;&nbsp;&nbsp;&nbsp; **Vm (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar volume [m3 / mol]

&nbsp;&nbsp;&nbsp;&nbsp; **x (1darray) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition [-]

&nbsp;&nbsp;&nbsp;&nbsp; **N (int, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Enskog approximation order (must be >= 2 for thermal effects)

&nbsp;&nbsp;&nbsp;&nbsp; **formulation (str, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  The NET formulation for which to compute the conductivity matrix.

&nbsp;&nbsp;&nbsp;&nbsp; **frame_of_reference (str, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  The frame of reference ('CoM', 'CoN', 'CoV', or 'solvent'). Default: 'CoM'.

&nbsp;&nbsp;&nbsp;&nbsp; **use_thermal_conductivity (callable, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  External thermal conductivity model. Assumed to have the signatureuse_thermal_conductivity(T, Vm, x), returning the thermal conductivity in [W / m K]. Defaults to None. If no model is supplied, KineticGas is used to compute thermal conductivity.  

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **ndarray :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  The conductivity matrix, contents will vary depending on the `formulation` kwarg. 

### `interdiffusion(self, T, Vm, x, N=None, use_independent=True, dependent_idx=None, frame_of_reference='CoN', use_binary=True, solvent_idx=None)`
Compute the interdiffusion coefficients [m^2 / s]. Default definition is

$$ J_i^{(n, n)} = - \sum_{j \neq l} D_{ij} \nabla n_j, \nabla T = \nabla p = F_k = 0 \forall k $$

where the flux, $J_i^{(n, n)}$ is on a molar basis, in the molar frame of reference, and $j \neq l$ is an
independent set of forces with $l=$ `dependent_idx`.
For fluxes in other frames of reference, use the `frame_of_reference` kwarg.
For the diffusion coefficients describing the fluxes' response to all forces (not independent) the definition
is:

$$ J_i^{(n, n)} = - \sum_j D_{ij} \nabla n_j, \nabla T = \nabla p = F_k = 0 \forall k $$

use the `use_independent` and `dependent_idx` kwargs to switch between these definitions.
See: Eq. (17-20) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature (K)

&nbsp;&nbsp;&nbsp;&nbsp; **Vm (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar volume (m3 / mol)

&nbsp;&nbsp;&nbsp;&nbsp; **x (array_like) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  composition (mole fractions)

&nbsp;&nbsp;&nbsp;&nbsp; **N (int, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Enskog approximation order. Default set on model initialisation.

&nbsp;&nbsp;&nbsp;&nbsp; **use_binary (bool, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  If the mixture is binary, and an independent set of fluxes and forces is considered, i.e.`use_independent=True`, the diffusion coefficients will be exactly equal with opposite sign. Setting `use_binary=True` will in that case only return the coefficient describing the independent flux-force relation. 

&nbsp;&nbsp;&nbsp;&nbsp; **use_independent (bool, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Return diffusion coefficients for independent set of forces.

&nbsp;&nbsp;&nbsp;&nbsp; **dependent_idx (int, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Index of the dependent molar density gradient (only if `use_independent=True`, the default behaviour).Defaults to last component, except when `frame_of_reference='solvent'`, in which case default is equal to `solvent_idx`. 

&nbsp;&nbsp;&nbsp;&nbsp; **frame_of_reference (str, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Which frame of reference the diffusion coefficients apply to. Defaultis `'CoN'`. Can be `'CoN'` (molar FoR), `'CoM'` (barycentric FoR), `'solvent'` (solvent FoR), `'zarate'` (See Memo on definitions of the diffusion coefficient), `'zarate_x'` ($D^{(x)}$ as defined by Ortiz de Zárate, doi 10.1140/epje/i2019-11803-2) or `'zarate_w'` ($D^{(w)}$ as defined by Ortiz de Zárate). 

&nbsp;&nbsp;&nbsp;&nbsp; **solvent_idx (int, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Index of component identified as solvent (only when using `frame_of_reference='solvent'`) 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **(ndarray or float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Diffusion coefficients, shape varies based on options and number of components. Unit [m^2 / s] 

### `interdiffusion_general(self, T, Vm, x, N=None)`
Compute the 'Kinetic CoM diffusion coefficients', defined by

$$ J_i^{(n, m)} = - \sum_j D_{ij} \nabla n_j, \nabla T = \nabla p = F_k = 0 \forall k $$

**For end-users, see `interdiffusion`**
See Eq. (19) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature (K)

&nbsp;&nbsp;&nbsp;&nbsp; **Vm (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar volume (m3 / mol)

&nbsp;&nbsp;&nbsp;&nbsp; **x (array_like) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  composition (mole fractions)

&nbsp;&nbsp;&nbsp;&nbsp; **N (int, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Enskog approximation order 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **(2D array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Array of the (not independent) $D^{(K, m)}$ diffusion coefficients. Unit [m^2 / s] 

### `resistivity_matrix(self, T, Vm, x, N=2, formulation='T-psi', frame_of_reference='CoM', use_thermal_conductivity=None)`
Compute the resistivity matrix $R = L^{-1}$, for use in NET calculations. The Flux/Force formulation used in the NET
model is selected using the `formulation` kwarg. Currently implemented formulations are:
----------------------------------------------------------------------------------------------

`'T-psi'`:

----------------------------------------------------

$$ J_q = L_{qq} \nabla (1 / T) - \sum_{i=1}^{N_c-1}(1 / T) L_{qi} \nabla_T \Psi_i $$

$$ J_i = L_{iq} \nabla (1 / T) - \sum_{j=1}^{N_c-1}(1 / T) L_{ij} \nabla_T \Psi_j $$

Where

$$ \Psi_i = \mu_i - \mu_{N_c} $$

and $N_c$ denotes the number of components. The last component is used as the dependent component.
The fluxes in this formulation are on a *mass basis*, and the formulation is implemented in the centre of
mass frame of reference. The formulation is only implemented for ideal gases.
----------------------------------------------------------------------------------------------
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature [K]

&nbsp;&nbsp;&nbsp;&nbsp; **Vm (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar volume [m3 / mol]

&nbsp;&nbsp;&nbsp;&nbsp; **x (1darray) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition [-]

&nbsp;&nbsp;&nbsp;&nbsp; **N (int, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Enskog approximation order (must be >= 2 for thermal effects)

&nbsp;&nbsp;&nbsp;&nbsp; **formulation (str, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  The NET formulation for which to compute the conductivity matrix.

&nbsp;&nbsp;&nbsp;&nbsp; **frame_of_reference (str, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  The frame of reference ('CoM', 'CoN', 'CoV', or 'solvent'). Default: 'CoM'.

&nbsp;&nbsp;&nbsp;&nbsp; **use_thermal_conductivity (callable, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  External thermal conductivity model. Assumed to have the signaturethermal_conductivity(T, Vm, x), returning the thermal conductivity in [W / m K]. Defaults to None. If no model is supplied, KineticGas is used to compute the thermal conductivity.  

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **ndarray :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  The resistivity matrix, contents will vary depending on the `formulation` kwarg. 

### `soret_coefficient(self, T, Vm, x, N=None)`
Compute the Soret coefficients, $S_{T, ij}$ defined as
$S_{T, ij} = \alpha_{ij} / T$
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature [K]

&nbsp;&nbsp;&nbsp;&nbsp; **Vm (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar volume [m3 / mol]

&nbsp;&nbsp;&nbsp;&nbsp; **x (array_like) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition

&nbsp;&nbsp;&nbsp;&nbsp; **N (int, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Enskog approximation order (>= 2) 

### `thermal_conductivity(self, T, Vm, x, N=None)`
Compute the thermal conductivity, $\lambda$. For models initialized with `is_idealgas=True`, the thermal
conductivity is not a function of density (i.e. $d \lambda / d V_m = 0$).
See Eq. (13) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature [K]

&nbsp;&nbsp;&nbsp;&nbsp; **Vm (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar volume [m3 / mol]

&nbsp;&nbsp;&nbsp;&nbsp; **x (array_like) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition [-]

&nbsp;&nbsp;&nbsp;&nbsp; **N (int, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Enskog approximation order (>= 2) 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **(float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  The thermal conductivity of the mixture. 

### `thermal_diffusion_coeff(self, T, Vm, x, N=None, use_independent=False, dependent_idx=None, frame_of_reference='CoN', solvent_idx=None)`
Compute thermal diffusion coefficients, $D_{T,i}$ [mol / m^2 s]
Default definition is

$$ J_i^{(n, n)} = D_{T,i} \nabla \ln T - \sum_j D_{ij} \nabla n_j, \nabla p = F_k = 0 \forall k $$

where the flux, $J_i^{(n, n)}$ is on a molar basis, in the molar frame of reference. For fluxes in other frames
of reference, use the 'frame_of_reference' kwarg. For the diffusion coefficients corresponding to an
independent set of forces, defined by

$$ J_i^{(n, n)} = D_{T,i} \nabla \ln T - \sum_{j \neq l} D_{ij} \nabla n_j, \nabla p = F_k = 0 \forall k $$

where $l$ is the index of the dependent molar density gradient, use the 'use_independent' and 'dependent_idx' kwargs.
See Eq. (23) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature [K]

&nbsp;&nbsp;&nbsp;&nbsp; **Vm (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar volume [m^3 / mol]

&nbsp;&nbsp;&nbsp;&nbsp; **x (array_like) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Mole fractions [-]

&nbsp;&nbsp;&nbsp;&nbsp; **N (int, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Enskog approximation order (>=2). Defaults to 2.

&nbsp;&nbsp;&nbsp;&nbsp; **use_independent (bool, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Return diffusion coefficients for independent set of forces.

&nbsp;&nbsp;&nbsp;&nbsp; **dependent_idx (int, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Index of the dependent molar density gradient (only if use_dependent=True).Defaults to last component. 

&nbsp;&nbsp;&nbsp;&nbsp; **frame_of_reference (str, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  What frame of reference the coefficients apply to. Valid options are`'CoM'` (centre of mass / barycentric), `'CoN'` (centre of moles), `'CoV'` (centre of volume) `'solvent'` (together with `solvent_idx`) or `'zarate'`, for the coefficients as defined by Ortiz de Zarate (doi 10.1140/epje/i2019-11803-2). 

&nbsp;&nbsp;&nbsp;&nbsp; **solvent_idx (int, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Index of component identified as solvent (only when using `frame_of_reference='solvent'`) 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **(1D array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Thermal diffusion coefficients. Unit [mol m^2 / s] 

### `thermal_diffusion_factor(self, T, Vm, x, N=None)`
Compute the thermal diffusion factors $\alpha_{ij}$, defined by

$$ \alpha_{ij} = k_{T, i} - k_{T, j} $$

where $k_{T,i}$ are the thermal diffusion ratios. This definition implies that

$$ \nabla \ln (n_i / n_j) = - \alpha_{ij} \nabla \ln T $$

when the mass fluxes vanish. The thermal diffusion factors are independent of the frame of reference.
See Eq. (29) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature [K]

&nbsp;&nbsp;&nbsp;&nbsp; **Vm (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar volume [m3 / mol]

&nbsp;&nbsp;&nbsp;&nbsp; **x (array_like) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition [-]

&nbsp;&nbsp;&nbsp;&nbsp; **N (int, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Enskog approximation order (>= 2) 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **(2D array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  The thermal diffusion factors of the mixture. Dimensionless. 

### `thermal_diffusion_ratio(self, T, Vm, x, N=None)`
Calculate the "independent" thermal diffusion ratios, $k_{T, i}$ defined by

$$ J_i^{(n, n)} = - \sum_{j \neq l} D_{ij}^{(I, n)} ( \nabla n_j + n_j k_{T, j} D_{T, j}^{(I, n)} \nabla \ln T ) $$

and

$$ \sum_i x_i \sum_j [\delta_{ij} + (4 \pi / 3) \sigma_{ij}^3 n_j M_{ij} \chi_{ij} - (n_j k_{T,j} / k_B T) ( d \mu_i / n_j )_{T,n_{l \neq j}}] = 0 $$

This definition implies that

$$ \nabla n_j = -n_j k_{T,j} \nabla \ln T\ \forall j $$

when all mass fluxes vanish. For models initialised with `is_idealgas=True`, the second equation is replaced with

$$ \sum_i x_i k_{T,i} = 1 $$

The thermal diffusion ratios are independent of the frame of reference.
See Eq. (26-27) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature [K]

&nbsp;&nbsp;&nbsp;&nbsp; **Vm (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar volume [m3 / mol]

&nbsp;&nbsp;&nbsp;&nbsp; **x (array_like) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition [-]

&nbsp;&nbsp;&nbsp;&nbsp; **N (int, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Enskog approximation order (>= 2) 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **(1D array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  The thermal diffusion ratio of each component. Unit Dimensionless. 

### `viscosity(self, T, Vm, x, N=None)`
Compute the shear viscosity, $\eta$. For models initialized with `is_idealgas=True`, the shear viscosity
is not a function of density (i.e. $d \eta / d V_m = 0). See Eq. (12) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature [K]

&nbsp;&nbsp;&nbsp;&nbsp; **Vm (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar volume [m3 / mol]

&nbsp;&nbsp;&nbsp;&nbsp; **x (array_like) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition [-]

&nbsp;&nbsp;&nbsp;&nbsp; **N (int, optional) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Enskog approximation order 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **(float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  The shear viscosity of the mixture. 

## Tp-property interfaces

Computing properties as a function of temperature and pressure. Simply forwards calls to the appropriate TV-property methods after computing molar volume using the internal equation of state object.

### Table of contents
  * [Tp-property interfaces](#tp-property-interfaces)
    * [interdiffusion_tp](#interdiffusion_tpself-t-p-x-nnone-use_independenttrue-dependent_idxnone-frame_of_referencecon-use_binarytrue-solvent_idxnone)
    * [thermal_coductivity_tp](#thermal_coductivity_tpself-t-p-x-nnone)
    * [thermal_diffusion_coeff_tp](#thermal_diffusion_coeff_tpself-t-p-x-nnone-use_independentfalse-dependent_idxnone-frame_of_referencecon-solvent_idxnone)
    * [thermal_diffusion_factor_tp](#thermal_diffusion_factor_tpself-t-p-x-nnone)
    * [viscosity_tp](#viscosity_tpself-t-p-x-nnone)


### `interdiffusion_tp(self, T, p, x, N=None, use_independent=True, dependent_idx=None, frame_of_reference='CoN', use_binary=True, solvent_idx=None)`
Compute molar volume using the internal equation of state (`self.eos`), assuming vapour, and pass the call to
`self.interdiffusion`. See `self.interdiffusion` for documentation.
 

### `thermal_coductivity_tp(self, T, p, x, N=None)`
Compute molar volume using the internal equation of state (`self.eos`), assuming vapour, and pass the call to
`self.thermal_conductivity`. See `self.thermal_conductivity` for documentation.
 

### `thermal_diffusion_coeff_tp(self, T, p, x, N=None, use_independent=False, dependent_idx=None, frame_of_reference='CoN', solvent_idx=None)`
Compute molar volume using the internal equation of state (`self.eos`), assuming vapour, and pass the call to
`self.thermal_diffusion_coeff`. See `self.thermal_diffusion_coeff` for documentation.
 

### `thermal_diffusion_factor_tp(self, T, p, x, N=None)`
Compute molar volume using the internal equation of state (`self.eos`), assuming vapour, and pass the call to
`self.thermal_diffusion_factor`. See `self.thermal_diffusion_factor` for documentation.
 

### `viscosity_tp(self, T, p, x, N=None)`
Compute molar volume using the internal equation of state (`self.eos`), assuming vapour, and pass the call to
`self.viscosity`. See `self.viscosity` for documentation.
 

## Frame of Reference transformations

Generate matrices for Frame of Reference transformations. See the supportingmaterial of [Revised Enskog Theory for Mie fluids](https://doi.org/10.1063/5.0149865) for details.

### Table of contents
  * [Frame of Reference transformations](#frame-of-reference-transformations)
    * [get_com_2_con_matr](#get_com_2_con_matrself-x)
    * [get_com_2_cov_matr](#get_com_2_cov_matrself-t-vm-x)
    * [get_com_2_for_matr](#get_com_2_for_matrself-t-vm-x-for-**kwargs)
    * [get_com_2_solv_matr](#get_com_2_solv_matrself-x-solvent_idx)
    * [get_solv_2_solv_matr](#get_solv_2_solv_matrself-x-prev_solv_idx-new_solv_idx)
    * [get_zarate_W_matr](#get_zarate_w_matrself-x-dependent_idx)
    * [get_zarate_X_matr](#get_zarate_x_matrself-x-dependent_idx)


### `get_com_2_con_matr(self, x)`
Get transformation matrix from centre of mass (CoM) to centre of moles (CoN).
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **x (array_like) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition [-] 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **(2d array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Transformation matrix $\Psi^{n \leftmapsto m}$ 

### `get_com_2_cov_matr(self, T, Vm, x)`
Get centre of mass (CoM) to centre of volume (CoV) transformation matrix
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature [K]

&nbsp;&nbsp;&nbsp;&nbsp; **Vm (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar volume [m3 / mol]

&nbsp;&nbsp;&nbsp;&nbsp; **x (array_like) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition [-] 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **2d array :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  The transformation matrix $\Psi^{V \leftmapsto m}$. 

### `get_com_2_for_matr(self, T, Vm, x, FoR, **kwargs)`
Dispatcher to get a specific 'change of frame of reference' matrix for transformation from
centre of mass to 'FoR'.
Returns the appropriate matrix for the transformations derived in Appendix A of ... to transform a flux
from the centre of mass frame of reference to the 'FoR' frame of reference by the transformation

$$ J^{(n, FoR)} = \psi @ J^{(n, m)} $$

where $\psi$ is the matrix returned by this method, and $J$ is the vector of (all) molar fluxes, with the
subscript indicating the frame of reference.
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature [K]

&nbsp;&nbsp;&nbsp;&nbsp; **Vm (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar volume [m3 / mol]

&nbsp;&nbsp;&nbsp;&nbsp; **x (array_like) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition [-] 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **(2Darray) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  The $N$ x $N$ transformation matrix to transform the fluxes, with $N$ being the number ofcomponents.  

### `get_com_2_solv_matr(self, x, solvent_idx)`
Get transformation matrix from centre of mass (CoM) to solvent (solvent) frame of reference
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **x (array_like) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition [-]

&nbsp;&nbsp;&nbsp;&nbsp; **solvent_idx (int) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  The component index of the solvent. 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **2darray :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  The transformation matrix $\Psi^{n_k \leftmapsto m}$, where $k$ is the `solvent_idx`. 

### `get_solv_2_solv_matr(self, x, prev_solv_idx, new_solv_idx)`
Get solvent-to-solvent frame of reference transformation matrix
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **x (array_like) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition [-]

&nbsp;&nbsp;&nbsp;&nbsp; **prev_solv_idx (int) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Component index of the old (current) solvent

&nbsp;&nbsp;&nbsp;&nbsp; **new_solv_idx (int) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Component index of the new solvent 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **2d array :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  The transformation matrix $\Psi^{n_k \leftmapsto n_l}$, where $k$ is `new_solv_idx` and $l$ is `prev_solv_idx`. 

### `get_zarate_W_matr(self, x, dependent_idx)`
Compute the matrix $W$ as defined by Zárate. See (Definition of frame-invariant thermodiffusion and Soret coefficients for ternary mixtures)
and memo on diffusion coefficient definitions.
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **x (array_like) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition [-]

&nbsp;&nbsp;&nbsp;&nbsp; **dependent_idx (int) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Index of the dependent species

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **2d array :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  The transformation matrix $W$ 

### `get_zarate_X_matr(self, x, dependent_idx)`
Compute the matrix $X$ as defined by Zárate. See (Definition of frame-invariant thermodiffusion and Soret coefficients for ternary mixtures)
and memo on diffusion coefficient definitions.
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **x (array_like) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition [-]

&nbsp;&nbsp;&nbsp;&nbsp; **dependent_idx (int) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Index of the dependent species

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **2d array :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  The transformation matrix $X$ 

## Interfaces to C++ methods

Lightweight wrappers for the most commonly used C++ side methods.

### Table of contents
  * [Interfaces to C++ methods](#interfaces-to-c++-methods)
    * [get_collision_diameters](#get_collision_diametersself-particle_density-t-x)
    * [get_conductivity_matrix](#get_conductivity_matrixself-particle_density-t-mole_fracs-nnone)
    * [get_conductivity_vector](#get_conductivity_vectorself-particle_density-t-mole_fracs-n)
    * [get_diffusion_vector](#get_diffusion_vectorself-particle_density-t-mole_fracs-nnone)
    * [get_rdf](#get_rdfself-particle_density-t-x)


### `get_collision_diameters(self, particle_density, T, x)`
Compute collision diameters given by Eq. (40) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
*Note* Returns zeros for models initialised with is_idealgas=True.
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **particle_density (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Particle density (not molar!) [1 / m3]

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature [K]

&nbsp;&nbsp;&nbsp;&nbsp; **x (list[float]) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition [-] 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **2d array :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Collision diameters [m], indexed by component pair. 

### `get_conductivity_matrix(self, particle_density, T, mole_fracs, N=None)`
Compute the elements of the matrix corresponding to the set of equations that must be solved for the
thermal response function sonine polynomial expansion coefficients:
Eq. (6) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **particle_density (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Particle density (not molar!) [1 / m3]

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature [K]

&nbsp;&nbsp;&nbsp;&nbsp; **mole_fracs (list[float]) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition [-]

&nbsp;&nbsp;&nbsp;&nbsp; **N (Optional, int) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Enskog approximation order. 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **2Darray :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  ($N N_c$ x $N N_c$) matrix, where $N$ is the Enskog approximation order and $N_c$ isthe number of components.  

### `get_conductivity_vector(self, particle_density, T, mole_fracs, N)`
Compute the right-hand side vector to the set of equations that must be solved for the
thermal response function Sonine polynomial expansion coefficients:
Eq. (6) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **particle_density (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Particle density (not molar!) [1 / m3]

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature [K]

&nbsp;&nbsp;&nbsp;&nbsp; **mole_fracs (list<float>) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition [-]

&nbsp;&nbsp;&nbsp;&nbsp; **N (Optional, int) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Enskog approximation order. 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **(1Darray) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  ($N N_c$,) vector, where $N$ is the Enskog approximation order and $N_c$ isthe number of components.  

### `get_diffusion_vector(self, particle_density, T, mole_fracs, N=None)`
Compute the right-hand side vector to the set of equations that must be solved for the
diffusive response function Sonine polynomial expansion coefficients.
Eq. (10) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **particle_density (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Particle density (not molar!) [1 / m3]

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature [K]

&nbsp;&nbsp;&nbsp;&nbsp; **mole_fracs (list[float]) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition [-]

&nbsp;&nbsp;&nbsp;&nbsp; **N (Optional, int) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Enskog approximation order. 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **(1Darray) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  ($N N_c^2$,) vector, where $N$ is the Enskog approximation order and $N_c$ isthe number of components.  

### `get_rdf(self, particle_density, T, x)`
Compute the radial distribution function at contact
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **particle_density (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Particle density (not molar!) [1 / m3]

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature [K]

&nbsp;&nbsp;&nbsp;&nbsp; **x (list[float]) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition [-] 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **2d array :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  RDF at contact, indexed by component pair. 

## Utility methods

Methods for various helper computations

### Table of contents
  * [Utility methods](#utility-methods)
    * [check_valid_composition](#check_valid_compositionself-x)
    * [compute_cond_vector](#compute_cond_vectorself-particle_density-t-mole_fracs-nnone)
    * [compute_diffusion_coeff_vector](#compute_diffusion_coeff_vectorself-particle_density-t-mole_fracs-nnone)
    * [compute_dth_vector](#compute_dth_vectorself-particle_density-t-mole_fracs-nnone)
    * [compute_visc_vector](#compute_visc_vectorself-t-particle_density-mole_fracs-nnone)
    * [get_Eij](#get_eijself-vm-t-x)
    * [get_P_factors](#get_p_factorsself-vm-t-x)
    * [reshape_diffusion_coeff_vector](#reshape_diffusion_coeff_vectorself-d)


### `check_valid_composition(self, x)`
Check that enough mole fractions are supplied for the initialised model. Also check that they sum to unity.
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **x (array_like) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition Raises 

&nbsp;&nbsp;&nbsp;&nbsp; **IndexError :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  If wrong number of mole fractions is supplied.

&nbsp;&nbsp;&nbsp;&nbsp; **RuntimeWarning :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  If mole fractions do not sum to unity. 

### `compute_cond_vector(self, particle_density, T, mole_fracs, N=None)`
Compute the thermal response function Sonine polynomial expansion coefficients by solving the set of equations

$$\Lambda \ell = \lambda$$

Corresponding to Eq. (6) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
Where $\Lambda$ is the matrix returned by the c++ method `get_conductivity_matrix`, and $\lambda$ is the vector
returned by the c++ method `get_conductivity_vector`.
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **particle_density (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Particle density (not molar!) [1 / m3]

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature [K]

&nbsp;&nbsp;&nbsp;&nbsp; **mole_fracs (list[float]) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition [-]

&nbsp;&nbsp;&nbsp;&nbsp; **N (Optional, int) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Enskog approximation order. 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **1Darray :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  ($N N_c$,) vector, where $N$ is the Enskog approximation order and $N_c$ isthe number of components. The vector is ordered as 

&nbsp;&nbsp;&nbsp;&nbsp; **l[:** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; N_c] = [$\ell_1^{(0)}$, $\ell_2^{(0)}$, ..., $\ell_{N_c}^{(0)}$]

&nbsp;&nbsp;&nbsp;&nbsp; **l[N_c :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  2 * N_c] = [$\ell_1^{(1)}$, $\ell_2^{(1)}$, ..., $\ell_{N_c}^{(1)}$]... etc ... Where subscripts indicate component indices, and superscripts are Enskog approximation order summation indices.  

### `compute_diffusion_coeff_vector(self, particle_density, T, mole_fracs, N=None)`
Compute the diffusive response function Sonine polynomial expansion coefficients by solving the set of equations

$$D d = \delta$$

Corresponding to Eq. (10) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
Where $D$ is the matrix returned by the c++ method `get_diffusion_matrix`, and $\delta$ is the vector
returned by the c++ method `get_diffusion_vector`.
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **particle_density (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Particle density (not molar!) [1 / m3]

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature [K]

&nbsp;&nbsp;&nbsp;&nbsp; **mole_fracs (list[float]) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition [-]

&nbsp;&nbsp;&nbsp;&nbsp; **N (Optional, int) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Enskog approximation order. 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **1Darray :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  ($N N_c^2$,) vector, where $N$ is the Enskog approximation order and $N_c$ isthe number of components, containing the diffusive response function sonine polynomial expansion coefficients ($d_{i, j}^{(q)}$). See `reshape_diffusive_coeff_vector` for help on practical usage.  

### `compute_dth_vector(self, particle_density, T, mole_fracs, N=None)`
Compute the coefficients $d_i^{(J = 0)}$, by solving the set of equations

$$\sum_j d_{i,j}^{(0)} d_j^{\vec{J} = 0} = \ell_i^{(0)}$$,

i.e. Eq. (15) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **particle_density (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Particle density (not molar!) [1 / m3]

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature [K]

&nbsp;&nbsp;&nbsp;&nbsp; **mole_fracs (list<float>) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition [-]

&nbsp;&nbsp;&nbsp;&nbsp; **N (Optional, int) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Enskog approximation order. 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **1Darray :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Vector of the $d_i^{(J = 0)}$ coefficients. 

### `compute_visc_vector(self, T, particle_density, mole_fracs, N=None)`
Compute the viscous response function Sonine polynomial expansion coefficients by solving the set of equations

$$\Beta b = \beta$$

Corresponding to Eq. (8) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
Where $\Beta$ is the matrix returned by the c++ method `get_viscosity_matrix`, and $\beta$ is the vector
returned by the c++ method `get_viscosity_vector`.
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **particle_density (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Particle density (not molar!) [1 / m3]

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature [K]

&nbsp;&nbsp;&nbsp;&nbsp; **mole_fracs (list<float>) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition [-]

&nbsp;&nbsp;&nbsp;&nbsp; **N (Optional, int) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Enskog approximation order. 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **1D array :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  ($N N_c$,) vector, where $N$ is the Enskog approximation order and $N_c$ isthe number of components, ordered as 

&nbsp;&nbsp;&nbsp;&nbsp; **b[:** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; N_c] = [b_1^{(0)}, b_2^{(0)}, ..., b_{N_c}^{(0)}]

&nbsp;&nbsp;&nbsp;&nbsp; **b[N_c:** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  2 * N_c] = [b_1^{(1)}, b_2^{(1)}, ..., b_{N_c}^{(1)}]... etc ... where subscripts denote component indices and superscripts denote Enskog approximation summation indices.  

### `get_Eij(self, Vm, T, x)`
Compute the factors

$$ ( n_i / k_B T ) (d \mu_i / d n_j)_{T, n_{k \neq j}}, $$

where $n_i$ is the molar density of species $i$.
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **Vm (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar volume [m3 / mol]

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature [K]

&nbsp;&nbsp;&nbsp;&nbsp; **x (array_like) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **(2D array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  The factors E[i][j] = $ ( n_i / k_B T ) (d \mu_i / d n_j)_{T, n_{k \neq j}}$, where $n_i$is the molar density of species $i$. Unit [1 / mol]  

### `get_P_factors(self, Vm, T, x)`
Compute the factors $\Xi_i = \sum_j E_{ij}$, where $E_{ij}$ are the factors computed by `get_Eij`.
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **Vm (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar volume [m3 / mol]

&nbsp;&nbsp;&nbsp;&nbsp; **T (float) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Temperature [K]

&nbsp;&nbsp;&nbsp;&nbsp; **x (array_like) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  Molar composition 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **(1D array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  The factors $\Xi_i$, Unit [1 / mol] 

### `reshape_diffusion_coeff_vector(self, d)`
The vector returned by `compute_diffusion_coeff_vector` contains the diffusive response function sonine
polynomial expansion coefficients (eg. $d_{i, j}^{(q)}$ ).
To more easily access the correct coefficients, this method reshapes the vector to a matrix indiced
as `d[i][q][j]` where i and j are component indices, and q refferes to the approximation order summation index.
 

#### Args:

&nbsp;&nbsp;&nbsp;&nbsp; **d (1D array) :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  The array returned by `get_diffusion_coeff_vector` 

#### Returns:

&nbsp;&nbsp;&nbsp;&nbsp; **3D array :** 

&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  The matrix of $d_{i, j}^{(q)}$ coefficients ordered as `d[i][q][j]` 

