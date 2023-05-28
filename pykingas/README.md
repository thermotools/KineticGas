<!--- 
Generated at: 2023-05-28T10:21:01.818266
This is an auto-generated file, generated using the script at KineticGas/docs/join_docs.py
The file is created by joining the contents of the files
    KineticGas/docs/markdown/
        header.md
        toc_pypi.md
        cite_acknowl_licence.md
        pypi_structure.md
        getting_started_py.md
        fluid_identifiers.md
--->
# KineticGas
KineticGas is an implementation of Revised Enskog Theory (RET) for spherical potentials. The most notable of which is the implementation of RET-Mie, the Revised Enskog Theory for Mie fluids. 

The package is implemented mostly in C++ to handle the numerical computations involved in evaluating the collision integrals and the radial distribution function at contact for the target fluids, with the possibility of setting up multithreading at compile time.

KineticGas can be used to predict diffusion coefficients, thermal diffusion coefficients, viscosities and thermal conductivities in gas mixtures, and is reliable over a large range of temperatures and pressures. The package also contains an extensive database of fluid parameters collected from the open literature.

![](https://github.com/thermotools/KineticGas/blob/main/docs/figures/all.gif?raw=true)


## Table of contents
   * [Getting started](#getting-started-in-python)
     * [Initializing a model](#Initializing-a-model)
     * [Making predictions](#Making-predictions)
   * [Advanced usage]()
     * [Modifying and adding fluids]()
   * [Structure](#Structure)
   * [Fluid indentifiers](#fluid-identifiers)

## Please cite

KineticGas has been developed throughout a series of two works. If you are referencing the package, please cite the works

   * [Revised Enskog theory for Mie fluids: Prediction of diffusion coefficients, thermal diffusion coefficients, viscosities and thermal conductivities]() (Vegard G. Jervell and Øivind Wilhelmsen, 2023)
   * [The Kinetic Gas theory of Mie fluids](https://ntnuopen.ntnu.no/ntnu-xmlui/handle/11250/3029213) (Vegard G. Jervell, 2022)

## Acknowledgments and sources
This implementation of the Revised Enskog solutions is build upon the work presented by M. López de Haro, E. D. G. Cohen, and J. Kincaid in the series of papers *The Enskog Theory for multicomponent mixtures I - IV*, J. Chem. Phys. (1983 - 1987) ([I](https://doi.org/10.1063/1.444985), [II](https://doi.org/10.1063/1.446388), [III](https://doi.org/10.1063/1.446463), [IV](https://doi.org/10.1063/1.452243)).

The implementation utilises the explicit summational expressions for the square bracket integrals published by Tompson, Tipton and Loyalka in *Chapman–Enskog solutions to arbitrary order in Sonine polynomials I - IV* (Physica A, E. J. Mech. - B) 2007-2009 ([I](https://doi.org/10.1016/j.physa.2006.12.001), [II](https://doi.org/10.1016/j.euromechflu.2008.09.002), [III](https://doi.org/10.1016/j.euromechflu.2008.12.002), [IV](https://doi.org/10.1016/j.euromechflu.2009.05.002)).

The work by T. Lafitte, A. Apostolakou, C. Avendaño, A. Galindo, C. Adjiman, E. Müller and G. Jackson, *Accurate statistical associating fluid theory for chain molecules formed from Mie segments* [J. Chem. Phys. 2013](https://doi.org/10.1063/1.4819786) is also of great importance to this implementation.

## Licence

The KineticGas package is distributed as free software under the MIT licence.



# Getting started: In Python

In addition to this explanation, some examples may be found in the [pyExamples directory](https://github.com/thermotools/KineticGas_private/tree/main/pyExamples).

## Initializing a model

The available models are `HardSphere` - The RET for Hard Spheres, `MieKinGas` - The RET-Mie. They are initialized by passing the appropriate component identifiers to the class constructors.

```Python
from pykingas.HardSphere import HardSphere
from pykingas.MieKinGas import MieKinGas

mie = MieKinGas('CO2,C1') # RET-Mie for CO2/CH4 mixture
hs = HardSphere('AR,KR,XE') # RET-HS for Ar/Kr/He mixture
```

The component identifiers are equivalent to the file names in the `pykingas/fluids` directory, and are consistent with the identifiers used by `ThermoPack`. A list of all available fluids and their identifiers can be found in the [Fluid identifiers](#fluid-identifiers) section.

### Note on pure components

*When doing computations for a single component, two mole fractions must be supplied.*

Internally pure components are treated as binary mixtures of equivalent species, such that a model initialized with e.g. `MieKinGas('H2')` will treat pure hydrogen as a mixture of "Hydrogen with hydrogen". This allows computation of the self-diffusion coefficient through the normal `interdiffusion` method, but carries the caveat mentioned above.

Properties are not dependent on the supplied mole fractions, but it has been found that for numerical stability, the choice `x = [0.5, 0.5` is best.

This may be changed in future versions, such that no mole fraction needs to be supplied when working with pure fluids.

### Specifying parameters

If we wish to pass specific parameters to the models, this is done through various keyword arguments, as
```Python
# Continued 
mie = MieKinGas('LJF,LJF', mole_weights=[5, 10], sigma=[2.5e-10, 3e-10], eps_div_k=[150, 200], la=[6, 7], lr=[12, 13])
```
the `mole_weights` argument sets the molar masses of the components, the `sigma` argument sets the mie-potential $\sigma$-parameters (in m), the `eps_div_k` argument sets the mie-potential $\epsilon$-parameters, the `la` argument sets the attractive exponents ($\lambda_a$), and the `lr` argument sets the repulsive exponents ($\lambda_r$).

Classes will only accept keyword arguments that are relevant to them, i.e.
```Python
hs = HardSphere('LJF,LJF', eps_div_k=[100, 200]) # Throws an error
```
will throw an error.

To specify the parameters for only one component, and use default parameters for another, set the parameter for the components that are to use default values to `None`, as
```Python
# Continued 
mie = MieKinGas('AR,KR', la=[None, 7], lr=[None, 14]) # Uses the default values for Ar, and specified values for Kr
mie = MieKinGas('AR,KR', la=[6, None], lr=[None, 14]) # Uses default la for Kr, and default lr for Ar.
```

For isotopic mixtures, one can specify masses in the same way

```Python
from pykingas.MieKinGas import MieKinGas
mie = MieKinGas('CH4,CH4,CH4,CH4', mole_weights=[16, 17, 18, 19]) # Isotopic mixture of 1-, 2-, 3-, and 4 times deuterised methane
```

### The Equation of State

KineticGas uses an Equation of State (EoS) internally to compute the derivatives of chemical potential with respect to molar density. Additionally, the `tp-inteface` methods for predicting transport coefficients use the EoS to compute molar volume at a given T, p, x. This each models stores its own equation of state in the `self.eos` attribute. By default, this is a `ThermoPack` equation of state object, which can be specified using the `use_eos` kwarg upon initialization, as

```Python
from pykingas.MieKinGas import MieKinGas
from thermopack.cubic import cubic

comps = 'AR,H2O' # The components we wish to model
eos = cubic(comps, 'SRK') # Soave-Redlich-Kwong EoS for Argon-water mixture
mie = MieKinGas(comps, use_eos=eos)
```

This can be useful if the components to be modeled do not have parameters for the default eos (`thermopack.saftvrmie` for `MieKinGas`), or if one wishes to use some other eos. 

In the latter case, the only requirement is that the EoS object implements a method with signature equivalent to `thermopack`'s `chemical_potential_tv`. If the `tp-interface` is to be used, the object must also implement a method with signature equivalent to `thermopack`'s `specific_volume`.

### Properties at infinite dilution

Properties at infinite dilution can be of interest. Note that at infinite dilution, viscosity, thermal conductivity, and the thermal diffusion factor are independent of density, while the diffusion coefficient and thermal diffusion coefficient are inversely proportional to the density. To initialize a model where the species have negligible covolume (i.e. the radial distribution function is uniformly equal to one), set the kwarg `is_idealgas=True`, as

```Python
from pykingas.MieKinGas import MieKinGas
mie = MieKinGas('H2', is_idealgas=True) # Properties of hydrogen at infinite dilution
```

## Making predictions

In addition to the methods here, a Tp-interface exists for the same methods, consisting of the methods `thermal_conductivity_tp`, `viscosity_tp`, `interdiffusion_tp`, `theramal_diffusion_coeff_tp` and `thermal_diffusion_factor_tp`. These methods are only wrappers for ease of use, that use the internal equation of state of the object (`self.eos`) to compute the molar volume at given (T, p, x) (assuming vapour phase), and passes the call to the methods documented here. Those methods have signatures equivalent to these, but with molar volume swapped out for pressure.

Please note that the Enskog solutions are explicit in density (not pressure), such that when making predictions as a function of pressure, an accurate equation of state is required to translate from a (T, V, n) state to a (T, p, n) state.

### Thermal conductivity

Thermal conductivities are predicted with the method `thermal_conductivity(self, T, Vm, x, N=None)`, where `T` is the temperature, `Vm` is the molar volume, `x` is the molar composition and `N` is the Enskog approximation order.

Example:

```Python
from pykingas.MieKinGas import MieKinGas

kin = MieKinGas('O2,N2,CO2,C1') # Mixture of air with carbon dioxide and methane, modeled with RET-Mie
T = 800 # Kelvin
Vm = 66.5 # cubic meter per mole, approximately equivalent to a pressure of 1 bar
x = [0.05, 0.25, 0.5, 0.2] # Molar composition

cond = kin.thermal_conductivity(T, Vm, x, N=2) # Thermal conductivity [W / m K]
```

### Shear viscosity

Shear viscosities are predicted with the method `viscosity(self, T, Vm, x, N=None)`, where `T` is the temperature, `Vm` is the molar volume, `x` is the molar composition and `N` is the Enskog approximation order.

Example:

```Python
from pykingas.MieKinGas import MieKinGas

kin = MieKinGas('O2,N2,CO2,C1') # Mixture of air with carbon dioxide and methane, modeled with RET-Mie
T = 800 # Kelvin
Vm = 66.5 # cubic meter per mole, approximately equivalent to a pressure of 1 bar
x = [0.05, 0.25, 0.5, 0.2] # Molar composition

visc = kin.viscosity(T, Vm, x, N=2) # Shear viscosity [Pa s]
```

### Diffusion coefficients

Diffusion coefficients may be defined in many different ways, and depend upon the frame of reference (FoR). For a more in-depth discussion on this see the supporting information of [Revised Enskog Theory for Mie fluids: Prediction of diffusion coefficients, thermal diffusion coefficients, viscosities and thermal conductivities.]()

The interface to all diffusion coefficients is the method `interdiffusion(self, T, Vm, x, N)`, where `T` is the temperature, `Vm` is the molar volume, `x` is the molar composition and `N` is the Enskog approximation order.

The default definition of the diffusion coefficient is 

$$J_i^{(n)} = - \sum_{i \neq l} D_{ij} \nabla n_j$$

where $J_i$ is the molar flux of species $i$ in the centre of moles (CoN) FoR, and $i \neq l$ are the independent *molar density* gradients. $l$ is by default the last component in the mixture, such that for a binary system, this reduces to

$$J_1 = - D \nabla n_1$$

The common Fickean diffusion coefficient. The diffusion coefficients are then computed as

```python
from pykingas.MieKinGas import MieKinGas

kin = MieKinGas('AR,KR') # RET-Mie for a mixture of argon and krypton
T = 300 # Kelvin
Vm = 25 # cubic meter per mole, approximately equivalent to a pressure of 1 bar
x = [0.3, 0.7] # Molar composition

D = kin.interdiffusion(T, Vm, x, N=2) # Binary diffusion coefficient [m^2 / s]
```

Note: For binary mixtures, if the kwarg `use_binary=True` and `use_independent=True` (default behaviour), only a single diffusion coefficient is returned (not an array).

#### Variations of the diffusion coefficient

To compute diffusion coefficients in other frames of reference, use the `frame_of_reference` kwarg, the valid options are `'CoN'` (centre of moles, default), `'CoM'` (centre of mass / barycentric), `'CoV'` (centre of volume), and `'solvent'`, in combination with the `solvent_idx` kwarg.

Example:

```Python
from pykingas.MieKinGas import MieKinGas

kin = MieKinGas('AR,KR') # RET-Mie for a mixture of argon and krypton
T = 300 # Kelvin
Vm = 25 # cubic meter per mole, approximately equivalent to a pressure of 1 bar
x = [0.3, 0.7] # Molar composition

D_CoN = kin.interdiffusion(T, Vm, x, N=2, frame_of_reference='CoN') # Diffusion coefficient in the CoN FoR
D_CoM = kin.interdiffusion(T, Vm, x, N=2, frame_of_reference='CoM') # Diffusion coefficient in the CoM FoR (barycentric)
D_CoV = kin.interdiffusion(T, Vm, x, N=2, frame_of_reference='CoV') # Diffusion coefficient in the CoV FoR
D_solv_Ar = kin.interdiffusion(T, Vm, x, N=2, frame_of_reference='solvent', solvent_idx=0) # Diffusion coefficient in the solvent FoR, with Argon as the solvent
D_solv_Kr = kin.interdiffusion(T, Vm, x, N=2, frame_of_reference='solvent', solvent_idx=1) # Diffusion coefficient in the solvent FoR, with Krypton as the solvent
```

When using the `solvent` FoR, the *dependent* molar density gradient is by default set to be the solvent.

To explicitly set the dependent molar density gradient (default is the last component), use the `dependent_idx` kwarg, as

```python
# Continued

D_1 = kin.interdiffusion(T, Vm, x, N=2, dependent_idx=0) # Diffusion coefficeint in the CoN FoR, with \nabla n_{Ar} as the dependent gradient
D_2 = kin.interdiffusion(T, Vm, x, N=2, dependent_idx=1) # Diffusion coefficeint in the CoN FoR, with \nabla n_{Kr} as the dependent gradient
```

The `dependent_idx`, the specifies the value of $l$ in the equation

$$J_i^{(n)} = - \sum_{i \neq l} D_{ij} \nabla n_j$$

defining the diffusion coefficient. The two diffusion coefficients computed above would thus correspond to the diffusion coefficients

$$J_1^{(n)} = D_1 \nabla n_2 $$
$$J_2^{(n)} = D_1 \nabla n_2 $$

and

$$J_1^{(n)} = D_2 \nabla n_1 $$
$$J_2^{(n)} = D_2 \nabla n_1 $$

where the superscript $^(n)$ denotes that the fluxes are in the centre of moles frame of reference.

To compute diffusion coefficients corresponding to a dependent set of fluxes and forces, defined by

$$J_i^{(FoR)} = - \sum_j D_{ij} \nabla n_j,$$

set the kwarg `use_independent=False`, as

```Python
# Continued

D = kin.interdiffusion(T, Vm, x, N=2, use_independent=False) # Dependent diffusion coefficients in the CoN FoR
```

For the current system this corresponds to the coefficients of the equation

$$J_1^{(n)} = - D[0, 0] \nabla n_1 - D[0, 1] \nabla n_2$$

and 

$$J_2^{(n)} = - D[1, 0] \nabla n_1 - D[1, 1] \nabla n_2.$$

where `D[i, j]` are the elements of the matrix returned by `kin.interdiffusion(T, Vm, x, N=2, use_independent=False)`.

The `frame_of_reference` kwarg works as normal when `use_independet=False`.

### Thermal diffusion

Thermal diffusion is characterised by several common coefficients, the thermal diffusion coefficients $D_{T,i}^{(FoR)}$, the thermal diffusion factor $\alpha_{ij}$, the thermal diffusion ratios $k_{T, i}$ and the Soret coefficients $S_{T,i}$.

Of these, the thermal diffusion coefficients, $D_{T,i}^{(FoR)}$, carry the same ambiguity as the diffusion coefficients in their dependency on the frame of reference (FoR) and choice of dependent gradient.

#### The Thermal diffusion factors

The thermal diffusion factor gives the ratio

$$\nabla \ln (x_i / x_j) = - \alpha_{ij} \nabla \ln T$$

in the absence of mass fluxes, and can be directly related to the Onsager phenomenological coefficients. They are computed as

```Python
from pykingas.MieKinGas import MieKinGas

kin = MieKinGas('C1,C3,CO2') # RET-Mie for a mixture of methane, propane and CO2
T = 300 # Kelvin
Vm = 25 # cubic meter per mole, approximately equivalent to a pressure of 1 bar
x = [0.3, 0.6, 0.1] # Molar composition

alpha = kin.thermal_diffusion_factor(T, Vm, x, N=2) # Thermal diffusion factors [dimensionless]
```

#### The thermal diffusion ratios

The thermal diffusion ratios satisfy the relation

$$\nabla n_i = - k_{T,i} \nabla \ln T$$

in the absence of mass fluxes, and can be directly related to the Onsager phenomenological coefficients. They are computed as

```Python
# Continued 
kT = kin.thermal_diffusion_ratio(T, Vm, x, N=2) # Thermal diffusion ratios [dimensionless]
```

#### The thermal diffusion coefficients

The thermal diffusion coefficients are by default defined by

$$J_i^{(n)} = D_{T, i} \nabla \ln T - \sum_{j \neq l} D_{ij} \nabla n_j,$$

where $J_i^{(n)}$ is the molar flux of species $i$ in the centre of moles (CoN) FoR, $\nabla n_j$ is the *molar density* gradient of component $j$, and $l$ is the index of the dependent gradient. This is computed by

```Python
from pykingas.MieKinGas import MieKinGas

kin = MieKinGas('C1,O2,CO2') # RET-Mie for a mixture of methane, oxygen and CO2
T = 300 # Kelvin
Vm = 25 # cubic meter per mole, approximately equivalent to a pressure of 1 bar
x = [0.3, 0.6, 0.1] # Molar composition

DT = kin.thermal_diffusion_coeff(T, Vm, x, N=2) # Thermal diffusion coefficients in the CoN FoR [mol / m s]
```

#### Variations of the thermal diffusion coefficients
For other frames of reference, use the `frame_of_reference` kwarg, with options equivalent to those for `interdiffusion`, that is: `'CoN'` (centre of moles, default), `'CoM'` (centre of mass / barycentric), `'CoV'` (centre of volume), and `'solvent'`, in combination with the `solvent_idx` kwarg.

Example:

```Python
# Continued
DT_CoN = kin.thermal_diffusion_coeff(T, Vm, x, N=2, frame_of_reference='CoN') # Thermal diffusion coefficient in the CoN FoR
DT_CoM = kin.thermal_diffusion_coeff(T, Vm, x, N=2, frame_of_reference='CoM') # Thermal diffusion coefficient in the CoM FoR (barycentric)
DT_CoV = kin.thermal_diffusion_coeff(T, Vm, x, N=2, frame_of_reference='CoV') # Thermal diffusion coefficient in the CoV FoR
DT_solv_C1 = kin.thermal_diffusion_coeff(T, Vm, x, N=2, frame_of_reference='solvent', solvent_idx=0) # Thermal diffusion coefficient in the solvent FoR, with methane as the solvent
DT_solv_C3 = kin.thermal_diffusion_coeff(T, Vm, x, N=2, frame_of_reference='solvent', solvent_idx=1) # Thermal diffusion coefficient in the solvent FoR, with propane as the solvent
DT_solv_CO2 = kin.thermal_diffusion_coeff(T, Vm, x, N=2, frame_of_reference='solvent', solvent_idx=2) # Thermal diffusion coefficient in the solvent FoR, with CO2 as the solvent
```

To explicitly select the dependent molar gradient (default is the last component), use the `dependent_idx` kwarg, equivalently to `interdiffusion`. 

Example:

```Python
# Continued
DT = kin.thermal_diffusion_coeff(T, Vm, x, N=2, dependent_idx=0) # Thermal diffusion coefficient in the CoN FoR, with \nabla n_{C1} as the dependent gradient
D = kin.interdiffusion(T, Vm, x, N=2, dependent_idx=0) # Diffusion coefficient in the CoN FoR with \nabla n_{C1} as the dependent gradient
```

This gives the coefficients corresponding to the flux equations

$$J_{C1} = D_{T}[0] \nabla \ln T - D[0, 1] \nabla n_{O2} - D[0, 2] \nabla n_{CO2}, $$

$$J_{O2} = D_{T}[1] \nabla \ln T - D[1, 1] \nabla n_{O2} - D[1, 2] \nabla n_{CO2}, $$

$$J_{CO2} = D_{T}[2] \nabla \ln T - D[2, 1] \nabla n_{O2} - D[2, 2] \nabla n_{CO2}. $$

To compute coefficients corresponding to flux equation with all forces and fluxes (not an independent set), set the kwarg `use_independent=False`, as

```Python
# Continued
DT = kin.thermal_diffusion_coeff(T, Vm, x, N=2, use_independent=False) # Thermal diffusion coefficient in the CoN FoR, with all gradients
D = kin.interdiffusion(T, Vm, x, N=2, use_independent=False) # Diffusion coefficient in the CoN FoR with all gradients
```

This gives the coefficients corresponding to the flux equations

$$J_{C1} = D_{T}[0] \nabla \ln T - D[0, 0] \nabla n_{C1} - D[0, 1] \nabla n_{O2} - D[0, 2] \nabla n_{CO2}, $$

$$J_{O2} = D_{T}[1] \nabla \ln T - D[1, 0] \nabla n_{C1} - D[1, 1] \nabla n_{O2} - D[1, 2] \nabla n_{CO2}, $$

$$J_{CO2} = D_{T}[2] \nabla \ln T - D[2, 0] \nabla n_{C1} - D[2, 1] \nabla n_{O2} - D[2, 2] \nabla n_{CO2}. $$

The `frame_of_reference` kwarg works as normal when setting `use_independent=False`.

# Fluid identifiers

| Fluid name | Fluid identifyer | CAS |
| ---------- | ---------------- | --- |
| Argon | AR | 7440-37-1 |
| Methane | C1 | 74-82-8 |
| Ethane | C2 | 74-84-0 |
| Propane | C3 | 74-98-6 |
| Carbon dioxide | CO2 | 124-38-9 |
| Deuterium | D2 | 7782-39-0 |
| Hydrogen | H2 | 1333-74-0 |
| Water | H2O | 7732-18-5 |
| Helium-4 | HE | 7440-59-7 |
| Krypton | KR | 7439-90-9 |
| Lennard-jones_fluid | LJF |  |
| Nitrogen | N2 | 7727-37-9 |
| N-decane | NC10 | 124-18-5 |
| N-pentadecane | NC15 | 629-62-9 |
| N-eicosane | NC20 | 112-95-8 |
| N-docosane | NC22 | 629-97-0 |
| N-butane | NC4 | 106-97-8 |
| N-pentan | NC5 | 109-66-0 |
| N-hexane | NC6 | 110-54-3 |
| N-heptane | NC7 | 142-82-5 |
| N-octane | NC8 | 111-65-9 |
| N-nonane | NC9 | 111-84-2 |
| Neon | NE | 7440-01-9 |
| Ortho-hydrogen | O-H2 | 1333-74-0 |
| Oxygen | O2 | 7782-44-7 |
| Para-hydrogen | P-H2 | 1333-74-0 |
| Xenon | XE | 7440-63-3 |


