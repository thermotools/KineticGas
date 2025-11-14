<!--- 
Generated at: 2025-09-12T17:43:51.253469
This is an auto-generated file, generated using the script at KineticGas/docs/join_docs.py
The file is created by joining the contents of the files
    KineticGas/docs/markdown/
        readme_parts/header.md
        readme_parts/toc_github.md
        metapages/cite_acknowl_licence.md
        vCurrent/source_build.md
        vCurrent/getting_started_py.md
        vCurrent/getting_started_cpp.md
        vCurrent/advanced.md
        vCurrent/structure.md
        vCurrent/fluid_identifiers.md
--->
# KineticGas

KineticGas is an implementation of Revised Enskog Theory (RET) for spherical potentials. The most notable of which is the implementation of RET-Mie, the Revised Enskog Theory for Mie fluids. 

The package is implemented mostly in C++ to handle the numerical computations involved in evaluating the collision integrals and the radial distribution function at contact for the target fluids, with the possibility of setting up multithreading at compile time.

KineticGas can be used to predict diffusion coefficients, thermal diffusion coefficients, viscosities and thermal conductivities in gas mixtures, and is reliable over a large range of temperatures and pressures. The package also contains an extensive database of fluid parameters collected from the open literature.

# [KineticGas homepage](https://thermotools.github.io/KineticGas)
The full documentation, with installation- and getting started-guides can be found on the [KineticGas homepage](https://thermotools.github.io/KineticGas).
This readme is only intended to provide a minimal introduction, and may be out-of-sync with the `pykingas` version currently
on `PyPI`.

![](https://thermotools.github.io/KineticGas/v2.0.0/graphics/all.gif?raw=true)


## Table of contents
   * [Installing KineticGas](#Installing-KineticGas)
   * [Getting started](#Getting-started:-In-Python)
     * [Python](#getting-started-in-python)
       * [Initializing a model](#Initializing-a-model)
       * [Making predictions](#Making-predictions)
     * [C++](#getting-started-in-c)
   * [Advanced usage](#advanced-usage)
     * [Modifying and adding fluids](#modifying-and-adding-fluids)
     * [Implementing new potentials](#implementing-new-potentials)
       * [The C++ side](#implementing-the-c-side)
       * [The Python side](#implementing-the-python-side)
   * [Program structure](#structure)
   * [File system](#file-system)
   * [Fluid indentifiers](#fluid-identifiers)

# Please Cite

KineticGas has been developed throughout several works. If you are referencing the package, please cite one or more of the associated works

* General usage
   * [Revised Enskog theory for Mie fluids: Prediction of diffusion coefficients, thermal diffusion coefficients, viscosities and thermal conductivities](https://doi.org/10.1063/5.0149865) (Vegard G. Jervell and Øivind Wilhelmsen, J. Chem. Phys. 2023)
   * [Predicting viscosities and thermal conductivities from dilute gas to dense liquid: Deriving fundamental transfer lengths for momentum and energy exchange in revised Enskog theory](https://pubs.aip.org/aip/jcp/article/161/23/234106/3325824/Predicting-viscosities-and-thermal-conductivities) (V. G. Jervell and Ø. Wilhelmsen, J. Chem. Phys. 2024)
   * [The Kinetic Gas theory of Mie fluids](https://ntnuopen.ntnu.no/ntnu-xmlui/handle/11250/3029213) (Vegard G. Jervell, 2022)
* Connection to Non-Equilibrium thermodynamics (Onsager coefficients)
   * [The influence of thermal diffusion on water migration through a porous insulation material](https://doi.org/10.1016/j.ijheatmasstransfer.2024.125576) (V. G. Jervell, M. Aa. Gjennestad, T. T. Trinh, Ø. Wilhelmsen, Int. J. Heat Mass Transfer, 2024)
* Transfer Lengths, the EWCA model
  * [Predicting viscosities and thermal conductivities from dilute gas to dense liquid: Deriving fundamental transfer lengths for momentum and energy exchange in revised Enskog theory](https://pubs.aip.org/aip/jcp/article/161/23/234106/3325824/Predicting-viscosities-and-thermal-conductivities) (V. G. Jervell and Ø. Wilhelmsen, J. Chem. Phys. 2024)
* The Lennard-Jones spline fluid
  * [Viscosity, thermal conductivity and self-diffusion coefficient of the Lennard Jones spline fluid: Evaluation of theories for a short-ranged potential](https://doi.org/10.1016/j.fluid.2025.114584) (J. S. Løken, V. G. Jervell, M. Hammer, B. Hafskjold, T. T. Trinh and Ø. Wilhelmsen, 2025)
  * [Revised Enskog theory and extended corresponding states models for the transport properties of the Lennard-Jones/spline fluid](https://ntnuopen.ntnu.no/ntnu-xmlui/handle/11250/3182899) (Johannes S. Løken, 2025)

Cite this repository as

```latex
@article{kineticgas_repo,
  title={{ThermoTools: KineticGas}},
  author={Vegard Gjeldvik Jervell and Johannes Salomonsen L{\o}ken},
  year={2025},
  howpublished={github.com/thermotools/kineticgas}
}
```

## Acknowledgments and sources
This implementation of the Revised Enskog solutions is build upon the work presented by M. López de Haro, E. D. G. Cohen, and J. Kincaid in the series of papers *The Enskog Theory for multicomponent mixtures I - IV*, J. Chem. Phys. (1983 - 1987) ([I](https://doi.org/10.1063/1.444985), [II](https://doi.org/10.1063/1.446388), [III](https://doi.org/10.1063/1.446463), [IV](https://doi.org/10.1063/1.452243)).

The implementation utilises the explicit summational expressions for the square bracket integrals published by Tompson, Tipton and Loyalka in *Chapman–Enskog solutions to arbitrary order in Sonine polynomials I - IV* (Physica A, E. J. Mech. - B) 2007-2009 ([I](https://doi.org/10.1016/j.physa.2006.12.001), [II](https://doi.org/10.1016/j.euromechflu.2008.09.002), [III](https://doi.org/10.1016/j.euromechflu.2008.12.002), [IV](https://doi.org/10.1016/j.euromechflu.2009.05.002)).

The work by T. Lafitte, A. Apostolakou, C. Avendaño, A. Galindo, C. Adjiman, E. Müller and G. Jackson, *Accurate statistical associating fluid theory for chain molecules formed from Mie segments* [J. Chem. Phys. 2013](https://doi.org/10.1063/1.4819786) is also of great importance to this implementation.

## Licence

The KineticGas package is distributed as free software under the MIT licence.

# Installing KineticGas

KineticGas is available on PyPi as the [`pykingas`](https://pypi.org/project/pykingas/) package, for python versions 3.8-3.11, compiled for MacOS running on Apple Silicon, Linux and Windows.

In addition, wheels versions of `KineticGas > 2.0.0` for macOS, Linux and Windows, as well as wheels for the latest version on GitHub can be downloaded [here](https://github.com/thermotools/KineticGas/releases). Instructions for installing with `pip` directly from a downloaded wheel are provided at the linked page.

For MacOS running on Intel, or other operating systems, KineticGas must currently be built from source or installed from one of the distributed wheels linked above.

A KineticGas C++ library is available, and can be built using `cmake` and `make`.

- [Python - `pykingas`](#python---pykingas)
  - [Dependencies](#dependencies)
  - [Building from source](#building-from-source)
    - [First Try](#first-try)
      - [Short explanation](#short-explanation)
    - [When something goes wrong](#when-something-goes-wrong)
- [C++](#c)
  - [Building and installing](#building-and-installing)
  - [Fluid file search path](#fluid-file-search-path)
  - [Linking to the KineticGas library](#linking-to-the-kineticgas-library)

# Python - `pykingas`

## Dependencies

The Python package dependencies are listed in the `pyproject.toml` file in the root directory of the package.

To compile the binary that is called from the python wrapper, [pybind11](https://pybind11.readthedocs.io/en/stable/) is required. `pybind11` is included in `cpp/external` as a git submodule, so cloning the `KineticGas` repository should provide you with the files you need.

A standalone C++ module, that works without the python wrapper is currently under development. See the branch `pure_cpp/` for the most up-to-date version there.


## Building from source

Python wheels for the latest version of KineticGas on `main` are built for macOS and Windows using `cibuildwheels`, and distributed [here](https://github.com/thermotools/KineticGas/releases/tag/Latest-beta).

A build system using `cmake` and `make` is set up to support Mac, Linux and Windows.

### First Try
If all goes well, running

```
git clone https://github.com/thermotools/KineticGas.git
cd KineticGas
git submodule update --init --recursive
mkdir build
cd build
cmake ..
make install
pip install ..
```

make sure to activate a virtual environment first if you want to avoid doing system-level installs.

#### Short explanation

The dynamic library `libpykingas` will be built and installed to the `pykingas` directory, additionally, the `fluids` directory containing the fluid parameter database is copied into the `pykingas` directory.

### When something goes wrong

*Note:* The build system has been changed relatively recently, and is less tested than the build system that was used in the `2.0.0` release. If you encounter issues, please don't hesitate to post an issue on github so that we can improve robustness.

* Warning that thermopack is not installed
  * The easiest way to obtain the `ThermoPack` dynamic library (which `KineticGas` needs) is likely to download the appropriate zip file [here](https://github.com/thermotools/thermopack/releases), unzip it, and set the environment variable `THERMOPACK_DIR` to the resulting directory (where `thermopack-config.cmake` is located).
    * On Linux and macOS: `export THERMOPACK_DIR=/path/to/thermopack-<system>/`
    * On Windows (powershell): `$THERMOPACK_DIR = C:\path\to\thermopack-<system>\thermopack-<system>`
    * To check that it is set correctly: `ls ${THERMOPACK_DIR}` should give a list of files including `thermopack-config.cmake`.
  * The `KineticGas` library has a dependency on the `ThermoPack` C++ wrapper. If you have not installed thermopack, the build system will generate a target from the `thermopack` submodule. Running `make install` should build and install this target, re-running `cmake ..` after building and installing `thermopack` should then give output telling you that `thermopack` has been found and is installed.
  * If you have installed thermopack, run `export THERMOPACK_DIR=<path/to/thermopack>`, to help `cmake` find your installation.


# C++

The KineticGas C++ library is built using `cmake` and `make`. All dependencies are included as git submodules under `cpp/external`, and should be properly retrieved when you clone the `KineticGas` repository and run `git submodule update --init --recursive`. 

*Note*: `KineticGas` depends on [`ThermoPack`](https://thermotools.github.io/thermopack/). If an installation of `ThermoPack` is not found, the build system will attempt to compile it as part of the build process. If you already have an installation of `ThermoPack`, setting the environment variable `THERMOPACK_DIR` to the root directory of `ThermoPack` (where `thermopack-config.cmake` is found), that installation of `ThermoPack` will be used istead of re-compiling. You can also download a binary distribution of ThermoPack at the [ThermoPack repository](https://github.com/thermotools/thermopack/releases).

## Building and installing

If all goes well, you should be able to build the `KineticGas` C++ library by running

```bash
git clone https://github.com/thermotools/KineticGas.git
cd KineticGas
git submodule update --init --recursive
mkdir build
cd build
cmake -Dpurecpp=ON -Dpylib=OFF ..
make install
```

This will provide you with the `lib/libkineticgas.[so/dylib/dll]` dynamic library, and the minimal example program `build/run_kineticgas`, which is built from the source file at `cpp/run_kineticgas.cpp`.

## Fluid file search path

By default, `KineticGas` will search for fluid files at the relative path `../fluids` (relative to the location of the `libkineticgas` dynamic library). This default search path can be changed by building with 
```bash
cmake -DFLUID_DIR=<path/to/fluids> ..
```
where supplying a relative path will result in the library searching for fluid files in the path relative to it's location (`KineticGas/lib`). Supplying absolute paths is also supported. To check
or change where your compiled `KineticGas` library is searching for fluid files, use the `[get/set]_fluid_dir` functions with signatures
```C++
// In utils.h
void set_fluid_dir(const std::string path); // supports both absolute and relative paths (relative to dynamic library location).
std::string get_fluid_dir(); // Current search path for fluid files
```

## Linking to the KineticGas library

An example program with a `CMakeLists.txt` demonstrating how you can include and link to the `KineticGas` library once it is installed is found in `KineticGas/cppExamples`. 

In short terms: Setting the environment variable `KINETICGAS_DIR` to the top-level directory of the KineticGas package (where `kineticgas-config.cmake` is found), should allow `cmake` to find the library using `find_library(KINETICGAS)`. Some convenience variables are set once the library is found:

* `KINETICGAS_ROOT` : Path to root directory of the package
* `KINETICGAS_INSTALLED` : `TRUE` if the dynamic library is found in the correct install location, `FALSE` otherwise
* `KINETICGAS_LIB` : Path to the `libkineticgas` dynamic library
* `KINETICGAS_INCLUDE` : List of include paths needed to include the kineticgas headers and dependencies
* `kineticgas` : Exported target, linking to this target should automatically add the appropriate directories to your include path. 

# Getting started - In Python

In addition to this explanation, some examples may be found in the [pyExamples directory](https://github.com/thermotools/KineticGas_private/tree/main/pyExamples).

## Initialising a model

The available models are `HardSphere` - The RET for Hard Spheres, `MieKinGas` - The RET-Mie. They are initialised by passing the appropriate component identifiers to the class constructors.

```Python
from pykingas.HardSphere import HardSphere
from pykingas.MieKinGas import MieKinGas

mie = MieKinGas('CO2,C1') # RET-Mie for CO2/CH4 mixture
hs = HardSphere('AR,KR,XE') # RET-HS for Ar/Kr/He mixture
```

The component identifiers are equivalent to the file names in the `pykingas/fluids` directory, and are consistent with the identifiers used by `ThermoPack`. A list of all available fluids and their identifiers can be found in the [Fluid identifiers](https://thermotools.github.io/KineticGas/vcurrent/fluid_identifiers.html) section.

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

## Working in reduced units

When working in reduced (Lennard-Jones) units

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
Vm = 0.0665 # cubic meter per mole, approximately equivalent to a pressure of 1 bar
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
Vm = 0.0665 # cubic meter per mole, approximately equivalent to a pressure of 1 bar
x = [0.05, 0.25, 0.5, 0.2] # Molar composition

visc = kin.viscosity(T, Vm, x, N=2) # Shear viscosity [Pa s]
```

### Diffusion coefficients

Diffusion coefficients may be defined in many different ways, and depend upon the frame of reference (FoR). For a more in-depth 
discussion on this see the supporting information of [Revised Enskog Theory for Mie fluids: Prediction of diffusion coefficients, 
thermal diffusion coefficients, viscosities and thermal conductivities.](https://pubs.aip.org/aip/jcp/article/158/22/224101/2895227/Revised-Enskog-theory-for-Mie-fluids-Prediction-of) For more details on the definitions available in the KineticGas package, see the [memo on definitions of the diffusion 
coefficient.](/KineticGas/memo/diffusion/diffusion_definitions.pdf)

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
Vm = 0.025 # cubic meter per mole, approximately equivalent to a pressure of 1 bar
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
Vm = 0.025 # cubic meter per mole, approximately equivalent to a pressure of 1 bar
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

Of these, the thermal diffusion coefficients, $D_{T,i}^{(FoR)}$, carry the same ambiguity as the diffusion coefficients 
in their dependency on the frame of reference (FoR) and choice of dependent gradient. For more details on the definitions 
available in the KineticGas package, see the [memo on definitions of the diffusion coefficient.](/KineticGas/memo/diffusion/diffusion_definitions.pdf)

#### The Thermal diffusion factors

The thermal diffusion factor gives the ratio

$$\nabla \ln (x_i / x_j) = - \alpha_{ij} \nabla \ln T$$

in the absence of mass fluxes, and can be directly related to the Onsager phenomenological coefficients. They are computed as

```Python
from pykingas.MieKinGas import MieKinGas

kin = MieKinGas('C1,C3,CO2') # RET-Mie for a mixture of methane, propane and CO2
T = 300 # Kelvin
Vm = 0.025 # cubic meter per mole, approximately equivalent to a pressure of 1 bar
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
Vm = 0.025 # cubic meter per mole, approximately equivalent to a pressure of 1 bar
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

# Getting started - In C++

For instructions on building the `KineticGas` C++ library, see the [installation guide](source_build.html#c).

A basic example showing initialization of a model is found in [`cppExamples/basic.cpp`](), the `cppExamples` directory also contains a `CMakeLists.txt` showing how to obtain the required headers for the `KineticGas` library, as well as link the library to your program.

## Initializing a model

To initialize a model, `include` the appropriate header file, and specify the components to model with a comma separated string, as
```C
#include "MieKinGas.h"

int main(){
    MieKinGas mie("HE,NE"); // Mixture of helium and neon.
}
```
The component identifiers used are equivalent to the file names of the [fluid files](https://thermotools.github.io/KineticGas/vcurrent/https://github.com/thermotools/KineticGas/tree/main/fluids), and are summarised [here](fluid_identifiers.html)

## Computing properties

The interfaces for property calculations are more or less equivalent to those used in Python. The major differences you should be aware of are
* Diffusion coefficients are returned as an `Eigen::MatrixXd`
  * *Note*: The optional `dependent_idx` argument to `interdiffusion` supports python-style negative indexing (i.e. `-1` is the last component).
* Vector properties (e.g. thermal diffusion coefficients) are returned as an `Eigen::VectorXd`
* Frames of reference are specified with the `FrameOfReference` enum, found in `utils.h`. Valid values are
  * `FrameOfReference::CoN` - Center of moles
  * `CoM` - Center of mass (barycentric)
  * `CoV` - Center of volume
  * `solvent` - Solvent, solvent index is the `dependent_idx`, which defaults to the last component.
  * `zarate`, `zarate_x`, and `zarate_w` - See the [memo](/KineticGas/memo/diffusion/diffusion_definitions.pdf)
  * See the python docs and the [memo](/KineticGas/memo/diffusion/diffusion_definitions.pdf) for more details on definitions of the diffusion coefficients.

## Selecting transfer length models

Use the methods
* `void KineticGas::set_transfer_length_model(int model_id)` - Set the transfer length model
* `std::pair<int, std::string> KineticGas::get_transfer_length_model()` - Return the current transfer length model key (`int`) and description (`string`)
* `std::map<int, std::string> KineticGas::get_valid_transfer_length_models()` - Get a map of valid models with descriptions

In addition, the enum `TransferLengthModel` in `utils.h` may be useful if you don't like remembering keys. The enum is used everywhere internally, and it is heavily
recommended to use it instead of manually specifying , in case keys for different models are changed in the future

# Advanced usage

* [Selecting transfer length models](#selecting-transfer-length-models)
* [Modifying and adding fluids](#modifying-and-adding-fluids)
* [Implementing new potentials](#implementing-new-potentials)

## Selecting transfer length models

For the computation of transfer lengths, several models can be used. All classes inheritting from `py_KineticGas` support the `get_valid_tl_models()` method,
which returns a `dict` with key-description pairs indicating the available transfer length models. Use the methods `get_tl_model()` and `set_tl_model(key)`
to see the active transfer length model, and to select another model.

## Modifying and adding fluids

All fluid parameters are accessed via the `.json` files in the `pykingas/fluids` directory. The structure of the files in the `pykingas/fluids` directory is

```json
<fluid_id.json>
{
    "ident": "<fluid identifier (optional)>",
    "formula": "<chemical formula (optional)>",
    "cas_number": "<optional>",
    "name": "<fluid name (optional)>",
    "aliases": [
            "<optional alias 1>",
            "<optional alias 2>"
      ],
      "mol_weight": <molar mass [g / mol]>,
      "<Potential identifier>" : {
        "default" : {
            "<some parameter>" : <value>,
            "<parameter 2" : <value>,
            "<parameter 3>" : <value>,
            etc...
            "bib_reference" : "<link to article or other reference to parameter set>"
        }
        "<alternative parameter set>" : {
            "<some parameter>" : <value>,
            "<parameter 2" : <value>,
            "<parameter 3>" : <value>,
            etc...
            "bib_reference" : "<link to article or other reference to parameter set>"
        }
      }
}
```

The currently supported `"<Potential identifier>"`'s are `"Mie"` (for RET-Mie) and `"HardSphere"` (for Hard sphere). Check the files in `pykingas/fluids` to see what fields are required for the different parameter sets. 

Other than the potential parameters, only the `"mol_weight"` field is strictly required. Filling in the other fields is recommended for consistency with existing code, in case it at some point becomes desirable to use them.

The identifier used for a fluid in `KineticGas` is equivalent to the name of the corresponding `<name>.json` file.

## Managing fluid file search path (only for C++)

By default, `KineticGas` will search for fluid files at the relative path `../fluids` (relative to the location of the `libkineticgas` dynamic library). This default search path can be changed by building with 
```bash
cmake -DFLUID_DIR=<path/to/fluids> ..
```
where supplying a relative path will result in the library searching for fluid files in the path relative to it's location (`KineticGas/lib`). Supplying absolute paths is also supported. To check
or change where your compiled `KineticGas` library is searching for fluid files, use the `[get/set]_fluid_dir` functions with signatures
```C++
// In utils.h
void set_fluid_dir(const std::string path); // supports both absolute and relative paths (relative to dynamic library location).
std::string get_fluid_dir(); // Current search path for fluid files
```

## Implementing new potentials

Functionality making it simple to implement new potentials is at the core of `KineticGas`. Broadly speaking, implementing a new potential consist of four steps: 

* Writing a class that inherits (directly or indirectly) from the `KineticGas` class on the C++ side
* Exposing the C++ class in `cpp/bindings.cpp`
* Writing a "mirror" class on the python side that inherits (directly or indirectly) from the `py_KineticGas` class on the python side.
* Adding appropriate parameter sets to the `pykingas/fluids` files.

### Implementing the C++ side

All classes that inherit from `KineticGas` must implement the methods `omega`, which returns the collision integrals, the method `model_rdf`, which returns the radial distribution function at contact, and the method `get_collision_diameters`, which returns the collision diameters. 

Out of these, the `omega` method is implemented in the  `Spherical` class which instead requires that inheritting classes implement the methods `potential`, `potential_derivative_r` and `potential_dblderivative_rr`, corresponding to the pair potential, and its first and second derivative wrt. distance. 

The options for implementing a new potential are then

 * Inherit `KineticGas`
   * Implement `omega` (The collision integrals)
   * Implement `model_rdf` (The radial distribution function at contact)
   * Implement `get_collision_diameters` (The collision diameters)
 * Inherit `Spherical`
   * Implement `potential` (The pair-potential)
   * Implement `potential_derivative_r` (Derivative of the pair-potential)
   * Implement `potential_dblderivative_rr` (Second derivative of the pair-potential)
   * Implement `model_rdf` (The radial distribution function at contact)
   * Implement `get_collision_diameters` (The collision diameters)

### Implementing the Python side

The Python-side class mirroring a C++ class has two responsibilities: Fetch the appropriate parameters from the `pykingas/fluids/*.json` files, initialize the `self.cpp_kingas` object and initialize the `self.eos` object (typically a `ThermoPack` eos object). The constructor should accept (at least) a string containing the fluid identifiers of a mixture.

The `py_KineticGas` constructor accepts the `comps` argument, which is a string of comma-separated fluid identifiers, fetches the corresponding `.json`-files, and stores them in the `self.fluids` attribute. The inherriting class needs only to call the `py_KineticGas` constructor, retrieve the appropriate parameters, and pass them to the constructor of the corresponding C++ class. A minimal example is:

```Python
class MyNewPotential(py_KineticGas)
    def __init__(self, comps):
        super().__init__(comps) # super() initializes self.mole_weights
        self.fluids = [self.fluids[i]['<paramter identifier>']["default"] for i in range(self.ncomps)]
        self.cpp_kingas = cpp_MyNewPotential(self.mole_weights, self.fluids['param 1'], self.fluids['param 2'], '... etc')
        self.eos = <Some ThermoPack EoS>(comps)
```

# Structure

See the [structure docs](https://github.com/thermotools/KineticGas/blob/main/docs/structure/structure.pdf) for more information.

The primary responsibilities of the python-side and C++ side of the package are

 * Python-side
   * KineticGas parent class
     * Compute transport coefficients using Sonine polynomial expansion coefficients, RDF at contact and collision diameter by C++ model, and thermodynamic factors supplied by ThermoPack model
   * Inheriting classes
     * Read parameters from fluid database
     * Initialize corresponding C++ model
     * Initialize corresponding ThermoPack model

 * C++ Side
   * KineticGas (abstract class)
     * Derived classes implement collision integrals, RDF at contact and collision diameter.
     * Evaluate square bracket integrals, using collision integrals implemented in derived classes
     * Build matrices to compute Sonine polynomial expansion coefficients using square bracket integrals and RDF at contact implemented in derived classes
   * Spherical (abstract class)
     * Numerical solvers for evaluating collision integrals
     * Derived classes must implement interaction potential with first and second derivatives.
   * MieKinGas (concrete class)
     * Implements interaction potential - such that collision integrals can be evaluated by methods in Spherical
     * Implements RDF at contact
     * Implements collision diameter
 
# File system

`cpp/` : The C++ source code and headers for `KineticGas`

`cpp/Integration/` : The C++ source code and headers for the integration module used to evaluate the collision integrals.

`pyExamples/` : Example files for doing computations in python

`cppExamples/`: Example files for C++

`pykingas/` : Python source code for the `pykingas` package

`pykingas/tests/` : Python-side test suite

`pykingas/fluids/` : Fluid parameter database

`Dockerfiles/` : (Not in use, should be made up to date)

`docs/` : Documentation

# Fluid identifiers

*Note* : Many of these fluid parameters have been pulled directly from the [ThermoPack](https://github.com/thermotools/thermopack) fluid database for SAFT-VR Mie parameters. In the cases where SAFT-VR Mie uses segment numbers $>1$ to describe the fluids, the parameter sets cannot be expected to be suitable for use with RET-Mie.

| Fluid name | Fluid identifier | CAS |
| ---------- |------------------| --- |
| Argon | AR               | 7440-37-1 |
| Methane | C1               | 74-82-8 |
| Ethane | C2               | 74-84-0 |
| Propane | C3               | 74-98-6 |
| Carbon dioxide | CO2              | 124-38-9 |
| Deuterium | D2               | 7782-39-0 |
| Hydrogen | H2               | 1333-74-0 |
| Water | H2O              | 7732-18-5 |
| Helium-4 | HE               | 7440-59-7 |
| Krypton | KR               | 7439-90-9 |
| Lennard-jones_fluid | LJF              |  |
| Nitrogen | N2               | 7727-37-9 |
| N-decane | NC10             | 124-18-5 |
| N-pentadecane | NC15             | 629-62-9 |
| N-eicosane | NC20             | 112-95-8 |
| N-docosane | NC22             | 629-97-0 |
| N-butane | NC4              | 106-97-8 |
| N-pentan | NC5              | 109-66-0 |
| N-hexane | NC6              | 110-54-3 |
| N-heptane | NC7              | 142-82-5 |
| N-octane | NC8              | 111-65-9 |
| N-nonane | NC9              | 111-84-2 |
| Neon | NE               | 7440-01-9 |
| Ortho-hydrogen | O-H2             | 1333-74-0 |
| Oxygen | O2               | 7782-44-7 |
| Para-hydrogen | P-H2             | 1333-74-0 |
| Xenon | XE               | 7440-63-3 |


