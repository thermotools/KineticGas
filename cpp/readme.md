# Kineticgas module
Contains the abstract `KineticGas` class. This class defines functions that solve the neccesary equations to produce 
the vectors and matrices required to determine diffusion coefficients and conductivities of dilute gas mixtures for an 
arbitrary potential.

The abstract class `Spherical` contains various utility methods that are used for evaluating the collision integrals,
transfer lengths, etc. for arbitrary spherical potentials.

Model classes such as `MieKinGas` implement a pair potential and, if they are to be used at non-infinite dilution,
the function `model_rdf`, which returns the RDF at contact. Note that specifically `MieKinGas` also overrides the 
collision integrals (`omega`) in `Spherical` to use a correlation when applicable.

* Various "multiparameter" potentials (typically used for fitting to ab initio data) are found in `multiparam.h`
* Feynman-Hibbs corrections: See `extensions.h` (the `FH_Corrected` template class)
* `Quantum` is an abstract class implementing quantum-mechanical calculations
  * Any class implementing a spherical potential can inherit from `Quantum` instead of `Spherical` to get access to this stuff

The `KineticGas` implementation is split up a bit:
* KineticGas.cpp: "base-functionality", which basically consists of obtaining the Sonine polynomial expansion coefficients
* KineticGas_properties.cpp: The actual transport-property methods
* KineticGas_mthr.cpp: Multithreading-related stuff (pre-calculation of expensive stuff)

The `Spherical` class is also split in two files:
* Spherical.cpp: All the base-algorithms for working with spherical potentials
* transfer_lengths.cpp: Transfer length models

We also have `eos_interface.h`, which contains some abstract templated stuff that lets us cleanly wrap an equation of state implemented some other place.
This is used when wrapping thermopack, and also if we use an equation of state implemented in python.

# Utility modules

## Global parameters and utility structs

The header `global_params.h` contains natural constants.

The header `utils.h` contains various `enum`s and `struct`s that are used throughout the code base, as well as the `[get/set]_fluid_dir` functions to 
manage the fetching of fluid files from the fluid database.

## Factorial module
Defines the types `Frac`, `Product` and `Fac`, representing exact fractions, products and factorials respectively. This
is required to avoid overflow and float truncation issues when evaluating fractions containing factorials that largely
cancel, but where the numerator and denominator cannot be evaluated separately. This is specifically the case in the
`A_*` and `B_*` methods of `KineticGas`.

In short, a factorial (`Fac`) is treated as a single integer (upper value) until the `.eval()` method is called. Similarly, 
a `Product` is treated as an array of integers and a singular double. A fraction consists of two `Product`s, the 
numerator and the denominator. The `*` and `/` operators are overloaded by letting `Fac` and `Product` be implicitly 
converted to `Frac` (i.e. a fraction with a denominator of `Product(1)`). Thus you get the slightly unintuitive 
arithmatic of 

`Product * Product => Frac`

`Fac * Fac => Frac`

`Product * Fac => Frac`

This is quite simply done to greatly reduce the number of neccesary operator overloads.

When the `Frac::eval()` method is called, integers in the numerator and denominator that cancel are set to unity 
before the two are evaluated. This effectivly allows the evaluation of fractions such as `Fac(32) / Fac(28)` without 
getting overflow issues that would appear without this module.

The module also contains the function `Product ipow(int base, int power)`, which evaluates a power as a `Product`,
such that we can do stuff like `(ipow(10, 500) / ipow(10, 499)).eval() == 10`, without ever worrying about overflow. 

## Integration module
Contains various numerical algorithms. The only non-trivial of which may be the function `integrate2d`.
and corresponding `mesh2d`, which are used to compute collision integrals.
Mathematical description of the algorithm used by `integrate2d` function can be found in
[The Kinetic Gas theory of Mie fluids](https://ntnuopen.ntnu.no/ntnu-xmlui/handle/11250/3029213).