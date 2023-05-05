# Factorial module
Defines the types `Frac`, `Product` and `Fac`, representing exact fractions, products and factorials respectively. This is required to avoid overflow and float truncation issues when evaluating fractions containing factorials that largely cancel, but where the numerator and denominator cannot be evaluated separately. This is specifically the case in the `A_*` and `B_*` methods of `KineticGas`.

In short, a factorial (`Fac`) is treated as an array of integers up until the `.eval()` method is called. Similarly, a `Product` is treated as an array of integers and a singular double. A fraction consists of two `Product`s, the numerator and the denominator. The `*` and `/` operators are overloaded by letting `Fac` and `Product` be implicitly converted to `Frac` (i.e. a fraction with a denominator of `Product(1)`). Thus you get the slightly unintuitive arithmatic of 

`Product * Product => Frac`

`Fac * Fac => Frac`

`Product * Fac => Frac`

This is quite simply done to greatly reduce the number of neccesary operator overloads.

When the `Frac::eval()` method is called, integers in the numerator and denominator that cancel are set to unity before the two are evaluated. This effectivly allows the evaluation of fractions such as `Fac(32) / Fac(28)` without getting overflow issues that would appear without this module.

# Kineticgas module
Contains the abstract `KineticGas` class. This class defines functions that solve the neccesary equations to produce the vectors and matrices required to determine diffusion coefficients and conductivities of dilute gas mixtures for an arbitrary potential.

`KineticGas` is inherited by `HardSphere` and the abstract class `Spherical`. `HardSphere` implements the analytic collision integrals for hard sphere mixtures. `Spherical` implements the solvers neccesary to evaluate the collision integrals for an arbitrary spherical potential, and is in turn inherited by `MieKinGas` and `PseudoHardSphere`, which implement the respective potential models with first and second derivatives. 

The model classes (`HardSphere`, `PseudoHardSphere` and `MieKinGas`) are initialized for a given binary mixture, by supplying the particle/potential parameters. The parameter database (`pykingas/XX.json`) is accessible through the python wrapper. Matrix inversion and evaluation of the transport coefficients is done on the python-side, because I didn't want to spend time learning to use `LAPACK`.

# Integration module
Functional module containing the primary functions `integrate2d()` and `mesh2d()`, used for adaptive refinement of an integration mesh. Module is used by `Spherical` to evaluate collision integrals. Mathematical description of the module can be found in [theory.pdf](https://github.com/vegardjervell/Kineticgas/blob/main/theory.pdf).