# Kineticgas module
Contains the abstract `KineticGas` class. This class defines functions that solve the neccesary equations to produce 
the vectors and matrices required to determine diffusion coefficients and conductivities of dilute gas mixtures for an 
arbitrary potential.

The abstract class `Spherical` contains various utility methods that are used for evaluating the collision integrals,
collision diameters, etc. for arbitrary spherical potentials.

Model classes such as `MieKinGas` implement a pair potential and, if they are to be used at non-infinite dilution,
the function `model_rdf`, which returns the RDF at contact. Note that specifically `MieKinGas` also overrides the 
collision integrals (`omega`) in `Spherical` to use a correlation when applicable.

# Utility modules

## Factorial module
Defines the types `Frac`, `Product` and `Fac`, representing exact fractions, products and factorials respectively. This
is required to avoid overflow and float truncation issues when evaluating fractions containing factorials that largely
cancel, but where the numerator and denominator cannot be evaluated separately. This is specifically the case in the
`A_*` and `B_*` methods of `KineticGas`.

In short, a factorial (`Fac`) is treated as an array of integers up until the `.eval()` method is called. Similarly, 
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

## Integration module
Functional module containing the primary functions `integrate2d()` and `mesh2d()`, used for adaptive refinement of an 
integration mesh. Module is used by `Spherical` to evaluate collision integrals. Mathematical description of the 
module can be found in [The Kinetic Gas theory of Mie fluids](https://ntnuopen.ntnu.no/ntnu-xmlui/handle/11250/3029213).

The purpose of the integration algorithm is, in short, to evaluate double integrals in which the integrand has low 
curvatures in some regions, and high curvatures in others, with a minimal number of function evaluations.