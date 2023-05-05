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
 

Stuff is illustrated here as well:

![](https://github.com/thermotools/KineticGas/blob/main/docs/structure/kineticgas_classes.pdf)

![](https://github.com/thermotools/KineticGas/blob/main/docs/structure/who_does_what.pdf)

# File system

`cpp/` : The C++ source code and headers for `KineticGas`

`cpp/Integration/` : The C++ source code and headers for the integration module used to evaluate the collision integrals.

`pyExamples` : Example files for doing computations

`pykingas/` : Python source code for the package

`pykingas/tests/` : Tests that are run after compiling

`pykingas/fluids/` : Fluid parameter database

`Dockerfiles/` : (Not in use, should be made up to date)

`docs/` : Documentation