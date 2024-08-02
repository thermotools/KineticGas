# KineticGas

KineticGas is an implementation of Revised Enskog Theory (RET) for spherical potentials. The most notable of which is the implementation of RET-Mie, the Revised Enskog Theory for Mie fluids. 

The package is implemented mostly in C++ to handle the numerical computations involved in evaluating the collision integrals and the radial distribution function at contact for the target fluids, with the possibility of setting up multithreading at compile time.

KineticGas can be used to predict diffusion coefficients, thermal diffusion coefficients, viscosities and thermal conductivities in gas mixtures, and is reliable over a large range of temperatures and pressures. The package also contains an extensive database of fluid parameters collected from the open literature.

# [KineticGas homepage](https://thermotools.github.io/KineticGas)
The full documentation, with installation- and getting started-guides can be found on the [KineticGas homepage](https://thermotools.github.io/KineticGas).
This readme is only intended to provide a minimal introduction, and may be out-of-sync with the `pykingas` version currently
on `PyPI`.

![](https://thermotools.github.io/KineticGas/v2.0.0/graphics/all.gif?raw=true)
