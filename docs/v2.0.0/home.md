---
version: 2.0.0
sidebar_file: sidebar_2.0.0.md
layout: home
title: ThermoPack v2.0.0
permalink: /v2.0.0/home.html
---

KineticGas version 2.0.0 was released on May 12th 2023, and is the first release of KineticGas with a fully functional 
implementation of RET-Mie to be made available  on the Python packaging index (PyPI).

KineticGas is an implementation of Revised Enskog Theory (RET) for spherical potentials. The most notable of which is the implementation of RET-Mie, the Revised Enskog Theory for Mie fluids. 

The package is implemented mostly in C++ to handle the numerical computations involved in evaluating the collision integrals and the radial distribution function at contact for the target fluids, with the possibility of setting up multithreading at compile time.

KineticGas can be used to predict diffusion coefficients, thermal diffusion coefficients, viscosities and thermal conductivities in gas mixtures, and is reliable over a large range of temperatures and pressures. The package also contains an extensive database of fluid parameters collected from the open literature.

![](https://github.com/thermotools/KineticGas/blob/main/docs/figures/all.gif?raw=true)