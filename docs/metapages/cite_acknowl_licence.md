---
layout: default
version: 
title: Please Cite
permalink: /please_cite.html
---

KineticGas has been developed throughout several works. If you are referencing the package, please cite the works

* General usage
   * [Revised Enskog theory for Mie fluids: Prediction of diffusion coefficients, thermal diffusion coefficients, viscosities and thermal conductivities](https://doi.org/10.1063/5.0149865) (Vegard G. Jervell and Øivind Wilhelmsen, 2023)
   * [The Kinetic Gas theory of Mie fluids](https://ntnuopen.ntnu.no/ntnu-xmlui/handle/11250/3029213) (Vegard G. Jervell, 2022)
* Connection to Non-Equilibrium thermodynamics (Onsager coefficients)
   * [The influence of thermal diffusion on water migration through a porous insulation material](10.1016/j.ijheatmasstransfer.2024.125576) (V. G. Jervell, M. Aa. Gjennestad, T. T. Trinh, Ø. Wilhelmsen, 2024)

## Acknowledgments and sources
This implementation of the Revised Enskog solutions is build upon the work presented by M. López de Haro, E. D. G. Cohen, and J. Kincaid in the series of papers *The Enskog Theory for multicomponent mixtures I - IV*, J. Chem. Phys. (1983 - 1987) ([I](https://doi.org/10.1063/1.444985), [II](https://doi.org/10.1063/1.446388), [III](https://doi.org/10.1063/1.446463), [IV](https://doi.org/10.1063/1.452243)).

The implementation utilises the explicit summational expressions for the square bracket integrals published by Tompson, Tipton and Loyalka in *Chapman–Enskog solutions to arbitrary order in Sonine polynomials I - IV* (Physica A, E. J. Mech. - B) 2007-2009 ([I](https://doi.org/10.1016/j.physa.2006.12.001), [II](https://doi.org/10.1016/j.euromechflu.2008.09.002), [III](https://doi.org/10.1016/j.euromechflu.2008.12.002), [IV](https://doi.org/10.1016/j.euromechflu.2009.05.002)).

The work by T. Lafitte, A. Apostolakou, C. Avendaño, A. Galindo, C. Adjiman, E. Müller and G. Jackson, *Accurate statistical associating fluid theory for chain molecules formed from Mie segments* [J. Chem. Phys. 2013](https://doi.org/10.1063/1.4819786) is also of great importance to this implementation.

## Licence

The KineticGas package is distributed as free software under the MIT licence.