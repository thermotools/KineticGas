---
layout: default
version: 
title: References
permalink: /please_cite.html
---

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
* Quantum mechanical methods and Feynman-Hibbs corrections
  * [The limits of Feynman–Hibbs corrections in capturing quantum-nuclear contributions to thermophysical properties](https://doi.org/10.1063/5.0295049) (V. G. Jervell and Ø. Wilhelmsen, J. Chem. Phys. 2025)
* Ab initio reference potentials
  * Argon, Neon, Hydrogen and Helium: [The limits of Feynman–Hibbs corrections in capturing quantum-nuclear contributions to thermophysical properties](https://doi.org/10.1063/5.0295049) (V. G. Jervell and Ø. Wilhelmsen, J. Chem. Phys. 2025)

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