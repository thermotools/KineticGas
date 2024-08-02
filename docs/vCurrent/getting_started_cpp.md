---
layout: default
version: 
title: Getting started - In C++
permalink: /vcurrent/getting_started_cpp.html
---

## Getting started: In C++

For instructions on building the `KineticGas` C++ library, see the [installation guide](installation#c).

A basic example showing initialization of a model is found in [`cppExamples/basic.cpp`](), the `cppExamples` directory also contains a `CMakeLists.txt` showing how to obtain the required headers for the `KineticGas` library, as well as link the library to your program.

## Initializing a model

To initialize a model, `include` the appropriate header file, and specify the components to model with a comma separated string, as
```C
#include "MieKinGas.h"

int main(){
    MieKinGas mie("HE,NE") // Mixture of helium and neon.
}
```
The component identifiers used are equivalent to the file names of the [fluid files](https://github.com/thermotools/KineticGas/tree/main/fluids), and are summarised [here](fluid_identifiers.html)

## Computing properties

The interfaces for property calculations are more or less equivalent to those used in Python. The major differences you should be aware of are
* Diffusion coefficients are returned as an `Eigen::MatrixXd`
* Frames of reference are specified with the `FrameOfReference` enum, found in `KineticGas.h`. Valid values are
  * `FrameOfReference::CoN` - Center of moles
  * `CoM` - Center of mass (barycentric)
  * `CoV` - Center of volume
  * `solvent` - Solvent, solvent index is the `dependent_idx`, which defaults to the last component.
  * `zarate`, `zarate_x`, and `zarate_w` - See the [memo](/KineticGas/memo/diffusion_definitions.pdf)
  * See the python docs and the [memo](/KineticGas/memo/diffusion_definitions.pdf) for more details on definitions of the diffusion coefficients.
  * *Note*: The optional `dependent_idx` argument to `interdiffusion` supports python-style negative indexing.

