---
layout: default
version: 
title: Advanced usage
permalink: /vcurrent/advanced.html
---

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