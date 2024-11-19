---
layout: default
version: 
title: Installing KineticGas
permalink: /vcurrent/source_build.html
---

KineticGas is available on PyPi as the [`pykingas`](https://pypi.org/project/pykingas/) package, for python versions 3.8-3.11, compiled for MacOS running on Apple Silicon, Linux and Windows.

In addition, wheels versions of `KineticGas > 2.0.0` for macOS, Linux and Windows, as well as wheels for the latest version on GitHub can be downloaded [here](https://github.com/thermotools/KineticGas/releases). Instructions for installing with `pip` directly from a downloaded wheel are provided at the linked page.

For MacOS running on Intel, or other operating systems, KineticGas must currently be built from source or installed from one of the distributed wheels linked above.

A KineticGas C++ library is available, and can be built using `cmake` and `make`.

- [Python - `pykingas`](#python---pykingas)
  - [Dependencies](#dependencies)
  - [Building from source](#building-from-source)
    - [First Try](#first-try)
      - [Short explanation](#short-explanation)
    - [When something goes wrong](#when-something-goes-wrong)
- [C++](#c)
  - [Building and installing](#building-and-installing)
  - [Fluid file search path](#fluid-file-search-path)
  - [Linking to the KineticGas library](#linking-to-the-kineticgas-library)

# Python - `pykingas`

## Dependencies

The Python package dependencies are listed in the `pyproject.toml` file in the root directory of the package.

To compile the binary that is called from the python wrapper, [pybind11](https://pybind11.readthedocs.io/en/stable/) is required. `pybind11` is included in `cpp/external` as a git submodule, so cloning the `KineticGas` repository should provide you with the files you need.

A standalone C++ module, that works without the python wrapper is currently under development. See the branch `pure_cpp/` for the most up-to-date version there.


## Building from source

Python wheels for the latest version of KineticGas on `main` are built for macOS and Windows using `cibuildwheels`, and distributed [here](https://github.com/thermotools/KineticGas/releases/tag/Latest-beta).

A build system using `cmake` and `make` is set up to support Mac, Linux and Windows.

### First Try
If all goes well, running

```
git clone https://github.com/thermotools/KineticGas.git
cd KineticGas
git submodule update --init --recursive
mkdir build
cd build
cmake ..
make install
pip install ..
```

make sure to activate a virtual environment first if you want to avoid doing system-level installs.

#### Short explanation

The dynamic library `libpykingas` will be built and installed to the `pykingas` directory, additionally, the `fluids` directory containing the fluid parameter database is copied into the `pykingas` directory.

### When something goes wrong

*Note:* The build system has been changed relatively recently, and is less tested than the build system that was used in the `2.0.0` release. If you encounter issues, please don't hesitate to post an issue on github so that we can improve robustness.

* Warning that thermopack is not installed
  * The easiest way to obtain the `ThermoPack` dynamic library (which `KineticGas` needs) is likely to download the appropriate zip file [here](https://github.com/thermotools/thermopack/releases), unzip it, and set the environment variable `THERMOPACK_DIR` to the resulting directory (where `thermopack-config.cmake` is located).
    * On Linux and macOS: `export THERMOPACK_DIR=/path/to/thermopack-<system>/`
    * On Windows (powershell): `$THERMOPACK_DIR = C:\path\to\thermopack-<system>\thermopack-<system>`
    * To check that it is set correctly: `ls ${THERMOPACK_DIR}` should give a list of files including `thermopack-config.cmake`.
  * The `KineticGas` library has a dependency on the `ThermoPack` C++ wrapper. If you have not installed thermopack, the build system will generate a target from the `thermopack` submodule. Running `make install` should build and install this target, re-running `cmake ..` after building and installing `thermopack` should then give output telling you that `thermopack` has been found and is installed.
  * If you have installed thermopack, run `export THERMOPACK_DIR=<path/to/thermopack>`, to help `cmake` find your installation.


# C++

The KineticGas C++ library is built using `cmake` and `make`. All dependencies are included as git submodules under `cpp/external`, and should be properly retrieved when you clone the `KineticGas` repository and run `git submodule update --init --recursive`. 

*Note*: `KineticGas` depends on [`ThermoPack`](https://thermotools.github.io/thermopack/). If an installation of `ThermoPack` is not found, the build system will attempt to compile it as part of the build process. If you already have an installation of `ThermoPack`, setting the environment variable `THERMOPACK_DIR` to the root directory of `ThermoPack` (where `thermopack-config.cmake` is found), that installation of `ThermoPack` will be used istead of re-compiling. You can also download a binary distribution of ThermoPack at the [ThermoPack repository](https://github.com/thermotools/thermopack/releases).

## Building and installing

If all goes well, you should be able to build the `KineticGas` C++ library by running

```bash
git clone https://github.com/thermotools/KineticGas.git
cd KineticGas
git submodule update --init --recursive
mkdir build
cd build
cmake -Dpurecpp=ON -Dpylib=OFF ..
make install
```

This will provide you with the `lib/libkineticgas.[so/dylib/dll]` dynamic library, and the minimal example program `build/run_kineticgas`, which is built from the source file at `cpp/run_kineticgas.cpp`.

## Fluid file search path

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

## Linking to the KineticGas library

An example program with a `CMakeLists.txt` demonstrating how you can include and link to the `KineticGas` library once it is installed is found in `KineticGas/cppExamples`. 

In short terms: Setting the environment variable `KINETICGAS_DIR` to the top-level directory of the KineticGas package (where `kineticgas-config.cmake` is found), should allow `cmake` to find the library using `find_library(KINETICGAS)`. Some convenience variables are set once the library is found:

* `KINETICGAS_ROOT` : Path to root directory of the package
* `KINETICGAS_INSTALLED` : `TRUE` if the dynamic library is found in the correct install location, `FALSE` otherwise
* `KINETICGAS_LIB` : Path to the `libkineticgas` dynamic library
* `KINETICGAS_INCLUDE` : List of include paths needed to include the kineticgas headers and dependencies
* `kineticgas` : Exported target, linking to this target should automatically add the appropriate directories to your include path. 