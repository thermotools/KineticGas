---
layout: default
version: 
title: Installing KineticGas
permalink: /vcurrent/source_build.html
---

KineticGas is available on PyPi as the [`pykingas`](https://pypi.org/project/pykingas/) package, for python versions 3.8-3.11, compiled for MacOS running on Apple Silicon, Linux and Windows.

For MacOS running on Intel, or other operating systems, KineticGas must currently be built from source.

## Dependencies

The Python package dependencies are listed in the `setup.py` file in the root directory of the package.

To compile the binary that is called from the python wrapper, [pybind11](https://pybind11.readthedocs.io/en/stable/) is required.

A standalone C++ module, that works without the python wrapper is currently under development. See branches under `pure_cpp/` for the most up-to-date version there.


## Building from source

A build system using `cmake` and `make` is set up to support Mac, Linux and Windows. For Mac machines running on intel chips, one compiler flag must be modified.

### First Try - For Mac/Linux
If all goes well, running

```
bash cpp/build.sh
pip install .
```

From the top level directory should provide you with an installation of the `KineticGas` python package `pykingas`.

For Mac's running on an intel chip, the compiler flag `-arch arm64` which is set in `cpp/CMakeLists.txt` must be removed or changed to `-arch x86_64`.

#### Short explanation

The `bash` script `cpp/build_kingas.sh` uses `cmake` and `make` to compile the binary that is called from the python module. Then it moves the binary to the `pykingas` directory.

### When something goes wrong

 * The variable `PYBIND11_ROOT`, set in `cpp/CMakeLists.txt` must contain the path to the root directory of your [`pybind11`](https://github.com/pybind/pybind11) installation.
   * If you don't have `pybind11`:
     * Run `git clone https://github.com/pybind/pybind11.git`
     * Set `PYBIND11_ROOT` in `cpp/CMakeLists.txt` to the resulting directory.
 * The system arcitecture to compile for, and the python version, are specified in `cpp/CMakeLists.txt`, modify these as needed.
 * The bash script `cpp/build.sh` sets the environment variables `CC` and `CXX`, these may also need to be modified for your system.
 * The python installation to build against can be specified with
   * `bash cpp/build.sh -DPYTHON_EXECUTABLE=<path/to/python>`
   * Where `<path/to/python>` can (usually) be replaced by `$(which python)`.
   * Alternatively, add the line `set(PYBIND11_PYTHON_VERSION 3)` to the top of the file `cpp/CMakeLists.txt`
   * Or: add the line `set(PYTHON_EXECUTABLE "<path/to/python>"`
 * If `cmake` starts looping infinitely:
   * You may be getting a message of the type :
```
You have changed variables that require your cache to be deleted.
Configure will be re-run and you may have to reset some variables.
The following variables have changed:
CMAKE_C_COMPILER= /Library/Developer/CommandLineTools/usr/bin/cc
CMAKE_CXX_COMPILER= /Library/Developer/CommandLineTools/usr/bin/c++
CMAKE_CXX_COMPILER= /Library/Developer/CommandLineTools/usr/bin/c++
```
   * To fix this: Ensure that there are no `set(CMAKE_CXX_COMPILER <path/to/compiler>)` or `set(CMAKE_C_COMPILER <path/to/compiler>)` statements in `cpp/CMakeLists.txt`. If you need to specify a compiler, do so by using the `export CC=<path/to/compiler>` and `export CXX=<path/to/compiler>` statements in `cpp/build.sh`.
     * It appears that everything works fine as long as the environment variables `CC` and `CXX` match the variables `CMAKE_C_COMPILER` and `CMAKE_CXX_COMPILER`.
     * **NOTE** : You may need to delete the file `cpp/release/CMakeCache.txt` for changes to take effect.
 * If you get an error message when the file `bindings.cpp` is compiling, that originates from the `pybind11` headers:
   * You are likely getting an error of the type
```
 error: address of overloaded function '<some_func>' does not match required type 'pybind11::overload_cast<some stuff>'
```
and 
```
error: static_assert failed due to requirement 'detail::integral_constant<bool, false>::value' "pybind11::overload_cast<...> requires compiling in C++14 mode"
```
   * and you are likely using `clang` compiled for the C++-11 standard. (The compiler located in `/usr/...` on Mac is likely `clang`, even though it is called `gcc`)
   * To fix the issue: 
     * Install `gcc` with [homebrew](https://formulae.brew.sh/formula/gcc).
     * Locate the compilers you've installed (`which gcc-13` and `which g++-13` should work if you installed gcc version 13.x.x)
     * Set the environment variables `CC` and `CXX` in `cpp/build.sh` to the path to these compilers by modifying the `export` statements. For example
```
export CC=/opt/homebrew/bin/gcc-13
export CXX=/opt/homebrew/bin/g++-13
```
   * **NOTE**: You may need to delete the file `cpp/release/CMakeCache.txt` for changes to take effect.
 * If none of the above works, please feel free to leave an issue.

### For Windows

Running `cmake` from the `cpp` directory should produce an MSVC solution file. Building this solution should generate the file `KineticGas_r.cp<python-version>-win_amd64.pyd` which will be displayed as a "python extension module". Copy this file to the `pykingas` directory, and run `pip install .` from the top-level directory (where `setup.py`) is found.
