---
layout: default
version: 
title: Installing KineticGas
permalink: /vcurrent/source_build.html
---

KineticGas is available on PyPi as the [`pykingas`](https://pypi.org/project/pykingas/) package, for python versions 3.8-3.11, compiled for MacOS running on Apple Silicon, Linux and Windows.

In addition, wheels versions of `KineticGas > 2.0.0` for macOS, Linux and Windows can be downloaded [here](https://github.com/thermotools/KineticGas/releases). Instructions for installing with `pip` directly from a downloaded wheel are provided at the linked page.

For MacOS running on Intel, or other operating systems, KineticGas must currently be built from source or installed from one of the distributed wheels linked above.

## Dependencies

The Python package dependencies are listed in the `pyproject.toml` file in the root directory of the package.

To compile the binary that is called from the python wrapper, [pybind11](https://pybind11.readthedocs.io/en/stable/) is required. `pybind11` is included in `cpp/external` as a git submodule, so cloning the `KineticGas` repository should provide you with the files you need.

A standalone C++ module, that works without the python wrapper is currently under development. See the branch `pure_cpp/` for the most up-to-date version there.


## Building from source

Python wheels for the latest version of KineticGas on `main` are built for macOS and Windows using `cibuildwheels`, and distributed [here](https://github.com/thermotools/KineticGas/releases).

A build system using `cmake` and `make` is set up to support Mac, Linux and Windows.

### First Try
If all goes well, running

```
git clone https://github.com/thermotools/KineticGas.git
cd KineticGas
mkdir build
cd build
cmake ..
make install
pip install ..
```

make sure to activate a virtual environment first if you want to avoid doing system-level installs.

#### Short explanation

The `bash` script `cpp/build_kingas.sh` uses `cmake` and `make` to compile the binary that is called from the python module. Then it moves the binary to the `pykingas` directory.

### When something goes wrong

*Note:* The build system has been changed relatively recently, and is less tested than the build system that was used in the `2.0.0` release. If you encounter issues, please don't hesitate to post an issue on github so that we can improve robustness. Also, the old build system should still work fine. So if you are having trouble, a workaround may be to download the build files in the `v2.0.0` tagged version on github and use those.

* Error when importing `pykingas`: If you get an error of the type
```
ImportError: dlopen(/.../venv/lib/python3.11/site-packages/pykingas/libpykingas.cpython-311-darwin.so, 0x0002): tried: '/.../venv/lib/python3.11/site-packages/pykingas/libpykingas.cpython-311-darwin.so' (mach-o file, but is an incompatible architecture (have (x86_64), need (arm64e)))
```
 * set the environment variables `CC` and `CXX` with
```
export CC=/opt/homebrew/bin/gcc-13
export CXX=/opt/homebrew/bin/g++-13
```
  * This can help force compilation for `arm64`
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
  * See also: [This stackoverflow question](https://stackoverflow.com/questions/73758291/is-there-a-way-to-specify-the-c-standard-of-clangd-without-recompiling-it) for info
* If none of the above works, please feel free to leave an issue.
