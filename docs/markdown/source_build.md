# Installing KineticGas

An up-to-date release of the KineticGas package is coming to PyPi as [pykingas](https://pypi.org/project/pykingas/#files), at the moment KineticGas must be built from source.

## Building from source

The current build system is set up to support Mac/Linux. If you create a build system for Windows, please leave a PR.

### First Try
If all goes well, running

```
bash cpp/build.sh
pip install .
```

From the top level directory should provide you with an installation of the `KineticGas` python package `pykingas`.

#### Short explanation

The `bash` script `cpp/build_kingas.sh` uses `cmake` and `make` to compile the binary that is called from the python module. Then it moves the binary to the `pykingas` directory.

### When something goes wrong

 * The variable `PYBIND11_ROOT`, set in `cpp/CMakeLists.txt` must contain the path to the root directory of your [`pybind11`](https://github.com/pybind/pybind11) installation.
   * If you don't have `pybind11`:
     * Run `git clone https://github.com/pybind/pybind11.git`
     * Set `PYBIND11_ROOT` in `cpp/CMakeLists.txt` to the resulting directory.
 * The system arcitecture to compile for, and the python version, are specified in `cpp/CMakeLists.txt`, modify these as needed.
 * The bash script `cpp/build.sh` sets the environment variables `CC` and `CXX`, these may also need to be modified for your system.


