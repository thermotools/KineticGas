# Installing KineticGas

KineticGas is available on PyPi as the [`pykingas`](https://pypi.org/project/pykingas/) package, compiled for MacOS running on Apple Silicon, for python versions 3.8-3.11.

For other operating systems, KineticGas must currently be built from source.

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
 * The python installation to build against can be specified with
   * `bash cpp/build.sh -DPYTHON_EXECUTABLE=<path/to/python>`
   * Where `<path/to/python>` can (usually) be replaced by `$(which python)`.
 * If none of the above works, please feel free to leave an issue.


