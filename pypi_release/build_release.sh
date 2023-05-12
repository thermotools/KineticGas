#!/bin/bash
# Script to build wheels for pypi
# Note: Uses $(which pythonX.Y) to get the path to various python interpreters
set -e

pykingas_version="2.0.0"
os_tag="macosx_11_0_arm64"
py_versions=("3.8" "3.9" "3.10" "3.11") # Need to build different wheels for every python version, because pybind11.
wheelhouse_dir="pypi_release/wheelhouse"
for py_version in "${py_versions[@]}"; do

  echo "Building wheel for python version : ${py_version}"

  [ -f "pykingas/KineticGas_r.so" ] && rm pykingas/KineticGas_r.so
  [ -d "cpp/release" ] && rm -rf cpp/release

  bash cpp/build.sh -DPYTHON_EXECUTABLE="$(which "python${py_version}")"

  python -m pip wheel --wheel-dir=${wheelhouse_dir} .
  mv pypi_release/wheelhouse/pykingas-${pykingas_version}-py3-none-any.whl pypi_release/dist/pykingas-${pykingas_version}-cp"${py_version//.}"-none-${os_tag}.whl

done

# echo "Finished Building : Uploading test to testpypi"
# twine upload -r testpypi pypi_release/dist/*

exit 0

# Code to upload to PyPi :
twine upload pypi_release/dist/*

exit 0

#Code to Delocate wheels
# It appears that there is nothing to delocate here (at least according to delocate-listdeps)

cd pypi_release/wheelhouse
delocate-wheel -w fixed_wheels -v pykingas-${pykingas_version}-py3-none-any.whl
cd fixed_wheels
mv pykingas-${pykingas_version}-py3-none-any.whl pykingas-${pykingas_version}-cp${py_version//.}-none-macosx_11_0_arm64.whl
cd ..
rm pykingas-${pykingas_version}-py3-none-any.whl # Preventing unintentionally overwriting the fixed wheel