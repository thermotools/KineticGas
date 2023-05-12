#!/bin/bash

set -e

pykingas_version="2.0.5"
py_versions=("3.8" "3.9" "3.10" "3.11")

for py_version in "${py_versions[@]}"; do

  source venv"${py_version//.}"/bin/activate

  echo "Testing install for Python version $(python --version)"
  python -m pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple pykingas==${pykingas_version}
  python pypi_release/test_import.py

done