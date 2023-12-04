# KineticGas test suite

This test module is intended to be run using `pytest`. To run the whole test suite, use
```bash
pip install pytest
cd tests/
pytest
```
to run a single test module, use
```bash
pytest <my_test_module.py>
```
To run a specific test function, use
```bash
pytest <my_test_module.py>::<my_function>
```
For more info on `pytest` see: https://docs.pytest.org/en/7.4.x/how-to/usage.html

## Adding new models
When implementing a new model, the model should be added to the `models` list in `tools.py`. Once that is done, it will
be piped through all existing tests.

Note: The current tests only check for physical consistency of the computed coefficients, for example that for an ideal
gas, $\sum x_i k_{T,i} = 1$, and that the order of components does not effect output, etc. This means that a model
can give completely nonsensical results (e.g. coefficients with the wrong units) and still pass all tests. However, 
given that you have an idea what you are doing, this test module helps you check that implemented methods satisfy
certain physical consistency checks.

## Adding new tests

All tests should use the
```python
from tools import models
@pytest.mark.parametrize('model', models)
```
decorator to run the test for all KineticGas models.

Tests should also run for at least two Enskog approximation orders, preferably one odd and one even. This helps us catch
potential indexing errors in methods that access the Sonine polynomial expansion coefficients.

Current test modules are:
* `test_binary_component_order.py` - Checks that order of components in binary mixtures does not affect output
* `test_summation_consistency.py` - Checks that various things that should sum to 1 or 0 do that.
* `test_binary_limits.py` - Checks that ternary coefficients reduce to the appropriate binary coefficients in the binary limit.
* `test_flux_transforms.py` - Checks that the frame of reference transformations give fluxes that appropriately sum to zero.