Unittest package for different parts of pykingas

All tests are run through the `__main__` file of `pykingas`. Requirements for a test package for new features are:

Must contain a function with a signature equivalent to `run_tests(do_plot=False, do_print=False)`. Said function must return a single value, with 0 signifying successful tests.

Failed tests should return a four-digit error code. The first digit should identify the file (ie. be different from all other test packages), the second digit should identify the test function that failed. The final two digits are used to identify where in the test function failure happened.

Current codes are:

1*** - Factorial module unittests
2*** - MieKingGs module unittests - Tests that the integrals in Spherical compute (does not check that they are correct)
3*** - Collision integral unittests - Tests that numerical collision integrals are accurate
4*** - Python side integration tests - Checks that output values for transport coefficients are equal to control values