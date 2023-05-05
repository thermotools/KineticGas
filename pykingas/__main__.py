'''
Author : Vegard Gjeldvik Jervell
Purpose : File to run pykingas as a module with the command 'python -m pykingas <args>'

Args are separated by a single space and can be supplied in any order. Valid args are:
-test : run test suite
    -print : Print information during testing
    -plot : Plot information during testing
    -silent : Suppress output during testing
    -force : Force plotting of results, even if test is successful
-debug : Use the debug build
'''

import sys, shutil, os, atexit

args = sys.argv

if '-test' in args:
    '''
        Test suite for pykingas package
        Runs the tests placed in the list 'test_pkgs'
        Each element is a function that accepts the keyword arguments 'do_plot' and 'do_print'
        The function should return a single value, 0 for successfull tests, not 0 for failed tests.
        Setting the kwargs to 'True' will yield additional information about the tests.
        Note that not all tests display any extra information.
    '''
        
    from pykingas.tests import mie_unittests, collision_integral_unittests, py_integration_tests, factorial_unittests
    from pykingas import bcolors

    print(f'{bcolors.HEADER}Testing from', __file__, f'{bcolors.ENDC}')

    test_pkgs = [factorial_unittests.run_tests, mie_unittests.run_tests, 
                collision_integral_unittests.run_tests, py_integration_tests.run_test]

    r = 0
    for test in test_pkgs:
        r = test(do_print=('-print' in args), do_plot=('-plot' in args))
        if r != 0:
            exit(r)
    exit(r)
