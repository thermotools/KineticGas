'''
Author: Vegard Gjeldvik Jervell
Purpose: Exposes functions and classes to a user-friendly interface

Note : Calling 'python <myscript.py> <arg>' where <arg> is either '-d', '-debug' or 'debug' will
        result in the debug build of the cpp-module being used.
        All methods and classes exported in the cpp-module are directly exposed through this interface, such
        that they can be imported as 'from pykingas import <thing_exposed_in_bindings.cpp>'.
'''

import sys, os
from contextlib import contextmanager

@contextmanager
def suppress_stdout(do_suppress=False):
    if do_suppress is True:
        with open(os.devnull, "w") as devnull:
            old_stdout = sys.stdout
            sys.stdout = devnull
            try:
                yield
            finally:
                sys.stdout = old_stdout
    else:
        yield

class bcolors: # For fancy (readable) printing during tests
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

args = sys.argv

if '-debug' in args or 'debug' in args or '-d' in args:
    from pykingas import KineticGas_d
    __cpp_Module__ = KineticGas_d
else:
    from pykingas import KineticGas_r
    __cpp_Module__ = KineticGas_r

# Expose everything in the __cpp_Module__ through the pykingas module
for _attr in dir(__cpp_Module__):
    if _attr[:2] != '__': #Exclude macros
        setattr(sys.modules[__name__], _attr, getattr(__cpp_Module__, _attr))

from pykingas import py_KineticGas
KineticGas = py_KineticGas.py_KineticGas
