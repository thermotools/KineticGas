'''
Author: Vegard Gjeldvik Jervell
Purpose: Exposes functions and classes to a user-friendly interface

Note : Calling 'python <myscript.py> <arg>' where <arg> is either '-d', '-debug' or 'debug' will
        result in the debug build of the cpp-module being used.
        All methods and classes exported in the cpp-module are directly exposed through this interface, such
        that they can be imported as 'from pykingas import <thing_exposed_in_bindings.cpp>'.
'''

import sys

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
from . import libpykingas

