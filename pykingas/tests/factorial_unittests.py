from pykingas import factorial_tests, bcolors
import sys

def run_tests(do_print=False, do_plot=False):
    r = factorial_tests()
    if r != 0:
        r += 1000
        print(f'{bcolors.FAIL}Factorial test failed with exit code :', r, f'{bcolors.ENDC}')
    else:
        print(f'{bcolors.OKGREEN}Factorial test was successful!{bcolors.ENDC}')
    return r

if __name__ == '__main__':
    run_tests(do_print=('-print' in sys.argv), do_plot=('-plot' in sys.argv))