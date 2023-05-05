from cgi import test
from pykingas import MieKinGas, HardSphere, PseudoHardSphere, bcolors, suppress_stdout
import matplotlib.pyplot as plt
import numpy as np
import sys, warnings

# Exit codes 3***
# Second digit identifies function, third digit identifies test

hardsphere = HardSphere.HardSphere('AR,C1')
p_hardsphere = PseudoHardSphere.PseudoHardSphere('AR,C1')

T = 300

def test_w_vs_HS(do_plot=False, do_print=False):
    # Return value upon failing a test is 200 + 10 * r + l
    # So a failure at r = 3, l = 2 gives exit code 232
    r_list = [1, 2, 3]
    l_list = [1, 2, 3]
    rgrid, lgrid = np.meshgrid(r_list, l_list)
    numeric = np.empty_like(rgrid, float)
    analytic = np.empty_like(rgrid, float)
    print(f'{bcolors.OKCYAN}Computing W-integrals for HS potential ...', f'{bcolors.ENDC}')
    r, v = 0, 0
    for ri in range(len(r_list)):
        for li in range(len(l_list)):
            numeric[ri, li] = p_hardsphere.cpp_kingas.w_integral(0, 0, T, l_list[li], r_list[ri])
            analytic[ri, li] = hardsphere.cpp_kingas.w_integral(0, 0, T, l_list[li], r_list[ri])

            if abs((numeric[ri, li] / analytic[ri, li]) - 1) > 2.5e-2:
                r = 000 + 10 * rgrid[ri, li] + lgrid[ri, li]
                v = 'For r = ' + str(rgrid[ri, li]) + ', l = ' + str(lgrid[ri, li]) + '\n' \
                    'Numeric HS dimetionless collision integral is : ' + str(numeric[ri, li]) + '\n' \
                    'Analytic HS dimetionless collision integral is : ' + str(analytic[ri, li])
                break
        if r != 0:
            break

    if do_plot is True and (r != 0 or '-force' in sys.argv):
        plot_w_vs_HS(rgrid, lgrid, numeric, analytic)

    return r, v

def plot_w_vs_HS(rgrid, lgrid, numeric, analytic):

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    ax.scatter(lgrid, rgrid, 100 * (numeric - analytic) / analytic)
    ax.set_xlabel(r'$r$ [-]')
    ax.set_ylabel(r'$\ell$ [-]')
    ax.set_zlabel(r'$\Delta_{HS}W_{r,\ell} / W^{HS}_{r,\ell}$ [%]')
    ax.set_title('Relative deviation between numeric and analytic\ndimentionless collision integrals (%)')
    plt.show()

def test_w_vs_high_res_mie(do_plot=False, do_print=False):
    # Return value upon failing a test is 200 + 10 * r + l
    # So a failure at r = 3, l = 2 gives exit code 232
    kin = MieKinGas.MieKinGas('HE,NE', N=1)
    r_list = [1, 2, 3]
    l_list = [1, 2, 3]
    rgrid, lgrid = np.meshgrid(r_list, l_list)
    numeric = np.empty_like(rgrid, float)
    high_res = None
    print(f'{bcolors.WARNING}Collision integral resolution test for Mie potentials is deactivated!')
    return 0, 0
    print(f'{bcolors.OKCYAN}Computing W-integrals for Mie potential ...', f'{bcolors.ENDC}')
    r, v = 0, 0
    for ri in range(len(r_list)):
        for li in range(len(l_list)):
            numeric[ri, li] = kin.cpp_kingas.w_integral(0, 0, T, l_list[li], r_list[ri])

            if abs((numeric[ri, li] / high_res[ri][li]) - 1) > 2.5e-2:
                r = 100 + 20 * rgrid[ri, li] + lgrid[ri, li]
                v = 'For r = ' + str(rgrid[ri, li]) + ', l = ' + str(lgrid[ri, li]) + '\n' \
                    'Numeric Mie dimetionless collision integral is : ' + str(numeric[ri, li]) + '\n' \
                    'High-res. Mie dimetionless collision integral is : ' + str(high_res[ri][li])
                break
        if r != 0:
            break

    if do_plot is True and (r != 0 or '-force' in sys.argv):
        plot_w_vs_HS(rgrid, lgrid, numeric, high_res)

    return r, v

def run_tests(do_plot=False, do_print=False):
    '''
        Submodule for testing collision integrals
        Each test in 'tests' must accept two arguments: 'do_plot' and 'do_print' and return two values
        The first is the exit status of the test (0 for successfull, !0 otherwise)
        The second value is some information about the test that failed
    '''
    tests = [test_w_vs_HS, test_w_vs_high_res_mie]
    r = 0
    for t in tests:
        with suppress_stdout('-silent' in sys.argv):
            r, v = t(do_plot, do_print)
        if r != 0:
            r += 3000
            if do_print:
                print(r, v)
            print(f'{bcolors.FAIL}Collision integral tests failed with exit code', r, f'{bcolors.ENDC}')
            break
    if r == 0:
        print(f'{bcolors.OKGREEN}Collision integral tests were successful!{bcolors.ENDC}')
    return r

if __name__ == '__main__':
    do_plot, do_print = False, False
    if '-plot' in sys.argv:
        do_plot = True
    if '-print' in sys.argv:
        do_print = True
    
    run_tests(do_plot=True, do_print=True)