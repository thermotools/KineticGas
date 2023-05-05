"""
Integration test module for the pykingas package

For all potential models a test value for the thermal conductivity, viscosity, diffusion matrix and thermal diffusion
-coefficients, -factors, and ratios have been computed and added to the control_value_dict of the corresponding test
function. Note that the control values are computed at high density, to ensure that the influence of the rdf is
significant.

The test simply consists of computing the transport coefficients at the same conditions at which the control values
were computed, and checking that they do not deviate from each other. This is a very simple but effective way of
testing that some change to the code has not effected the output of the package, and that all computations run.
"""

from pykingas import bcolors
from pykingas.MieKinGas import MieKinGas
from pykingas.HardSphere import HardSphere

# Exit codes 4***
# Second digit identifies potential mode
# Final digit identifies transport parameter test

FLT_EPS = 1e-8

def diffusion_test(model, model_id):
    control_value_dict = {'HS' : [[ 0.00040355,  0.00017183,  0.        ],
                                  [ 0.00011033,  0.00033149,  0.        ],
                                  [-0.00051388, -0.00050332,  0.        ]],
                          'Mie' : [[ 0.00052809,  0.00028933,  0.        ],
                                   [ 0.00018595,  0.00040727,  0.        ],
                                   [-0.00071404, -0.0006966,   0.        ]]}
    T = 300
    x = [0.3, 0.2, 0.5]
    Vm = 0.24
    D = model.interdiffusion(T, Vm, x, N=2, use_binary=False)
    if max(abs(D - control_value_dict[model_id]).flatten()) > FLT_EPS:
        return 1, D - control_value_dict[model_id]
    return 0, 0


def conductivity_test(model, model_id):
    control_value_dict = {'HS': 0.005790548371162743,
                          'Mie': 0.0064809028677221955}
    T = 300
    x = [0.3, 0.2, 0.5]
    Vm = 0.24
    k = model.thermal_conductivity(T, Vm, x, N=2)
    if abs(k - control_value_dict[model_id]) > FLT_EPS:
        return 2, k - control_value_dict[model_id]
    return 0, 0

def viscosity_test(model, model_id):
    control_value_dict = {'HS': 2.4497519210969225e-05,
                          'Mie': 2.57665247077739e-05}
    T = 300
    x = [0.3, 0.2, 0.5]
    Vm = 0.24

    eta = model.viscosity(T, Vm, x, N=2)
    if abs(eta - control_value_dict[model_id]) > FLT_EPS:
        return 3, eta - control_value_dict[model_id]
    return 0, 0

def thermal_diffusion_coeff_test(model, model_id):
    control_value_dict = {'HS': [-0.00042014,  0.0002648,   0.00015534],
                          'Mie': [-0.00047644,  0.00025725,  0.00021919]}
    T = 300
    x = [0.3, 0.2, 0.5]
    Vm = 0.24

    DT = model.thermal_diffusion_coeff(T, Vm, x, N=2)
    if max(abs(DT - control_value_dict[model_id])) > FLT_EPS:
        return 4, DT - control_value_dict[model_id]
    return 0, 0

def thermal_diffusion_ratio_test(model, model_id):
    control_value_dict = {'HS': [2.58654748e+00,  1.12416881e+00, -2.04536190e-03],
                          'Mie': [2.63136806e+00,  1.05616360e+00, -1.05074199e-03]}
    T = 300
    x = [0.3, 0.2, 0.5]
    Vm = 0.24

    kT = model.thermal_diffusion_ratio(T, Vm, x, N=2)

    if max(abs(kT - control_value_dict[model_id])) > FLT_EPS:
        return 5, kT - control_value_dict[model_id]
    return 0, 0

def thermal_diffusion_factor_test(model, model_id):
    control_value_dict = {'HS': [[ 0.,          1.46237868,  2.58859284],
                                 [-1.46237868,  0.,          1.12621417],
                                 [-2.58859284, -1.12621417,  0.        ]],
                          'Mie': [[ 0.,          1.57520447,  2.6324188],
                                   [-1.57520447,  0.,          1.05721434],
                                   [-2.6324188,  -1.05721434,  0.        ]]}
    T = 300
    x = [0.3, 0.2, 0.5]
    Vm = 0.24

    alpha_T = model.thermal_diffusion_factor(T, Vm, x, N=2)
    if max(abs(alpha_T - control_value_dict[model_id]).flatten()) > FLT_EPS:
        print(alpha_T)
        print(control_value_dict[model_id])
        return 6, alpha_T - control_value_dict[model_id]
    return 0, 0

def run_test(do_plot=False, do_print=False):
    # control values are precomputed values for the conditions in this file.
    # If they change, that means something's up...

    tests = [diffusion_test, conductivity_test, viscosity_test, thermal_diffusion_coeff_test,
             thermal_diffusion_ratio_test, thermal_diffusion_factor_test]
    comps = 'AR,KR,HE'
    models = [HardSphere(comps), MieKinGas(comps)]
    ids = ['HS', 'Mie']
    r = 0
    for i, (model, ident) in enumerate(zip(models, ids)):
        for test in tests:
            r, v = test(model, ident)
            if r != 0:
                r += 4000 + i * 100
                if do_print:
                    print(r, v)
                print(f'{bcolors.FAIL}Python-side integration tests failed with exit code', r, f'{bcolors.ENDC}')
                break
    if r == 0:
        print(f'{bcolors.OKGREEN}Python-side integration tests were successful!{bcolors.ENDC}')

    return r

if __name__ == '__main__':
    run_test()




