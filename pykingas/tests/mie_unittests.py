from pykingas import MieKinGas, bcolors, suppress_stdout
from scipy.constants import Boltzmann
import numpy as np
import sys

# Exit codes 2***
# Second digit identifies test function, third and fourth digit identifies specific test

kin = MieKinGas.MieKinGas('AR,C1')
sigma = kin.sigma_ij[0, 0]

T, g, b = 300, 2, 7.5 * sigma

func = lambda r: kin.cpp_kingas.theta_integrand(0, 0, T, r, g, b)
theta = lambda R_min, N: kin.cpp_kingas.theta(0, 0, T, R_min, g, b, N)

R = kin.cpp_kingas.get_R(0, 0, T, g, b)

FLTEPS = 1e-12

def theta_lim():
    b_list = np.linspace(0, 150, 30) * sigma
    t_list = np.empty_like(b_list)
    for i, bi in enumerate(b_list):
        R = kin.cpp_kingas.get_R(0, 0, T, g, bi)
        t_list[i] = kin.cpp_kingas.theta(0, 0, T, g, bi)
    
    if any(np.isnan(t_list)) or any(np.isinf(t_list)):
        return 200, t_list
    elif any(np.isinf(t_list)):
        return 201, t_list
    elif abs(t_list[-1] - t_list[-2]) > FLTEPS:
        return 202, abs(t_list[-1] - t_list[-2])
    elif abs(t_list[-1] - np.pi / 2) > FLTEPS:
        return 203, abs(t_list[-1] - np.pi / 2)

    g_list = np.linspace(1e-5, 25, 30)
    for i, gi in enumerate(g_list):
        R = kin.cpp_kingas.get_R(0, 0, T, gi, b)
        t_list[i] = kin.cpp_kingas.theta(0, 0, T, gi, b)

    if any(np.isnan(t_list)):
        return 210, t_list
    elif any(np.isinf(t_list)):
        return 211, t_list
    elif abs(t_list[-1] - t_list[-2]) > 1e-6:
        return 212, abs(t_list[-1] - t_list[-2])

    else:
        return 0, 0

def chi_lim():
    b_list = np.linspace(0, 150, 30) * sigma
    chi_list = np.empty_like(b_list)
    for i, bi in enumerate(b_list):
        chi_list[i] = kin.cpp_kingas.chi(0, 0, T, g, bi)

    if any(np.isnan(chi_list)):
        return 300, chi_list
    elif any(np.isinf(chi_list)):
        return 301, chi_list
    elif abs(chi_list[-1] - chi_list[-2]) > FLTEPS:
        return 302, chi_list[-1] - chi_list[-2]
    elif abs(chi_list[-1]) > FLTEPS:
        return 303, chi_list[-1]

    g_list = np.linspace(1e-5, 10, 30)
    for i, gi in enumerate(g_list):
        chi_list[i] = kin.cpp_kingas.chi(0, 0, T, gi, b)

    if any(np.isnan(chi_list)):
        return 310, chi_list
    elif any(np.isinf(chi_list)):
        return 311, chi_list
    elif abs(chi_list[-1] - chi_list[-2]) > 1e-7:
        return 312, chi_list[-1] - chi_list[-2]

    else:
        return 0, 0

def rt_to_xy(r, t):
    x = r * np.cos(t)
    y = r * np.sin(t)
    return x, y

def xy_to_rt(x, y):
    r = np.sqrt(x**2 + y**2)
    t = np.arccos(x / r)
    return r, t

def vec_len(vec):
    return np.sqrt(np.sum(vec**2))

def normalize_vec(vec):
    return vec / vec_len(vec)

def total_energy(r, t, g):
    return kin.cpp_kingas.potential(0, 0, r) * kin.m0[0][0] / (np.prod(kin.mole_weights)) + 0.5 * vec_len(g)**2

def potential_energy(r, t):
    return kin.cpp_kingas.potential(0, 0, r) * kin.m0[0][0] / (np.prod(kin.mole_weights))

def get_path(T, g, b, y0=5):
    g = g * np.sqrt(2 * Boltzmann * T * kin.m0[0][0] / np.prod(kin.mole_weights))
    y0 = y0 * sigma
    b = b * sigma
    g = np.array([0, - g])  # Rett nedover
    x = b
    y = y0
    r0, t = xy_to_rt(x, y)
    r = r0

    x_list = [x]
    y_list = [y]
    g_list = [vec_len(g)]
    E_list = [total_energy(r, t, g)]
    F = kin.cpp_kingas.potential_derivative_r(0, 0, r) * kin.m0[0][0] / (np.prod(kin.mole_weights))
    F_vec = - F * normalize_vec(np.array([x, y]))
    a = F_vec

    dt = - 0.1 * (sigma / g[1]) # 10 tidssteg for å bevege seg 1 sigma
    i = 0
    while r <= r0:
        pos = np.array([x_list[i], y_list[i]])  # Posisjon
        r, t = xy_to_rt(pos[0], pos[1])
        if (np.sum(g * normalize_vec(pos)) < 0 # Partikkelen er på vei mot potensialet
                and (E_list[0] - potential_energy(r, t) < 0.05 * abs(E_list[0])) # Potensiell energi er veldig stor
                and vec_len(g) * dt < 5e-2 * sigma): # Tidssteg er veldig lite
            g = g - 2 * normalize_vec(pos) * np.dot(g, normalize_vec(pos)) # Behandle potensialet som en hard kule (speil g-vektor om pos-vektor)

        pos += g * dt  # Ny posisjon
        r, t = xy_to_rt(pos[0], pos[1])

        if potential_energy(r, t) > E_list[0]: #Sørger for energibevaring
            dt *= 0.5 # Reduser tidssteg og beregn forflytning på nytt
        else:
            g = g + a * dt  # Ny hastighet
            g = normalize_vec(g) * np.sqrt(2 * (E_list[0] - potential_energy(r, t))) # Korrigerer for energibevaring
            dt = 0.01 * (sigma / vec_len(g)) # 2 tidssteg for å bevege seg 1 sigma
            x_list.append(pos[0])
            y_list.append(pos[1])
            g_list.append(np.sqrt(np.sum(g**2)))
            E_list.append(total_energy(r, t, g))

            F = kin.cpp_kingas.potential_derivative_r(0, 0, r) * kin.m0[0][0] / (np.prod(kin.mole_weights))
            F_vec = - F * normalize_vec(np.array(pos))
            a = F_vec
            i += 1

            if i > 800 and np.dot(g, pos) < 0:
                break

    return np.array(x_list) / sigma, np.array(y_list) / sigma

def get_chi_from_path(x, y):
    g_in = np.array([x[1], y[1]]) - np.array([x[0], y[0]])
    g_out = np.array([x[-1], y[-1]]) - np.array([x[-2], y[-2]])
    chi = np.arccos(np.dot(g_in, g_out) / (vec_len(g_in) * vec_len(g_out)))
    if g_out[0] < 0:
        chi = - chi
    return chi

import matplotlib.pyplot as plt

def collision():
    g_list = [1.5]
    b_list = [0.9]
    failed = False
    rval = 0
    for bi in b_list:
        for gi in g_list:
            x, y = get_path(T, gi, bi)
            chi_path = get_chi_from_path(x, y)
            chi_pred = kin.cpp_kingas.chi(0, 0, T, gi, bi * sigma)
            if abs(chi_path - chi_pred) > 0.06:
                failed = True
                rval = chi_path - chi_pred
                break

        if failed is True:
            break

    if failed is True:
        return 130, rval
    return 0, 0

def run_tests(do_print=False, do_plot=False):
    '''
        Submodule for testing functions used in computing collision integrals
        Each test in 'tests' must accept two arguments: 'do_plot' and 'do_print' and return two values
        The first is the exit status of the test (0 for successfull, !0 otherwise)
        The second value is some information about the test that failed
    '''
    tests = [theta_lim, chi_lim, collision]
    if do_plot:
        print('Plotting of mie tests is not implemented!')
    r = 0
    for t in tests:
        with suppress_stdout('-silent' in sys.argv):
            r, v = t()
        if r != 0:
            r += 2000
            if do_print:
                print(r, v)
            print(f'{bcolors.FAIL}Mie tests failed with exit code :', r, f'{bcolors.ENDC}')
            break
    if r == 0:
        print(f'{bcolors.OKGREEN}Mie tests were successful!{bcolors.ENDC}')
    return r

if __name__ == '__main__':
    run_tests(do_print=('-print' in sys.argv), do_plot=('-plot' in sys.argv))