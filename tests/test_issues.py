"""
Test that various criteria that have been adressed in previous issues or PR's are still fulfilled.

Each test here should have a reference to the related issue or PR.
"""
from pykingas.MieKinGas import MieKinGas
from scipy.constants import Avogadro, Boltzmann
import numpy as np
from tools import check_eq, check_eq_lst, check_eq_arr, check_eq_rel
import pytest

@pytest.mark.parametrize('lr', [20, 30, 40, 50])
def test_very_repulse(lr):
    """
    Issue #25 (https://github.com/thermotools/KineticGas/issues/25)
    Check that viscosity and thermal conductivity run for highly repulsive systems.
    """
    m = 10 # Molar mass (g / mol)
    mie = MieKinGas('LJF', mole_weights=[m, m], lr=[lr, lr])
    unt = mie.get_reducing_units()

    test_visc_vals = {20: 0.2158860328573334, 
                      30: 0.21695389791954947, 
                      40: 0.21751138069457196, 
                      50: 0.2178548098137264}

    test_cond_vals = {20: 1.0020395956063672, 
                      30: 0.9989005623538744, 
                      40: 0.9973621421137406, 
                      50: 0.9964959608177728}

    T_red = 2.0
    rho_red = 0.1
    T = T_red * unt.T
    rho = rho_red * unt.rho

    visc = mie.viscosity(T, 1 / rho, [0.5, 0.5], N=2) / unt.visc
    cond = mie.thermal_conductivity(T, 1 / rho, [0.5, 0.5], N=2) / unt.tcond

    assert check_eq(visc, test_visc_vals[lr], tol=1e-4)
    assert check_eq(cond, test_cond_vals[lr], tol=1e-4)

@pytest.mark.parametrize('m', [10, 30])
@pytest.mark.parametrize('sigma', [2.9e-10, 3.5e-10])
@pytest.mark.parametrize('eps_div_k', [100, 220])
@pytest.mark.parametrize('la', [6, 8])
@pytest.mark.parametrize('lr', [12, 20])
def test_singlecomp_binary(m, sigma, eps_div_k, la, lr):
    """
    Test that initializing a single component, and a binary mixture of identical species gives the same output.
    Note: Only works if the "binary mixture of identical species" recieves the mole fractions [0.5, 0.5], because
            for single components, mole fraction inputs are converted to [0.5, 0.5] in order to support supplying [1]
            for single components.
    See: PR #27
    """
    kin1 = MieKinGas('LJF', mole_weights=[m, m], sigma=[sigma, sigma],
                     eps_div_k=[eps_div_k, eps_div_k], la=[la, la], lr=[lr, lr])
    kin2 = MieKinGas('AR,KR', mole_weights=[m, m],
                     sigma=[sigma, sigma], eps_div_k=[eps_div_k, eps_div_k], la=[la, la], lr=[lr, lr])

    print(kin1.sigma)
    print(kin2.sigma)
    print(kin1.epsilon)
    print(kin2.epsilon)

    T = 300
    Vm = 3e-4
    p = 10e5
    x1 = 0.3
    
    D1 = kin1.interdiffusion(T, Vm, [x1, 1 - x1], N=2)
    D2 = kin2.interdiffusion(T, Vm, [x1, 1 - x1], N=2)
    assert check_eq_rel(D1, D2)

    l1 = kin1.viscosity(T, Vm, [x1, 1 - x1], N=2)
    l2 = kin2.viscosity(T, Vm, [x1, 1 - x1], N=2)
    assert check_eq_rel(l1, l2)

    l1 = kin1.viscosity_tp(T, p, [x1, 1 - x1], N=2)
    l2 = kin2.viscosity_tp(T, p, [x1, 1 - x1], N=2)
    assert check_eq_rel(l1, l2)

    k1 = sum(kin1.thermal_conductivity(T, Vm, [x1, 1 - x1], N=2, contributions='td').values()) # Exclude internal (Cp-dependent) contribution
    k2 = sum(kin2.thermal_conductivity(T, Vm, [x1, 1 - x1], N=2, contributions='td').values()) # Exclude internal (Cp-dependent) contribution
    assert check_eq_rel(k1, k2)

    k1 = sum(kin1.thermal_conductivity_tp(T, p, [x1, 1 - x1], N=2, contributions='td').values()) # Exclude internal (Cp-dependent) contribution
    k2 = sum(kin2.thermal_conductivity_tp(T, p, [x1, 1 - x1], N=2, contributions='td').values()) # Exclude internal (Cp-dependent) contribution
    assert check_eq_rel(k1, k2)

    DT1 = kin1.thermal_diffusion_coeff(T, Vm, [0.5, 0.5], N=2)
    DT2 = kin2.thermal_diffusion_coeff(T, Vm, [0.5, 0.5], N=2)
    assert check_eq_arr(DT1, DT2)

if __name__ == '__main__':
    cond_dat, visc_dat = {}, {}
    for lr in [20, 30, 40, 50]:
        visc_dat[lr], cond_dat[lr] = test_very_repulse(lr)
    
    print(visc_dat)
    print(cond_dat)