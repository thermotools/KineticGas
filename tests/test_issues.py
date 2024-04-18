from pykingas.MieKinGas import MieKinGas
from scipy.constants import Avogadro, Boltzmann
import numpy as np
from tools import check_eq, check_eq_lst, check_eq_arr
import pytest

@pytest.mark.parametrize('lr', [20, 30, 40, 50])
def test_very_repulse(lr):
    """
    Issue #25 (https://github.com/thermotools/KineticGas/issues/25)
    Check that viscosity and thermal conductivity run for highly repulsive systems.
    """
    m = 10
    mie = MieKinGas('LJF', mole_weights=[m, m], lr=[lr, lr])

    sigma = mie.sigma[0][0]
    eps = mie.epsilon[0]
    eps_div_k = eps / Boltzmann

    test_visc_vals = {20 : 0.2136582471368343,
                        30 : 0.21469802292722762,
                        40 : 0.21527111014341763,
                        50 : 0.21563637353012097}

    test_cond_vals = {20 : 0.7991821474987899,
                        30 : 0.8022966936925214,
                        40 : 0.8039439619090943,
                        50 : 0.8049669080299001}

    T_red = 2.0
    rho_red = 0.1
    T = T_red * eps_div_k
    rho = rho_red * Avogadro * sigma**3

    visc_unit = np.sqrt(eps * (m * 1e-3 / Avogadro)) / sigma**2
    cond_unit = Boltzmann * np.sqrt(eps/ (m * 1e-3 / Avogadro)) / sigma**2

    visc = mie.viscosity(T, 1 / rho, [0.5, 0.5], N=2) / visc_unit
    cond = mie.thermal_conductivity(T, 1 / rho, [0.5, 0.5], N=2) / cond_unit

    assert check_eq(visc, test_visc_vals[lr])
    assert check_eq(cond, test_cond_vals[lr])

@pytest.mark.parametrize('m', [10, 30])
@pytest.mark.parametrize('sigma', [2.9e-10, 3.5e-10])
@pytest.mark.parametrize('eps_div_k', [100, 220])
@pytest.mark.parametrize('la', [6, 8])
@pytest.mark.parametrize('lr', [12, 20])
def test_singlecomp_binary(m, sigma, eps_div_k, la, lr):
    """
    Test that initializing a single component, and a binary mixture of identical species gives the same output.
    See: PR #27
    """
    kin1 = MieKinGas('AR', mole_weights=[m * Avogadro * 1e3, m * Avogadro * 1e3], sigma=[sigma, sigma],
                     eps_div_k=[eps_div_k, eps_div_k], la=[la, la], lr=[lr, lr])
    kin2 = MieKinGas('AR,KR', mole_weights=[m * Avogadro * 1e3, m * Avogadro * 1e3],
                     sigma=[sigma, sigma], eps_div_k=[eps_div_k, eps_div_k], la=[la, la], lr=[lr, lr])

    T = 300
    Vm = 3e-4
    p = 10e5
    for x1 in (0.2, 0.6, 0.9999):
        D1 = kin1.interdiffusion(T, Vm, [x1, 1 - x1], N=2)
        D2 = kin2.interdiffusion(T, Vm, [x1, 1 - x1], N=2)
        assert check_eq(D1, D2)

        l1 = kin1.viscosity(T, Vm, [x1, 1 - x1], N=2)
        l2 = kin2.viscosity(T, Vm, [x1, 1 - x1], N=2)
        assert check_eq(l1, l2)

        l1 = kin1.viscosity_tp(T, p, [x1, 1 - x1], N=2)
        l2 = kin2.viscosity_tp(T, p, [x1, 1 - x1], N=2)
        assert check_eq(l1, l2)

        k1 = kin1.thermal_conductivity(T, Vm, [x1, 1 - x1], N=2)
        k2 = kin2.thermal_conductivity(T, Vm, [x1, 1 - x1], N=2)
        assert check_eq(k1, k2)

        k1 = kin1.thermal_conductivity_tp(T, p, [x1, 1 - x1], N=2)
        k2 = kin2.thermal_conductivity_tp(T, p, [x1, 1 - x1], N=2)
        assert check_eq(k1, k2)

        a1 = kin1.thermal_diffusion_factor(T, Vm, [x1, 1 - x1], N=2)
        a2 = kin2.thermal_diffusion_factor(T, Vm, [x1, 1 - x1], N=2)
        assert check_eq_arr(a1, a2)

        DT1 = kin1.thermal_diffusion_coeff(T, Vm, [x1, 1 - x1], N=2)
        DT2 = kin2.thermal_diffusion_coeff(T, Vm, [x1, 1 - x1], N=2)
        assert check_eq_arr(DT1, DT2)

        kT1 = kin1.thermal_diffusion_ratio(T, Vm, [x1, 1 - x1], N=2)
        kT2 = kin2.thermal_diffusion_ratio(T, Vm, [x1, 1 - x1], N=2)
        assert check_eq_arr(kT1, kT2)