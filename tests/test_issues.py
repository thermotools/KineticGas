from pykingas.MieKinGas import MieKinGas
from scipy.constants import Avogadro, Boltzmann
import numpy as np
from tools import check_eq
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