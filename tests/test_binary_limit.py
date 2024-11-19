"""
Module to test that ternary mixtures reduce to the correct binary mixtures in the binary limit.
"""
import copy
import numpy as np
from tools import models
import pytest
import matplotlib.pyplot as plt

def get_ternary_x(x_bin, x3):
    return np.array([(1 - x3) * x_bin[0], (1 - x3) * x_bin[1], x3])

@pytest.mark.parametrize('model', models)
@pytest.mark.parametrize('N', [1, 2])
def test_viscosity(model, N, silent=True):
    T = 300
    Vm = 1 / 10
    x_bin = [0.3, 0.7]
    binary = model('HE,C1')
    ternary = model('HE,C1,O2')

    eta_b = binary.viscosity(T, Vm, x_bin, N=N)
    eta_t = ternary.viscosity(T, Vm, get_ternary_x(x_bin, 0.5), N=N)
    diff = abs(eta_t - eta_b)
    x3_lst = np.logspace(-1, -6, 10)

    if silent is False:
        plt.plot(x3_lst, [ternary.viscosity(T, Vm, get_ternary_x(x_bin, x3), N=N) for x3 in x3_lst])
        plt.plot(x3_lst, eta_b * np.ones_like(x3_lst))
        plt.show()

    for x3 in x3_lst:
        eta_t = ternary.viscosity(T, Vm, get_ternary_x(x_bin, x3), N=N)
        new_diff = abs(eta_t - eta_b)
        assert new_diff < diff
        diff = new_diff

    assert (diff / eta_b) < 1e-3

@pytest.mark.parametrize('model', models)
@pytest.mark.parametrize('N', [2, 3])
def test_thermal_conductivity(model, N, silent=True):
    T = 300
    Vm = 1 / 10
    x_bin = [0.3, 0.7]
    binary = model('HE,C1')
    ternary = model('HE,C1,O2')

    cond_b = binary.thermal_conductivity(T, Vm, x_bin, N=N)
    cond_t = ternary.thermal_conductivity(T, Vm, get_ternary_x(x_bin, 0.5), N=N)
    diff = abs(cond_t - cond_b)
    x3_lst = np.logspace(-1, -6, 10)

    if silent is False:
        plt.plot(x3_lst, [ternary.thermal_conductivity(T, Vm, get_ternary_x(x_bin, x3), N=N) for x3 in x3_lst])
        plt.plot(x3_lst, cond_b * np.ones_like(x3_lst))
        plt.show()

    for x3 in x3_lst:
        cond_t = ternary.thermal_conductivity(T, Vm, get_ternary_x(x_bin, x3), N=N)
        new_diff = abs(cond_t - cond_b)
        assert new_diff < diff
        diff = new_diff

@pytest.mark.parametrize('model', models)
@pytest.mark.parametrize('N', [1, 2])
def test_diffusion(model, N, silent=True):
    T = 300
    Vm = 1 / 10
    x_bin = [0.3, 0.7]
    binary = model('HE,C1')
    ternary = model('HE,C1,O2')

    D_b = binary.interdiffusion(T, Vm, x_bin, N=N, dependent_idx=1)
    D_t = ternary.interdiffusion(T, Vm, get_ternary_x(x_bin, 0.5), N=N, dependent_idx=1)
    diff = abs(D_t[0][0] - D_b)
    x3_lst = np.logspace(-1, -5, 10)

    if silent is False:
        plt.plot(x3_lst, [ternary.interdiffusion(T, Vm, get_ternary_x(x_bin, x3), N=N, dependent_idx=1)[0][0] for x3 in x3_lst])
        plt.plot(x3_lst, D_b * np.ones_like(x3_lst))
        plt.show()

    for x3 in x3_lst:
        D_t = ternary.interdiffusion(T, Vm, get_ternary_x(x_bin, x3), N=N, dependent_idx=1)
        new_diff = abs(D_t[0][0] - D_b)
        assert new_diff < diff
        diff = new_diff

    assert (diff / D_b) < 1e-3


@pytest.mark.parametrize('model', models)
@pytest.mark.parametrize('N', [2, 3])
def test_th_diffusion(model, N, silent=True):
    T = 300
    Vm = 1 / 10
    x_bin = [0.3, 0.7]
    binary = model('HE,C1')
    ternary = model('HE,C1,O2')

    DT_b = binary.thermal_diffusion_coeff(T, Vm, x_bin, N=N, dependent_idx=1)
    DT_t = ternary.thermal_diffusion_coeff(T, Vm, get_ternary_x(x_bin, 0.5), N=N, dependent_idx=1)
    diff = abs(DT_t[:-1] - DT_b)
    x3_lst = np.logspace(-1, -5, 10)

    if silent is False:
        plt.plot(x3_lst, [ternary.thermal_diffusion_coeff(T, Vm, get_ternary_x(x_bin, x3), N=N, dependent_idx=1)[:-1] for x3 in x3_lst])
        plt.plot(x3_lst, DT_b * np.ones_like(x3_lst))
        plt.show()

    for x3 in x3_lst:
        DT_t = ternary.thermal_diffusion_coeff(T, Vm, get_ternary_x(x_bin, x3), N=N, dependent_idx=1)
        new_diff = abs(DT_t[:-1] - DT_b)
        assert all(new_diff < diff)
        diff = copy.deepcopy(new_diff)

    assert all((diff / DT_b) < 1e-3)

@pytest.mark.parametrize('model', models)
@pytest.mark.parametrize('N', [2, 3])
def test_soret(model, N, silent=True):
    T = 300
    Vm = 1 / 10
    x_bin = [0.3, 0.7]
    bin_12 = model('H2O,O2', is_idealgas=True)
    ternary = model('H2O,O2,N2', is_idealgas=True)

    x_bin = [0.6, 0.4]
    ST_1b2 = bin_12.soret_coefficient(T, Vm, x_bin, N=N, use_zarate=True)[0]
    ST_1t, ST_2t = ternary.soret_coefficient(T, Vm, get_ternary_x(x_bin, 0.5), N=N, use_zarate=True)
    diff = abs((ST_1t - ST_2t) - ST_1b2)
    x3_lst = np.logspace(-1, -5, 10)

    print(ST_1t - ST_2t, np.diff([ST_1t, ST_2t]), ST_1b2)

    if silent is False:
        plt.plot(x3_lst, [- np.diff(ternary.soret_coefficient(T, Vm, get_ternary_x(x_bin, x3), N=N, use_zarate=True)) for x3 in x3_lst])
        plt.plot(x3_lst, ST_1b2 * np.ones_like(x3_lst))
        plt.xscale('log')
        plt.show()

    for x3 in x3_lst:
        ST_1t, ST_2t = ternary.soret_coefficient(T, Vm, get_ternary_x(x_bin, x3), N=N, use_zarate=True)
        new_diff = abs((ST_1t - ST_2t) - ST_1b2)
        assert new_diff < diff
        diff = new_diff
    
    bin_23 = model('O2,N2', is_idealgas=True)
    def get_ternary_x1(x1):
        return [x1, x_bin[0] * (1 - x1), x_bin[1] * (1 - x1)]

    ST_2b3 = bin_23.soret_coefficient(T, Vm, x_bin, N=N, use_zarate=True)[0]
    _, ST_2t = ternary.soret_coefficient(T, Vm, get_ternary_x1(0.11), N=N, use_zarate=True)
    diff = abs(ST_2t - ST_2b3)
    x1_lst = np.logspace(-1, -5, 10)
    for x1 in x1_lst:
        _, ST_2t = ternary.soret_coefficient(T, Vm, get_ternary_x1(x1), N=N, use_zarate=True)
        new_diff = abs(ST_2t - ST_2b3)
        assert new_diff < diff
        diff = new_diff

    bin_13 = model('H2O,N2', is_idealgas=True)
    def get_ternary_x2(x2):
        return [x_bin[0] * (1 - x2), x2, x_bin[1] * (1 - x2)]
    
    ST_1b3 = bin_13.soret_coefficient(T, Vm, x_bin, N=N, use_zarate=True)[0]
    ST_1t, _ = ternary.soret_coefficient(T, Vm, get_ternary_x2(0.5), N=N, use_zarate=True)
    diff = abs(ST_1t - ST_1b3)
    x2_lst = np.logspace(-1, -5, 10)
    for x2 in x2_lst:
        ST_1t, _ = ternary.soret_coefficient(T, Vm, get_ternary_x2(x2), N=N, use_zarate=True)
        new_diff = abs(ST_1t - ST_1b3)
        assert new_diff < diff
        diff = new_diff


if __name__ == "__main__":
    test_soret(models[0], 2, silent=True)