"""
Module to ensure that swapping component order in binary mixtures does not change output
Essentially computes a bunch of properties, flips the component order, and recomputes, to check that output is the same.
"""
from tools import check_eq, check_eq_arr, check_eq_lst, Components, report, models
import pytest

@pytest.mark.parametrize('model', models)
@pytest.mark.parametrize('N', [1, 2, 3])
@pytest.mark.parametrize('FoR', ['CoN', 'CoM', 'CoV', 'solvent', 'zarate', 'zarate_x', 'zarate_w'])
def test_diffusion(model, N, FoR, silent=True):
    T = 300
    Vm = 1 / 20
    x_lst = [0.2, 0.8]

    comps = Components('AR,KR', x_lst)
    D12_lst = [None, None]
    for i in range(2):
        x = comps.get_x(i)
        kin = model(comps.get_comps(i))
        D12_lst[i] = kin.interdiffusion(T, Vm, x, N=N, frame_of_reference=FoR, dependent_idx=i, solvent_idx=i,
                                        use_independent=True, use_binary=True)

    assert check_eq_lst(D12_lst)
    if silent is False:
        report(check_eq_lst(D12_lst), [model, N, FoR, D12_lst])

@pytest.mark.parametrize('model', models)
@pytest.mark.parametrize('N', [2, 3])
@pytest.mark.parametrize('FoR', ['CoN', 'CoM', 'CoV', 'solvent', 'zarate'])
def test_th_diffusion(model, N, FoR, silent=True):
    T = 300
    Vm = 1 / 20
    x_lst = [0.2, 0.8]

    comps = Components('AR,KR', x_lst)
    DT_lst = [None, None]

    for i in range(2):
        x = comps.get_x(i)
        kin = model(comps.get_comps(i))
        DT = kin.thermal_diffusion_coeff(T, Vm, x, N=N, frame_of_reference=FoR, dependent_idx=i, solvent_idx=i,
                                                 use_independent=True)
        if FoR != 'zarate':
            DT_lst[i] = DT[i]
        else:
            DT_lst[i] = DT[0]

    assert check_eq_lst(DT_lst)
    if silent is False:
        report(check_eq_lst(DT_lst), [model, N, FoR, DT_lst])

@pytest.mark.parametrize('model', models)
@pytest.mark.parametrize('N', [2, 3])
def test_conductivity(model, N, silent=True):
    T = 300
    Vm = 1 / 20
    x_lst = [0.2, 0.8]

    comps = Components('AR,KR', x_lst)
    cond_lst = [None, None]
    for i in range(2):
        x = comps.get_x(i)
        kin = model(comps.get_comps(i))
        cond_lst[i] = kin.thermal_conductivity(T, Vm, x, N=N)

    assert check_eq_lst(cond_lst)
    if silent is False:
        report(check_eq_lst(cond_lst), [model, N, cond_lst])

@pytest.mark.parametrize('model', models)
@pytest.mark.parametrize('N', [2, 3])
def test_viscosity(model, N, silent=True):
    T = 300
    Vm = 1 / 20
    x_lst = [0.2, 0.8]

    comps = Components('AR,KR', x_lst)
    visc_lst = [None, None]

    for i in range(2):
        x = comps.get_x(i)
        kin = model(comps.get_comps(i))
        visc_lst[i] = kin.viscosity(T, Vm, x, N=N)

    assert check_eq_lst(visc_lst)
    if silent is False:
        report(check_eq_lst(visc_lst), [model, N, visc_lst])

@pytest.mark.parametrize('model', models)
@pytest.mark.parametrize('N', [2, 3])
def test_thermal_diffusion_factor(model, N, silent=True):
    T = 300
    Vm = 1 / 20
    x_lst = [0.2, 0.8]

    comps = Components('AR,KR', x_lst)
    alpha_lst = [None, None]
    for i in range(2):
        x = comps.get_x(i)
        kin = model(comps.get_comps(i))
        alpha_lst[i] = kin.thermal_diffusion_factor(T, Vm, x, N=N)[i][1 - i]

    assert check_eq_lst(alpha_lst)
    if silent is False:
        report(check_eq_lst(alpha_lst), [model, N, alpha_lst])

@pytest.mark.parametrize('model', models)
@pytest.mark.parametrize('N', [2, 3])
def test_th_diffusion_ratio(model, N, silent=True):
    T = 300
    Vm = 1 / 20
    x_lst = [0.2, 0.8]

    comps = Components('AR,KR', x_lst)
    kT_lst = [None, None]
    for i in range(2):
        x = comps.get_x(i)
        kin = model(comps.get_comps(i))
        kT_lst[i] = comps.rollback_lst(kin.thermal_diffusion_ratio(T, Vm, x, N=N), i)

    assert check_eq_arr(kT_lst[0], kT_lst[1])
    if silent is False:
        report(check_eq_lst(kT_lst), [model, N, kT_lst])

@pytest.mark.parametrize('model', models)
@pytest.mark.parametrize('N', [2, 3])
def _test_soret(model, N, silent=True):
    T = 300
    Vm = 1 / 20
    x_lst = [0.2, 0.8]

    comps = Components('AR,KR', x_lst)
    ST_lst = [None, None]
    for i in range(2):
        x = comps.get_x(i)
        kin = model(comps.get_comps(i))
        ST_lst[i] = kin.soret_coefficient(T, Vm, x, N=N, use_zarate=True)[0]

    assert check_eq_lst(ST_lst)
    if silent is False:
        report(check_eq_lst(ST_lst), [model, N, ST_lst])