import pandas as pd
import numpy as np
import os, ast
from scipy.constants import Avogadro
from tools import models, check_eq_arr, check_eq_rel
from pykingas.py_KineticGas import py_KineticGas
from pykingas.LJSpline import LJSpline
import pytest

@pytest.mark.parametrize('__model', models)
def test_diffusion(__model, overwrite=False):
    data = pd.read_csv(f'{os.path.dirname(__file__)}/data/diffusion_{__model.__name__}.csv')
    data['x'] = data['x'].apply(lambda x: ast.literal_eval(x))
    data['D'] = data['D'].apply(lambda x: ast.literal_eval(x))
    D_new = [None for _ in range(len(data['x']))]
    mixtures = set(data['mix'])

    is_unchanged = True
    for mix in mixtures:
        mixdata = data[data['mix'] == mix]
        model = __model(mix)
        ncomps = len(mix.split(','))
        unt = model.get_reducing_units()
        for i, T, rho_red, x, dep in zip(mixdata.index, mixdata['T'], mixdata['rho*'], mixdata['x'], mixdata['dependent']):
            dep = int(dep) if (not np.isnan(dep)) else np.nan
            rho = rho_red * unt.rho
            if ncomps == 1:
                D = model.interdiffusion(T, 1 / rho, x, N=2, use_binary=False)
            else:
                D = model.interdiffusion(T, 1 / rho, x, N=2, dependent_idx=dep, use_binary=False)
            
            if (overwrite is False) and (check_eq_arr(D, np.array(data['D'][i])) is False):
                is_unchanged = False
                print('Failure at : ', mix, T, rho_red, x, dep, rho, ': \n\t', data['D'][i], '\n\t', D.tolist(), '\n\t', (np.array(data['D'][i]) - D).tolist())
            else:
                D_new[i] = D.tolist()
        print(f'Finished {mix}')
    
    if overwrite is True:
        data['D'] = pd.Series(D_new)
        data.to_csv(f'{os.path.dirname(__file__)}/data/diffusion_{__model.__name__}.csv')
    else:
        assert is_unchanged

def __test_scalar(__model, method, overwrite=False):
    data = pd.read_csv(f'{os.path.dirname(__file__)}/data/{method.__name__}_{__model.__name__}.csv')
    data['x'] = data['x'].apply(lambda x: ast.literal_eval(x))
    val_new = [np.nan for _ in range(len(data))]
    mixtures = set(data['mix'])

    is_unchanged = True
    for mix in mixtures:
        mixdata = data[data['mix'] == mix]
        model = __model(mix)
        unt = model.get_reducing_units()
        for i, T, rho_red, x, oldval in zip(mixdata.index, mixdata['T'], mixdata['rho*'], mixdata['x'], mixdata['value']):
            rho = rho_red * unt.rho
            val = method(model, T, 1 / rho, x, N=2)
            
            if (overwrite is False) and (check_eq_rel(val, oldval) is False):
                is_unchanged = False
                print('Failure at : ', mix, T, rho_red, x, rho, ': \n\t', oldval, ' => ', val, ', Relative error : ', (val - oldval) / oldval)
            else:
                val_new[i] = val
        print(f'Finished {mix}')

    if overwrite is True:
        data['value'] = pd.Series(val_new)
        data.to_csv(f'{os.path.dirname(__file__)}/data/{method.__name__}_{__model.__name__}.csv')
    else:
        assert is_unchanged

@pytest.mark.parametrize('__model', models)
def test_conductivity(__model, overwrite=False):
    method = py_KineticGas.thermal_conductivity
    __test_scalar(__model, method, overwrite=overwrite)

@pytest.mark.parametrize('__model', models)
def test_viscosity(__model, overwrite=False):
    method = py_KineticGas.viscosity
    __test_scalar(__model, method, overwrite=overwrite)

## LJ/Spline Tests ###

def test_ljs():
    sig = 3.42e-10
    eps_div_k = 124.0
    mw = 40.0
    ljs = LJSpline()
    data = pd.read_csv(f'{os.path.dirname(__file__)}/data/transport_properties_ljs.csv')
    unchanged = True
    for t, r, cond_old, visc_old, Drho_old in zip(data['T'], data['rho'], data['cond'], data['visc'], data['Drho']):
        T_SI = t*eps_div_k
        V_SI = Avogadro * sig**3 / r
        cond_new = ljs.thermal_conductivity(T_SI, V_SI, [0.5,0.5])
        visc_new = ljs.viscosity(T_SI, V_SI, [0.5,0.5])
        Drho_new = ljs.selfdiffusion(T_SI, V_SI) / V_SI 
        if not (check_eq_rel(cond_new, cond_old) and check_eq_rel(visc_new, visc_old) and check_eq_rel(Drho_new, Drho_old)):
            unchanged = False
            print('Failure at : T =', t, ' rho =', r)
            print(cond_new, cond_old)
            print(visc_new, visc_old)
            print(Drho_new, Drho_old)
    assert unchanged

#############################################################################################################################
#############      METHODS BELOW HAVE BEEN USED TO GENERATE THE TEST GRIDS BUT SHOULD NOT BE NEEDED IN THE FUTURE ###########
#############      TO RE-GENERATE TEST-GRIDS, CALL THE TEST FUNCTIONS WITH `overwrite=True`                       ###########
#############################################################################################################################

def __compute_scalar(__model, method):
    data = {'mix' : [], 'T' : [], 'rho*' : [], 'x' : [], 'value' : []}

    def add_to_data(mix, T, rho_red, x, val):
        data['mix'].append(mix)
        data['T'].append(T)
        data['rho*'].append(rho_red)
        data['x'].append(x)
        data['value'].append(val)

    for mix in mixtures:
        model = __model(mix)
        ncomps = len(mix.split(','))
        unt = model.get_reducing_units()
        for T in T_lst:
            for rho_red in rho_red_lst:
                rho = rho_red * unt.rho
                if ncomps == 1:
                    x = [1.]
                    
                    val = method(model, T, 1 / rho, x, N=2)
                    add_to_data(mix, T, rho_red, x, val)
                elif ncomps == 2:
                    for x1 in x1_lst:
                        x = [x1, 1 - x1]
                        val = method(model, T, 1 / rho, x, N=2)
                        add_to_data(mix, T, rho_red, x, val)
                else:
                    for x1 in x1_lst:
                        for r in x2_x3_ratios:
                            x = [x1, r * (1 - x1) / (r + 1), (1 - x1) / (r + 1)]
                            val = method(model, T, 1 / rho, x, N=2)
                            add_to_data(mix, T, rho_red, x, val)
        print(f'Finished {mix}')
    
    df = pd.DataFrame(data)
    df.to_csv(f'{os.path.dirname(__file__)}/data/{method.__name__}_{__model.__name__}.csv')

def __compute_tcond(model):
    method = py_KineticGas.thermal_conductivity
    __compute_scalar(model, method)

def __compute_visc(model):
    method = py_KineticGas.viscosity
    __compute_scalar(model, method)


def __compute_diffusion(__model):
    data = {'mix' : [], 'T' : [], 'rho*' : [], 'x' : [], 'D' : [], 'dependent' : []}

    def add_to_data(mix, T, rho_red, x, D, dep):
        data['mix'].append(mix)
        data['T'].append(T)
        data['rho*'].append(rho_red)
        data['x'].append(x)
        data['D'].append(D.tolist())
        data['dependent'].append(dep)

    for mix in mixtures:
        model = __model(mix)
        ncomps = len(mix.split(','))
        unt = model.get_reducing_units()
        for T in T_lst:
            for rho_red in rho_red_lst:
                rho = rho_red * unt.rho
                if ncomps == 1:
                    x = [1.]
                    D = model.interdiffusion(T, 1 / rho, x, N=2, use_binary=False)
                    add_to_data(mix, T, rho_red, x, D, np.nan)
                elif ncomps == 2:
                    for x1 in x1_lst:
                        x = [x1, 1 - x1]
                        for dep in range(2):
                            D = model.interdiffusion(T, 1 / rho, x, N=2, dependent_idx=dep, use_binary=False)
                            add_to_data(mix, T, rho_red, x, D, dep)
                else:
                    for x1 in x1_lst:
                        for r in x2_x3_ratios:
                            x = [x1, r * (1 - x1) / (r + 1), (1 - x1) / (r + 1)]
                            for dep in range(3):
                                D = model.interdiffusion(T, 1 / rho, x, N=2, dependent_idx=dep)
                                add_to_data(mix, T, rho_red, x, D, dep)
        print(f'Finished {mix}')
    
    df = pd.DataFrame(data)
    df.to_csv(f'{os.path.dirname(__file__)}/data/diffusion_{__model.__name__}.csv')


if __name__ == '__main__':

    # test_conductivity(models[1], overwrite=False)
    # test_viscosity(models[1], overwrite=False)
    # test_diffusion(models[1], overwrite=False)

    mixtures = ['AR', 'H2,O2', 'CO2,C1,N2']
    T_lst = np.linspace(150, 800, 5)
    rho_red_lst = np.linspace(0.01, 0.7, 10)
    x1_lst = np.linspace(1e-3, 1 - 1e-3, 5)
    x2_x3_ratios = np.linspace(0.1, 0.9, 3)

    __compute_diffusion(models[0])
    print('Finished Diffusion')
    __compute_tcond(models[0])
    print('Finished Tcond')
    __compute_visc(models[0])
    print('Finished Visc')
    
    exit(0)