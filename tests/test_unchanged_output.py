import pandas as pd
import numpy as np
import os, ast
from argparse import ArgumentParser
from scipy.constants import Avogadro
from tools import models, check_eq_arr, check_eq_rel
from pykingas.py_KineticGas import py_KineticGas
from pykingas.LJSpline import LJSpline
import pytest
import time

@pytest.mark.parametrize('__model', models)
def test_diffusion(__model, overwrite=False, fail_fast=True):
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
                ostr = f'Failure at : {mix}, {T}, {rho_red}, {x}, {dep}, {rho}\n\t{data['D'][i]}\n\t{D.tolist()}\n\t{(np.array(data['D'][i]) - D).tolist()}'
                if fail_fast is True:
                    assert is_unchanged, ostr
                else:
                    print(ostr)
            else:
                D_new[i] = D.tolist()
        print(f'Finished {mix}')
    
    if overwrite is True:
        data['D'] = pd.Series(D_new)
        data.to_csv(f'{os.path.dirname(__file__)}/data/diffusion_{__model.__name__}.csv')
    else:
        assert is_unchanged

def __test_scalar(__model, method, overwrite=False, fail_fast=True):
    data = pd.read_csv(f'{os.path.dirname(__file__)}/data/{method.__name__}_{__model.__name__}.csv')
    data['x'] = data['x'].apply(lambda x: eval(x))
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
                ostr = f'Failure at : {mix}, {T}, {rho_red}, {x}, {rho}\n\tChange: {oldval} => {val}, Relative error : {(val - oldval) / oldval}'
                if fail_fast:
                    assert False, ostr
                else:
                    rel_err = (val - oldval) / oldval
                    if abs(rel_err) > 5e-4: 
                        print(ostr)
            else:
                val_new[i] = val
        print(f'Finished {mix}')

    if overwrite is True:
        data['value'] = pd.Series(val_new)
        data.to_csv(f'{os.path.dirname(__file__)}/data/{method.__name__}_{__model.__name__}.csv')
    else:
        assert is_unchanged

@pytest.mark.parametrize('__model', models)
def test_conductivity(__model, overwrite=False, fail_fast=True):
    method = py_KineticGas.thermal_conductivity
    __test_scalar(__model, method, overwrite=overwrite, fail_fast=fail_fast)

@pytest.mark.parametrize('__model', models)
def test_viscosity(__model, overwrite=False, fail_fast=True):
    method = py_KineticGas.viscosity
    __test_scalar(__model, method, overwrite=overwrite, fail_fast=fail_fast)

## LJ/Spline Tests ###

def test_ljs(overwrite=False):
    sig = 3.42e-10
    eps_div_k = 124.0
    mw = 40.0
    ljs = LJSpline()
    data = pd.read_csv(f'{os.path.dirname(__file__)}/data/transport_properties_ljs.csv')
    unchanged = True
    for i, t, r, cond_old, visc_old, Drho_old in zip(data.index, data['T'], data['rho'], data['cond'], data['visc'], data['Drho']):
        T_SI = t*eps_div_k
        V_SI = Avogadro * sig**3 / r
        cond_new = ljs.thermal_conductivity(T_SI, V_SI, [0.5,0.5])
        visc_new = ljs.viscosity(T_SI, V_SI, [0.5,0.5])
        Drho_new = ljs.selfdiffusion(T_SI, V_SI) / V_SI 
        if overwrite is True:
            data.loc[i, 'cond'] = cond_new
            data.loc[i, 'visc'] = visc_new
            data.loc[i, 'Drho'] = Drho_new
            print(f'Changes at ({t:.1f}, {r:.3f}) : {(cond_new - cond_old) / cond_old:.1e}, {(visc_new - visc_old) / visc_old:.1e}, {(Drho_new - Drho_old) / Drho_old:.1e}')
            continue
        if not (check_eq_rel(cond_new, cond_old, tol=1e-6) and check_eq_rel(visc_new, visc_old, tol=1e-6) and check_eq_rel(Drho_new, Drho_old, tol=1e-6)):
            unchanged = False
            print('Failure at : T =', t, ' rho =', r)
            print(f'\tCond : {cond_old} => {cond_new} ({abs((cond_new - cond_old) / cond_old)})')
            print(f'\tVisc : {visc_old} => {visc_new} ({abs((visc_new - visc_old) / visc_old)})')
            print(f'\tDiff : {Drho_old} => {Drho_new} ({abs((Drho_new - Drho_old) / Drho_old)})')
    
    if overwrite is True:
        data.to_csv(f'{os.path.dirname(__file__)}/data/transport_properties_ljs.csv')
    assert unchanged

#############################################################################################################################
#############      METHODS BELOW HAVE BEEN USED TO GENERATE THE TEST GRIDS BUT SHOULD NOT BE NEEDED IN THE FUTURE ###########
#############      TO RE-GENERATE TEST-GRIDS, RUN THIS FILE (see bottom or run with --help)                       ###########
#############################################################################################################################

def __compute_scalar(__model, method):
    data = {'mix' : [], 'T' : [], 'rho*' : [], 'x' : [], 'value' : []}

    def add_to_data(mix, T, rho_red, x, val):
        data['mix'].append(mix)
        data['T'].append(T)
        data['rho*'].append(rho_red)
        data['x'].append(x)
        data['value'].append(val)
    
    def to_floats(x):
        return [float(xi) for xi in x]

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
                        x = to_floats([x1, 1 - x1]) # Prevent np.float64, because that will mess up our stored file
                        val = method(model, T, 1 / rho, x, N=2)
                        add_to_data(mix, T, rho_red, x, val)
                else:
                    for x1 in x1_lst:
                        for r in x2_x3_ratios:
                            x = to_floats([x1, r * (1 - x1) / (r + 1), (1 - x1) / (r + 1)]) # Prevent np.float64, because that will mess up our stored file
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

    parser = ArgumentParser()
    parser.add_argument('-c', '--cond', action='store_true')
    parser.add_argument('-v', '--visc', action='store_true')
    parser.add_argument('-d', '--diff', action='store_true')
    parser.add_argument('--ljs', action='store_true')
    parser.add_argument('-f', '--fail-fast', action='store_true')
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()

    test = not args.write
    USING_MODEL = models[0]
    print(f'Using model: {USING_MODEL.__name__}')
    time.sleep(1)
    if test:
        print(f'Running tests for model: {USING_MODEL.__name__}')
        if args.cond:
            test_conductivity(USING_MODEL, overwrite=False, fail_fast=args.fail_fast)
        if args.visc:
            test_viscosity(USING_MODEL, overwrite=False, fail_fast=args.fail_fast)
        if args.diff:
            test_diffusion(USING_MODEL, overwrite=False, fail_fast=args.fail_fast)
        if args.ljs:
            test_ljs(overwrite=False)

        exit(0)

    if args.write:
        print(f'Starting overwrite of results for model: {USING_MODEL.__name__}')

        mixtures = ['AR', 'H2,O2', 'CO2,C1,N2']
        T_lst = np.linspace(150, 800, 5)
        rho_red_lst = np.linspace(0.01, 0.7, 10)
        x1_lst = np.linspace(1e-3, 1 - 1e-3, 5)
        x2_x3_ratios = np.linspace(0.1, 0.9, 3)

        if args.diff:
            print('Start overwrite diffusion ...')
            __compute_diffusion(USING_MODEL)
            print('Finished Diffusion')
        if args.cond:
            print('Start overwrite thermal conductivity ...')
            __compute_tcond(USING_MODEL)
            print('Finished Tcond')
        if args.visc:
            print('Start overwrite viscosity ...')
            __compute_visc(USING_MODEL)
            print('Finished Visc')
        if args.ljs:
            print('Start overwrite LJs ...')
            test_ljs(overwrite=True)
            print('Finished LJs')
    
    exit(0)