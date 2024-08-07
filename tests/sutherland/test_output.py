import numpy as np
import pandas as pd
from pykingas.QuantumMie import QuantumMie
from scipy.constants import Avogadro
import os
import pytest


rdf_Tlst = np.linspace(0.8, 3.5, 10)
rdf_rholst = np.linspace(1e-6, 0.7, 10)

singledata = pd.read_csv(f'{os.path.dirname(__file__)}/../data/rdf_single.csv')
binarydata = pd.read_csv(f'{os.path.dirname(__file__)}/../data/rdf_binary.csv')

def gen_rdf_single():
    out = {'T' : [], 'rho' : [], 'g0' : [], 'g1' : [], 'g2' : [], 'g2c' : []}

    kin = QuantumMie('H2', FH_orders=2)
    unt = kin.get_reducing_units()
    for Tr in rdf_Tlst:
        for rho_r in rdf_rholst:
            out['T'].append(Tr)
            out['rho'].append(rho_r)
            g0, g1, g2, g2c = kin.get_rdf_terms(rho_r * unt.rho * Avogadro, Tr * unt.T, [0.3, 0.7])
            out['g0'].append(g0[0][0])
            out['g1'].append(g1[0][0])
            out['g2'].append(g2[0][0])
            out['g2c'].append(g2c[0][0])
    
    pd.DataFrame(out).to_csv(f'{os.path.dirname(__file__)}/../data/rdf_single.csv')

def gen_rdf_binary():
    out = {'T' : [], 'rho' : []
           , 'g0_11' : [], 'g1_11' : [], 'g2_11' : [], 'g2c_11' : []
           , 'g0_12' : [], 'g1_12' : [], 'g2_12' : [], 'g2c_12' : []
           , 'g0_22' : [], 'g1_22' : [], 'g2_22' : [], 'g2c_22' : []}

    kin = QuantumMie('H2,NE', FH_orders=2)
    unt = kin.get_reducing_units()
    for Tr in rdf_Tlst:
        for rho_r in rdf_rholst:
            out['T'].append(Tr)
            out['rho'].append(rho_r)
            g0, g1, g2, g2c = kin.get_rdf_terms(rho_r * unt.rho * Avogadro, Tr * unt.T, [0.3, 0.7])
            out['g0_11'].append(g0[0][0])
            out['g1_11'].append(g1[0][0])
            out['g2_11'].append(g2[0][0])
            out['g2c_11'].append(g2c[0][0])
            out['g0_12'].append(g0[0][1])
            out['g1_12'].append(g1[0][1])
            out['g2_12'].append(g2[0][1])
            out['g2c_12'].append(g2c[0][1])
            out['g0_22'].append(g0[1][1])
            out['g1_22'].append(g1[1][1])
            out['g2_22'].append(g2[1][1])
            out['g2c_22'].append(g2c[1][1])
    
    pd.DataFrame(out).to_csv(f'{os.path.dirname(__file__)}/../data/rdf_binary.csv')

@pytest.mark.parametrize('rho_r', rdf_rholst)
def test_rdf_single(rho_r):
    kin = QuantumMie('H2', FH_orders=2)
    unt = kin.get_reducing_units()
    
    for Tr in sorted(set(singledata['T'])):
        data = singledata[(abs(singledata['T'] - Tr) < 1e-2) & (abs(singledata['rho'] - rho_r) < 1e-3)]
        g0, g1, g2, g2c = kin.get_rdf_terms(rho_r * unt.rho * Avogadro, Tr * unt.T, [0.3, 0.7])
        assert abs(g0[0][0] - data['g0'].iloc[0]) < 1e-10
        assert abs(g1[0][0] - data['g1'].iloc[0]) < 1e-10
        assert abs(g2[0][0] - data['g2'].iloc[0]) < 1e-10
        assert abs(g2c[0][0] - data['g2c'].iloc[0]) < 1e-10

@pytest.mark.parametrize('rho_r', rdf_rholst)
def test_rdf_binary(rho_r):
    kin = QuantumMie('H2,NE', FH_orders=2)
    unt = kin.get_reducing_units()
    
    for Tr in sorted(set(binarydata['T'])):
        data = binarydata[(abs(binarydata['T'] - Tr) < 1e-2) & (abs(binarydata['rho'] - rho_r) < 1e-3)]
        g0, g1, g2, g2c = kin.get_rdf_terms(rho_r * unt.rho * Avogadro, Tr * unt.T, [0.3, 0.7])
        assert abs(g0[0][0] - data['g0_11'].iloc[0]) < 1e-10
        assert abs(g1[0][0] - data['g1_11'].iloc[0]) < 1e-10
        assert abs(g2[0][0] - data['g2_11'].iloc[0]) < 1e-10
        assert abs(g2c[0][0] - data['g2c_11'].iloc[0]) < 1e-10
        assert abs(g0[0][1] - data['g0_12'].iloc[0]) < 1e-10
        assert abs(g1[0][1] - data['g1_12'].iloc[0]) < 1e-10
        assert abs(g2[0][1] - data['g2_12'].iloc[0]) < 1e-10
        assert abs(g2c[0][1] - data['g2c_12'].iloc[0]) < 1e-10
        assert abs(g0[1][1] - data['g0_22'].iloc[0]) < 1e-10
        assert abs(g1[1][1] - data['g1_22'].iloc[0]) < 1e-10
        assert abs(g2[1][1] - data['g2_22'].iloc[0]) < 1e-10
        assert abs(g2c[1][1] - data['g2c_22'].iloc[0]) < 1e-10
        assert g0[0][1] == g0[1][0]
        assert g1[0][1] == g1[1][0]
        assert g2[0][1] == g2[1][0]
        assert g2c[0][1] == g2c[1][0]

# gen_rdf_single()
# gen_rdf_binary()
for rho_r in rdf_rholst:    
    test_rdf_binary(rho_r)
    test_rdf_single(rho_r)