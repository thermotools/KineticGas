import copy
import numpy as np
from pykingas.MieKinGas import MieKinGas
from pykingas.HardSphere import HardSphere
from pykingas.Sutherland import S_MieKinGas
from pykingas.QuantumMie import QuantumMie

FLTEPS = 1e-10

models = [MieKinGas, HardSphere, S_MieKinGas, QuantumMie]

def check_eq(a, b, eps=FLTEPS):
    if abs(a - b) > eps:
        return False
    return True

def check_eq_arr(a, b):
    if any(abs(a.flatten() - b.flatten()) > FLTEPS):
        return False
    return True

def check_eq_lst(lst):
    if any(abs(lst - lst[0]) > FLTEPS):
        return False
    return True

def report(passed, info):
    if passed is False:
        passed = 'FAILED'
    else:
        passed = 'passed'

    print(f'{passed} : ' + ', '.join(str(i) for i in info))

class Components:

    def __init__(self, comps, x):
        self.comps = comps.split(',')
        self.x = copy.deepcopy(x)
        self.ncomps = len(self.comps)

    def get_comps(self, roll):
        return ','.join([self.comps[i - roll] for i in range(self.ncomps)])

    def get_x(self, roll):
        return np.array([self.x[i - roll] for i in range(self.ncomps)])

    def rollback_lst(self, lst, rolled):
        return np.array([lst[i + rolled if (i + rolled < self.ncomps) else i + rolled - self.ncomps] for i in range(self.ncomps)])