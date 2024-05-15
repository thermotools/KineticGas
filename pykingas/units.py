"""
The dumb struct Units, holder for the reducing scale of common variables, for easy construction
from a given model.

Usage:

    model = MieKinGas('AR')
    unt = model.get_reducing_units()
    T_red = np.linspace(1, 3)
    T = T_red * unt.T # T in kelvin
    rho = 0.3 * unt.rho # density in mol / m3
    D = model.interdiffusion(T, 1 / rho, [1]) # D in m2 / s
    D_red = D / unt.D

    etc.
"""
from scipy.constants import Boltzmann, Avogadro
import numpy as np

class Units:
    def __init__(self, m_unit, L_unit, T_unit):
        self.T = T_unit
        self.E = self.T * Boltzmann
        self.L = L_unit
        self.m = m_unit
        self.V = self.L ** 3
        self._time_unit = np.sqrt(self.m * self.L ** 2 / self.E)
        self.rho = 1 / (Avogadro * self.V)
        self.D = L_unit ** 2 / self._time_unit
        self.p = self.E / self.V
        self.visc = self.p * self._time_unit
        self.tcond = self.E / (self._time_unit * self.T * self.L)
