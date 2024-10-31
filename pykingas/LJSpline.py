'''
Author: Johannes Salomonsen Løken
Purpose: Wrapper for the LJ/Spline class.
'''

import numpy as np
from numpy import pi
from scipy.constants import Avogadro, Boltzmann as kB
from scipy.optimize import root
from .libpykingas import cpp_LJSpline
from pykingas.py_KineticGas import py_KineticGas

class LJSpline(py_KineticGas):

    def __init__(self, sig, eps_div_k, mole_weights):
        if not type(mole_weights) is float:
            raise TypeError("MW must be a float!")
        mole_weights = [mole_weights,mole_weights]
        super().__init__("AR", mole_weights=mole_weights, is_idealgas=False, is_single_component=True)
        if not type(sig) is float:
            raise TypeError("Sigma must be a float!")
        if not type(eps_div_k) is float:
            raise TypeError("Epsilon_div_k must be a float!")    
        self.sigma = np.array([[sig,sig],[sig,sig]])
        self.epsilon = kB*np.array([[eps_div_k,eps_div_k],[eps_div_k,eps_div_k]])
        self.cpp_kingas = cpp_LJSpline(self.mole_weights, self.sigma, self.epsilon,self.is_idealgas,self._is_singlecomp)