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

    def __init__(self, sig, eps_div_k, mole_weight, N = 3, is_ideal = False):
        """
        Constructor
        Initializes a LJSpline object. Can only be used for pure components. 
        &&
        Args:
            sig (float): sigma-parameter of LJ/spline potential [m]
            eps_div_k (float): ratio of epsilon-parameter of LJ/spline and Boltzmann's constant [-]
            mole_weight (float): molar mass [g/mol]
        """
        if not type(mole_weight) is float:
            raise TypeError("MW must be a float!")
        mole_weight = [mole_weight,mole_weight]
        super().__init__("AR", mole_weight=mole_weight, N = N, is_idealgas=is_ideal, is_single_component=True)
        if not type(sig) is float:
            raise TypeError("Sigma must be a float!")
        if not type(eps_div_k) is float:
            raise TypeError("Epsilon_div_k must be a float!")
        if not type(mole_weight) is float:
            raise TypeError("Mole_weight must be a float!")        
        self.sigma = np.array([[sig,sig],[sig,sig]])
        self.epsilon = kB*np.array([[eps_div_k,eps_div_k],[eps_div_k,eps_div_k]])
        self.cpp_kingas = cpp_LJSpline(self.mole_weight, self.sigma, self.epsilon,self.is_idealgas,self._is_singlecomp)
    def get_omega_star(self, T, r, l, i = 0, j = 0):
        return self.cpp_kingas.omega_star(i,j,T,r,l)
    def get_omega_star_approx(self, T, r, l, i = 0, j = 0):
        return self.cpp_kingas.omega_star_approx(i,j,T,r,l)
    def get_omega(self, T, r, l, i = 0, j = 0):
        return self.cpp_kingas.omega(i, j, l, r, T)