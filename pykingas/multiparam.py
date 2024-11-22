from pykingas.py_KineticGas import py_KineticGas
from pykingas.Quantum import Quantum
from .libpykingas import cpp_ModTangToennis, cpp_TangToennisParam, cpp_HFD_B2, cpp_Patowski
from scipy.constants import Boltzmann, Avogadro
import numpy as np
from scipy.integrate import quad
from scipy.optimize import root

class MultiParam(py_KineticGas):

    def __init__(self, comps):
        """Constructor
        Initialize modified Tang-Toennies potential
        &&
        Args:
            comps (str) : Single component identifier
            parameter_ref (str, optional) : Identifier for parameter set to use
        """
        super().__init__(comps, is_idealgas=True)

    def potential(self, r):
        """Utility
        Evaluate potential
        &&
        Args:
            r (float) : Distance (m)
        
        Returns:
            float : Potential (J)
        """
        return self.cpp_kingas.potential(0, 0, r)

    def potential_r(self, r):
        """Utility
        Evaluate potential derivative wrt. distance
        &&
        Args:
            r (float) : Distance (m)
        
        Returns:
            float : Potential derivative wrt. distance (N)
        """
        return self.cpp_kingas.potential_r(0, 0, r)

    def potential_rr(self, r):
        """Utility
        Evaluate potential second derivative wrt. distance
        &&
        Args:
            r (float) : Distance (m)
        
        Returns:
            float : Potential second derivative wrt. distance (N / m)
        """        
        return self.cpp_kingas.potential_rr(0, 0, r)

    def second_virial(self, T):
        """Utility
        Compute second virial coefficient
        &&
        Args:
            T (float) : Temperature (K)
        
        Returns: 
            float : Second virial coefficient
        """
        integrand = lambda r_aa: (1 - np.exp(- self.potential(r_aa * 1e-10) / (Boltzmann * T))) * 4 * np.pi * r_aa ** 2
        r1, r2 = self.Re * 1e10, 1.5 * self.Re * 1e10
        B0 = quad(integrand, 0, r1)[0]
        B1 = quad(integrand, r1, r2)[0]
        B2 = quad(integrand, r2, np.inf)[0]
        return 0.5 * (B0 + B1 + B2) * Avogadro * 1e-30

    def vdw_alpha(self):
        """Utility
        Get the dimensionless Van der Waals alpha-parameter
        &&
        Returns:
            float : alpha (-)
        """
        integrand = lambda r_aa: self.potential(r_aa * 1e-10) * r_aa ** 2
        E = quad(integrand, self.sigma * 1e10, 1.5 * self.sigma * 1e10)[0]
        E += quad(integrand, 1.5 * self.sigma * 1e10, np.inf)[0]
        E *= 1e-30
        return - E / ((self.eps_div_k * Boltzmann) * self.sigma ** 3)

class ModTangToennies(MultiParam):
    
    def __init__(self, comps, parameter_ref='default'):
        """Constructor
        Initialize modified Tang-Toennies potential
        &&
        Args:
            comps (str) : Single component identifier
            parameter_ref (str, optional) : Identifier for parameter set to use
        """
        super().__init__(comps, is_idealgas=True)

        potential = self.fluids[0]['ModTangToennis'][parameter_ref]
        param = cpp_TangToennisParam(potential['A_div_k'], potential['b'],
                                     potential['A_tilde_div_k'], potential['a'],
                                     potential['a_tilde'], potential['eps_div_k'], potential['Re'],
                                     potential['sigma'], potential['C']
                                     )
        self.eps_div_k = potential['eps_div_k']
        self.sigma = potential['sigma']
        self.Re = potential['Re'] * 1e-9
        self.cpp_kingas = cpp_ModTangToennis(param, self.mole_weights, np.ones((2, 2)) * self.sigma, self.is_idealgas)

class HFD_B2(Quantum):

    def __init__(self, comps):
        super().__init__(comps)
        self.cpp_kingas = cpp_HFD_B2(comps)

class Patowski(MultiParam, Quantum):

    def __init__(self, comps):
        super().__init__(comps)
        self.cpp_kingas = cpp_Patowski(comps)
        self.param = self.cpp_kingas.get_param()