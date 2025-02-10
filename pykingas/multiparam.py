from pykingas.py_KineticGas import py_KineticGas
from pykingas.Quantum import Quantum, FH_Corrected
from .libpykingas import cpp_ModTangToennis, cpp_TangToennisParam, cpp_HFD_B2, cpp_FH_HFD_B2, cpp_Patowski, cpp_PatowskiFH
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
    
    def get_r_min(self, i, j):
        return self.cpp_kingas.get_r_min(i, j)

class ModTangToennies(MultiParam):
    
    def __init__(self, comps, parameter_ref='default'):
        """Constructor
        Initialize modified Tang-Toennies potential
        &&
        Args:
            comps (str) : Single component identifier
            parameter_ref (str, optional) : Identifier for parameter set to use
        """
        super().__init__(comps)

        potential = self.fluids[0]['ModTangToennis'][parameter_ref]
        self.eps_div_k = potential['eps_div_k']
        self.sigma = potential['sigma']
        self.Re = potential['Re'] * potential['L_unit']
        self.cpp_kingas = cpp_ModTangToennis(comps, True, 'default') # param, self.mole_weights, self.is_idealgas)

class HFD_B2(MultiParam, Quantum):

    def __init__(self, comps, quantum_active=True):
        super().__init__(comps)
        self.cpp_kingas = cpp_HFD_B2(comps)
        self.cpp_kingas.set_quantum_active(quantum_active)

class FH_HFD_B2(MultiParam, FH_Corrected):

    def __init__(self, comps, FH_order=1):
        super().__init__(comps)
        self.cpp_kingas = cpp_FH_HFD_B2(comps, FH_order)
        self.cpp_kingas.set_quantum_active(False)


class Patowski(MultiParam, Quantum):

    def __init__(self, comps, quantum_active=True):
        super().__init__(comps)
        self.cpp_kingas = cpp_Patowski(comps)
        self.param = self.cpp_kingas.get_param()
        self.cpp_kingas.set_quantum_active(quantum_active)

# class PatowskiFH1(MultiParam, Quantum):
# 
#     def __init__(self, comps):
#         super().__init__(comps)
#         self.cpp_kingas = cpp_PatowskiFH1(comps)
#         self.param = self.cpp_kingas.get_param()
#     
#     def potential(self, r, T):
#         return self.cpp_kingas.potential(0, 0, r, T)

class PatowskiFH(MultiParam, FH_Corrected):

    def __init__(self, comps, FH_order=1):
        super().__init__(comps)
        self.cpp_kingas = cpp_PatowskiFH(comps, FH_order)
        # if FH_order == 1:
        #     self.cpp_kingas = cpp_PatowskiFH1(comps)
        # elif FH_order == 2:
        #     self.cpp_kingas = cpp_PatowskiFH2(comps)
        # elif FH_order == 3:
        #     self.cpp_kingas = cpp_PatowskiFH3(comps)
        # else:
        #     raise IndexError(f'FH order {FH_order} not available in Python!')
        
        self.cpp_kingas.set_quantum_active(False)
    
    def potential(self, r, T):
        return self.cpp_kingas.potential(0, 0, r, T)
    
    def potential_r(self, r, T):
        return self.cpp_kingas.potential_r(0, 0, r, T)
    
    def potential_rr(self, r, T):
        return self.cpp_kingas.potential_rr(0, 0, r, T)
     
    def set_FH_order(self, FH_order):
        self.cpp_kingas.set_FH_order(FH_order)