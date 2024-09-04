from .py_KineticGas import py_KineticGas
from .libpykingas import cpp_ModTangToennis, cpp_TangToennisParam, cpp_AT_TangToennies
from .units import Units
from scipy.constants import Boltzmann, Avogadro
import numpy as np
from scipy.integrate import quad
from scipy.optimize import root

class ModTangToennis(py_KineticGas):

    def __init__(self, comps, parameter_ref='default'):
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

    def potential(self, r):
        return self.cpp_kingas.potential(0, 0, r)

    def potential_r(self, r):
        return self.cpp_kingas.potential_r(0, 0, r)

    def potential_rr(self, r):
        return self.cpp_kingas.potential_rr(0, 0, r)

    def second_virial(self, T):
        integrand = lambda r_aa: (1 - np.exp(- self.potential(r_aa * 1e-10) / (Boltzmann * T))) * 4 * np.pi * r_aa ** 2
        r1, r2 = self.Re * 1e10, 1.5 * self.Re * 1e10
        B0 = quad(integrand, 0, r1)[0]
        B1 = quad(integrand, r1, r2)[0]
        B2 = quad(integrand, r2, np.inf)[0]
        return 0.5 * (B0 + B1 + B2) * Avogadro * 1e-30

    def vdw_alpha(self):
        integrand = lambda r_aa: self.potential(r_aa * 1e-10) * r_aa ** 2
        E = quad(integrand, self.sigma * 1e10, 1.5 * self.sigma * 1e10)[0]
        E += quad(integrand, 1.5 * self.sigma * 1e10, np.inf)[0]
        E *= 1e-30
        return - E / ((self.eps_div_k * Boltzmann) * self.sigma ** 3)

    def get_eff_la(self, lr):
        alpha = self.vdw_alpha()
        mie_alpha = lambda la: mie_C(la, lr) * ((1 / (la - 3)) - (1 / (lr - 3)))
        sol = root(lambda la : mie_alpha(la) - alpha, x0=np.array([6.0]))
        return sol.x[0]

    def get_eff_lr(self, la):
        alpha = self.vdw_alpha()
        mie_alpha = lambda lr: mie_C(la, lr) * ((1 / (la - 3)) - (1 / (lr - 3)))
        sol = root(lambda lr : mie_alpha(lr) - alpha, x0=np.array([15.0]))
        return sol.x[0]

def mie_C(la, lr):
    return (lr / (lr - la)) * (lr / la)**(la / (lr - la))


class AT_TangToennies(py_KineticGas):

    def __init__(self, comps, parameter_ref='default'):
        super().__init__(comps, is_idealgas=True)

        self.cpp_kingas = cpp_AT_TangToennies(comps, self.is_idealgas, parameter_ref)

        param = self.cpp_kingas.get_param()
        self.eps_div_k = param.eps_div_k
        self.sigma = param.sigma
        self.Re = param.Re
    
    def potential(self, r, rho):
        return self.cpp_kingas.potential(0, 0, r, rho * Avogadro)

    def potential_r(self, r, rho):
        return self.cpp_kingas.potential_r(0, 0, r, rho * Avogadro)

    def potential_rr(self, r, rho):
        return self.cpp_kingas.potential_rr(0, 0, r, rho * Avogadro)

    def get_reducing_units(self):
        """Utility
        Get reducing units for this model, as a `Units` struct. See `units.py`.
        &&
        Args:
            comp_idx (int, optional) : Which component to use for reducing units, defaults to first component

        Returns:
            Units : Struct holding the reducing units
        """
        return Units(self.mole_weights[0], self.sigma, self.eps_div_k)