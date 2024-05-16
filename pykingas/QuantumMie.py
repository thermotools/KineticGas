from pykingas import MieType, cpp_QuantumMie
from pykingas.MieKinGas import MieKinGas
from thermopack.saftvrqmie import saftvrqmie
import numpy as np
from warnings import warn
import copy
from scipy.optimize import root
from scipy.constants import Boltzmann, Avogadro

class QuantumMie(MieType.MieType):

    def __init__(self, comps,
                 mole_weights=None, sigma=None, eps_div_k=None,
                 la=None, lr=None, lij=0, kij=0,
                 N=4, FH_orders=None, is_idealgas=False,
                 parameter_ref='default', use_eos=None):
        """
        If parameters are explicitly supplied, these will be used instead of those in the database

        Args:
            comps (str): Comma-separated list of components
            FH_orders (list[int] or int) : Feynman-Hibbs correction orders (0 = Standard Mie potential,
                                                                     1 = 1st order correction,
                                                                     2 = 2nd order correction)
        """
        if isinstance(FH_orders, int):
            FH_orders = [FH_orders for _ in range(len(comps.split(',')))]

        if FH_orders is None:
            potentials = ['q-Mie' for _ in range(len(comps.split(',')))]
        else:
            potentials = [f'Mie-FH{fh}' if (fh is not None) else 'q-Mie' for fh in FH_orders]

        super().__init__(comps, potentials,
                            mole_weights=mole_weights, sigma=sigma,
                            eps_div_k=eps_div_k, la=la, lr=lr, lij=lij, kij=kij,
                            N=N, parameter_ref=parameter_ref, is_idealgas=is_idealgas)

        fh_orders_db = np.array([self.fluids[i]['FH_order'] for i in range(self.ncomps)])
        self.__FH_orders = fh_orders_db
        self.cpp_kingas = cpp_QuantumMie(self.mole_weights, self.sigma_ij, self.epsilon_ij, self.la, self.lr, 
                                            self.__FH_orders, is_idealgas, self._is_singlecomp)

        if self.is_idealgas is False:
            if use_eos is None:
                self.eos = saftvrqmie()
                if parameter_ref == 'default':
                    self.eos.init(comps)
                else:
                    self.eos.init(comps, parameter_reference=parameter_ref)
            else:
                self.eos = use_eos

    def get_fh_order(self, ci=None):
        return self.__FH_orders[ci] if (ci is not None) else copy.deepcopy(self.__FH_orders)

    def get_sigma_eff(self, T):
        """Utility
        Compute effective sigma parameters

        Args:
            T (float) : Temperature [K]
        Returns:
            2d array : Effective sigma parameters [m]
        """
        return self.cpp_kingas.get_sigma_eff(T)

    def get_epsilon_eff(self, T):
        """Utility
        Compute effective epsilon parameter

        Args:
            T (float) : Temperature [K]
        Returns:
            2d array : Effective epsilon parameters [J]
        """
        return self.cpp_kingas.get_epsilon_eff(T)

    def get_sigma_min(self, T):
        """Utility
        Compute position of the potential minimum

        Args:
            T (float) : Temperature [K]
        Returns:
            2d array : Position of potential minimum [m]
        """
        return self.cpp_kingas.get_sigma_min(T)

    def get_C(self):
        """Utility
        Get the Mie-potential C prefactors

        Returns:
            2d array : prefactors [-]
        """
        return self.cpp_kingas.C

    def potential(self, i, j, r, T):
        """Utility
        Evaluate the interaction potential between types i and j at distance r

        Args:
            i, j (int) : Component indices
            r (float) : Distance [m]
            T (float) : Temperature [K]
        Returns:
            float : Interaction potential [J]
        """
        return self.cpp_kingas.potential(i, j, r, T)

    def potential_r(self, i, j, r, T):
        """Utility
        Evaluate the force between types i and j at distance r

        Args:
            i, j (int) : Component indices
            r (float) : Distance [m]
            T (float) : Temperature [K]
        Returns:
            float : Interaction force [N]
        """
        return self.cpp_kingas.potential_derivative_r(i, j, r, T)

    def potential_rr(self, i, j, r, T):
        """Utility
        Second derivative of potential between types i and j at distance r

        Args:
            i, j (int) : Component indices
            r (float) : Distance [m]
            T (float) : Temperature [K]
        Returns:
            float : Second derivative of potential [N / m]
        """
        return self.cpp_kingas.potential_dblderivative_rr(i, j, r, T)

    def get_effective_mie_model(self, T, fit_la=False):
        """Utility
        Fit parameters of an effective Mie potential that has the same root, well depth, position of minimum, and
        first and second derivatives at the root. If `fit_la=False` (default), only gives equal first derivative at root.

        Args:
             T (float) : Temperature [K]
             fit_la (bool) : If false, set la=6, and only solve for lr. If True, fit both la and lr.
        Returns:
            MieKinGas : An initialised model with the effective parameters.
        """

        sigma = np.diag(self.get_sigma_eff(T))
        eps = np.diag(self.get_epsilon_eff(T))
        r_min = np.diag(self.get_sigma_min(T))
        d = np.array([self.potential_r(i, i, sigma[i], T) for i in range(self.ncomps)])

        if fit_la is True:
            d2 = [self.potential_rr(i, i, sigma[i], T) for i in range(self.ncomps)]
            A = sigma / r_min
            B = - eps / (d * sigma)
            C = - (eps / (d * sigma)) - (d2 * eps / d ** 2)

            sol = root(lambda la: A ** la + B * la + C, x0=np.array([6.0 for _ in range(self.ncomps)]))
            lambda_a = sol.x
        else:
            lambda_a = [6 for _ in range(self.ncomps)]

        lambda_r = - (d * sigma / eps) * (sigma / r_min) ** lambda_a
        comps = ','.join(self.comps)
        mw = self.mole_weights * Avogadro * 1e3
        mie = MieKinGas(comps, mole_weights=mw, sigma=sigma, eps_div_k=eps / Boltzmann,
                        la=lambda_a, lr=lambda_r, use_eos=self.eos)
        return mie
