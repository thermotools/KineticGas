from pykingas import MieType, cpp_QuantumMie
from thermopack.saftvrqmie import saftvrqmie
import numpy as np
from warnings import warn

class QuantumMie(MieType.MieType):

    def __init__(self, comps,
                 mole_weights=None, sigma=None, eps_div_k=None,
                 la=None, lr=None, lij=0, kij=0,
                 N=4, FH_order=None, is_idealgas=False,
                 parameter_ref='default', use_eos=None):
        '''
        :param comps (str): Comma-separated list of components

        If parameters are explicitly supplied, these will be used instead of those in the database
        :param FH_order (int) : Feynman-Hibbs correction order (0 = Standard Mie potential, 
                                                                1 = 1st order correction, 
                                                                2 = 2nd order correction)
        '''

        super().__init__(comps, 'q-Mie',
                            mole_weights=mole_weights, sigma=sigma,
                            eps_div_k=eps_div_k, la=la, lr=lr, lij=lij, kij=kij,
                            N=N, parameter_ref=parameter_ref)

        fh_orders = np.array([self.fluids[i]['FH_order'] for i in range(self.ncomps)])
        if FH_order is None:
            FH_order = self.fluids[0]['FH_order']
        elif any(FH_order != fh_orders):
            warn('You have explicitly supplied a different FH_order that that found in the fluid parameter file!\n'
                 'Are you sure this is correct?', category=ResourceWarning, stacklevel=2)
        for i in range(self.ncomps):
            if any(fh_orders != fh_orders[i]):
                warn('There were components using different FH-orders! Are you sure this is correct?',
                     category=ResourceWarning, stacklevel=2)

        self.__FH_order = FH_order
        self.cpp_kingas = cpp_QuantumMie(self.mole_weights, self.sigma_ij, self.epsilon_ij, self.la, self.lr, 
                                            self.__FH_order, is_idealgas)
        if self.is_idealgas is False:
            if use_eos is None:
                self.eos = saftvrqmie()
                if parameter_ref == 'default':
                    self.eos.init(comps)
                else:
                    self.eos.init(comps, parameter_reference=parameter_ref)
            else:
                self.eos = use_eos

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
