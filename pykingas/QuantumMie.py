from pykingas import MieType, cpp_QuantumMie
from thermopack.saftvrqmie import saftvrqmie
import numpy as np
from warnings import warn

class QuantumMie(MieType.MieType):

    def __init__(self, comps,
                 mole_weights=None, sigma=None, eps_div_k=None,
                 la=None, lr=None, lij=0, kij=0,
                 N=4, FH_order=None, is_idealgas=False,
                 parameter_ref='default'):
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
        self.eos = saftvrqmie()
        self.eos.init(comps, parameter_reference=parameter_ref)