'''
Author: Vegard Gjeldvik Jervell
Purpose: Wrapper for the MieKinGas class. Calls the MieType initializer with the appropriate parameter set identifier.
'''
from . import MieType
from .libpykingas import cpp_MieKinGas
from thermopack.saftvrmie import saftvrmie
from thermopack.saft import saft

class MieKinGas(MieType.MieType):

    def __init__(self, comps,
                 mole_weights=None, sigma=None, eps_div_k=None,
                 la=None, lr=None, lij=0, kij=0,
                 N=4, is_idealgas=False, use_eos=None,
                 parameter_ref='default', use_default_eos_param=False):
        """Constructor
        If parameters are explicitly supplied through optional arguments, these will be used instead of those in the database.
        To supply specific parameters for only some components, give `None` for the components that should use the database
        value
        &&
        Args:
            comps (str) : Comma-separated list of components
            mole_weights (1D array) : Molar weights [g/mol]
            sigma (1D array) : sigma-parameters [m]
            eps_div_k (1D array) : epsilon parameter / Boltzmann constant [-]
            la, lr (1D array) : attractive and repulsive exponent of the pure components [-]
            lij (float) : Mixing parameter for sigma (lij > 0 => smaller sigma_12, lij < 0 => larger sigma_12)
            kij (float) : Mixing parameter for epsilon (kij > 0 => favours mixing, kij < 0 => favours separation)
            use_eos (thermopack eos object, optional) : EoS to use (initialized), defaults to `saftvrmie`
            use_default_eos_param (bool) : If `False` (default), ensure that the EoS and RET-model use the same parameters
                                            (if applicable). If `False`, do not forward specified parameters to the EoS.
        """
        super().__init__(comps, 'Mie',
                    mole_weights=mole_weights, sigma=sigma,
                    eps_div_k=eps_div_k, la=la, lr=lr, lij=lij, kij=kij,
                    N=N, is_idealgas=is_idealgas,
                    parameter_ref=parameter_ref)

        self.__update_cpp_kingas_param__()
        if self.is_idealgas is False:
            if use_eos is None:
                self.eos = saftvrmie()
                if parameter_ref == 'default':
                    self.eos.init(comps)
                else:
                    self.eos.init(comps, parameter_reference=parameter_ref)
            else:
                self.eos = use_eos
            
            if isinstance(self.eos, saft) and (use_default_eos_param is False):
                self.__update_eos_param__()
        
    def __update_cpp_kingas_param__(self):
        """Internal
        See MieType
        """
        self.cpp_kingas = cpp_MieKinGas(self.mole_weights, self.sigma_ij, self.epsilon_ij, self.la, self.lr,
                                        self.is_idealgas)
    
    def set_sigma(self, sigma, update_eos=True):
        """Utility
        See MieType
        """
        super().set_sigma(sigma, update_eos=update_eos)
        self.__update_cpp_kingas_param__()
    
    def set_eps_div_k(self, eps_div_k, update_eos=True):
        """Utility
        See MieType
        """
        super().set_eps_div_k(eps_div_k, update_eos=update_eos)
        self.__update_cpp_kingas_param__()
    
    def set_la(self, la, update_eos=True):
        """Utility
        See MieType
        """
        super().set_la(la, update_eos=update_eos)
        self.__update_cpp_kingas_param__()
    
    def set_lr(self, lr, update_eos=True):
        """Utility
        See MieType
        """
        super().set_lr(lr, update_eos=update_eos)
        self.__update_cpp_kingas_param__()
            
    
