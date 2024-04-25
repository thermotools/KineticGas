'''
Author: Vegard Gjeldvik Jervell
Purpose: Wrapper for the MieKinGas class. Calls the MieType initializer with the appropriate parameter set identifier.
'''
from pykingas import cpp_MieKinGas, MieType
from thermopack.saftvrmie import saftvrmie

class MieKinGas(MieType.MieType):

    def __init__(self, comps,
                 mole_weights=None, sigma=None, eps_div_k=None,
                 la=None, lr=None, lij=0, kij=0,
                 N=4, is_idealgas=False, use_eos=None,
                 parameter_ref='default'):
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
            use_eos : (thermopack eos object, optional) EoS to use (initialized), defaults to `saftvrmie`
        """
        super().__init__(comps, 'Mie',
                    mole_weights=mole_weights, sigma=sigma,
                    eps_div_k=eps_div_k, la=la, lr=lr, lij=lij, kij=kij,
                    N=N, is_idealgas=is_idealgas,
                    parameter_ref=parameter_ref)

        self.cpp_kingas = cpp_MieKinGas(self.mole_weights, self.sigma_ij, self.epsilon_ij, self.la, self.lr, self.is_idealgas)
        if self.is_idealgas is False:
            if use_eos is None:
                self.eos = saftvrmie()
                if parameter_ref == 'default':
                    self.eos.init(comps)
                else:
                    self.eos.init(comps, parameter_reference=parameter_ref)
            else:
                self.eos = use_eos

    def get_vdw_alpha(self):
        lr, la = self.lr[0][0], self.la[0][0]
        C = (lr / (lr - la)) * (lr / la) ** (la / (lr - la))
        print('mie : ', lr, la, C)
        return C * ((1 / (la - 3)) - (1 / (lr - 3)))
