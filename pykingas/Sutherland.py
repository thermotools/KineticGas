from pykingas.py_KineticGas import py_KineticGas
from pykingas import cpp_Sutherland
import numpy as np
from scipy.constants import Boltzmann
from warnings import warn
from thermopack.saftvrmie import saftvrmie

class Sutherland(py_KineticGas):

    def __init__(self, mole_weights, sigma, eps_div_k, C, lambdas, N=3, is_idealgas=False):
        """
        Note: When initializing for a specific set of components, the inheriting class should call the py_KineticGas
        constructor first, in order to retrieve the component parameters and initialize various attributes, then
        compute the coefficients and exponents required to call the Sutherland constructor.
        """

        # Assuming that the fluids attribute will be set only if the py_KineticGas constructor has already been called.
        if not hasattr(self, 'fluids'):
            ncomps = len(mole_weights)
            comps = ','.join(['PSEUDO' for _ in range(ncomps)])
            super().__init__(comps, mole_weights=mole_weights, N=N, is_idealgas=is_idealgas)

        assert np.shape(sigma) == (self.ncomps, self.ncomps)
        assert np.shape(eps_div_k) == (self.ncomps, self.ncomps)
        assert np.shape(C)[0] == np.shape(lambdas)[0]
        assert np.shape(C)[1:] == (self.ncomps, self.ncomps)
        assert np.shape(lambdas)[1:] == (self.ncomps, self.ncomps)

        self.C = C
        self.lambdas = lambdas
        self.sigma = sigma
        self.eps_div_k = eps_div_k

        self.cpp_kingas = cpp_Sutherland(mole_weights, sigma, eps_div_k * Boltzmann, C, lambdas, is_idealgas)
    
    def get_sigma_eff(self):
        return self.cpp_kingas.get_sigma_eff()
    
    def get_sigma_min(self):
        return self.cpp_kingas.get_sigma_min()

    def get_epsilon_eff(self):
        return self.cpp_kingas.get_epsilon_eff()

class S_MieKinGas(Sutherland):

    def __init__(self, comps,
                 mole_weights=None, sigma=None, eps_div_k=None,
                 la=None, lr=None, lij=0, kij=0,
                 N=4, is_idealgas=False, use_eos=None,
                 parameter_ref='default'):
        
        """Constructor
        If optional parameters are supplied, these are used instead of the parameters found in the database. To supply specific parameters for only some components, give `None` for the components that should use the database
        value
        &&
        Args:
            comps (str) : Comma-separated list of components
            mole_weights (optional, 1D array) : Molar masses [g/mol]
            sigma (optional, 1D array) : hard-sphere diameters [m]
            eps_div_k (optional, 1D array) : epsilon parameter / Boltzmann constant [-]
            la, lr (optional, 1D array) : attractive and repulsive exponent of the pure components [-]
            lij (optional, float) : Mixing parameter for sigma (lij > 0 => smaller sigma_12, lij < 0 => larger sigma_12)
            kij (optional, float) : Mixing parameter for epsilon (kij > 0 => favours mixing, kij < 0 => favours separation)
        """
        py_KineticGas.__init__(self, comps, mole_weights=mole_weights, N=N, is_idealgas=is_idealgas)

        if self.is_idealgas is False:
            if use_eos is None:
                self.eos = saftvrmie()
                if parameter_ref == 'default':
                    self.eos.init(comps)
                else:
                    self.eos.init(comps, parameter_reference=parameter_ref)
            else:
                self.eos = use_eos

        self.lij = lij
        self.kij = kij
        potential = 'Mie'

        try:
            self.fluids = [self.fluids[i][potential][parameter_ref] for i in range(self.ncomps)]
        except KeyError:
            for i in range(self.ncomps):
                if parameter_ref not in self.fluids[i][potential]:
                    warn('Missing parameter_ref ' + parameter_ref + ' for component ' + self.fluids[i]['ident'],
                         stacklevel=2)
            raise KeyError('Missing parameters ' + parameter_ref + ' for compontents ' + comps)

        if eps_div_k is None:
            eps_div_k = [self.fluids[i]['eps_div_k'] for i in range(self.ncomps)]
        elif None in eps_div_k:
            for i in range(self.ncomps):
                if eps_div_k[i] is None:
                    eps_div_k[i] = self.fluids[i]['eps_div_k']
        elif self._is_singlecomp is True:
            eps_div_k = [eps_div_k[0] for _ in range(2)]
        eps_div_k = np.array(eps_div_k)
        assert eps_div_k.shape == (self.ncomps,)
        self.epsilon_ij = self.get_epsilon_matrix(eps_div_k, kij)
        self.epsilon = np.diag(self.epsilon_ij)

        if la is None:
            la = [self.fluids[i]['lambda_a'] for i in range(self.ncomps)]
        elif None in la:
            for i in range(self.ncomps):
                if la[i] is None:
                    la[i] = self.fluids[i]['lambda_a']
        elif self._is_singlecomp is True:
            la = [la[0] for _ in range(2)]

        la = np.array(la)
        assert la.shape == (self.ncomps,)
        self.la = self.get_lambda_matrix(la, 0)

        if lr is None:
            lr = [self.fluids[i]['lambda_r'] for i in range(self.ncomps)]
        elif None in lr:
            for i in range(self.ncomps):
                if lr[i] is None:
                    lr[i] = self.fluids[i]['lambda_r']
        elif self._is_singlecomp is True:
            lr = [lr[0] for _ in range(2)]

        lr = np.array(lr)
        assert lr.shape == (self.ncomps,)
        self.lr = self.get_lambda_matrix(lr, self.lij)

        if sigma is None:
            sigma = [self.fluids[i]['sigma'] for i in range(self.ncomps)]
        elif None in sigma:
            for i in range(self.ncomps):
                if sigma[i] is None:
                    sigma[i] = self.fluids[i]['sigma']
        elif self._is_singlecomp is True:
            sigma = [sigma[0] for _ in range(2)]

        sigma = np.array(sigma)
        assert sigma.shape == (self.ncomps,)
        self.sigma_ij = self.get_sigma_matrix(sigma)
        self.sigma = self.sigma_ij

        self.C = (self.lr / (self.lr - self.la)) * (self.lr / self.la) ** (self.la / (self.lr - self.la))

        super().__init__(self.mole_weights, self.sigma, self.epsilon_ij / Boltzmann, [self.C, - self.C], [self.lr, self.la], is_idealgas=is_idealgas)
    
    def get_epsilon_matrix(self, eps_div_k, kij):
        """Utility
        Compute matrix of well-depths, given well depth of each component
        Warning: Use of mixing parameters is not thouroughly tested.
        &&
        Args:
            eps_div_k (1d array) : Well depth parameter of each component
            kij (2d array) : Not in use, internal parameter `self.kij` is used for mixing.

        Returns:
            2d array : Well depth for each interaction pair.
        """
        # Apply mixing rules
        epsilon = np.array(eps_div_k) * Boltzmann
        return (np.ones((self.ncomps, self.ncomps)) - self.kij * (np.ones((self.ncomps, self.ncomps)) - np.identity(self.ncomps))) * np.sqrt(
            epsilon * np.vstack(epsilon))  # Only apply mixing parameter kij to the off-diagonals

    def get_sigma_matrix(self, sigma):
        """Utility
        Compute interaction parameter $sigma$ for each particle pair, applying mixing parameters given by `self.lij`.
        Warning: Use of mixing parameters is not thouroughly tested.
        &&
        Args:
            sigma (1D array) : sigma-parameters [m]

        Retunrs:
            2d array : N x N matrix of sigma parameters, where sigma_ij = 0.5 * (sigma_i + sigma_j), if self.lij = 0.
                        Warning: Use of mixing parameters is not thouroughly tested.

        """

        sigma_ij = (np.ones((self.ncomps, self.ncomps)) - self.lij * (np.ones((self.ncomps, self.ncomps)) - np.identity(self.ncomps)))\
                * 0.5 * np.sum(np.meshgrid(sigma, np.vstack(sigma)), axis=0)  # Only apply mixing parameter lij to the off-diagonals

        return sigma_ij

    def get_lambda_matrix(self, lambdas, lij):
        r"""Utility
        Compute pair-interaction $\lambda_r$ parameters, apply mixing parameter.
        &&
        Args:
            lambdas (1d array) : Repulsive exponents for each pure-component interaction potential
            lij (1d array) : Mixing parameters

        Returns:
            2d array : Repulsive exponent for each pair-interaction.
        """
        l = np.array(lambdas)
        return 3 + (1 - lij) * np.sqrt((l - 3) * np.vstack(l - 3))
    
    def get_C_matrix(self):
        return (self.lr / (self.lr - self.la)) * (self.lr / self.la) ** (self.la / (self.lr - self.la))

class PureSutherland(Sutherland):

    def __init__(self, mole_weights, sigma, eps_div_k, C, lambdas, N=3, is_idealgas=False):
        """Constructor
        Construct a Sutherland model for a pure component
        &&
        Args:
            mole_weights (float) : Molar mass [g / mol]
            sigma (float) : Size scale [m]
            eps_div_k (float) : Energy scale divided by Boltzmanns constant [K]
            C (Iterable[float]) : Sutherland coefficients
            lambdas (Iterable[float]) : Sutherland exponents
            N (int, optional) : Default Enskog approximation order, defaults to 3.
            is_idealgas (bool, optional) : Whether model is at infinite dilution, defaults to false.
        """
        mole_weights = np.ones((2, 2)) * mole_weights
        sigma = np.ones((2, 2)) * sigma
        eps_div_k = np.ones((2, 2)) * eps_div_k
        C = np.ones((2, 2)) * C
        lambdas = np.ones((2, 2)) * lambdas
        super().__init__(mole_weights, sigma, eps_div_k, C, lambdas, N=N, is_idealgas=is_idealgas)