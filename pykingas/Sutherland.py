from pykingas.py_KineticGas import py_KineticGas
from .libpykingas import cpp_Sutherland, cpp_ExtSutherland
import numpy as np
from scipy.constants import Boltzmann
from warnings import warn
from thermopack.saftvrmie import saftvrmie
from thermopack.saftvrss import saftvrss
from pykingas.units import Units

class Sutherland(py_KineticGas):

    def __init__(self, mole_weights, sigma, eps_div_k, C, lambdas, N=2, is_idealgas=False, use_eos=None, singlecomp=False):
        """
        Note: When initializing for a specific set of components, the inheriting class should call the py_KineticGas
        constructor first, in order to retrieve the component parameters and initialize various attributes, then
        compute the coefficients and exponents required to call the Sutherland constructor.
        """

        # Assuming that the fluids attribute will be set only if the py_KineticGas constructor has already been called.
        if not hasattr(self, 'fluids'):
            comps = ','.join(['LJF' for _ in mole_weights])
            super().__init__(comps, mole_weights=mole_weights, N=N, is_idealgas=is_idealgas)

        if np.shape(sigma) != (self.ncomps, self.ncomps): raise ValueError(f'Incorrect sigma array shape : {np.shape(sigma)}, expected {(self.ncomps, self.ncomps)}')
        if np.shape(eps_div_k) != (self.ncomps, self.ncomps): raise ValueError(f'Incorrect eps_div_k array shape : {np.shape(sigma)}, expected {(self.ncomps, self.ncomps)}')
        if np.shape(C) != np.shape(lambdas): raise ValueError(f'Shapes C and lambda ({np.shape(C)}, {np.shape(lambdas)}) do not match!')

        self.C = C
        self.lambdas = lambdas
        self.sigma = sigma
        self.eps_div_k = eps_div_k

        self.cpp_kingas = cpp_Sutherland(self.mole_weights, self.sigma, self.eps_div_k * Boltzmann, self.C, self.lambdas,
                                         is_idealgas, self._is_singlecomp)
        if (use_eos is None) and (self.eos is None):
            self.eos = saftvrss(','.join(self.comps), init_from_db='SAFT-VR-MIE')
            for i in range(self.ncomps):
                for j in range(self.ncomps):
                    self.eos.set_pair_potential_params(i + 1, j + 1, self.C[:, i, j], self.lambdas[:, i, j], self.sigma[i][j],
                                                       self.eps_div_k[i][j])
        self.cpp_kingas.set_eos(self.eos)

    @staticmethod
    def init_single(mole_weights, sigma, eps_div_k, C, lambdas, N=2, is_idealgas=False):
        """Constructor
        Initialize a pure-fluid Sutherland-Sum
        &&
        Args:
            mole_weights (float) : Mole weight of species (g / mol)
            sigma (float) : Size parameter (m)
            eps_div_k (float) : Energy parameter divided by Boltzmann constant (K)
            C (1d array) : Coefficient of each term
            lambdas (1d array) : Exponent of each term
            N (int, optional) : Default Enskog approximation order (default 2)
            is_idealgas (bool, optional) : Whether model is strictly at infinite dilution (default False)

        Returns:
             Sutherland : Initialized single-component model.
        """
        if len(lambdas) != len(C):
            raise ValueError(f"Number of coefficients ({len(C)}) does not match number of exponents ({len(lambdas)}!")

        eos = saftvrss('LJF', init_from_db='SAFT-VR-MIE')
        eos.set_pair_potential_params(1, 1, C, lambdas, sigma, eps_div_k)

        mole_weights = [mole_weights]
        sigma = Sutherland.mix_sigma([sigma, sigma])
        eps_div_k = Sutherland.mix_epsilon([eps_div_k, eps_div_k])
        C = Sutherland.mix_C_array([C, C])
        lambdas = Sutherland.mix_lambda_array([lambdas, lambdas])

        return Sutherland(mole_weights, sigma, eps_div_k, C, lambdas, N=N, is_idealgas=is_idealgas, use_eos=eos)

    @staticmethod
    def init_single_reduced(C, lambdas, N=2, is_idealgas=False):
        """Constructor
        Initialize a pure-fluid Sutherland-Sum, with default values for sigma, eps_div_k and mass
        &&
        Args:
            C (1d array) : Coefficient of each term
            lambdas (1d array) : Exponent of each term
            N (int, optional) : Default Enskog approximation order (default 2)
            is_idealgas (bool, optional) : Whether model is strictly at infinite dilution (default False)

        Returns:
             Sutherland : Initialized single-component model.
        """
        assert len(lambdas) == len(C)
        m = 10
        sigma = 3e-10
        eps_div_k = 100
        return Sutherland.init_single(m, sigma, eps_div_k, C, lambdas, N=N, is_idealgas=is_idealgas)

    @staticmethod
    def mix_sigma(sigma_vec, lij=None):
        """Utility
        Generate matrix of sigma parameters for each component pair from vector of pure component sigma
        &&
        Args:
            sigma_vec (1d array) : Pure component sigma (m)
            lij (2d array) : Not yet in use

        Returns:
            2d array : Sigma for all interactions
        """
        if lij is not None:
            raise NotImplementedError('Mixing parameters not implemented for Sutherland!')

        return 0.5 * (np.vstack([sigma_vec for _ in range(len(sigma_vec))]) + np.reshape(sigma_vec, (len(sigma_vec), 1)))

    @staticmethod
    def mix_epsilon(eps_vec, kij=None):
        """Utility
        Generate matrix of well-depth parameters for each component pair from vector of pure component well-depths
        &&
        Args:
            eps_vec (1d array) : Pure component well-depth (K)
            kij (2d array) : Not yet in use

        Returns:
            2d array : Well depth for all interactions
        """
        if kij is not None:
            raise NotImplementedError('Mixing parameters not implemented for Sutherland!')

        return np.sqrt(np.vstack([eps_vec for _ in range(len(eps_vec))]) * np.reshape(eps_vec, (len(eps_vec), 1)))

    @staticmethod
    def mix_C_array(C_arr, mix_rules=None, mix_params=None):
        """Utility
        Generate matrix of potential coefficients for each component pair from vector of pure component coefficients.
        Currently uses a geometric combining rule
        &&
        Args:
            C_arr (2d array) : Pure component coefficients, organised as C_arr[comp_idx][term_idx]
            mix_rules (None) : Not yet in use
            mix_params (None) : Not yet in use

        Returns:
            3d array : Coefficients for all interactions, of shape (nterms, ncomps, ncomos), where nterms is the largest
                        number of terms for any pure component.
        """
        if (mix_rules is not None) or (mix_params is not None):
            raise NotImplementedError('Mixing params not implemented for Sutherland')

        nterms = max([len(c_vec) for c_vec in C_arr])
        ncomps = len(C_arr)
        C_arr_in = np.zeros((ncomps, nterms))
        for i in range(ncomps):
            for j in range(len(C_arr[i])):
                C_arr_in[i][j] = C_arr[i][j]
        C_mixed = np.zeros((nterms, ncomps, ncomps))
        for term_i in range(nterms):
            for ci in range(ncomps):
                for cj in range(ncomps):
                    if np.sign(C_arr_in[ci][term_i]) != np.sign(C_arr_in[cj][term_i]):
                        warn('Mixing terms with opposite signs!', RuntimeWarning, stacklevel=2)
                    C_mixed[term_i][ci][cj] = np.sign(C_arr_in[ci][term_i]) * np.sqrt(C_arr_in[ci][term_i] * C_arr_in[cj][term_i])

        return C_mixed

    @staticmethod
    def mix_lambda_array(lambda_arr):
        """Utility
        Generate matrix of potential exponents for each component pair from vector of pure component exponents.
        Currently uses a geometric combining rule
        &&
        Args:
            lambda_arr (2d array) : Pure component coefficients, organised as C_arr[comp_idx][term_idx]

        Returns:
            3d array : Exponents for all interactions, of shape (nterms, ncomps, ncomos), where nterms is the largest
                        number of terms for any pure component.
        """
        return Sutherland.mix_C_array(lambda_arr)

    def get_reducing_units(self, comp_idx=0):
        """Utility
        Get reducing units for this model, as a `Units` struct. See `units.py`.
        &&
        Args:
            comp_idx (int, optional) : Which component to use for reducing units

        Returns:
            Units : Struct holding the reducing units
        """
        return Units(self.mole_weights[comp_idx], self.sigma[comp_idx][comp_idx], self.eps_div_k[comp_idx][comp_idx])

    @staticmethod
    def extract_components(array, i=None, j=None):
        """Utility
        """
        if (i is not None) and (j is not None):
            return array[i][j]
        elif i is not None:
            return array[i][i]
        return array

    def get_sigma_eff(self, i=None, j=None):
        r"""Potential
        Get effective size parameter (root of potential). If only one index is supplied, that pure component value
        is returned. If no indexes are supplied, a 2d array of all component pairs is returned.
        &&
        Args:
            i (int, optional) : First component
            j (int, optional) : Second component

        Returns:
            float or 2d array : Root of potential for the given interaction pair, or for all interaction pairs.
        """
        return self.extract_components(self.cpp_kingas.get_sigma_eff(), i, j)

    def get_sigma_min(self, i=None, j=None):
        r"""Potential
        Get location of the potential minimum. If only one index is supplied, that pure component value
        is returned. If no indexes are supplied, a 2d array of all component pairs is returned.
        &&
        Args:
            i (int, optional) : First component
            j (int, optional) : Second component

        Returns:
            float or 2d array : Location of the potential minimum for the given interaction pair, or for all interaction pairs.
        """
        return self.extract_components(self.cpp_kingas.get_sigma_min(), i, j)

    def get_epsilon_eff(self, i=None, j=None):
        r"""Potential
        Get potential well depth. If only one index is supplied, that pure component value
        is returned. If no indexes are supplied, a 2d array of all component pairs is returned.
        &&
        Args:
            i (int, optional) : First component
            j (int, optional) : Second component

        Returns:
            float or 2d array : Potential well depth for the given interaction pair, or for all interaction pairs.
        """
        return self.extract_components(self.cpp_kingas.get_epsilon_eff(), i, j)

    def get_vdw_alpha(self, i=None, j=None):
        r"""Potential
        Get the dimensionless van der Waals $\alpha$-parameter, defined as

        $$ \alpha = - \epsilon_{e}^{-1} \sigma_e^{-3} \int_{\sigma_e}^{\infty} \phi_{ij}(r) r^2 dr $$

        where $\sigma_e$ and $\epsilon_e$ are the potential root and well depth.
        &&
        Args:
            i (int) : Component 1
            j (int) : Component 2

        Returns:
            float : Dimensionless van der Waals $\alpha$-parameter
        """
        return self.extract_components(self.cpp_kingas.get_vdw_alpha(), i, j)

    def get_BH_diameters(self, T):
        return self.cpp_kingas.get_BH_diameters(T)
    
    def get_rdf_terms(self, rho, T, x):
        return self.cpp_kingas.get_rdf_terms(rho, T, x)

    def potential(self, i, j, r):
        r"""Potential
        Evaluate potential at the position r
        &&
        Args:
            i (int) : Component 1
            j (int) : Component 2
            r (float) : Position (m)

        Returns:
            float : Value of pair-potential (J)
        """
        return self.cpp_kingas.potential(i, j, r)

    def potential_r(self, i, j, r):
        r"""Potential
        Evaluate potential derivative at the position r
        &&
        Args:
            i (int) : Component 1
            j (int) : Component 2
            r (float) : Position (m)

        Returns:
            float : Value of pair-potential derivative (N)
        """
        return self.cpp_kingas.potential_derivative_r(i, j, r)

    def potential_rr(self, i, j, r):
        r"""Potential
        Evaluate potential second derivative at the position r
        &&
        Args:
            i (int) : Component 1
            j (int) : Component 2
            r (float) : Position (m)

        Returns:
            float : Value of pair-potential second derivative (N/m)
        """
        return self.cpp_kingas.potential_dblderivative_rr(i, j, r)


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
                eos = saftvrmie()
                if parameter_ref == 'default':
                    eos.init(comps)
                else:
                    eos.init(comps, parameter_reference=parameter_ref)
            else:
                eos = use_eos
        else:
            eos = None

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

        super().__init__(self.mole_weights, self.sigma, self.epsilon_ij / Boltzmann, [self.C, - self.C], [self.lr, self.la],
                         is_idealgas=is_idealgas, use_eos=eos)
    
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

class ExtSutherland(py_KineticGas):

    def __init__(self, mole_weights, sigma, eps_div_k, C, lambdas, beta_exp, rho_exp, N=2, is_idealgas=False, use_eos=None, singlecomp=False):
        if not hasattr(self, 'fluids'):
            comps = ','.join(['LJF' for _ in mole_weights])
            super().__init__(comps, mole_weights=mole_weights, N=N, is_idealgas=is_idealgas)

        if np.shape(sigma) != (self.ncomps, self.ncomps): raise ValueError(f'Incorrect sigma array shape : {np.shape(sigma)}, expected {(self.ncomps, self.ncomps)}')
        if np.shape(eps_div_k) != (self.ncomps, self.ncomps): raise ValueError(f'Incorrect eps_div_k array shape : {np.shape(sigma)}, expected {(self.ncomps, self.ncomps)}')
        if np.shape(C) != np.shape(lambdas): raise ValueError(f'Shapes C and lambda ({np.shape(C)}, {np.shape(lambdas)}) do not match!')

        self.C = C
        self.lambdas = lambdas
        self.sigma = sigma
        self.eps_div_k = eps_div_k
        self._is_singlecomp = singlecomp

        self.cpp_kingas = cpp_ExtSutherland(self.mole_weights, self.sigma, self.eps_div_k * Boltzmann, self.C, self.lambdas,
                                            beta_exp, rho_exp, is_idealgas, self._is_singlecomp)
        print(f"Use eos : {use_eos}, self.eos : {self.eos}")
        if (use_eos is None) and (self.eos is None):
            self.eos = saftvrss(','.join(self.comps), init_from_db='SAFT-VR-MIE')
            for i in range(self.ncomps):
                for j in range(self.ncomps):
                    self.eos.set_pair_potential_params(i + 1, j + 1, self.C[:, i, j], self.lambdas[:, i, j], self.sigma[i][j],
                                                       self.eps_div_k[i][j])
        else:
            self.eos = use_eos

        self.cpp_kingas.set_eos(self.eos)

    @staticmethod
    def init_single(mole_weights, sigma, eps_div_k, C, lambdas, beta_exp, rho_exp, N=2, is_idealgas=False):
        """Constructor
        Initialize a pure-fluid Sutherland-Sum
        &&
        Args:
            mole_weights (float) : Mole weight of species (g / mol)
            sigma (float) : Size parameter (m)
            eps_div_k (float) : Energy parameter divided by Boltzmann constant (K)
            C (1d array) : Coefficient of each term
            lambdas (1d array) : Exponent of each term
            beta_exp (1d array) : Inverse temperature exponent of each term
            rho_exp (1d array) : Density exponent of each term
            N (int, optional) : Default Enskog approximation order (default 2)
            is_idealgas (bool, optional) : Whether model is strictly at infinite dilution (default False)

        Returns:
             Sutherland : Initialized single-component model.
        """
        if len(lambdas) != len(C):
            raise ValueError(f"Number of coefficients ({len(C)}) does not match number of exponents ({len(lambdas)}!")

        eos = saftvrss('LJF', init_from_db='SAFT-VR-MIE')
        eos.set_pair_potential_params(1, 1, C, lambdas, sigma, eps_div_k)

        mole_weights = [mole_weights]
        sigma = ExtSutherland.mix_sigma([sigma, sigma])
        eps_div_k = ExtSutherland.mix_epsilon([eps_div_k, eps_div_k])
        C = ExtSutherland.mix_C_array([C, C])
        lambdas = ExtSutherland.mix_lambda_array([lambdas, lambdas])
        beta_exp = ExtSutherland.mix_lambda_array([beta_exp, beta_exp])
        rho_exp = ExtSutherland.mix_lambda_array([rho_exp, rho_exp])

        return ExtSutherland(mole_weights, sigma, eps_div_k, C, lambdas, beta_exp, rho_exp, N=N, is_idealgas=is_idealgas, use_eos=eos, singlecomp=True)

    @staticmethod
    def init_single_reduced(C, lambdas, beta_exp, rho_exp, N=2, is_idealgas=False):
        """Constructor
        Initialize a pure-fluid Sutherland-Sum, with default values for sigma, eps_div_k and mass
        &&
        Args:
            C (1d array) : Coefficient of each term
            lambdas (1d array) : Exponent of each term
            beta_exp (1d array) : Inverse temperature exponent of each term
            rho_exp (1d array) : Density exponent of each term
            N (int, optional) : Default Enskog approximation order (default 2)
            is_idealgas (bool, optional) : Whether model is strictly at infinite dilution (default False)

        Returns:
             Sutherland : Initialized single-component model.
        """
        assert len(lambdas) == len(C)
        m = 10
        sigma = 3e-10
        eps_div_k = 100
        return ExtSutherland.init_single(m, sigma, eps_div_k, C, lambdas, beta_exp, rho_exp, N=N, is_idealgas=is_idealgas)

    @staticmethod
    def mix_sigma(sigma_vec, lij=None):
        """Utility
        Generate matrix of sigma parameters for each component pair from vector of pure component sigma
        &&
        Args:
            sigma_vec (1d array) : Pure component sigma (m)
            lij (2d array) : Not yet in use

        Returns:
            2d array : Sigma for all interactions
        """
        if lij is not None:
            raise NotImplementedError('Mixing parameters not implemented for Sutherland!')

        return 0.5 * (np.vstack([sigma_vec for _ in range(len(sigma_vec))]) + np.reshape(sigma_vec, (len(sigma_vec), 1)))

    @staticmethod
    def mix_epsilon(eps_vec, kij=None):
        """Utility
        Generate matrix of well-depth parameters for each component pair from vector of pure component well-depths
        &&
        Args:
            eps_vec (1d array) : Pure component well-depth (K)
            kij (2d array) : Not yet in use

        Returns:
            2d array : Well depth for all interactions
        """
        if kij is not None:
            raise NotImplementedError('Mixing parameters not implemented for Sutherland!')

        return np.sqrt(np.vstack([eps_vec for _ in range(len(eps_vec))]) * np.reshape(eps_vec, (len(eps_vec), 1)))

    @staticmethod
    def mix_C_array(C_arr, mix_rules=None, mix_params=None):
        """Utility
        Generate matrix of potential coefficients for each component pair from vector of pure component coefficients.
        Currently uses a geometric combining rule
        &&
        Args:
            C_arr (2d array) : Pure component coefficients, organised as C_arr[comp_idx][term_idx]
            mix_rules (None) : Not yet in use
            mix_params (None) : Not yet in use

        Returns:
            3d array : Coefficients for all interactions, of shape (nterms, ncomps, ncomos), where nterms is the largest
                        number of terms for any pure component.
        """
        if (mix_rules is not None) or (mix_params is not None):
            raise NotImplementedError('Mixing params not implemented for Sutherland')

        nterms = max([len(c_vec) for c_vec in C_arr])
        ncomps = len(C_arr)
        C_arr_in = np.zeros((ncomps, nterms))
        for i in range(ncomps):
            for j in range(len(C_arr[i])):
                C_arr_in[i][j] = C_arr[i][j]
        C_mixed = np.zeros((nterms, ncomps, ncomps))
        for term_i in range(nterms):
            for ci in range(ncomps):
                for cj in range(ncomps):
                    if np.sign(C_arr_in[ci][term_i]) != np.sign(C_arr_in[cj][term_i]):
                        warn('Mixing terms with opposite signs!', RuntimeWarning, stacklevel=2)
                    C_mixed[term_i][ci][cj] = np.sign(C_arr_in[ci][term_i]) * np.sqrt(C_arr_in[ci][term_i] * C_arr_in[cj][term_i])

        return C_mixed

    @staticmethod
    def mix_lambda_array(lambda_arr):
        """Utility
        Generate matrix of potential exponents for each component pair from vector of pure component exponents.
        Currently uses a geometric combining rule
        &&
        Args:
            lambda_arr (2d array) : Pure component coefficients, organised as C_arr[comp_idx][term_idx]

        Returns:
            3d array : Exponents for all interactions, of shape (nterms, ncomps, ncomos), where nterms is the largest
                        number of terms for any pure component.
        """
        return Sutherland.mix_C_array(lambda_arr)

    def get_reducing_units(self, comp_idx=0):
        """Utility
        Get reducing units for this model, as a `Units` struct. See `units.py`.
        &&
        Args:
            comp_idx (int, optional) : Which component to use for reducing units

        Returns:
            Units : Struct holding the reducing units
        """
        return Units(self.mole_weights[comp_idx], self.sigma[comp_idx][comp_idx], self.eps_div_k[comp_idx][comp_idx])

    @staticmethod
    def extract_components(array, i=None, j=None):
        """Utility
        """
        if (i is not None) and (j is not None):
            return array[i][j]
        elif i is not None:
            return array[i][i]
        return array

    def get_sigma_eff(self, i=None, j=None):
        r"""Potential
        Get effective size parameter (root of potential). If only one index is supplied, that pure component value
        is returned. If no indexes are supplied, a 2d array of all component pairs is returned.
        &&
        Args:
            i (int, optional) : First component
            j (int, optional) : Second component

        Returns:
            float or 2d array : Root of potential for the given interaction pair, or for all interaction pairs.
        """
        return self.extract_components(self.cpp_kingas.get_sigma_eff(), i, j)

    def get_sigma_min(self, i=None, j=None):
        r"""Potential
        Get location of the potential minimum. If only one index is supplied, that pure component value
        is returned. If no indexes are supplied, a 2d array of all component pairs is returned.
        &&
        Args:
            i (int, optional) : First component
            j (int, optional) : Second component

        Returns:
            float or 2d array : Location of the potential minimum for the given interaction pair, or for all interaction pairs.
        """
        return self.extract_components(self.cpp_kingas.get_sigma_min(), i, j)

    def get_epsilon_eff(self, i=None, j=None):
        r"""Potential
        Get potential well depth. If only one index is supplied, that pure component value
        is returned. If no indexes are supplied, a 2d array of all component pairs is returned.
        &&
        Args:
            i (int, optional) : First component
            j (int, optional) : Second component

        Returns:
            float or 2d array : Potential well depth for the given interaction pair, or for all interaction pairs.
        """
        return self.extract_components(self.cpp_kingas.get_epsilon_eff(), i, j)

    def get_vdw_alpha(self, i=None, j=None):
        r"""Potential
        Get the dimensionless van der Waals $\alpha$-parameter, defined as

        $$ \alpha = - \epsilon_{e}^{-1} \sigma_e^{-3} \int_{\sigma_e}^{\infty} \phi_{ij}(r) r^2 dr $$

        where $\sigma_e$ and $\epsilon_e$ are the potential root and well depth.
        &&
        Args:
            i (int) : Component 1
            j (int) : Component 2

        Returns:
            float : Dimensionless van der Waals $\alpha$-parameter
        """
        return self.extract_components(self.cpp_kingas.get_vdw_alpha(), i, j)

    def get_BH_diameters(self, T):
        return self.cpp_kingas.get_BH_diameters(T)

    def get_rdf_terms(self, rho, T, x):
        return self.cpp_kingas.get_rdf_terms(rho, T, x)

    def potential(self, i, j, r):
        r"""Potential
        Evaluate potential at the position r
        &&
        Args:
            i (int) : Component 1
            j (int) : Component 2
            r (float) : Position (m)

        Returns:
            float : Value of pair-potential (J)
        """
        return self.cpp_kingas.potential(i, j, r)

    def potential_r(self, i, j, r):
        r"""Potential
        Evaluate potential derivative at the position r
        &&
        Args:
            i (int) : Component 1
            j (int) : Component 2
            r (float) : Position (m)

        Returns:
            float : Value of pair-potential derivative (N)
        """
        return self.cpp_kingas.potential_derivative_r(i, j, r)

    def potential_rr(self, i, j, r):
        r"""Potential
        Evaluate potential second derivative at the position r
        &&
        Args:
            i (int) : Component 1
            j (int) : Component 2
            r (float) : Position (m)

        Returns:
            float : Value of pair-potential second derivative (N/m)
        """
        return self.cpp_kingas.potential_dblderivative_rr(i, j, r)


class AT_Sutherland(ExtSutherland):

    def __init__(self, mole_weights, sigma, eps_div_k, C, lambdas, beta_exp, rho_exp, N=2, is_idealgas=False, use_eos=None, singlecomp=False):
        super().__init__(mole_weights, sigma, eps_div_k, C, lambdas, beta_exp, rho_exp, N=N, is_idealgas=is_idealgas, use_eos=use_eos, singlecomp=singlecomp)

    @staticmethod
    def init_single(mole_weights, sigma, eps_div_k, C, lambdas, beta_exp, N=2, is_idealgas=False, at_correction=1, at_alpha=0.1):
        """Constructor
        Initialize a pure-fluid Sutherland-Sum
        &&
        Args:
            mole_weights (float) : Mole weight of species (g / mol)
            sigma (float) : Size parameter (m)
            eps_div_k (float) : Energy parameter divided by Boltzmann constant (K)
            C (1d array) : Coefficient of each term
            lambdas (1d array) : Exponent of each term
            beta_exp (1d array) : Inverse temperature exponent of each term
            rho_exp (1d array) : Density exponent of each term
            N (int, optional) : Default Enskog approximation order (default 2)
            is_idealgas (bool, optional) : Whether model is strictly at infinite dilution (default False)

        Returns:
             Sutherland : Initialized single-component model.
        """
        rho_exp = [0 for _ in C]
        nterms = len(C)
        if at_correction == 1:
            for k in range(nterms):
                Ci, li, bi = C[k], lambdas[k], beta_exp[k]
                C.append(- Ci * at_alpha)
                lambdas.append(li)
                beta_exp.append(bi)
                rho_exp.append(1)
        elif at_correction == 2:
            for k in range(nterms):
                Ci, li, bi = C[k], lambdas[k], beta_exp[k]
                if Ci < 0:
                    C.append(- Ci * at_alpha)
                    lambdas.append(li)
                    beta_exp.append(bi)
                    rho_exp.append(1)

        print(f'Init with:\n\tC : {C}\n\tlamb : {lambdas}\n\tbeta_exp : {beta_exp}\n\trho_exp : {rho_exp}')
        return ExtSutherland.init_single(mole_weights, sigma, eps_div_k, C, lambdas, beta_exp, rho_exp, N=N, is_idealgas=is_idealgas)

    @staticmethod
    def init_single_reduced(C, lambdas, beta_exp, N=2, is_idealgas=False, at_correction=1, at_alpha=0.1):
        """Constructor
        Initialize a pure-fluid Sutherland-Sum, with default values for sigma, eps_div_k and mass
        &&
        Args:
            C (1d array) : Coefficient of each term
            lambdas (1d array) : Exponent of each term
            beta_exp (1d array) : Inverse temperature exponent of each term
            rho_exp (1d array) : Density exponent of each term
            N (int, optional) : Default Enskog approximation order (default 2)
            is_idealgas (bool, optional) : Whether model is strictly at infinite dilution (default False)

        Returns:
             Sutherland : Initialized single-component model.
        """
        assert len(lambdas) == len(C)
        m = 10
        sigma = 3e-10
        eps_div_k = 100
        return AT_Sutherland.init_single(m, sigma, eps_div_k, C, lambdas, beta_exp, N=N, is_idealgas=is_idealgas, at_correction=at_correction, at_alpha=at_alpha)
