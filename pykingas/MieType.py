'''
Author : Vegard Gjeldvik Jervell
Contains : Parent class for all 'Mie-Type' potentials, including MieKinGas ...
'''
import numpy as np
from scipy.constants import Boltzmann, Avogadro
from pykingas.py_KineticGas import py_KineticGas
from warnings import warn
from pykingas.units import Units
from thermopack.saft import saft
import abc

class MieType(py_KineticGas):

    def __init__(self, comps, potential,
                 mole_weights=None, sigma=None, eps_div_k=None,
                 la=None, lr=None, lij=0, kij=0,
                 N=4, is_idealgas=False,
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
        super().__init__(comps, mole_weights=mole_weights, N=N, is_idealgas=is_idealgas)

        self.lij = lij
        self.kij = kij

        
        if isinstance(potential, str):
            potential = [potential for _ in range(self.ncomps)]
        elif self._is_singlecomp is True:
            potential = [potential[0], potential[0]]

        if 'q-Mie' in potential:
            for i, pot in enumerate(potential):
                if pot == 'q-Mie':
                    potential[i] = f'Mie-FH{self.fluids[i]["default_fh_order"]}'
        
        try:
            self.fluids = [self.fluids[i][potential[i]][parameter_ref] for i in range(self.ncomps)]
        except KeyError:
            for i in range(self.ncomps):
                if potential[i] not in self.fluids[i].keys():
                    raise KeyError(f'Component {self.comps[i]} does not have parameters for {potential[i]}.')
                if parameter_ref not in self.fluids[i][potential[i]]:
                    warn('Missing parameter_ref ' + parameter_ref + ' for component ' + self.fluids[i]['ident'],
                         stacklevel=2)
            raise KeyError('Missing parameters ' + parameter_ref + ' for compontents ' + comps)

        self.sigma, self.sigma_ij, self.epsilon, self.epsilon_ij, self.la, self.lr = [None] * 6
        MieType.set_eps_div_k(self, eps_div_k, update_eos=False)
        MieType.set_la(self, la, update_eos=False)
        MieType.set_lr(self, lr, update_eos=False)
        MieType.set_sigma(self, sigma, update_eos=False)
    
    def set_sigma(self, sigma, update_eos=True):
        """Utility
        Set the size parameter. Note: Running with `update_eos=False` will result in a model using
        different parameters for collision integrals and the RDF than for the equation of state.

        Args:
            sigma (list[float or None]) : Size parameter for each pure component. Use `None` for default values. Unit (m).
            update_eos (bool) : If True (default) also update the eos.
        """
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
        
        if update_eos is True:
            self.__update_eos_param__()
    
    def set_lr(self, lr, update_eos=True):
        """Utility
        Set the repulsive exponent. Note: Running with `update_eos=False` will result in a model using
        different parameters for collision integrals and the RDF than for the equation of state.

        Args:
            lr (list[float or None]) : Repulsive exponent for each pure component. Use `None` for default values.
            update_eos (bool) : If True (default) also update the eos.
        """
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
        
        if update_eos is True:
            self.__update_eos_param__()
    
    def set_la(self, la, update_eos=True):
        """Utility
        Set the attractive exponent. Note: Running with `update_eos=False` will result in a model using
        different parameters for collision integrals and the RDF than for the equation of state.

        Args:
            la (list[float or None]) : Attractive exponent for each pure component. Use `None` for default values.
            update_eos (bool) : If True (default) also update the eos.
        """
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
        
        if update_eos is True:
            self.__update_eos_param__()
    
    def set_eps_div_k(self, eps_div_k, update_eos=True):
        """Utility
        Set the well depth parameter. Note: Running with `update_eos=False` will result in a model using
        different parameters for collision integrals and the RDF than for the equation of state.

        Args:
            eps_div_k (list[float or None]) : Well depth parameter for each pure component. Use `None` for default values. Unit (K).
            update_eos (bool) : If True (default) also update the eos.
        """
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
        self.epsilon_ij = self.get_epsilon_matrix(eps_div_k, self.kij)
        self.epsilon = np.diag(self.epsilon_ij)
        
        if update_eos is True:
            self.__update_eos_param__()
    
    def __update_eos_param__(self):
        """Internal
        Update the EoS model to use the same parameters as this model is currently set with.

        Raises:
            AttributeError : If the EoS object does not support setting parameters
            Warning : If the EoS object is using segment numbers different from 1.
        """
        if (self.eos is None) or (self.is_idealgas is True):
            return

        if not isinstance(self.eos, saft):
            raise AttributeError(f"Model EoS object ({self.eos}) does not support setting potential parameters!")
        
        if self._is_singlecomp:
            set_ncomps = 1
        else:
            set_ncomps = self.ncomps
        
        for i in range(set_ncomps):
            segment_number = self.eos.get_pure_fluid_param(i + 1)[0]
            if segment_number != 1:
                warn(f"Current EoS uses segment number {segment_number} != 1 for component {self.comps[i]}", Warning, stacklevel=2)
            params = (segment_number, self.sigma[i][i], self.epsilon[i] / Boltzmann, self.la[i][i], self.lr[i][i])
            self.eos.set_pure_fluid_param(i + 1, *params)

    @abc.abstractmethod
    def __update_cpp_kingas_param__(self):
        """Internal
        Re-Initialize the C++ module with the current parameters in this model. Inheriting classes are responsible
        for initializing the correct model.
        """
        pass

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
        r"""Utility
        Compute interaction parameter $\sigma$ for each particle pair, applying mixing parameters given by `self.lij`.
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

    def get_BH_diameters(self, T):
        """Utility
        Compute Barker-Henderson diameters

        Args:
            T (float) : Temperature [K]

        Returns:
            2d array : Barker-Henderson Diameters, indexed by component pair [m]
        """
        return self.cpp_kingas.get_BH_diameters(T)

    def potential(self, i, j, r):
        """Utility
        Evaluate the interaction potential between types i and j at distance r

        Args:
            i, j (int) : Component indices
            r (float) : Distance [m]
        Returns:
            float : Interaction potential [J]
        """
        return self.cpp_kingas.potential(i, j, r)

    def potential_r(self, i, j, r):
        """Utility
        Evaluate the derivative of the interaction potential between types i and j at distance r

        Args:
            i, j (int) : Component indices
            r (float) : Distance [m]
        Returns:
            float : First derivative of interaction potential [N]
        """
        return self.cpp_kingas.potential_derivative_r(i, j, r)

    def potential_rr(self, i, j, r):
        """Utility
        Evaluate the second derivative of the interaction potential between types i and j at distance r

        Args:
            i, j (int) : Component indices
            r (float) : Distance [m]
        Returns:
            float : Second derivative of interaction potential [N / m]
        """
        return self.cpp_kingas.potential_dblderivative_rr(i, j, r)

    def saft_rdf(self, T, Vm, x, order=2, g2_correction=True):
        """cpp-interface
        Compute the radial distribution function at contact
        &&
        Args:
            T (float) : Temperature [K]
            Vm (float) : Molar volume [m3/mol]
            x (list[float]) : Molar composition [-]
            order (int) : Pertubation order
            g2_correction (bool) : Use correction factor for g2?

        Returns:
            2d array : RDF at contact, indexed by component pair.
        """
        particle_density = (1 / Vm) * Avogadro
        key = tuple((particle_density, T, tuple(x)))
        if key in self.computed_rdf.keys():
            return self.computed_rdf[key]

        rdf = self.cpp_kingas.saft_rdf(particle_density, T, x, order, g2_correction)
        self.computed_rdf[key] = rdf
        return rdf

    def get_reducing_units(self, comp_idx=0):
        """Utility
        Get reducing units for this model, as a `Units` struct. See `units.py`.
        &&
        Args:
            comp_idx (int, optional) : Which component to use for reducing units, defaults to first component

        Returns:
            Units : Struct holding the reducing units
        """
        return Units(self.mole_weights[comp_idx], self.sigma[comp_idx][comp_idx], self.epsilon_ij[comp_idx][comp_idx] / Boltzmann)