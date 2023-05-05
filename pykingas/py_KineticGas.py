'''
Author: Vegard Gjeldvik Jervell
Purpose: Parent class for python-wrappers to the KineticGas models. This class takes care of matrix inversion 
        and the final computations that yield transport coefficients, as well as reading fluid parameters from
        the 'fluids/XX.json' files.
'''

import numpy as np
import json
import scipy.linalg as lin
from scipy.constants import Boltzmann, Avogadro, pi
from scipy.linalg import block_diag
from scipy.integrate import quad
from pykingas import bcolors, suppress_stdout
import warnings, sys, time, atexit, os

FLT_EPS = 1e-12
__fluid_db_path__ = os.path.dirname(__file__) + '/fluids/'

kB = Boltzmann

def k_delta(i, j): # Kronecker delta
    if i == j:
        return 1
    return 0

class py_KineticGas:

    def __init__(self, comps, mole_weights=None, N=3, is_idealgas=False):
        '''
        :param comps (str): Comma-separated list of components, following Thermopack-convention
        :param mole_weights (1d array) : Mole weights [g/mol]. Will be used instead of database values if provided
        :param is_idealgas (bool) : If true, radial distribution function is unity, if false use radial distribution function of model
        :param parameter_ref (str) : Id for parameter set to use
        '''

        self._is_singlecomp = False
        self.default_N = N
        self.is_idealgas=is_idealgas
        self.computed_d_points = {} # dict of state points in which (d_1, d0, d1) have already been computed
        self.computed_a_points = {}  # dict of state points in which (a_1, a1) have already been computed
        self.computed_b_points = {}  # dict of state points in which (b_1, b1) have already been computed

        self.comps = comps.split(',')
        self.ncomps = len(self.comps)
        if self.ncomps == 1:
            self._is_singlecomp = True
            self.comps = [self.comps[0], self.comps[0]]
            self.ncomps = 2
        self.fluids = [json.load(open(__fluid_db_path__+c+'.json', 'r')) for c in self.comps]

        if mole_weights is None:
            mole_weights = np.array([self.fluids[i]['mol_weight'] for i in range(self.ncomps)])
        elif None in mole_weights:
            for i in range(self.ncomps):
                if mole_weights[i] is None:
                    mole_weights[i] = self.fluids[i]['mol_weight']
        elif len(mole_weights) > 1:
            mole_weights = np.array([mole_weights[i] for i in range(self.ncomps)])
        elif self._is_singlecomp is True:
            mole_weights = np.array([mole_weights[0] for i in range(self.ncomps)])

        self.mole_weights = np.array(mole_weights) * 1e-3 / Avogadro
        self.m = np.array([m for m in self.mole_weights])

        self.m0 = np.empty((self.ncomps, self.ncomps))
        self.M = np.empty_like(self.m0)
        for i in range(self.ncomps):
            for j in range(self.ncomps):
                self.m0[i][j] = self.m[i] + self.m[j]
                self.M[i][j] = self.m[i] / self.m0[i][j]

        self.cpp_kingas = None
        self.eos = None

    #####################################################
    #                   Utility                         #
    #####################################################
    def check_valid_composition(self, x):
        if abs(sum(x) - 1) > FLT_EPS:
            warnings.warn('Mole fractions do not sum to unity, sum(x) = '+str(sum(x)), stacklevel=2)
        elif len(x) != self.ncomps:
            raise IndexError(str(len(x)) +' mole fractions were supplied for a '+str(self.ncomps)+' component mixture!\n'
                            'Note: Single-component mixtures are treated as binary! (use x = [0.5, 0.5])')


    def get_Eij(self, Vm, T, x):
        """
        Compute the factors $ ( n_i / k_B T ) (d \mu_i / d n_j)_{T, n_{k \neq j}}$, where $n_i$ is the molar density
        of species $i$.

        :param Vm (float) : Molar volume [m3 / mol]
        :param T (float) : Temperature [K]
        :param x (list<float>) : Molar composition

        :return (2D array) : The factors E[i][j] = $ ( n_i / k_B T ) (d \mu_i / d n_j)_{T, n_{k \neq j}}$, where $n_i$
                                is the molar density of species $i$.
        """
        if self.is_idealgas is True:
            Vm = 1e6
        if self._is_singlecomp is True:
            x = [0.5, 0.5]
            _, dmudn_pure = self.eos.chemical_potential_tv(T, Vm, [0.5], dmudn=True)
            dmudrho_pure = Vm * dmudn_pure
            dmudrho = np.zeros((2, 2))
            rho = 1 / Vm
            dmudrho[0, 0] = dmudrho[1, 1] = dmudrho_pure + Avogadro * Boltzmann * T / rho
            dmudrho[0, 1] = dmudrho[1, 0] = dmudrho_pure - Avogadro * Boltzmann * T / rho
            dmudrho /= Avogadro
        else:
            _, dmudn = self.eos.chemical_potential_tv(T, Vm, x, dmudn=True)
            dmudrho = Vm * dmudn / Avogadro

        n = np.array(x) / Vm
        Eij = np.empty_like(dmudrho)
        for i in range(self.ncomps):
            for j in range(self.ncomps):
                Eij[i][j] = (n[i] / (kB * T)) * dmudrho[i][j]

        return Eij

    def get_P_factors(self, Vm, T, x):
        """
        Compute the factors $\Xi_i = \sum_j E_{ij}$, where $E_{ij}$ are the factors computed by `get_Eij`.

        :param Vm (float) : Molar volume [m3 / mol]
        :param T (float) : Temperature [K]
        :param x (list<float>) : Molar composition

        :return (1D array) : The factors $\Xi_i$
        """
        E = self.get_Eij(Vm, T, x)
        P = np.zeros(self.ncomps)
        for i in range(self.ncomps):
            P[i] = sum(E[:, i])
        return P

    #####################################################
    #          Transport property computations          #
    #####################################################

    def interdiffusion(self, T, Vm, x, N=None,
                       use_independent=True, dependent_idx=None,
                        frame_of_reference='CoN', use_binary=True,
                        solvent_idx=None):
        '''
        Compute the interdiffusion coefficients. Default definition is

        $J_i^{(n, n)} = - \sum_{j \neq l} D_{ij} \nabla n_j,$ $\nabla T = \nabla p = F_k = 0$ \forall $k$

        where the flux, $J_i^{(n, n)}$ is on a molar basis, in the molar frame of reference, and $j \neq l$ is an
        independent set of forces with l=dependent_idx.
        For fluxes in other frames of reference, use the `frame_of_reference` kwarg.
        For the diffusion coefficients describing the fluxes' response to all forces (not independent) the definition
        is:

        $J_i^{(n, n)} = - \sum_j D_{ij} \nabla n_j,$ $\nabla T = \nabla p = F_k = 0$ \forall $k$

        use the `use_independent` and `dependent_idx` kwargs to switch between these definitions.

        :param T: Temperature (K)
        :param Vm: Molar volume (m3 / mol)
        :param x: composition (mole fractions)
        :param N: (int) Enskog approximation order
        :param use_binary: (bool) If the mixture is binary, and an independent set of fluxes and forces is considered, i.e.
            `use_independent=True`, the diffusion coefficients will be exactly equal with opposite sign. Setting
            `use_binary=True` will in that case only return the coefficient describing the independent flux-force relation.
        :param use_independent: (bool) Return diffusion coefficients for independent set of forces.
        :param dependent_idx: (int) Index of the dependent molar density gradient (only if use_dependent=True).
            Defaults to last component, except when `frame_of_reference='solvent'`, in which case default is equal
            to kwargs['solvent_idx'].
        :param frame_of_reference (str) : Which frame of reference the diffusion coefficients apply to. Can be
            'CoN' (molar FoR), 'CoM' (barycentric FoR) or 'solvent' (solvent FoR). See the 'solvent_idx' kwarg for
            information on selecting the solvent index.
        :param solvent_idx (int) : Index of component identified as solvent (only when using `frame_of_reference='solvent'`)
        :return: array of diffusion coefficients
        '''
        D = self.interdiffusion_general(T, Vm, x, N=N)
        # psi = Transformation matrix from 'centre of mass' to 'frame_of_reference'
        # get_com_2_for_matr() dispatches the call to specific functions for different frames of reference.
        psi = self.get_com_2_for_matr(T, Vm, x, frame_of_reference, solvent_idx=solvent_idx)
        D = psi @ D
        if use_independent is True:
            if dependent_idx is None:
                if frame_of_reference == 'solvent':
                    dependent_idx = solvent_idx
                else:
                    dependent_idx = self.ncomps - 1
            P = self.get_P_factors(Vm, T, x)
            for i in range(self.ncomps):
                for j in range(self.ncomps):
                    D[i, j] -= (P[j] / P[dependent_idx]) * D[i, dependent_idx]
            if use_binary is True and self.ncomps == 2:
                independent_idx = 1 - dependent_idx
                return D[independent_idx][independent_idx] # Return the independent fluxes response to the independent force
        return D

    def interdiffusion_general(self, T, Vm, x, N=None):
        """
        Compute the 'Kinetic CoM diffusion coefficients', defined by

        $J_i^{(n, m)} = - \sum_j D_{ij} \nabla n_j,$ $\nabla T = \nabla p = F_k = 0$ \forall $k$

        For end-users, see: self.interdiffusion()
        :param T: Temperature (K)
        :param Vm: Molar volume (m3 / mol)
        :param x: composition (mole fractions)
        :param N: (int) Enskog approximation order

        :return (2D array): Array of the (not independent) $D^{(K, m)}$ diffusion coefficients
        """
        x = np.array(x)
        if N is None:
            N = self.default_N

        self.check_valid_composition(x)
        particle_density = Avogadro / Vm
        d = self.compute_diffusion_coeff_vector(particle_density, T, x, N=N)
        d = self.reshape_diffusion_coeff_vector(d)
        E = self.get_Eij(Vm, T, x)
        Dij = np.zeros((self.ncomps, self.ncomps))

        for i in range(self.ncomps):
            for j in range(self.ncomps):
                for k in range(self.ncomps):
                    Dij[i][j] += d[i][0][k] * E[k, j]
                Dij[i][j] *= x[i] / (2 * particle_density)
        return Dij

    def thermal_diffusion_coeff(self, T, Vm, x, N=None,
                                use_independent=False, dependent_idx=None,
                                frame_of_reference='CoN', solvent_idx=None):
        '''
        Compute thermal diffusion coefficients, D_{T,i} [mol / m^2 s]
        Default definition is

        $$J_i^{(n, n)} = D_{T,i} \nabla \ln T - \sum_j D_{ij} \nabla n_j, \nabla p = F_k = 0 \forall k$$

        where the flux, J_i^{(n, n)} is on a molar basis, in the molar frame of reference. For fluxes in other frames
        of reference, use the 'frame_of_reference' kwarg. For the diffusion coefficients corresponding to an
        independent set of forces, defined by

        $$J_i^{(n, n)} = D_{T,i} \nabla \ln T - \sum_{j \neq l} D_{ij} \nabla n_j, \nabla p = F_k = 0 \forall k$$

        where l is the index of the dependent molar density gradient, use the 'use_independent' and 'dependent_idx' kwargs.

        :param T: Temperature [K]
        :param Vm: Molar volume [m^3 / mol]
        :param x: Mole fractions [-]
        :param N: Enskog approximation order (>=2)
        :param use_independent: (bool) Return diffusion coefficients for independent set of forces.
        :param dependent_idx: (int) Index of the dependent molar density gradient (only if use_dependent=True).
                                Defaults to last component.
        :param frame_of_reference: What frame of reference the molar fluxes are measured in
        :param solvent_idx (int) : Index of component identified as solvent (only when using `frame_of_reference='solvent'`)
        :return: array of thermal diffusion coefficients
        '''
        if N is None:
            N = self.default_N
        self.check_valid_composition(x)

        if N < 2:
            warnings.warn('Thermal diffusion is a 2nd order phenomena, cannot be computed for N < 2 (got N = '
                          + str(N) + ')', RuntimeWarning, stacklevel=2)
            return np.full(self.ncomps, np.nan)

        particle_density = Avogadro / Vm
        d = self.compute_diffusion_coeff_vector(particle_density, T, x, N=N)
        d = self.reshape_diffusion_coeff_vector(d)
        a = self.compute_cond_vector(particle_density, T, x, N=N)
        E = self.get_Eij(Vm, T, x)
        P = self.get_P_factors(Vm, T, x)
        rdf = self.cpp_kingas.get_rdf(particle_density, T, x)
        cd = self.cpp_kingas.get_contact_diameters(particle_density, T, x)

        b = np.empty((self.ncomps, self.ncomps)) # Precomputing some factors that are used many places later
        for j in range(self.ncomps):
            for k in range(self.ncomps):
                b[j, k] = k_delta(j, k) + (4 * np.pi / 3) * particle_density * x[k] * cd[j][k] ** 3 * self.M[j, k] \
                          * rdf[j][k]

        DT = np.zeros(self.ncomps)
        for i in range(self.ncomps):
            outer_sum = 0
            for j in range(self.ncomps):
                inner_sum = 0
                for k in range(self.ncomps):
                    inner_sum += b[j, k]
                outer_sum += d[i][0][j] * x[j] * inner_sum
            DT[i] = (x[i] / 2) * (a[i] - outer_sum)
        psi = self.get_com_2_for_matr(T, Vm, x, frame_of_reference, solvent_idx=solvent_idx)
        DT = psi @ DT

        if use_independent is True:
            if dependent_idx is None:
                dependent_idx = self.ncomps - 1

            Dij = self.interdiffusion_general(T, Vm, x, N=N)
            Dij = psi @ Dij
            for i in range(self.ncomps):
                tmp = 0
                for k in range(self.ncomps):
                    for m in range(self.ncomps):
                        tmp -= particle_density * x[m] * b[m, k]

                DT[i] -= (Dij[i, dependent_idx] / P[dependent_idx]) * tmp

        return DT / Avogadro # Dividing by Avogadros number to convert unit from particle number to mole number

    def thermal_diffusion_ratio(self, T, Vm, x, N=None):
        """
        Calculate the "independent" thermal diffusion ratios, $k_{T, i}$ defined by

        $$J_i^{(n, n)} = - \sum_{j \neq l} D_{ij}^{(I, n)} ( \nabla n_j + n_j k_{T, j} D_{T, j}^{(I, n)} \nabla \ln T )$$

        and

        $$\sum_i x_i \sum_j [\delta_{ij} + (4 \pi / 3) \sigma_{ij}^3 n_j M_{ij} \chi_{ij} - (n_j k_{T,j} / k_B T) ( d \mu_i / n_j )_{T,n_{l \neq j}}] = 0$$

        This definition implies that

        $$\nabla n_j = -n_j k_{T,j} \nabla \ln T\ \forall j$$ when all mass fluxes vanish.

        Note that the thermal diffusion ratios are independent of the frame of reference.

        :param T (float) : Temperature [K]
        :param Vm (float) : Molar volume [m3 / mol]
        :param x (list<float>) : Molar composition [-]
        :param N (Optional, int) : Enskog approximation order (>= 2)
        :return (1D array) : The thermal diffusion ratio of each component.
        """
        if N is None:
            N = self.default_N

        if N < 2:
            warnings.warn('Thermal diffusion is a 2nd order phenomena, cannot be computed for N < 2 (got N = '
                          + str(N) + ')', RuntimeWarning, stacklevel=2)
            return np.full(self.ncomps, np.nan)

        particle_density = Avogadro / Vm
        DT = self.thermal_diffusion_coeff(T, Vm, x, N=N, frame_of_reference='CoM',
                                          use_independent=True, dependent_idx=self.ncomps - 1)
        Dij = self.interdiffusion(T, Vm, x, N=N, frame_of_reference='CoM',
                                  use_binary=False, use_independent=True, dependent_idx=self.ncomps - 1)
        rdf = self.cpp_kingas.get_rdf(particle_density, T, x)
        cd = self.cpp_kingas.get_contact_diameters(particle_density, T, x)
        P = self.get_P_factors(Vm, T, x)
        E = self.get_Eij(Vm, T, x)
        A = np.zeros((self.ncomps, self.ncomps))

        for i in range(self.ncomps - 1):
            for j in range(self.ncomps):
                A[i, j] = - Dij[i, j] * x[j] * (1 / Vm) # density in moles, because DT has unit moles (not particles)

        for i in range(self.ncomps):
            A[-1, i] = x[i] * P[i]

        DT[-1] = 0 # Overwriting the DT vector for the final element of the b-vector of A @ kT = b
        for i in range(self.ncomps):
            for j in range(self.ncomps):
                DT[-1] += x[i] * (k_delta(i, j) + (4 * np.pi / 3) * particle_density * x[j] * cd[i][j]**3 * self.M[i, j] * rdf[i][j])

        kT = np.linalg.solve(A, DT)

        return kT

    def thermal_diffusion_factor(self, T, Vm, x, N=None):
        """
        Compute the thermal diffusion factors $\alpha_{ij}$, defined by

        $$\alpha_{ij} = k_{T, i} - k_{T, j}$$

        where $k_{T,i}$ are the thermal diffusion ratios. This definition implies that

        $$\nabla \ln (n_i / n_j) = - \alpha_{ij} \nabla \ln T$$ when the mass fluxes vanish.

        The thermal diffusion factors are independent of the frame
        of reference.

        :param T (float) : Temperature [K]
        :param Vm (float) : Molar volume [m3 / mol]
        :param x (list<float>) : Molar composition [-]
        :param N (Optional, int) : Enskog approximation order (>= 2)

        :return (2D array) : The thermal diffusion factors of the mixture.
        """
        if N is None:
            N = self.default_N

        if N < 2:
            warnings.warn('Thermal diffusion is a 2nd order phenomena, cannot be computed for N < 2 (got N = '
                          + str(N) + ')', RuntimeWarning, stacklevel=2)
            return np.full((self.ncomps, self.ncomps), np.nan)

        kT = self.thermal_diffusion_ratio(T, Vm, x, N=N)
        alpha = np.empty((self.ncomps, self.ncomps))
        for i in range(self.ncomps):
            for j in range(self.ncomps):
                alpha[i, j] = kT[i] - kT[j]

        return alpha

    def soret_coefficient(self, T, Vm, x, N=None):
        """
        Compute the Soret coefficients, $S_{T, ij}$
        """
        if N < 2:
            warnings.warn('Thermal diffusion is a 2nd order phenomena, cannot be computed for N < 2 (got N = '
                          + str(N) + ')', RuntimeWarning, stacklevel=2)
            return np.full((self.ncomps, self.ncomps), np.nan)
        alpha = self.thermal_diffusion_factor(T, Vm, x, N=N)
        return alpha / T

    def thermal_conductivity(self, T, Vm, x, N=None):
        """
        Compute the thermal conductivity, $\lambda$. For models initialized with `is_idealgas=True`, the thermal
        conductivity is not a function of density (i.e. $d \lambda / d V_m = 0).

        :param T (float) : Temperature [K]
        :param Vm (float) : Molar volume [m3 / mol]
        :param x (list<float>) : Molar composition [-]
        :param N (Optional, int) : Enskog approximation order

        :return (float) : The thermal conductivity of the mixture.
        """
        if N is None:
            N = self.default_N
        self.check_valid_composition(x)
        particle_density = Avogadro / Vm
        a = self.compute_cond_vector(particle_density, T, x, N=N)
        rdf = self.cpp_kingas.get_rdf(particle_density, T, x)
        K = self.cpp_kingas.get_K_factors(particle_density, T, x)
        cd = self.cpp_kingas.get_contact_diameters(particle_density, T, x)
        d = self.compute_diffusion_coeff_vector(particle_density, T, x, N=N)
        d = self.reshape_diffusion_coeff_vector(d)

        lambda_dblprime = 0
        if self.is_idealgas is False: # lambda_dblprime is only nonzero when density corrections are present, and vanishes at infinite dilution
            for i in range(self.ncomps):
                for j in range(self.ncomps):
                    lambda_dblprime += particle_density ** 2 * np.sqrt(2 * pi * self.m[i] * self.m[j] * Boltzmann * T / (self.m[i] + self.m[j])) \
                                        * (x[i] * x[j]) / (self.m[i] + self.m[j]) * (cd[i][j] ** 4) * rdf[i][j]
            lambda_dblprime *= (4 * Boltzmann / 3)

        lambda_prime = 0
        dth = self.compute_dth_vector(particle_density, T, x, N=N)
        for i in range(self.ncomps):
            tmp = 0
            for k in range(self.ncomps):
                tmp += d[i, 1, k] * dth[k]
            lambda_prime += x[i] * K[i] * (a[self.ncomps + i] - tmp)
        lambda_prime *= (5 * Boltzmann / 4)

        cond = lambda_prime + lambda_dblprime
        return cond

    def viscosity(self, T, Vm, x, N=None):
        """
        Compute the shear viscosity, $\eta$. For models initialized with `is_idealgas=True`, the shear viscosity
        is not a function of density (i.e. $d \eta / d V_m = 0).

        :param T (float) : Temperature [K]
        :param Vm (float) : Molar volume [m3 / mol]
        :param x (list<float>) : Molar composition [-]
        :param N (Optional, int) : Enskog approximation order

        :return (float) : The shear viscosity of the mixture.
        """
        if N is None:
            N = self.default_N
        self.check_valid_composition(x)

        particle_density = Avogadro / Vm
        b = self.compute_visc_vector(T, particle_density, x, N=N)
        cd = self.cpp_kingas.get_contact_diameters(particle_density, T, x)
        rdf = self.cpp_kingas.get_rdf(particle_density, T, x)
        K_prime = self.cpp_kingas.get_K_prime_factors(particle_density, T, x)

        eta_prime = 0
        for i in range(self.ncomps):
            eta_prime += K_prime[i] * x[i] * b[i]
        eta_prime *= Boltzmann * T / 2

        eta_dblprime = 0
        if self.is_idealgas is False: # eta_dblprime is only nonzero when density corrections are present, and vanish at infinite dilution
            for i in range(self.ncomps):
                for j in range(self.ncomps):
                    eta_dblprime += np.sqrt(self.m[i] * self.m[j] / (self.m[i] + self.m[j])) * x[i] * x[j] * cd[i][j]**4 * rdf[i][j]
            eta_dblprime *= 4 * particle_density**2 * np.sqrt(2 * np.pi * Boltzmann * T) / 15

        return eta_prime + eta_dblprime

    def bulk_viscosity(self, T, Vm, x, N=None):
        """
        Not implemented
        """
        if N is None:
            N = self.default_N
        self.check_valid_composition(x)
        raise NotImplementedError("Bulk viscosity is not implemented yet. See 'Multicomp docs' for more info.")

    #####################################################
    #                    Tp-interface                   #
    # Compute properties as a function of (T, p, x) by  #
    # Using self.eos to compute molar volume            #
    #####################################################

    def interdiffusion_tp(self, T, p, x, N=None,
                            use_independent=True, dependent_idx=None,
                            frame_of_reference='CoN', use_binary=True,
                            solvent_idx=None
                          ):
        """
        Compute molar volume using the internal equation of state (`self.eos`), assuming vapour, and pass the call to
        `self.interdiffusion`. See `self.interdiffusion` for documentation.
        """
        Vm, = self.eos.specific_volume(T, p, x, self.eos.VAPPH) # Assuming vapour phase
        return self.interdiffusion(T, Vm, x, N=N, use_independent=use_independent,
                                   dependent_idx=dependent_idx, frame_of_reference=frame_of_reference,
                                   use_binary=use_binary, solvent_idx=solvent_idx)

    def thermal_diffusion_coeff_tp(self, T, p, x, N=None,
                                    use_independent=False, dependent_idx=None,
                                    frame_of_reference='CoN', solvent_idx=None
                                    ):
        """
        Compute molar volume using the internal equation of state (`self.eos`), assuming vapour, and pass the call to
        `self.thermal_diffusion_coeff`. See `self.thermal_diffusion_coeff` for documentation.
        """
        Vm, = self.eos.specific_volume(T, p, x, self.eos.VAPPH)  # Assuming vapour phase
        return self.thermal_diffusion_coeff(T, Vm, x, N=N, use_independent=use_independent,
                                   dependent_idx=dependent_idx, frame_of_reference=frame_of_reference,
                                   solvent_idx=solvent_idx)

    def thermal_diffusion_factor_tp(self, T, p, x, N=None):
        """
        Compute molar volume using the internal equation of state (`self.eos`), assuming vapour, and pass the call to
        `self.thermal_diffusion_factor`. See `self.thermal_diffusion_factor` for documentation.
        """
        Vm, = self.eos.specific_volume(T, p, x, self.eos.VAPPH)  # Assuming vapour phase
        return self.thermal_diffusion_factor(T, Vm, x, N=N)

    def thermal_coductivity_tp(self, T, p, x, N=None):
        """
        Compute molar volume using the internal equation of state (`self.eos`), assuming vapour, and pass the call to
        `self.thermal_conductivity`. See `self.thermal_conductivity` for documentation.
        """
        Vm, = self.eos.specific_volume(T, p, x, self.eos.VAPPH)  # Assuming vapour phase
        return self.thermal_conductivity(T, Vm, x, N=N)

    def viscosity_tp(self, T, p, x, N=None):
        """
        Compute molar volume using the internal equation of state (`self.eos`), assuming vapour, and pass the call to
        `self.viscosity`. See `self.viscosity` for documentation.
        """
        Vm, = self.eos.specific_volume(T, p, x, self.eos.VAPPH)  # Assuming vapour phase
        return self.viscosity(T, Vm, x, N=N)

    #####################################################
    #        Frame of Reference Transformations         #
    #####################################################

    def get_com_2_for_matr(self, T, Vm, x, FoR, **kwargs):
        '''
        Dispatcher to get a specific 'change of frame of reference' matrix for transformation from
        centre of mass to 'FoR'

        Returns the appropriate matrix for the transformations derived in Appendix A of ... to transform a flux
        from the centre of mass frame of reference to the 'FoR' frame of reference by the transformation

        $$J^{(n, FoR)} = \psi @ J^{(n, m)}$$

        where $\psi$ is the matrix returned by this method, and $J$ is the vector of (all) molar fluxes, with the
        subscript indicating the frame of reference.

        :param T (float) : Temperature [K]
        :param Vm (float) : Molar volume [m3 / mol]
        :param x (list<float>) : Molar composition [-]

        :return (2Darray) : The $N$ x $N$ transformation matrix to transform the fluxes, with $N$ being the number of
                            components.
        '''
        if FoR == 'CoN':
            return self.get_com_2_con_matr(x)
        elif FoR == 'solvent':
            return self.get_com_2_solv_matr(x, kwargs['solvent_idx'])
        elif FoR == 'CoM':
            return np.identity(self.ncomps)
        elif FoR == 'CoV':
            return self.get_com_2_cov_matr(T, Vm, x)
        else:
            print('Invalid frame of reference key :', FoR)
            print("Valid keys are 'CoN' (centre of moles), 'CoM' (centre of mass), 'solvent' (solvent),")
            print("and 'CoV' (centre of volume).")
            print("When using 'solvent', the solvent compontent index must be supplied via the 'solvent_idx' kwarg")
            raise KeyError('Invalid frame of reference key : ' + FoR)

    def get_com_2_con_matr(self, x):
        psi = np.identity(self.ncomps)
        w = x * self.m / sum(x * self.m)
        k = 0 # Arbitrary component (tested)
        for i in range(self.ncomps):
            for j in range(self.ncomps):
                psi[i, j] -= x[i] * (1 - w[j] * x[k] / (x[j] * w[k]))
        return psi

    def get_com_2_solv_matr(self, x, solvent_idx):
        if solvent_idx < 0:
            solvent_idx += self.ncomps
        psi = np.identity(self.ncomps)
        m = self.m * Avogadro * 1e3
        wt_frac = m * x / sum(m * x)
        for i in range(self.ncomps):
            for j in range(self.ncomps):
                psi[i, j] += x[i] * ((wt_frac[j] / x[j]) - (k_delta(j, solvent_idx) / x[solvent_idx]))
        return psi

    def get_com_2_cov_matr(self, T, Vm, x):
        p, = self.eos.pressure_tv(T, Vm, x)
        _, dvdn = self.eos.specific_volume(T, p, x, self.eos.VAPPH, dvdn=True)
        psi = np.identity(self.ncomps)
        m = self.m * Avogadro * 1e3
        wt_frac = m * x / sum(m * x)
        for i in range(self.ncomps):
            for j in range(self.ncomps):
                psi[i, j] += x[i] * ((wt_frac[j] / x[j]) - (dvdn[j] / Vm))
        return psi

    def get_solv_2_solv_matr(self, x, prev_solv_idx, new_solv_idx):
        psi = np.identity(self.ncomps)
        k = prev_solv_idx # Arbitrary component (tested)
        for i in range(self.ncomps):
            for j in range(self.ncomps):
                psi[i, j] -= x[i] * (k_delta(j, new_solv_idx) - k_delta(j, prev_solv_idx) * k_delta(k, new_solv_idx)) / \
                             x[new_solv_idx]
        return psi

    ######################################################
    #         Interface to top-level C++ methods         #
    ######################################################
    
    def get_conductivity_matrix(self, particle_density, T, mole_fracs, N=None):
        """
        Compute the elements of the matrix corresponding to the set of equations that must be solved for the
        thermal response function sonine polynomial expansion coefficients

        :param particle_density (float) : Particle density (not molar!) [1 / m3]
        :param T (float) : Temperature [K]
        :param mole_fracs (list<float>) : Molar composition [-]
        :param N (Optional, int) : Enskog approximation order.

        :return (2Darray) : ($N N_c$ x $N N_c$) matrix, where $N$ is the Enskog approximation order and $N_c$ is
                            the number of components.
        """
        if N is None:
            N = self.default_N
        self.check_valid_composition(mole_fracs)
        return self.cpp_kingas.get_conductivity_matrix(particle_density, T, mole_fracs, N)

    def get_conductivity_vector(self, particle_density, T, mole_fracs, N):
        """
        Compute the right-hand side vector to the set of equations that must be solved for the
        thermal response function Sonine polynomial expansion coefficients.

        :param particle_density (float) : Particle density (not molar!) [1 / m3]
        :param T (float) : Temperature [K]
        :param mole_fracs (list<float>) : Molar composition [-]
        :param N (Optional, int) : Enskog approximation order.

        :return (1Darray) : ($N N_c$,) vector, where $N$ is the Enskog approximation order and $N_c$ is
                            the number of components.
        """
        if N is None:
            N = self.default_N
        self.check_valid_composition(mole_fracs)
        return self.cpp_kingas.get_conductivity_vector(particle_density, T, mole_fracs, N)

    def get_diffusion_vector(self, particle_density, T, mole_fracs, N=None):
        """
        Compute the right-hand side vector to the set of equations that must be solved for the
        diffusive response function Sonine polynomial expansion coefficients.

        :param particle_density (float) : Particle density (not molar!) [1 / m3]
        :param T (float) : Temperature [K]
        :param mole_fracs (list<float>) : Molar composition [-]
        :param N (Optional, int) : Enskog approximation order.

        :return (1Darray) : ($N N_c^2$,) vector, where $N$ is the Enskog approximation order and $N_c$ is
                            the number of components.
        """
        if N is None:
            N = self.default_N
        self.check_valid_composition(mole_fracs)
        return np.array(self.cpp_kingas.get_diffusion_vector(particle_density, T, mole_fracs, N))

    #####################################################
    #          Evaluation of Matrix equations           #
    #####################################################

    def compute_diffusion_coeff_vector(self, particle_density, T, mole_fracs, N=None):
        """
        Compute the diffusive response function Sonine polynomial expansion coefficients by solving the set of equations

        $$D d = \delta$$

        Corresponding to Equations ... in ...
        Where $D$ is the matrix returned by the c++ method `get_diffusion_matrix`, and $\delta$ is the vector
        returned by the c++ method `get_diffusion_vector`.

        :param particle_density (float) : Particle density (not molar!) [1 / m3]
        :param T (float) : Temperature [K]
        :param mole_fracs (list<float>) : Molar composition [-]
        :param N (Optional, int) : Enskog approximation order.

        :return (1Darray) : ($N N_c^2$,) vector, where $N$ is the Enskog approximation order and $N_c$ is
                            the number of components.
        """
        if N is None:
            N = self.default_N
        self.check_valid_composition(mole_fracs)
        if (T, particle_density, tuple(mole_fracs), N) in self.computed_d_points.keys():
            return np.array(self.computed_d_points[(T, particle_density, tuple(mole_fracs), N)])

        diffusion_matr = self.cpp_kingas.get_diffusion_matrix(particle_density, T, mole_fracs, N)
        diffusion_vec = self.get_diffusion_vector(particle_density, T, mole_fracs, N=N)

        if any(np.isnan(np.array(diffusion_matr).flatten())):
            warnings.warn('Diffusion-matrix contained NAN elements!')
            d = np.array([np.nan for _ in diffusion_vec])
        else:
            d = lin.solve(diffusion_matr, diffusion_vec)

        self.computed_d_points[(T, particle_density, tuple(mole_fracs), N)] = tuple(d)
        return d

    def reshape_diffusion_coeff_vector(self, d):
        """
        The vector returned by `compute_diffusion_coeff_vector` contains coefficients ordered as in equation ...
        in ...  (eg. $d_{i, j}^{(q)}$ ). To more easily access the correct coefficients, this method reshapes the vector to a matrix indiced
        as `d[i][q][j]` where i and j are component indices, and q refferes to the approximation order summation index.

        :param d (1D array) : The array returned by `get_diffusion_coeff_vector`
        :return (3D array) : The matrix of $d_{i, j}^{(q)}$ coefficients ordered as `d[i][q][j]`
        """

        N = len(d) // (self.ncomps * self.ncomps)
        matr = np.zeros((self.ncomps, N, self.ncomps))
        for i in range(self.ncomps):
            for q in range(N):
                for k in range(self.ncomps):
                    matr[i, q, k] = d[N * self.ncomps * k + self.ncomps * q + i]
        return matr
    
    def compute_dth_vector(self, particle_density, T, mole_fracs, N=None):
        """
        Compute the coefficients $d_i^{(J = 0)}$, by solving the set of equations

        $$\sum_j d_{i,j}^{(0)} d_j^{\vec{J} = 0} = \ell_i^{(0)}$$

        :param particle_density (float) : Particle density (not molar!) [1 / m3]
        :param T (float) : Temperature [K]
        :param mole_fracs (list<float>) : Molar composition [-]
        :param N (Optional, int) : Enskog approximation order.

        :return (1Darray) : Vector of the $d_i^{(J = 0)}$ coefficients.
        """
        if N is None:
            N = self.default_N
        self.check_valid_composition(mole_fracs)
        a = self.compute_cond_vector(particle_density, T, mole_fracs, N=N)[:self.ncomps]
        d = self.compute_diffusion_coeff_vector(particle_density, T, mole_fracs, N=N)
        d = self.reshape_diffusion_coeff_vector(d)

        Dth = np.ones((self.ncomps, self.ncomps))
        for i in range(self.ncomps - 1): # last row is used for sum(d_th) = 0 condition
            for j in range(self.ncomps):
                Dth[i, j] = d[i][0][j]

        Dth[-1] *= max(abs(Dth).flatten()) # Conditioning the matrix
        a[-1] = 0 # sum(d_th) = 0 condition
        if any(np.isnan(Dth.flatten())):
            warnings.warn('Dth matrix contained NAN elements')
            dth = np.full_like(a, np.nan)
        else:
            dth = lin.solve(Dth, a)

        return dth

    def compute_cond_vector(self, particle_density, T, mole_fracs, N=None):
        """
        Compute the thermal response function Sonine polynomial expansion coefficients by solving the set of equations

        $$\Lambda \ell = \lambda$$

        Corresponding to Equations ... in ...
        Where $\Lambda$ is the matrix returned by the c++ method `get_conductivity_matrix`, and $\lambda$ is the vector
        returned by the c++ method `get_conductivity_vector`.

        :param particle_density (float) : Particle density (not molar!) [1 / m3]
        :param T (float) : Temperature [K]
        :param mole_fracs (list<float>) : Molar composition [-]
        :param N (Optional, int) : Enskog approximation order.

        :return (1Darray) : ($N N_c$,) vector, where $N$ is the Enskog approximation order and $N_c$ is
                            the number of components.
        """

        if N is None:
            N = self.default_N
        self.check_valid_composition(mole_fracs)
        if (T, particle_density, tuple(mole_fracs), N) in self.computed_a_points.keys():
            return np.array(self.computed_a_points[(T, particle_density, tuple(mole_fracs), N)])

        L = self.cpp_kingas.get_conductivity_matrix(particle_density, T, mole_fracs, N)
        L = np.array(L)
        L[0] *= max(L.flatten()) / abs(min(L[0, :2]))
        l = self.cpp_kingas.get_conductivity_vector(particle_density, T, mole_fracs, N)
        l = np.array(l)
        if any(np.isnan(np.array(L).flatten())):
            warnings.warn('A-matrix contained NAN elements!')
            a = np.array([np.nan for _ in l])
        else:
            a = lin.solve(L, l)

        self.computed_a_points[(T, particle_density, tuple(mole_fracs), N)] = tuple(a)
        return np.array(a)
    
    def compute_visc_vector(self, T, particle_density, mole_fracs, N=None):
        """
        Compute the viscous response function Sonine polynomial expansion coefficients by solving the set of equations

        $$\Beta b = \beta$$

        Corresponding to Equations ... in ...
        Where $\Beta$ is the matrix returned by the c++ method `get_viscosity_matrix`, and $\beta$ is the vector
        returned by the c++ method `get_viscosity_vector`.

        :param particle_density (float) : Particle density (not molar!) [1 / m3]
        :param T (float) : Temperature [K]
        :param mole_fracs (list<float>) : Molar composition [-]
        :param N (Optional, int) : Enskog approximation order.

        :return (1Darray) : ($N N_c$,) vector, where $N$ is the Enskog approximation order and $N_c$ is
                            the number of components.
        """
        if N is None:
            N = self.default_N
        self.check_valid_composition(mole_fracs)
        if (T, particle_density, tuple(mole_fracs), N) in self.computed_b_points.keys():
            return self.computed_b_points[(T, particle_density, tuple(mole_fracs), N)]
        
        B = self.cpp_kingas.get_viscosity_matrix(particle_density, T, mole_fracs, N)
        beta = self.cpp_kingas.get_viscosity_vector(particle_density, T, mole_fracs, N)

        if any(np.isnan(np.array(B).flatten())):
            warnings.warn('Viscosity matrix contained NAN elements!')
            b = np.array([np.nan for _ in beta])
        else:
            b = lin.solve(B, beta)

        self.computed_b_points[(T, particle_density, tuple(mole_fracs), N)] = tuple(b)
        return b