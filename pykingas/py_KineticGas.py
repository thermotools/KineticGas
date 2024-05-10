'''
Author: Vegard Gjeldvik Jervell
Purpose: Parent class for python-wrappers to the KineticGas models. This class takes care of matrix inversion 
        and the final computations that yield transport coefficients, as well as reading fluid parameters from
        the 'fluids/XX.json' files.
'''
import copy

import numpy as np
import json
import scipy.linalg as lin
from scipy.constants import Boltzmann, Avogadro, pi, gas_constant
import warnings, os

FLT_EPS = 1e-12
__fluid_db_path__ = os.path.dirname(__file__) + '/fluids/'

kB = Boltzmann

def k_delta(i, j): # Kronecker delta
    if i == j:
        return 1
    return 0

def compress_diffusion_matr(M, dependent_idx):
    """Utility
    Remove the dependent row and column from a diffusion matrix, returning an (N - 1 x N - 1) matrix.
    &&
    Args:
        M (array_like) : Diffusion matrix, shape (N, N)
        dependent_idx (int) : Index of the dependent species
    Returns:
        array : (N - 1) x (N - 1) array of independent diffusion coefficients, where N is the number of components.
    """
    M_indep = np.empty((len(M) - 1, len(M) - 1))
    M_i = 0
    for i in range(len(M)):
        if i == dependent_idx:
            continue
        M_j = 0
        for j in range(len(M)):
            if j == dependent_idx:
                continue
            M_indep[M_i][M_j] = M[i][j]
            M_j += 1
        M_i += 1
    return M_indep

class IdealGas:
    """
    Class to use for the eos attribute when is_idealgas == True.
    Also serves as an example of a minimal implementation of an eos-object that can
    be used if one wishes to supply a custom eos through the use_eos kwarg (See: MieKinGas)
    """
    def __init__(self, comps):
        self.ncomps = len(comps.split(','))
        self.VAPPH = 1 # Required to be compatible with a ThermoPack-style eos object. The value of the flag can be whatever you want

    def chemical_potential_tv(self, T, V, n, dmudn=True):
        """
        The chemical potential of an ideal gas mixture. Note: Only dmudn is actually used by the py_KineticGas class.
        Also note: We need to have V in the signature, even though it is not used, in order to be compatible with
        the signature of chemical_potential_tv in thermopack.
        &&
        Args:
            T (float) : Temperature [K]
            V (float) : Volume [m3]
            n (list[float]) : Mole numbers [mol]
            dmudn (bool) : Flag to activate derivative calculation
        Returns
            tuple : chemical potential [J / mol] and derivatives
        """
        n = np.array(n)

        if dmudn is True:
            dmudn = np.identity(self.ncomps) * gas_constant * T / n
            return 0, dmudn
        raise NotImplementedError('Only dmudn is implemented, because that is all we need for the pykingas package.')

    def specific_volume(self, T, p, n, phase, dvdn=False):
        """
        Compute molar volume for an ideal gas. Note that we must have `n` and `phase` in the signature in
        order to be compatible with the signature of equation of state objects from ThermoPack. Partial molar volumes
        must be implemented in order to use the solvent frame of reference.
        &&
        Args:
            T (float) : Temperature [K]
            p (float) : Pressure [Pa]
            n (list[float]) : mole numbers [mol]
            phase (int) : phase flag (see: ThermoPack)
            dvdn (bool) : Compute partial molar volumes

        Returns:
            (float,) : molar volume [m3 / mol]
        """
        v = gas_constant * T / p
        if dvdn is True:
            return v, np.array([v for _ in range(self.ncomps)])
        return v, # Because ThermoPack v2.1.0 returns everything as a tuple, we need to return a tuple.

    def pressure_tv(self, T, Vm, x):
        """
        Compute pressure for an ideal gas, we must take `x` as an argument to be compatible with the generic ThermoPack
        signature. This method is required for the solvent `frame_of_reference` diffusion and thermal diffusion coefficients.
        &&
        Args:
            T (float) : Temperature [K]
            Vm (float) : molar volume [m3 / mol]
            x (list[float]) : mole fractions [-]

        Returns:
            (tuple) : (Pressure,)
        """
        return gas_constant * T / Vm

class py_KineticGas:

    def __init__(self, comps, mole_weights=None, N=3, is_idealgas=False):
        """Constructor
        &&
        Args:
            comps (str): Comma-separated list of components, following ThermoPack-convention
            mole_weights (1d array) : Mole weights [g/mol]. Will be used instead of database values if provided
            is_idealgas (bool) : If true, radial distribution function is unity, if false use radial distribution function of model
                                In addition, several density-dependent factors are set to zero, to ensure consistency with
                                the ideal gas law / Gibbs-Duhem for an ideal gas.
        """
        self._is_singlecomp = False
        self.default_N = N
        self.is_idealgas = is_idealgas
        self.computed_d_points = {} # dict of state points in which (d_1, d0, d1) have already been computed
        self.computed_a_points = {}  # dict of state points in which (a_1, a1) have already been computed
        self.computed_b_points = {}  # dict of state points in which (b_1, b1) have already been computed
        self.computed_kT = {} # dict of state points in which thermal diffusion ratios have been computed
        self.computed_cd = {} # dict of state points in which the collision diameters have been computed
        self.computed_rdf = {}  # dict of state points in which the RDFs have been computed

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
        if self.is_idealgas is True:
            self.eos = IdealGas(comps)
        else:
            self.eos = None

    #####################################################
    #                   Utility                         #
    #####################################################
    def check_valid_composition(self, x):
        """Utility
        Check that enough mole fractions are supplied for the initialised model. Also check that they sum to unity.
        &&
        Args:
            x (array_like) : Molar composition

        Raises
            IndexError : If wrong number of mole fractions is supplied.
            RuntimeWarning : If mole fractions do not sum to unity.
        """
        if abs(sum(x) - 1) > FLT_EPS:
            warnings.warn('Mole fractions do not sum to unity, sum(x) = '+str(sum(x)), RuntimeWarning, stacklevel=2)
        elif len(x) != self.ncomps:
            raise IndexError(str(len(x)) +' mole fractions were supplied for a '+str(self.ncomps)+' component mixture!\n'
                            'Note: Single-component mixtures are treated as binary! (use x = [0.5, 0.5])')


    def get_Eij(self, Vm, T, x):
        r"""Utility
        Compute the factors

        $$ ( n_i / k_B T ) (d \mu_i / d n_j)_{T, n_{k \neq j}}, $$

        where $n_i$ is the molar density of species $i$.
        &&
        Args:
            Vm (float) : Molar volume [m3 / mol]
            T (float) : Temperature [K]
            x (array_like) : Molar composition

        Returns:
            (2D array) : The factors E[i][j] = $ ( n_i / k_B T ) (d \mu_i / d n_j)_{T, n_{k \neq j}}$, where $n_i$
                                is the molar density of species $i$. Unit [1 / mol]
        """
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
        r"""Utility
        Compute the factors $\Xi_i = \sum_j E_{ij}$, where $E_{ij}$ are the factors computed by `get_Eij`.
        &&
        Args:
            Vm (float) : Molar volume [m3 / mol]
            T (float) : Temperature [K]
            x (array_like) : Molar composition

        Returns:
            (1D array) : The factors $\Xi_i$, Unit [1 / mol]
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
        r"""TV-property
        Compute the interdiffusion coefficients [m^2 / s]. Default definition is

        $$ J_i^{(n, n)} = - \sum_{j \neq l} D_{ij} \nabla n_j, \nabla T = \nabla p = F_k = 0 \forall k $$

        where the flux, $J_i^{(n, n)}$ is on a molar basis, in the molar frame of reference, and $j \neq l$ is an
        independent set of forces with $l=$ `dependent_idx`.
        For fluxes in other frames of reference, use the `frame_of_reference` kwarg.
        For the diffusion coefficients describing the fluxes' response to all forces (not independent) the definition
        is:

        $$ J_i^{(n, n)} = - \sum_j D_{ij} \nabla n_j, \nabla T = \nabla p = F_k = 0 \forall k $$

        use the `use_independent` and `dependent_idx` kwargs to switch between these definitions.
        See: Eq. (17-20) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
        &&
        Args:
            T (float) : Temperature (K)
            Vm (float) : Molar volume (m3 / mol)
            x (array_like) : composition (mole fractions)
            N (int, optional) : Enskog approximation order. Default set on model initialisation.
            use_binary (bool, optional) : If the mixture is binary, and an independent set of fluxes and forces is considered, i.e.
                `use_independent=True`, the diffusion coefficients will be exactly equal with opposite sign. Setting
                `use_binary=True` will in that case only return the coefficient describing the independent flux-force relation.
            use_independent (bool, optional) : Return diffusion coefficients for independent set of forces.
            dependent_idx (int, optional) : Index of the dependent molar density gradient (only if `use_independent=True`, the default behaviour).
                Defaults to last component, except when `frame_of_reference='solvent'`, in which case default is equal
                to `solvent_idx`.
            frame_of_reference (str, optional) : Which frame of reference the diffusion coefficients apply to. Default
                is `'CoN'`. Can be `'CoN'` (molar FoR), `'CoM'` (barycentric FoR), `'solvent'` (solvent FoR), `'zarate'` (See Memo on
                definitions of the diffusion coefficient), `'zarate_x'` ($D^{(x)}$ as defined by Ortiz de Zárate, doi 10.1140/epje/i2019-11803-2)
                or `'zarate_w'` ($D^{(w)}$ as defined by Ortiz de Zárate).
            solvent_idx (int, optional) : Index of component identified as solvent (only when using `frame_of_reference='solvent'`)

        Returns:
            (ndarray or float) : Diffusion coefficients, shape varies based on options and number of components. Unit [m^2 / s]
        """
        if dependent_idx is None:
            dependent_idx = self.ncomps - 1
        while dependent_idx < 0:
            dependent_idx += self.ncomps

        D = self.interdiffusion_general(T, Vm, x, N=N)
        # psi = Transformation matrix from 'centre of mass' to 'frame_of_reference'
        # get_com_2_for_matr() dispatches the call to specific functions for different frames of reference.
        if frame_of_reference == 'zarate_x':
            D = self.interdiffusion(T, Vm, x, N=N, frame_of_reference='CoN', dependent_idx=dependent_idx, use_independent=True, use_binary=False)
            return compress_diffusion_matr(D, dependent_idx)
        elif frame_of_reference == 'zarate':
            X = self.get_zarate_X_matr(x, dependent_idx)
            D_x = self.interdiffusion(T, Vm, x, N=N, frame_of_reference='zarate_x', dependent_idx=dependent_idx, use_binary=False)
            return X @ D_x @ np.linalg.inv(X)
        elif frame_of_reference == 'zarate_w':
            W = self.get_zarate_W_matr(x, dependent_idx)
            D_z = self.interdiffusion(T, Vm, x, N=N, frame_of_reference='zarate', dependent_idx=dependent_idx, use_binary=False)
            return W @ D_z @ np.linalg.inv(W)

        psi = self.get_com_2_for_matr(T, Vm, x, frame_of_reference, solvent_idx=solvent_idx)
        D = psi @ D
        if use_independent is True:
            if dependent_idx is None:
                if frame_of_reference == 'solvent':
                    dependent_idx = solvent_idx
                else:
                    dependent_idx = self.ncomps - 1
            P = self.get_P_factors(Vm, T, x)
            D_dep = copy.deepcopy(D)
            for i in range(self.ncomps):
                for j in range(self.ncomps):
                    D[i, j] -= (P[j] / P[dependent_idx]) * D_dep[i, dependent_idx]
            if use_binary is True and self.ncomps == 2:
                independent_idx = 1 - dependent_idx
                return D[independent_idx][independent_idx] # Return the independent fluxes response to the independent force
        return D

    def interdiffusion_general(self, T, Vm, x, N=None):
        r"""TV-property
        Compute the 'Kinetic CoM diffusion coefficients', defined by

        $$ J_i^{(n, m)} = - \sum_j D_{ij} \nabla n_j, \nabla T = \nabla p = F_k = 0 \forall k $$

        **For end-users, see `interdiffusion`**
        See Eq. (19) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
        &&
        Args:
            T (float) : Temperature (K)
            Vm (float) : Molar volume (m3 / mol)
            x (array_like) : composition (mole fractions)
            N (int, optional) : Enskog approximation order

        Returns:
            (2D array) : Array of the (not independent) $D^{(K, m)}$ diffusion coefficients. Unit [m^2 / s]
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
        r"""TV-Property
        Compute thermal diffusion coefficients, $D_{T,i}$ [mol / m^2 s]
        Default definition is

        $$ J_i^{(n, n)} = D_{T,i} \nabla \ln T - \sum_j D_{ij} \nabla n_j, \nabla p = F_k = 0 \forall k $$

        where the flux, $J_i^{(n, n)}$ is on a molar basis, in the molar frame of reference. For fluxes in other frames
        of reference, use the 'frame_of_reference' kwarg. For the diffusion coefficients corresponding to an
        independent set of forces, defined by

        $$ J_i^{(n, n)} = D_{T,i} \nabla \ln T - \sum_{j \neq l} D_{ij} \nabla n_j, \nabla p = F_k = 0 \forall k $$

        where $l$ is the index of the dependent molar density gradient, use the 'use_independent' and 'dependent_idx' kwargs.
        See Eq. (23) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
        &&
        Args:
            T (float) : Temperature [K]
            Vm (float) : Molar volume [m^3 / mol]
            x (array_like) : Mole fractions [-]
            N (int, optional) : Enskog approximation order (>=2). Defaults to 2.
            use_independent (bool, optional) : Return diffusion coefficients for independent set of forces.
            dependent_idx (int, optional) : Index of the dependent molar density gradient (only if use_dependent=True).
                                Defaults to last component.
            frame_of_reference (str, optional) : What frame of reference the coefficients apply to. Valid options are
                                        `'CoM'` (centre of mass / barycentric), `'CoN'` (centre of moles), `'CoV'` (centre of volume)
                                        `'solvent'` (together with `solvent_idx`) or `'zarate'`, for the coefficients as
                                        defined by Ortiz de Zarate (doi 10.1140/epje/i2019-11803-2).
            solvent_idx (int, optional) : Index of component identified as solvent (only when using `frame_of_reference='solvent'`)

        Returns:
            (1D array) : Thermal diffusion coefficients. Unit [mol m^2 / s]
        """
        if N is None:
            N = self.default_N
        if dependent_idx is None:
            dependent_idx = self.ncomps - 1
        while dependent_idx < 0:
            dependent_idx += self.ncomps

        self.check_valid_composition(x)

        if N < 2:
            warnings.warn('Thermal diffusion is a 2nd order phenomena, cannot be computed for N < 2 (got N = '
                          + str(N) + ')', RuntimeWarning, stacklevel=2)
            return np.full(self.ncomps, np.nan)

        use_zarate = False
        if frame_of_reference.lower() == 'zarate':
            frame_of_reference = 'CoN'
            use_zarate = True
            use_independent = True

        particle_density = Avogadro / Vm
        d = self.compute_diffusion_coeff_vector(particle_density, T, x, N=N)
        d = self.reshape_diffusion_coeff_vector(d)
        a = self.compute_cond_vector(particle_density, T, x, N=N)
        P = self.get_P_factors(Vm, T, x)
        rdf = self.get_rdf(particle_density, T, x)
        cd = self.get_collision_diameters(particle_density, T, x)

        b = np.empty((self.ncomps, self.ncomps)) # Precomputing some factors that are used many places later
        if self.is_idealgas is True:
            b = np.identity(self.ncomps)
        else:
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
                        tmp += particle_density * x[m] * b[m, k]

                DT[i] += (Dij[i, dependent_idx] / P[dependent_idx]) * tmp

        DT /= Avogadro # Dividing by Avogadros number to convert unit from particle number to mole number
        if use_zarate is True:
            D = self.interdiffusion(T, Vm, x, N=N, use_independent=True, dependent_idx=dependent_idx, frame_of_reference='CoN',
                                    use_binary=False)
            DT_indep = np.empty(self.ncomps - 1)
            D_indep = compress_diffusion_matr(D, dependent_idx)
            DT_idx = 0
            x_factor = np.empty(self.ncomps - 1)

            for i in range(self.ncomps):
                if i == dependent_idx:
                    continue
                x_factor[DT_idx] = x[i]
                DT_indep[DT_idx] = DT[i]
                DT_idx += 1
            X = self.get_zarate_X_matr(x, dependent_idx)
            c = 1 / Vm
            DT = np.linalg.solve(- c * X, (DT_indep + c * D_indep @ x_factor) / T)

        return DT

    def thermal_diffusion_ratio(self, T, Vm, x, N=None):
        r"""TV-property
        Calculate the "independent" thermal diffusion ratios, $k_{T, i}$ defined by

        $$ J_i^{(n, n)} = - \sum_{j \neq l} D_{ij}^{(I, n)} ( \nabla n_j + n_j k_{T, j} D_{T, j}^{(I, n)} \nabla \ln T ) $$

        and

        $$ \sum_i x_i \sum_j [\delta_{ij} + (4 \pi / 3) \sigma_{ij}^3 n_j M_{ij} \chi_{ij} - (n_j k_{T,j} / k_B T) ( d \mu_i / n_j )_{T,n_{l \neq j}}] = 0 $$

        This definition implies that

        $$ \nabla n_j = -n_j k_{T,j} \nabla \ln T\ \forall j $$

        when all mass fluxes vanish. For models initialised with `is_idealgas=True`, the second equation is replaced with

        $$ \sum_i x_i k_{T,i} = 1 $$

        The thermal diffusion ratios are independent of the frame of reference.
        See Eq. (26-27) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
        &&
        Args:
            T (float) : Temperature [K]
            Vm (float) : Molar volume [m3 / mol]
            x (array_like) : Molar composition [-]
            N (int, optional) : Enskog approximation order (>= 2)

        Returns:
            (1D array) : The thermal diffusion ratio of each component. Unit Dimensionless.
        """
        if N is None:
            N = self.default_N

        key = tuple((T, Vm, tuple(x), N))
        if key in self.computed_kT.keys():
            return self.computed_kT[key]

        if N < 2:
            warnings.warn('Thermal diffusion is a 2nd order phenomena, cannot be computed for N < 2 (got N = '
                          + str(N) + ')', RuntimeWarning, stacklevel=2)
            return np.full(self.ncomps, np.nan)

        particle_density = Avogadro / Vm
        DT = self.thermal_diffusion_coeff(T, Vm, x, N=N, frame_of_reference='CoM',
                                          use_independent=True, dependent_idx=self.ncomps - 1)
        Dij = self.interdiffusion(T, Vm, x, N=N, frame_of_reference='CoM',
                                  use_binary=False, use_independent=True, dependent_idx=self.ncomps - 1)
        rdf = self.get_rdf(particle_density, T, x)
        cd = self.get_collision_diameters(particle_density, T, x)
        P = self.get_P_factors(Vm, T, x)
        A = np.zeros((self.ncomps, self.ncomps))

        for i in range(self.ncomps - 1):
            for j in range(self.ncomps):
                A[i, j] = - Dij[i, j] * x[j] * (1 / Vm) # density in moles, because DT has unit moles (not particles)

        for i in range(self.ncomps):
            A[-1, i] = x[i] * P[i]

        # Overwriting the DT vector for the final element of the b-vector of A @ kT = b
        if self.is_idealgas is True:
            DT[-1] = 1
        else:
            DT[-1] = 0
            for i in range(self.ncomps):
                for j in range(self.ncomps):
                    DT[-1] += x[i] * (k_delta(i, j) + (4 * np.pi / 3) * particle_density * x[j] * cd[i][j]**3 * self.M[i, j] * rdf[i][j])

        kT = np.linalg.solve(A, DT)

        self.computed_kT[key] = tuple(kT)
        return kT

    def thermal_diffusion_factor(self, T, Vm, x, N=None):
        r"""TV-property
        Compute the thermal diffusion factors $\alpha_{ij}$, defined by

        $$ \alpha_{ij} = k_{T, i} - k_{T, j} $$

        where $k_{T,i}$ are the thermal diffusion ratios. This definition implies that

        $$ \nabla \ln (n_i / n_j) = - \alpha_{ij} \nabla \ln T $$

        when the mass fluxes vanish. The thermal diffusion factors are independent of the frame of reference.
        See Eq. (29) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
        &&
        Args:
            T (float) : Temperature [K]
            Vm (float) : Molar volume [m3 / mol]
            x (array_like) : Molar composition [-]
            N (int, optional) : Enskog approximation order (>= 2)

        Returns:
            (2D array) : The thermal diffusion factors of the mixture. Dimensionless.
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

    def soret_coefficient(self, T, Vm, x, N=None, use_zarate=True, dependent_idx=-1):
        r"""TV-Property
        Compute the Soret coefficients, $S_{T, ij}$. If `use_zarate=False`, the Soret coefficient is defined by

        $$ S_{T, ij} = \alpha_{ij} / T $$

        where $\alpha_{ij}$ are the thermal diffusion factors. If `use_zarate=True`, uses the definition proposed by
        Ortiz de Zarate in (doi 10.1140/epje/i2019-11803-2), i.e.

        $$ X^{-1} D^{(x)} X (S_T) = (D_T) $$

        where $(S_T)$ and $(D_T)$ indicate the vectors of Soret coefficients and thermal diffusion coefficients.
        Or, following the notation in the memo on definitions of diffusion and thermal diffusion coefficients,

        $$ D^{(z)} (S_T) = (D_T). $$

        The Soret coefficients defined this way satisfy

        $$ X (S_T) \nabla T = - (\nabla x) $$

        and

        $$ W (S_T) \nabla T = - (\nabla w) $$

        which in a binary mixture reduces to the same as the definition used when `use_zarate=False`, i.e.

        $$ S_T = - \frac{\nabla x_1}{x_1 (1 - x_1) \nabla T}$$

        if species 2 is the dependent species.

        &&
        Args:
            T (float) : Temperature [K]
            Vm (float) : Molar volume [m3 / mol]
            x (array_like) : Molar composition
            N (int, optional) : Enskog approximation order (>= 2)
            use_zarate (bool, optional) : Use Ortiz de Zarate formulation of the Soret coefficient, as given in the
                                        method description and (doi 10.1140/epje/i2019-11803-2). Defaults to `True`.
            dependent_idx (int) : Only applicable when `use_zarate=True` (default behaviour). The index of the dependent
                                    species. Defaults to the last species.
        """
        if N is None:
            N = self.default_N

        if N < 2:
            warnings.warn('Thermal diffusion is a 2nd order phenomena, cannot be computed for N < 2 (got N = '
                          + str(N) + ')', RuntimeWarning, stacklevel=2)
            return np.full((self.ncomps, self.ncomps), np.nan)

        if use_zarate is True:
            D = self.interdiffusion(T, Vm, x, N=N, frame_of_reference='zarate', dependent_idx=dependent_idx, use_binary=False)
            DT = self.thermal_diffusion_coeff(T, Vm, x, N=N, frame_of_reference='zarate', dependent_idx=dependent_idx)
            return np.linalg.solve(D, DT)

        alpha = self.thermal_diffusion_factor(T, Vm, x, N=N)
        return alpha / T

    def thermal_conductivity(self, T, Vm, x, N=None, idealgas=None, include_internal=True, contributions='all', p=None):
        r"""TV-Property
        Compute the thermal conductivity, $\lambda$. For models initialized with `is_idealgas=True`, the thermal
        conductivity is not a function of density (i.e. $d \lambda / d V_m = 0$).
        See Eq. (13) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
        &&
        Args:
            T (float) : Temperature [K]
            Vm (float) : Molar volume [m3 / mol]
            x (array_like) : Molar composition [-]
            N (int, optional) : Enskog approximation order (>= 2)
            idealgas (bool, optional) : Return infinite dilution value? Defaults to model default (set on init).
            include_internal (bool, optional) : Include contribution from internal degrees of freedom, computed using
                                                Eucken equation (doi.org/10.6028/NIST.IR.8209). Defaults to True.
            contributions (str, optional) : Return only specific contributions, can be ('all', '(i)nternal',
                                            '(t)ranslational', '(d)ensity', or several of the above, such as 'tid' or 'td'.
                                            If several contributions are selected, these are returned in an array of
                                            contributions in the same order as indicated in the supplied flag.

        Returns:
            (float) : The thermal conductivity of the mixture.
        """
        if N is None:
            N = self.default_N
        if idealgas is None:
            idealgas = self.is_idealgas

        self.check_valid_composition(x)
        particle_density = Avogadro / Vm
        a = self.compute_cond_vector(particle_density, T, x, N=N)
        rdf = self.get_rdf(particle_density, T, x)
        K = self.cpp_kingas.get_K_factors(particle_density, T, x)
        cd = self.get_collision_diameters(particle_density, T, x)
        d = self.compute_diffusion_coeff_vector(particle_density, T, x, N=N)
        d = self.reshape_diffusion_coeff_vector(d)

        lambda_prime = 0
        dth = self.compute_dth_vector(particle_density, T, x, N=N)
        for i in range(self.ncomps):
            tmp = 0
            for k in range(self.ncomps):
                tmp += d[i, 1, k] * dth[k]
            lambda_prime += x[i] * K[i] * (a[self.ncomps + i] - tmp)
        lambda_prime *= (5 * Boltzmann / 4)

        lambda_int = 0
        if include_internal is True:
            f_int = 1.32e3
            eta_0 = self.viscosity(T, Vm, x, N=N, idealgas=True)
            Cp = 0
            M = sum(self.mole_weights * x) * Avogadro * 1e3
            for i in range(self.ncomps):
                _, Cpi_id = self.eos.idealenthalpysingle(T, 1 if self._is_singlecomp else i + 1, dhdt=True)
                # _, Cpi_id = self.eos.idealenthalpysingle(T, 1, dhdt=True)
                # Mi = self.mole_weights[i] * Avogadro * 1e3  # Mole weight in g / mol
                Cp += x[i] * Cpi_id

            Cp_factor = (Cp - 5 * gas_constant / 2) / M
            lambda_int = f_int * eta_0 * Cp_factor

        lambda_dblprime = 0
        if idealgas is False:  # lambda_dblprime is only nonzero when density corrections are present, and vanishes at infinite dilution
            for i in range(self.ncomps):
                for j in range(self.ncomps):
                    lambda_dblprime += particle_density ** 2 * np.sqrt(
                        2 * pi * self.m[i] * self.m[j] * Boltzmann * T / (self.m[i] + self.m[j])) \
                                       * (x[i] * x[j]) / (self.m[i] + self.m[j]) * (cd[i][j] ** 4) * rdf[i][j]
            lambda_dblprime *= (4 * Boltzmann / 3)

        cond = lambda_prime + lambda_int + lambda_dblprime
        if contributions == 'all':
            return cond

        contribs = {'t' : lambda_prime, 'i' : lambda_int, 'd' : lambda_dblprime}
        if len(contributions) > 1:
            return np.array([contribs[c] for c in contributions])
        else:
            return contribs[contributions]

    def viscosity(self, T, Vm, x, N=None, idealgas=None):
        r"""TV-Property
        Compute the shear viscosity, $\eta$. For models initialized with `is_idealgas=True`, the shear viscosity
        is not a function of density (i.e. $d \eta / d V_m = 0). See Eq. (12) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
        &&
        Args:
            T (float) : Temperature [K]
            Vm (float) : Molar volume [m3 / mol]
            x (array_like) : Molar composition [-]
            N (int, optional) : Enskog approximation order
            idealgas (bool, optional) : Use infinite dilution value? Defaults to model default value (set on init)

        Returns:
            (float) : The shear viscosity of the mixture.
        """
        if N is None:
            N = self.default_N
        if idealgas is None:
            idealgas = self.is_idealgas
        self.check_valid_composition(x)
        particle_density = Avogadro / Vm
        b = self.compute_visc_vector(T, particle_density, x, N=N)
        K_prime = self.cpp_kingas.get_K_prime_factors(particle_density, T, x)

        eta_prime = 0
        for i in range(self.ncomps):
            eta_prime += K_prime[i] * x[i] * b[i]
        eta_prime *= Boltzmann * T / 2

        eta_dblprime = 0
        if idealgas is False: # eta_dblprime is only nonzero when density corrections are present, and vanish at infinite dilution
            cd = self.get_collision_diameters(particle_density, T, x)
            rdf = self.get_rdf(particle_density, T, x)
            for i in range(self.ncomps):
                for j in range(self.ncomps):
                    eta_dblprime += np.sqrt(self.m[i] * self.m[j] / (self.m[i] + self.m[j])) * x[i] * x[j] * cd[i][j]**4 * rdf[i][j]
            eta_dblprime *= 4 * particle_density**2 * np.sqrt(2 * np.pi * Boltzmann * T) / 15

        return eta_prime + eta_dblprime

    def bulk_viscosity(self, T, Vm, x, N=None):
        """TV-property
        Not implemented

        Raises:
            NotImplementedError
        """
        if N is None:
            N = self.default_N
        self.check_valid_composition(x)
        raise NotImplementedError("Bulk viscosity is not implemented yet. See 'Multicomp docs' for more info.")

    def conductivity_matrix(self, T, Vm, x, N=2, formulation='T-psi', frame_of_reference='CoM', use_thermal_conductivity=None):
        r"""TV-Property
        Compute the conductivity matrix $L$, for use in NET calculations. The Flux/Force formulation used in the NET
        model is selected using the `formulation` kwarg. Currently implemented formulations are:
        ----------------------------------------------------------------------------------------------

        `'T-psi'`:

        ----------------------------------------------------

        $$ J_q = L_{qq} \nabla (1 / T) - \sum_{i=1}^{N_c-1}(1 / T) L_{qi} \nabla_T \Psi_{i; p, x} $$

        $$ J_i = L_{iq} \nabla (1 / T) - \sum_{j=1}^{N_c-1}(1 / T) L_{ij} \nabla_T \Psi_{j; p, x} $$

        Where

        $$ \Psi_i = \mu_i - \mu_{N_c} $$

        and $N_c$ denotes the number of components. The last component is used as the dependent component.
        The fluxes in this formulation are on a *mass basis*, and the formulation is implemented in the centre of
        mass frame of reference. The formulation is only implemented for ideal gases.
        ----------------------------------------------------------------------------------------------
        &&
        Args:
            T (float) : Temperature [K]
            Vm (float) : Molar volume [m3 / mol]
            x (1darray) : Molar composition [-]
            N (int, optional) : Enskog approximation order (must be >= 2 for thermal effects)
            formulation (str, optional) : The NET formulation for which to compute the conductivity matrix.
            frame_of_reference (str, optional) : The frame of reference ('CoM', 'CoN', 'CoV', or 'solvent'). Default: 'CoM'.
            use_thermal_conductivity (callable, optional) : External thermal conductivity model. Assumed to have the signature
                    use_thermal_conductivity(T, Vm, x), returning the thermal conductivity in [W / m K]. Defaults to None.
                    If no model is supplied, KineticGas is used to compute thermal conductivity.

        Returns:
            ndarray : The conductivity matrix, contents will vary depending on the `formulation` kwarg.
        """
        if formulation == 'T-psi':
            L = np.zeros((self.ncomps, self.ncomps))
            rho = Avogadro / Vm
            d = self.compute_diffusion_coeff_vector(rho, T, x, N=N)
            d = self.reshape_diffusion_coeff_vector(d)

            d0 = d[:, 0, :]
            for i in range(self.ncomps - 1):
                for j in range(self.ncomps - 1):
                    L[i + 1][j + 1] = d0[i][j] * x[i] * x[j] * self.m[i] * self.m[j] / (2 * Boltzmann)

            k_T = self.thermal_diffusion_ratio(T, Vm, x, N=N)
            for i in range(self.ncomps - 1):
                for j in range(self.ncomps - 1):
                    L[0, i + 1] -= L[i + 1, j + 1] * Boltzmann * T * (((1 - k_T[j]) / self.m[j]) - ((1 - k_T[-1]) / self.m[-1]))

                L[i + 1, 0] = L[0, i + 1]

            if use_thermal_conductivity is None:
                cond = self.thermal_conductivity(T, Vm, x, N=N)
            else:
                cond = use_thermal_conductivity(T, Vm, x)

            L[0, 0] = T ** 2 * cond
            for i in range(self.ncomps - 1):
                L[0, 0] -= Boltzmann * T * L[0, i + 1] * (((1 - k_T[i]) / self.m[i]) - ((1 - k_T[-1]) / self.m[-1]))

            return L

        raise KeyError(f'Invalid formulation : {formulation}')

    def resistivity_matrix(self, T, Vm, x, N=2, formulation='T-psi', frame_of_reference='CoM', use_thermal_conductivity=None):
        r"""TV-property
        Compute the resistivity matrix $R = L^{-1}$, for use in NET calculations. The Flux/Force formulation used in the NET
        model is selected using the `formulation` kwarg. Currently implemented formulations are:
        ----------------------------------------------------------------------------------------------

        `'T-psi'`:

        ----------------------------------------------------

        $$ J_q = L_{qq} \nabla (1 / T) - \sum_{i=1}^{N_c-1}(1 / T) L_{qi} \nabla_T \Psi_i $$

        $$ J_i = L_{iq} \nabla (1 / T) - \sum_{j=1}^{N_c-1}(1 / T) L_{ij} \nabla_T \Psi_j $$

        Where

        $$ \Psi_i = \mu_i - \mu_{N_c} $$

        and $N_c$ denotes the number of components. The last component is used as the dependent component.
        The fluxes in this formulation are on a *mass basis*, and the formulation is implemented in the centre of
        mass frame of reference. The formulation is only implemented for ideal gases.
        ----------------------------------------------------------------------------------------------
        &&
        Args:
            T (float) : Temperature [K]
            Vm (float) : Molar volume [m3 / mol]
            x (1darray) : Molar composition [-]
            N (int, optional) : Enskog approximation order (must be >= 2 for thermal effects)
            formulation (str, optional) : The NET formulation for which to compute the conductivity matrix.
            frame_of_reference (str, optional) : The frame of reference ('CoM', 'CoN', 'CoV', or 'solvent'). Default: 'CoM'.
            use_thermal_conductivity (callable, optional) : External thermal conductivity model. Assumed to have the signature
                    thermal_conductivity(T, Vm, x), returning the thermal conductivity in [W / m K]. Defaults to None.
                    If no model is supplied, KineticGas is used to compute the thermal conductivity.

        Returns:
            ndarray : The resistivity matrix, contents will vary depending on the `formulation` kwarg.
        """
        if formulation == 'T-psi':
            L = self.conductivity_matrix(T, Vm, x, N=N, formulation=formulation, frame_of_reference=frame_of_reference,
                                         use_thermal_conductivity=use_thermal_conductivity)
            return np.linalg.inv(L)


        raise KeyError(f'Invalid formulation : {formulation}')

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
        """Tp-property
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
        """Tp-property
        Compute molar volume using the internal equation of state (`self.eos`), assuming vapour, and pass the call to
        `self.thermal_diffusion_coeff`. See `self.thermal_diffusion_coeff` for documentation.
        """
        Vm, = self.eos.specific_volume(T, p, x, self.eos.VAPPH)  # Assuming vapour phase
        return self.thermal_diffusion_coeff(T, Vm, x, N=N, use_independent=use_independent,
                                   dependent_idx=dependent_idx, frame_of_reference=frame_of_reference,
                                   solvent_idx=solvent_idx)

    def thermal_diffusion_factor_tp(self, T, p, x, N=None):
        """Tp-Property
        Compute molar volume using the internal equation of state (`self.eos`), assuming vapour, and pass the call to
        `self.thermal_diffusion_factor`. See `self.thermal_diffusion_factor` for documentation.
        """
        Vm, = self.eos.specific_volume(T, p, x, self.eos.VAPPH)  # Assuming vapour phase
        return self.thermal_diffusion_factor(T, Vm, x, N=N)

    def thermal_coductivity_tp(self, T, p, x, N=None):
        """Tp-property
        Compute molar volume using the internal equation of state (`self.eos`), assuming vapour, and pass the call to
        `self.thermal_conductivity`. See `self.thermal_conductivity` for documentation.
        """
        Vm, = self.eos.specific_volume(T, p, x, self.eos.VAPPH)  # Assuming vapour phase
        return self.thermal_conductivity(T, Vm, x, N=N, p=p)

    def viscosity_tp(self, T, p, x, N=None):
        """Tp-property
        Compute molar volume using the internal equation of state (`self.eos`), assuming vapour, and pass the call to
        `self.viscosity`. See `self.viscosity` for documentation.
        """
        Vm, = self.eos.specific_volume(T, p, x, self.eos.VAPPH)  # Assuming vapour phase
        return self.viscosity(T, Vm, x, N=N)

    #####################################################
    #        Frame of Reference Transformations         #
    #####################################################

    def get_com_2_for_matr(self, T, Vm, x, FoR, **kwargs):
        r"""FOR-Transform
        Dispatcher to get a specific 'change of frame of reference' matrix for transformation from
        centre of mass to 'FoR'.
        Returns the appropriate matrix for the transformations derived in Appendix A of ... to transform a flux
        from the centre of mass frame of reference to the 'FoR' frame of reference by the transformation

        $$ J^{(n, FoR)} = \psi @ J^{(n, m)} $$

        where $\psi$ is the matrix returned by this method, and $J$ is the vector of (all) molar fluxes, with the
        subscript indicating the frame of reference.
        &&
        Args:
            T (float) : Temperature [K]
            Vm (float) : Molar volume [m3 / mol]
            x (array_like) : Molar composition [-]

        Returns:
            (2Darray) : The $N$ x $N$ transformation matrix to transform the fluxes, with $N$ being the number of
                            components.
        """
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
        r"""FOR-Transform
        Get transformation matrix from centre of mass (CoM) to centre of moles (CoN).
        &&
        Args:
            x (array_like) : Molar composition [-]

        Returns:
            (2d array) : Transformation matrix $\Psi^{n \leftmapsto m}$
        """
        psi = np.identity(self.ncomps)
        w = x * self.m / sum(x * self.m)
        for i in range(self.ncomps):
            for j in range(self.ncomps):
                psi[i, j] += x[i] * ((w[j] / x[j]) - 1) #  (1 - w[j] * x[k] / (x[j] * w[k]))
        return psi

    def get_com_2_solv_matr(self, x, solvent_idx):
        r"""FOR-Transform
        Get transformation matrix from centre of mass (CoM) to solvent (solvent) frame of reference
        &&
        Args:
            x (array_like) : Molar composition [-]
            solvent_idx (int) : The component index of the solvent.

        Returns:
            2darray : The transformation matrix $\Psi^{n_k \leftmapsto m}$, where $k$ is the `solvent_idx`.
        """
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
        r"""FOR-Transform
        Get centre of mass (CoM) to centre of volume (CoV) transformation matrix
        &&
        Args:
            T (float) : Temperature [K]
            Vm (float) : Molar volume [m3 / mol]
            x (array_like) : Molar composition [-]

        Returns:
            2d array : The transformation matrix $\Psi^{V \leftmapsto m}$.
        """
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
        r"""FOR-Transform
        Get solvent-to-solvent frame of reference transformation matrix
        &&
        Args:
            x (array_like) : Molar composition [-]
            prev_solv_idx (int) : Component index of the old (current) solvent
            new_solv_idx (int) : Component index of the new solvent

        Returns:
            2d array : The transformation matrix $\Psi^{n_k \leftmapsto n_l}$, where $k$ is `new_solv_idx` and $l$ is `prev_solv_idx`.
        """
        psi = np.identity(self.ncomps)
        k = prev_solv_idx # Arbitrary component (tested)
        for i in range(self.ncomps):
            for j in range(self.ncomps):
                psi[i, j] -= x[i] * (k_delta(j, new_solv_idx) - k_delta(j, prev_solv_idx) * k_delta(k, new_solv_idx)) / \
                             x[new_solv_idx]
        return psi

    def get_zarate_X_matr(self, x, dependent_idx):
        """FOR-Transform
        Compute the matrix $X$ as defined by Zárate. See (Definition of frame-invariant thermodiffusion and Soret coefficients for ternary mixtures)
        and memo on diffusion coefficient definitions.
        &&
        Args:
            x (array_like) : Molar composition [-]
            dependent_idx (int) : Index of the dependent species
        Returns:
            2d array : The transformation matrix $X$
        """
        while dependent_idx < 0:
            dependent_idx += self.ncomps

        X = np.zeros((self.ncomps - 1, self.ncomps - 1))
        i_idx = 0
        for i in range(self.ncomps):
            if i == dependent_idx:
                continue
            X[i_idx][i_idx] = x[i]
            j_idx = 0
            for j in range(self.ncomps):
                if j == dependent_idx:
                    continue
                X[i_idx][j_idx] -= x[i] * x[j]
                j_idx += 1
            i_idx += 1
        return X

    def get_zarate_W_matr(self, x, dependent_idx):
        """FOR-Transform
        Compute the matrix $W$ as defined by Zárate. See (Definition of frame-invariant thermodiffusion and Soret coefficients for ternary mixtures)
        and memo on diffusion coefficient definitions.
        &&
        Args:
            x (array_like) : Molar composition [-]
            dependent_idx (int) : Index of the dependent species
        Returns:
            2d array : The transformation matrix $W$
        """
        wt_fracs = x * self.mole_weights / sum(x * self.mole_weights)
        return self.get_zarate_X_matr(wt_fracs, dependent_idx)

    ######################################################
    #         Interface to top-level C++ methods         #
    ######################################################

    def get_conductivity_matrix(self, particle_density, T, mole_fracs, N=None):
        """cpp-interface
        Compute the elements of the matrix corresponding to the set of equations that must be solved for the
        thermal response function sonine polynomial expansion coefficients:
        Eq. (6) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
        &&
        Args:
            particle_density (float) : Particle density (not molar!) [1 / m3]
            T (float) : Temperature [K]
            mole_fracs (list[float]) : Molar composition [-]
            N (Optional, int) : Enskog approximation order.

        Returns:
            2Darray : ($N N_c$ x $N N_c$) matrix, where $N$ is the Enskog approximation order and $N_c$ is
                        the number of components.
        """
        if N is None:
            N = self.default_N
        self.check_valid_composition(mole_fracs)
        return self.cpp_kingas.get_conductivity_matrix(particle_density, T, mole_fracs, N)

    def get_conductivity_vector(self, particle_density, T, mole_fracs, N):
        """cpp-interface
        Compute the right-hand side vector to the set of equations that must be solved for the
        thermal response function Sonine polynomial expansion coefficients:
        Eq. (6) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
        &&
        Args:
            particle_density (float) : Particle density (not molar!) [1 / m3]
            T (float) : Temperature [K]
            mole_fracs (list<float>) : Molar composition [-]
            N (Optional, int) : Enskog approximation order.

        Returns:
            (1Darray) : ($N N_c$,) vector, where $N$ is the Enskog approximation order and $N_c$ is
                            the number of components.
        """
        if N is None:
            N = self.default_N
        self.check_valid_composition(mole_fracs)
        return self.cpp_kingas.get_conductivity_vector(particle_density, T, mole_fracs, N)

    def get_diffusion_vector(self, particle_density, T, mole_fracs, N=None):
        """cpp-interface
        Compute the right-hand side vector to the set of equations that must be solved for the
        diffusive response function Sonine polynomial expansion coefficients.
        Eq. (10) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
        &&
        Args:
            particle_density (float) : Particle density (not molar!) [1 / m3]
            T (float) : Temperature [K]
            mole_fracs (list[float]) : Molar composition [-]
            N (Optional, int) : Enskog approximation order.

        Returns:
            (1Darray) : ($N N_c^2$,) vector, where $N$ is the Enskog approximation order and $N_c$ is
                            the number of components.
        """
        if N is None:
            N = self.default_N
        self.check_valid_composition(mole_fracs)
        return np.array(self.cpp_kingas.get_diffusion_vector(particle_density, T, mole_fracs, N))

    def get_collision_diameters(self, particle_density, T, x):
        """cpp-interface
        Compute collision diameters given by Eq. (40) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
        *Note* Returns zeros for models initialised with is_idealgas=True.
        &&
        Args:
            particle_density (float) : Particle density (not molar!) [1 / m3]
            T (float) : Temperature [K]
            x (list[float]) : Molar composition [-]

        Returns:
            2d array : Collision diameters [m], indexed by component pair.
        """
        key = tuple((particle_density, T))
        if key in self.computed_cd.keys():
            return self.computed_cd[key]

        cd = self.cpp_kingas.get_collision_diameters(particle_density, T, x)
        self.computed_cd[key] = cd
        return cd

    def get_rdf(self, particle_density, T, x):
        """cpp-interface
        Compute the radial distribution function at contact
        &&
        Args:
            particle_density (float) : Particle density (not molar!) [1 / m3]
            T (float) : Temperature [K]
            x (list[float]) : Molar composition [-]

        Returns:
            2d array : RDF at contact, indexed by component pair.
        """
        key = tuple((particle_density, T, tuple(x)))
        if key in self.computed_rdf.keys():
            return self.computed_rdf[key]

        rdf = self.cpp_kingas.get_rdf(particle_density, T, x)
        self.computed_rdf[key] = rdf
        return rdf

    #####################################################
    #          Evaluation of Matrix equations           #
    #####################################################

    def compute_diffusion_coeff_vector(self, particle_density, T, mole_fracs, N=None):
        r"""Utility
        Compute the diffusive response function Sonine polynomial expansion coefficients by solving the set of equations

        $$D d = \delta$$

        Corresponding to Eq. (10) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
        Where $D$ is the matrix returned by the c++ method `get_diffusion_matrix`, and $\delta$ is the vector
        returned by the c++ method `get_diffusion_vector`.
        &&
        Args:
            particle_density (float) : Particle density (not molar!) [1 / m3]
            T (float) : Temperature [K]
            mole_fracs (list[float]) : Molar composition [-]
            N (Optional, int) : Enskog approximation order.

        Returns:
            1Darray : ($N N_c^2$,) vector, where $N$ is the Enskog approximation order and $N_c$ is
                            the number of components, containing the diffusive response function sonine polynomial
                            expansion coefficients ($d_{i, j}^{(q)}$). See `reshape_diffusive_coeff_vector` for help on practical usage.
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
        r"""Utility
        The vector returned by `compute_diffusion_coeff_vector` contains the diffusive response function sonine
        polynomial expansion coefficients (eg. $d_{i, j}^{(q)}$ ).
        To more easily access the correct coefficients, this method reshapes the vector to a matrix indiced
        as `d[i][q][j]` where i and j are component indices, and q refferes to the approximation order summation index.
        &&
        Args:
            d (1D array) : The array returned by `get_diffusion_coeff_vector`

        Returns:
            3D array : The matrix of $d_{i, j}^{(q)}$ coefficients ordered as `d[i][q][j]`
        """

        N = len(d) // (self.ncomps * self.ncomps)
        matr = np.zeros((self.ncomps, N, self.ncomps))
        for i in range(self.ncomps):
            for q in range(N):
                for k in range(self.ncomps):
                    matr[i, q, k] = d[N * self.ncomps * k + self.ncomps * q + i]
        return matr
    
    def compute_dth_vector(self, particle_density, T, mole_fracs, N=None):
        r"""Utility
        Compute the coefficients $d_i^{(J = 0)}$, by solving the set of equations

        $$\sum_j d_{i,j}^{(0)} d_j^{\vec{J} = 0} = \ell_i^{(0)}$$,

        i.e. Eq. (15) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
        &&
        Args:
            particle_density (float) : Particle density (not molar!) [1 / m3]
            T (float) : Temperature [K]
            mole_fracs (list<float>) : Molar composition [-]
            N (Optional, int) : Enskog approximation order.

        Returns:
            1Darray : Vector of the $d_i^{(J = 0)}$ coefficients.
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
        r"""Utility
        Compute the thermal response function Sonine polynomial expansion coefficients by solving the set of equations

        $$\Lambda \ell = \lambda$$

        Corresponding to Eq. (6) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
        Where $\Lambda$ is the matrix returned by the c++ method `get_conductivity_matrix`, and $\lambda$ is the vector
        returned by the c++ method `get_conductivity_vector`.
        &&
        Args:
            particle_density (float) : Particle density (not molar!) [1 / m3]
            T (float) : Temperature [K]
            mole_fracs (list[float]) : Molar composition [-]
            N (Optional, int) : Enskog approximation order.

        Returns:
            1Darray : ($N N_c$,) vector, where $N$ is the Enskog approximation order and $N_c$ is
                            the number of components. The vector is ordered as
                            l[:N_c] = [$\ell_1^{(0)}$, $\ell_2^{(0)}$, ..., $\ell_{N_c}^{(0)}$]
                            l[N_c : 2 * N_c] = [$\ell_1^{(1)}$, $\ell_2^{(1)}$, ..., $\ell_{N_c}^{(1)}$]
                            ... etc ...
                            Where subscripts indicate component indices, and superscripts are Enskog approximation
                            order summation indices.
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
        r"""Utility
        Compute the viscous response function Sonine polynomial expansion coefficients by solving the set of equations

        $$\Beta b = \beta$$

        Corresponding to Eq. (8) in RET for Mie fluids (https://doi.org/10.1063/5.0149865)
        Where $\Beta$ is the matrix returned by the c++ method `get_viscosity_matrix`, and $\beta$ is the vector
        returned by the c++ method `get_viscosity_vector`.
        &&
        Args:
            particle_density (float) : Particle density (not molar!) [1 / m3]
            T (float) : Temperature [K]
            mole_fracs (list<float>) : Molar composition [-]
            N (Optional, int) : Enskog approximation order.

        Returns:
            1D array : ($N N_c$,) vector, where $N$ is the Enskog approximation order and $N_c$ is
                            the number of components, ordered as
                            b[:N_c] = [b_1^{(0)}, b_2^{(0)}, ..., b_{N_c}^{(0)}]
                            b[N_c: 2 * N_c] = [b_1^{(1)}, b_2^{(1)}, ..., b_{N_c}^{(1)}]
                            ... etc ...
                            where subscripts denote component indices and superscripts denote Enskog approximation
                            summation indices.
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