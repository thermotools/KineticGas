'''
Author : Vegard Gjeldvik Jervell
Contains : Parent class for all 'Mie-Type' potentials, including MieKinGas, QuantumMie and LJSpline
            Implements mixing rules and initializer common for all
Usage : Accepts a list of parameter dicts as the second argument to the initializer, as well as parameters
        See sub-classes for examples.
'''
import numpy as np
from scipy.constants import Boltzmann, Avogadro
from scipy.integrate import quad
from pykingas import KineticGas
from copy import deepcopy
from warnings import warn

class MieType(KineticGas):

    def __init__(self, comps, potential,
                 mole_weights=None, sigma=None, eps_div_k=None,
                 la=None, lr=None, lij=0, kij=0,
                 N=4, is_idealgas=False,
                 parameter_ref='default'):
        '''
        :param comps (str): Comma-separated list of components

        If parameters are explicitly supplied, these will be used instead of those in the database
        :param mole_weights : (1D array) Molar weights [g/mol]
        :param sigma : (1D array) hard-sphere diameters [m]
        :param eps_div_k : (1D array) epsilon parameter / Boltzmann constant [-]
        :param la, lr : (1D array) attractive and repulsive exponent of the pure components [-]
        :param lij : (float) Mixing parameter for sigma (lij > 0 => smaller sigma_12, lij < 0 => larger sigma_12)
        :param kij : (float) Mixing parameter for epsilon (kij > 0 => favours mixing, kij < 0 => favours separation)
        :param use_db (bool) : Use precomputed database values for omega_integrals if available
        '''
        super().__init__(comps, mole_weights=mole_weights, N=N, is_idealgas=is_idealgas)

        self.lij = lij
        self.kij = kij

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

    def update_omega_db(self):
        self.__omega_db.update(self.cpp_kingas.omega_map)
        self.__omega_db.dump()

    def get_avg_R(self, T, x):
        v_bar_1 = np.sqrt(3 * Boltzmann * T / self.mole_weights[0]) * np.sqrt(
            self.m0 * self.M[0] * self.M[1] / (2 * Boltzmann * T))
        v_bar_2 = np.sqrt(3 * Boltzmann * T / self.mole_weights[1]) * np.sqrt(
            self.m0 * self.M[0] * self.M[1] / (2 * Boltzmann * T))

        v_min = min(v_bar_1 - v_bar_2, v_bar_2 - v_bar_1)
        v_max = v_bar_1 + v_bar_2
        v_bar_12 = (v_bar_1 + v_bar_2) / 2
        R_11 = quad(lambda b: self.cpp_kingas.get_R(0, 0, T, v_bar_1, b), 0, self.sigma_ij[0, 0])[0] / self.sigma_ij[0, 0]
        R_12 = quad(lambda b: self.cpp_kingas.get_R(0, 1, T, v_bar_12, b), 0, self.sigma_ij[0, 1])[0] / self.sigma_ij[
            0, 1]
        R_22 = quad(lambda b: self.cpp_kingas.get_R(1, 1, T, v_bar_2, b), 0, self.sigma_ij[1, 1])[0] / self.sigma_ij[1, 1]

        return x[0] ** 2 * R_11 + 2 * x[0] * x[1] * R_12 + x[1] ** 2 * R_22

    def get_density_correction(self, Vm, T, x):
        particle_density = Avogadro / Vm
        avg_R = self.get_avg_R(T, x)
        eta = (np.pi / 6) * particle_density * avg_R ** 3  # (x[0] * BH_sigma[0, 0]**3 + x[1] * BH_sigma[1, 1]**3)
        g = (1 - 0.5 * eta) / (1 - eta) ** 3
        return g * particle_density

    def get_epsilon_matrix(self, eps_div_k, kij):
        # Apply mixing rules
        epsilon = np.array(eps_div_k) * Boltzmann
        return (np.ones((self.ncomps, self.ncomps)) - self.kij * (np.ones((self.ncomps, self.ncomps)) - np.identity(self.ncomps))) * np.sqrt(
            epsilon * np.vstack(epsilon))  # Only apply mixing parameter kij to the off-diagonals

    def get_sigma_matrix(self, sigma):
        '''
        Apply mixing rules

        :param sigma: (1D array) sigma-parameters [m]
        :return: N x N matrix of sigma parameters, where sigma_ij = 0.5 * (sigma_i + sigma_j),
                such that the diagonal is the diameter of each component, and off-diagonals are the cross-collision distances.
        '''

        sigma_ij = (np.ones((self.ncomps, self.ncomps)) - self.lij * (np.ones((self.ncomps, self.ncomps)) - np.identity(self.ncomps)))\
                * 0.5 * np.sum(np.meshgrid(sigma, np.vstack(sigma)), axis=0)  # Only apply mixing parameter lij to the off-diagonals

        return sigma_ij

    def get_lambda_matrix(self, lambdas, lij):
        l = np.array(lambdas)
        return 3 + (1 - lij) * np.sqrt((l - 3) * np.vstack(l - 3))