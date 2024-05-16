#include "ModTangToennis.h"
#include <algorithm>

double partialfactorial(int start, int stop){
    double fac = 1.;
    start = (start == 0) ? 1 : start;
    for (int i = start; i <= stop; i++){
        fac *= i;
    }
    return fac;
}

ModTangToennis::ModTangToennis(TangToennisParam param, vector1d mole_weights, vector2d sigma, bool is_idealgas)
    : Spherical(mole_weights, sigma, is_idealgas, true), param{param}
{
    if (!is_idealgas) throw std::runtime_error("Modified Tang-Toennis only implemented for ideal gas!");
}

double ModTangToennis::potential(int i, int j, double r){
    r *= 1e9; // Using nm internally
    if (r < 0.4 * param.Re){
        return (param.A_tilde / r) * exp(- param.a_tilde * r) * BOLTZMANN;
    }
    double u = param.A * exp(param.a1 * r + param.a2 * pow(r, 2) + param.am1 * pow(r, -1) + param.am2 * pow(r, -2));
    double exp_prefactor = exp(- param.b * r);
    for (int n = 3; n <= 8; n++){
        double tmp = 0.;
        int k = 0;
        for (; k <= std::min(2 * n, 10); k++){
            tmp += pow(param.b * r, k) / partialfactorial(1, k);
        }
        for (; k <= 2 * n; k++){ // Prevent factorial overflow
            tmp += (pow(param.b * r, 10) / partialfactorial(1, 10)) * (pow(param.b * r, k - 10) / partialfactorial(11, k));
        }
        u -= (param.C[n - 3] / pow(r, 2 * n)) * (1 - exp_prefactor * tmp);
    }
    return u * BOLTZMANN;
}

double ModTangToennis::potential_derivative_r(int i, int j, double r){
    if (r * 1e9 < 0.4 * param.Re){
        r *= 1e9;
        return - BOLTZMANN * param.A_tilde * exp(- param.a_tilde * r) * (param.a_tilde * r + 1) / pow(r, 2);
    }
    double eps = 1e-6;
    double u1 = potential(i, j, r * (1. + eps));
    double um1 = potential(i, j, r * (1. - eps));
    return (u1 - um1) / (2. * eps * r);
}

double ModTangToennis::potential_dblderivative_rr(int i, int j, double r){
    if (r * 1e9 < 0.4 * param.Re){
        r *= 1e9;
        return BOLTZMANN * param.A_tilde * exp(- param.a_tilde * r) * (pow(param.a_tilde * r, 2) + 2. * param.a_tilde * r + 2) / pow(r, 3);
    }
    double eps = 1e-6;
    double u1 = potential(i, j, r * (1. + eps));
    double u0 = potential(i, j, r);
    double um1 = potential(i, j, r * (1. - eps));
    return (u1 - 2 * u0 + um1) / pow(2. * eps * r, 2);
}