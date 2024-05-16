/*
    Author : Vegard Gjeldvik Jervell

    Contains :

    Major base-functions for KineticGas. These functions are required for all implementations, regardless of what
    potential model they use. Contains summational expressions for the bracket integrals, and interfaces to
    get the matrices and vectors required to evaluate diffusion coefficients.
*/

#include "KineticGas.h"
#include <vector>
#include <algorithm>
#include <thread>
#include <functional>
#include <math.h>
#include <iostream>
#include "pybind11/pybind11.h"
#include <chrono>

#ifdef DEBUG
#define _LIBCPP_DEBUG 1
#endif

// --------------------------------------------------------------------------------------------------- //
// -------------------------------Constructor and helper functions------------------------------------ //

KineticGas::KineticGas(const std::vector<double> mole_weights, bool is_idealgas, bool is_singlecomp)
  : Ncomps{static_cast<unsigned long>(mole_weights.size())},
    is_idealgas{is_idealgas}, is_singlecomp{is_singlecomp}, m{mole_weights}
    {

    m0 = std::vector<std::vector<double>>(Ncomps, std::vector<double>(Ncomps));
    M = std::vector<std::vector<double>>(Ncomps, std::vector<double>(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j < Ncomps; j++){
            M[i][j] = (m[i] / (m[i] + m[j]));
            m0[i][j] = (m[i] + m[j]);
        }
    }
}

std::vector<double> KineticGas::get_wt_fracs(const std::vector<double> mole_fracs){
    std::vector<double> wt_fracs(Ncomps, 0.);
    double tmp{0.};
    for (int i = 0; i < Ncomps; i++){
        tmp += mole_fracs[i] * m[i];
    }
    for (int i = 0; i < Ncomps; i++){
        wt_fracs[i] = mole_fracs[i] * m[i] / tmp;
    }
    return wt_fracs;
}

// --------------------------------------------------------------------------------------------------- //
//             K-factors, neccesary for computations above infinite dilution                           //
// --------------------------------------------------------------------------------------------------- //

std::vector<double> KineticGas::get_K_factors(double rho, double T, const std::vector<double>& mole_fracs){

    if (is_idealgas) return std::vector<double>(Ncomps, 1.);

    std::vector<std::vector<double>> rdf = get_rdf(rho, T, mole_fracs);
    std::vector<double> K(Ncomps, 0.);
    std::vector<std::vector<double>> cd = get_collision_diameters(rho, T, mole_fracs);
    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j < Ncomps; j++){
            K[i] += mole_fracs[j] * pow(cd[i][j], 3) * M[i][j] * M[j][i] * rdf[i][j];
        }
        K[i] *= (24. * PI * rho / 15);
        K[i] += 1;
    }
    return K;
}

std::vector<double> KineticGas::get_K_prime_factors(double rho, double T, const std::vector<double>& mole_fracs){

    if (is_idealgas) return std::vector<double>(Ncomps, 1.);

    std::vector<std::vector<double>> rdf = get_rdf(rho, T, mole_fracs);
    std::vector<double> K_prime(Ncomps, 0.);
    std::vector<std::vector<double>> cd = get_collision_diameters(rho, T, mole_fracs);
    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j < Ncomps; j++){
            K_prime[i] += mole_fracs[j] * M[j][i] * pow(cd[i][j], 3.) * rdf[i][j];
        }
        K_prime[i] *= 8. * PI * rho / 15;
        K_prime[i] += 1.;
    }
    return K_prime;
}

// --------------------------------------------------------------------------------------------------- //
// ---------------------------------Matrix and vector generators-------------------------------------- //
    /*
        These are the methods used to generate the matrices and vectors of the linear sets of equations that determine
        the thermal-, diffusive- and viscous response functions' Sonine polynomial expansion coefficients. That is,
        to get the expansion coefficients of the thermal response function (denoted "\ell"), the set of equations
            \Lambda \ell = \lambda
        must be solved. The method 'get_conductivity_matrix' returns the left-hand side matrix \Lambda, while the
        method 'get_conductivity_vector' returns the right-hand side vector \lambda.
        The same is true for the 'get_<*>_matrix' and 'get_<*>_vector', with <*> = {diffusion, viscosity}.

        As of the current implementation, the matrix equations are solved on the python-side of the package.
    */

std::vector<std::vector<double>> KineticGas::get_conductivity_matrix(double rho, double T, const std::vector<double>& x, int N){
    std::vector<std::vector<double>> rdf = get_rdf(rho, T, x);
    std::vector<std::vector<double>> matr(N * Ncomps, std::vector<double>(N * Ncomps, 0.));
    std::vector<double> wt_fracs = get_wt_fracs(x);

    precompute_conductivity_omega(N, T); // Compute all the collision integrals required for this conductivity

    // Build omega_k row of the matrix
    for (int i = 0; i < Ncomps; i++){
            matr[0][i] = wt_fracs[i];
    }

    int row, col;
    // Build Lambda_1^(p>0) block
    for (int p = 1; p < N; p++){
        for (int q = 0; q < N; q++){
            for (int j = 0; j < Ncomps; j++){
                row = p;
                col = Ncomps * q + j;
                matr[row][col] = x[0] * x[j] * rdf[0][j] * H_ij(p, q, 0, j, T);
                if (j == 0){
                    for (int l = 0; l < Ncomps; l++){
                        matr[row][col] += x[0] * x[l] * rdf[0][l] * H_i(p, q, 0, l, T);
                    }
                }
                matr[row][col] *= 8. * sqrt(m[0] * m[j]) / (75. * pow(BOLTZMANN, 2) * T);
            }
        }
    }
    // Build Lambda_{(i > 1)} block
    for (int p = 0; p < N; p++){
        for (int q = 0; q < N; q++){
            for (int i = 1; i < Ncomps; i++){
                for (int j = 0; j < Ncomps; j++){
                    row = N + (Ncomps - 1) * p + (i - 1);
                    col = Ncomps * q + j;
                    matr[row][col] = x[i] * x[j] * rdf[i][j] * H_ij(p, q, i, j, T);
                    if (i == j){
                        for (int l = 0; l < Ncomps; l++){
                            matr[row][col] += x[i] * x[l] * rdf[i][l] * H_i(p, q, i, l, T);
                        }
                    }
                    matr[row][col] *= 8. * sqrt(m[i] * m[j]) / (75. * pow(BOLTZMANN, 2) * T);
                }
            }
        }
    }
    return matr;
}

std::vector<double> KineticGas::get_conductivity_vector(double rho, double T, const std::vector<double>& x, int N){
    std::vector<double> K = get_K_factors(rho, T, x);
    std::vector<double> l(Ncomps * N, 0.);
    l[1] = x[0] * K[0] * 4. / (5. * BOLTZMANN); // Lambda_1^{(p > 0} block
    for (int i = 1; i < Ncomps; i++){
        l[N + (Ncomps - 1) + (i - 1)] = x[i] * K[i] * 4. / (5. * BOLTZMANN); // Lambda_{i > 1} block
    }
    return l;
}


std::vector<std::vector<double>> KineticGas::get_diffusion_matrix(double rho, double T, const std::vector<double>& x, int N){
    std::vector<std::vector<double>> rdf = get_rdf(rho, T, x);
    std::vector<std::vector<double>> matr(N*pow(Ncomps, 2), std::vector<double>(N*pow(Ncomps, 2), 0.));
    std::vector<double> wt_fracs = get_wt_fracs(x);

    precompute_diffusion_omega(N, T); // Compute all the collision integrals required for this diffusion matrix

    for (int k = 0; k < Ncomps; k++){ // Build omega_k block of the matrix
        for (int i = 0; i < Ncomps; i++){
            matr[k][N * Ncomps * k + i] = wt_fracs[i];
        }
    }
    int row, col;
    // Build Lambda_1^(p>0) block
    for (int k = 0; k < Ncomps; k++){
        for (int p = 1; p < N; p++){
            for (int q = 0; q < N; q++){
                for (int j = 0; j < Ncomps; j++){
                    row = Ncomps + (N - 1) * k + (p - 1);
                    col = N * Ncomps * k + q * Ncomps + j;
                    matr[row][col] = x[0] * x[j] * rdf[0][j] * H_ij(p, q, 0, j, T);
                    if (j == 0){
                        for (int l = 0; l < Ncomps; l++){
                            matr[row][col] += x[0] * x[l] * rdf[0][l] * H_i(p, q, 0, l, T);
                        }
                    }
                    matr[row][col] *= 8. * sqrt(m[0] * m[j]) / (75. * pow(BOLTZMANN, 2) * T);
                }
            }
        }
    }
    // Build Lambda_{(i > 1)} block
    for (int k = 0; k < Ncomps; k++){
        for (int p = 0; p < N; p++){
            for (int q = 0; q < N; q++){
                for (int i = 1; i < Ncomps; i++){
                    for (int j = 0; j < Ncomps; j++){
                        row = N * Ncomps + N * (Ncomps - 1) * k + (Ncomps - 1) * p + (i - 1);
                        col = Ncomps * N * k + Ncomps * q + j;
                        matr[row][col] = x[i] * x[j] * rdf[i][j] * H_ij(p, q, i, j, T);
                        if (i == j){
                            for (int l = 0; l < Ncomps; l++){
                                matr[row][col] += x[i] * x[l] * rdf[i][l] * H_i(p, q, i, l, T);
                            }
                        }
                        matr[row][col] *= 8. * sqrt(m[i] * m[j]) / (75. * pow(BOLTZMANN, 2) * T);
                    }
                }
            }
        }
    }
    return matr;
}

std::vector<double> KineticGas::get_diffusion_vector(double rho, double T, const std::vector<double>& x, int N){
    std::vector<double> delta_vector(N * Ncomps * Ncomps, 0.);
    std::vector<double> wt_fracs = get_wt_fracs(x);

    for (int k = 0; k < Ncomps; k++){
        for (int i = 1; i < Ncomps; i++){
            delta_vector[Ncomps * N + N * (Ncomps - 1) * k + (i - 1)] = - 8. * wt_fracs[i] / (25. * BOLTZMANN);
            if (i == k){
                delta_vector[Ncomps * N + N * (Ncomps - 1) * k + (i - 1)] += 8. / (25. * BOLTZMANN);
            }
        }
        
    }
    return delta_vector;
}

std::vector<std::vector<double>> KineticGas::get_viscosity_matrix(double rho, double T, const std::vector<double>&x, int N){
    std::vector<std::vector<double>> rdf = get_rdf(rho, T, x);
    std::vector<std::vector<double>> viscosity_mat(Ncomps * N, std::vector<double>(Ncomps * N, 0.));

    precompute_viscosity_omega(N, T); // Compute all the collision integrals required for this viscosity

    for (int p = 0; p < N; p++){
        for (int i = 0; i < Ncomps; i++){
            for (int q = 0; q < N; q++){
                for (int j = 0; j < Ncomps; j++){
                    viscosity_mat[p * Ncomps + i][q * Ncomps + j] = x[i] * x[j] * rdf[i][j] * L_ij(p, q, i, j, T);
                    if (i == j){
                        for (int k = 0; k < Ncomps; k++){
                            viscosity_mat[p * Ncomps + i][q * Ncomps + j] += x[i] * x[k] * rdf[i][k] * L_i(p, q, i, k, T);
                        }
                    }
                    viscosity_mat[p * Ncomps + i][q * Ncomps + j] *= 2. / (5. * BOLTZMANN * T);
                }
            }
        }
    }
    return viscosity_mat;
}
std::vector<double> KineticGas::get_viscosity_vector(double rho, double T, const std::vector<double>& x, int N){
    std::vector<double> K_prime = get_K_prime_factors(rho, T, x);
    std::vector<double> viscosity_vec(N * Ncomps, 0.);
    // Elements with p != 0 are 0, therefore only iterate over i.
    for (int i = 0; i < Ncomps; i++){
        viscosity_vec[i] = 2. * x[i] * K_prime[i] / (BOLTZMANN * T);
    }
    return viscosity_vec;
}


// --------------------------------------------------------------------------------------------------- //
//                                   Linear combination weights                                        //
    /*
        The square bracket integrals required to generate the matrices corresponding to the sets of equations that
        must be solved to determine the response functions Sonine polynomial expansion coefficients can be expressed
        as linear combinations of the collision integrals (omega integrals).

        The methods in this section generate the weights of these linear combinations using the explicit summational
        expressions derived by Thompson, Tipton and Lloyalka in the series of papers
        Chapman–Enskog solutions to arbitrary order in Sonine polynomials (I - IV),
        European Journal of Mechanics - B/Fluids.
        doi :
        (I)   : https://doi.org/10.1016/j.physa.2006.12.001
        (II)  : https://doi.org/10.1016/j.euromechflu.2008.09.002
        (III) : https://doi.org/10.1016/j.euromechflu.2008.12.002
        (IV)  : https://doi.org/10.1016/j.euromechflu.2009.05.002
    */


// --------------------------------------------------------------------------------------------------- //
//                 A_pqrl factors, used for conductivity and diffusion                                 //
// --------------------------------------------------------------------------------------------------- //

double KineticGas::A(int p, int q, int r, int l) const {
    double value{0.0};
    int max_i = std::min(std::min(p, q), std::min(r, p + q + 1 - r));
    for (int i = l - 1; i <= max_i; i++){
        value += ((ipow(8, i) * Fac(p + q - 2 * i) * ipow(-1, l + r + i) * Fac(r + 1) * Fac(2 * (p + q + 2 - i)) * ipow(4, r)) /
                (Fac(p - i) * Fac(q - i) * Fac(l) * Fac(i + 1 - l) * Fac(r - i) * Fac(p + q + 1 - i - r) * Fac(2 * r + 2)
                    * Fac(p + q + 2 - i) * ipow(4, p + q + 1))) * ((i + 1 - l) * (p + q + 1 - i - r) - l * (r - i));
    }
    return value;
}

double KineticGas::A_prime(int p, int q, int r, int l, double tmp_M1, double tmp_M2) const {
    double F = (pow(tmp_M1, 2) + pow(tmp_M2, 2)) / (2 * tmp_M1 * tmp_M2);
    double G = (tmp_M1 - tmp_M2) / tmp_M2;

    int max_i = std::min(p, std::min(q, std::min(r, p + q + 1 - r)));
    int max_k;
    int max_w;

    Frac p1{1}, p2{1};

    double value{0.0};
    for (int i = l - 1; i <= max_i; i++ ){
        max_w = std::min(p, std::min(q, p + q + 1 - r)) - i;
        max_k = std::min(l, i);
        for (int k = l - 1; k <= max_k; k++){
            for (int w = 0; w <= max_w; w++){
                p1 = ((ipow(8, i) * Fac(p + q - 2 * i - w) * ipow(-1, r + i) * Fac(r + 1) * Fac(2 * (p + q + 2 - i - w)) * ipow(2, 2 * r) * pow(F, i - k) * pow(G, w)
                        * ((ipow(2, 2 * w - 1) * pow(tmp_M1, i) * pow(tmp_M2, p + q - i - w)) * 2)
                        * (tmp_M1 * (p + q + 1 - i - r - w) * delta(k, l) - tmp_M2 * (r - i) * delta(k, l - 1))
                        ));
                p2 = (Fac(p - i - w) * Fac(q - i - w) * Fac(r - i) * Fac(p + q + 1 - i - r - w) * Fac(2 * r + 2) * Fac(p + q + 2 - i - w) * ipow(4, p + q + 1) * Fac(k) * Fac(i - k) * Fac(w));
                value += p1 / p2; // NB: Gives Bus error if unless v1 and v2 are initialized before adding to value... pls help
            }
        }
    }

    return value;
}

// --------------------------------------------------------------------------------------------------- //
//                            B_pqrl factors, used for viscosity                                       //
// --------------------------------------------------------------------------------------------------- //

inline Frac poch(int z, int n){ // Pochhammer notation (z)_n, used in the B-factors
    return Fac(z + n - 1) / Fac(z - 1);
}

double KineticGas::B_prime(int p, int q, int r, int l, double M1, double M2) const {
    if (r == p + q + 3) return 0.0;
    double val{0.0};
    double inner{0.0};
    Frac num;
    Frac denom;
    double G = (M1 - M2) / M2;
    double F = (pow(M1, 2) + pow(M2, 2)) / (2.0 * M1 * M2);
    int w_max;
    int i_max = std::min(p, std::min(q, std::min(r, p + q + 2 - r)));
    for (int i = l - 2; i <= i_max; i++){
        inner = 0.0;
        w_max = std::min(p, std::min(q, p + q + 2 - r)) - i;
        for (int w = 0; w <= w_max; w++){
            num = poch(p + 1 - i - w, w) * poch(q + 1 - i - w, w) * poch(p + q + 3 - i - r - w, w)
                    * ipow(2, 2 * w - 2) * pow(G, w) * poch(p + q + 4 - i - w, w) * pow(M1 * M2, i)
                    * pow(M2, p + q - 2 * i - w) * pow(F, i + 2 - l) * pow(M1, 2) * 4;

            denom = Fac(w) * poch(p + q + 1 - 2 * i - w, w) * poch(2 * (p + q + 3 - i) - 2 * w + 1, w)
                    * poch(2 * (p + q + 3 - i) - w + 1, w) * Fac(l) * Fac(i + 2 - l);

            inner += (num / denom) * ((3.0 / 2.0) * pow(M2 / M1, 2) * l * (l - 1) * (r - i) * (r - i - 1)
                    - (2.0 / F) * (M2 / M1) * l * (i + 2 - l) * (r - i) * (p + q + 2 - i - r - w)
                    + pow(1.0 / F, 2) * (i + 1 - l) * (i + 2 - l)
                        * ((p + q + 1 - i - r - w) * (p + q + 2 - i - r - w)
                            - 0.5 * pow(M2 / M1, 2) * (r - i) * (r - i - 1)));
        }

        num = ipow(2, 2 * r) * ipow(8, i) * Fac(p + q - 2 * i) * ipow(-1, r + i)
                * (1 - delta(i, -1)) * Fac(r + 1) * Fac(2 * (p + q + 3 - i));
        denom = ipow(4, p + q + 2) * Fac(p - i) * Fac(q - i) * Fac(r - i)
                * Fac(p + q + 2 - i - r) * Fac(2 * r + 2) * Fac(p + q + 3 - i);

        val += (num / denom).eval() * inner;
    }
    return val;
}

double KineticGas::B_dblprime(int p, int q, int r, int l, double M1, double M2) const {
    if (r == p + q + 3) return 0.0;
    double val{0.0};
    Frac num;
    Frac denom;
    Frac prefactor;
    int i_max = std::min(p, std::min(q, std::min(r, p + q + 2 - r)));
    for (int i = l - 2; i <= i_max; i++){
        num = ipow(2, 2 * r) * ipow(8, i) * Fac(p + q - 2 * i) * ipow(-1, r + i)
                * (1 - delta(i, -1)) * Fac(r + 1) * Fac(2 * (p + q + 3 - i))
                * ipow(-1, l);
        denom = ipow(4, p + q + 2) * Fac(p - i) * Fac(q - i) * Fac(r - i) * Fac(p + q + 2 - i - r)
                * Fac(2 * r + 2) * Fac(p + q + 3 - i) * Fac(l) * Fac(i + 2 - l);

        prefactor = num / denom;
        val += prefactor * (i + 1 - l) * (i + 2 - l) * ((p + q + 1 - i - r) * (p + q + 2 - i - r)
                                                        - 0.5 * (r - i) * (r - i - 1));

        val += prefactor * (3.0 / 2.0) * (l - 1) * l * (r - i) * (r - i - 1);
        val += prefactor * (- 2) * l * (i + 2 - l) * (r - i) * (p + q + 2 - i - r);
    }
    val *= pow(M2, p + 1) * pow(M1, q + 1);
    return val;
}

// --------------------------------------------------------------------------------------------------- //
//                                    Square bracket integrals                                         //
    /*
        The square bracket integrals are computed as linear combinations of the collision integrals (omega integrals)
        as identified by Chapman and Cowling in "The mathematical theory of non-uniform gases (Cambridge University
        press, 1970).
    */


// --------------------------------------------------------------------------------------------------- //
//                      H-integrals, used for conductivity and diffusion                               //
// --------------------------------------------------------------------------------------------------- //

double KineticGas::H_ij(int p, int q, int i, int j, double T){
    double tmp_M1{M[i][j]}, tmp_M2{M[j][i]};
    double value{0.0};
    int max_l = std::min(p, q) + 1;
    int max_r;
    for (int l = 1; l <= max_l; l++){
        max_r = p + q + 2 - l;
        for (int r = l; r <= max_r; r++){
            value += A(p, q, r, l) * omega(i, j, l, r, T);
        }
    }
    value *= 8 * pow(tmp_M2, p + 0.5) * pow(tmp_M1, q + 0.5);
    return value;
}

double KineticGas::H_i(int p, int q, int i, int j, double T){
    double tmp_M1{M[i][j]}, tmp_M2{M[j][i]};
    double value{0.0};

    int max_l = std::min(p, q) + 1;
    int max_r;
    for (int l = 1; l <= max_l; l++){
        max_r = p + q + 2 - l;
        for (int r = l; r <= max_r; r++){
            value += A_prime(p, q, r, l, tmp_M1, tmp_M2) * omega(i, j, l, r, T);
        }
    }
    value *= 8;
    return value;
}

// --------------------------------------------------------------------------------------------------- //
//                               L-integrals, used for viscosity                                       //
// --------------------------------------------------------------------------------------------------- //

double KineticGas::L_ij(int p, int q, int i, int j, double T){
    double val{0.0};
    double M1{M[i][j]};
    double M2{M[j][i]};
    for (int l = 1; l <= std::min(p, q) + 2; l++){
        for (int r = l; r <= p + q + 4 - l; r++){
            val += B_dblprime(p, q, r, l, M1, M2) * omega(i, j, l, r, T);
        }
    }
    val *= 16.0 / 3.0;
    return val;
}
double KineticGas::L_i(int p, int q, int i, int j, double T){
    double val{0.0}, M1{M[i][j]}, M2{M[j][i]};
    for (int l = 1; l <= std::min(p, q) + 2; l++){
        for (int r = l; r <= p + q + 4 - l; r++){
            val += B_prime(p, q, r, l, M1, M2) * omega(i, j, l, r, T);
        }
    }
    val *= 16.0 / 3.0;
    return val;
}
