#include "QuantumMie.h"
#include <iostream>
#include <mutex>
#include <thread>

QuantumMie::QuantumMie(vector1d mole_weights, vector2d sigma, vector2d eps, vector2d la, vector2d lr, std::vector<int> FH_order, bool is_idealgas)
        : Sutherland(mole_weights, sigma, eps, 6, is_idealgas), FH_order{FH_order}, Q_factors(4, vector2d(Ncomps, vector1d(Ncomps, 0.)))
    {
        for (size_t i = 0; i < Ncomps; i++){
            for (size_t j = i; j < Ncomps; j++){
                C[0][i][j] = C[0][j][i] = (lr[i][j] / (lr[i][j] - la[i][j]))
                                                * pow(lr[i][j] / la[i][j], (la[i][j] / (lr[i][j] - la[i][j])));
                C[1][i][j] = C[1][j][i] = - C[0][i][j];

                lambda[0][i][j] = lambda[0][j][i] = lr[i][j];
                lambda[1][i][j] = lambda[1][j][i] = la[i][j];
                lambda[2][i][j] = lambda[2][j][i] = lr[i][j] + 2.;
                lambda[3][i][j] = lambda[3][j][i] = la[i][j] + 2.;
                lambda[4][i][j] = lambda[4][j][i] = lr[i][j] + 4.;
                lambda[5][i][j] = lambda[5][j][i] = la[i][j] + 4.;

                double mu1_inv{0.}, mu2_inv{0.};
                if (FH_order[i] > 0) mu1_inv += 1. / m[i];
                if (FH_order[j] > 0) mu1_inv += 1. / m[j];
                if (FH_order[i] > 1) mu2_inv += 1. / m[i];
                if (FH_order[j] > 1) mu2_inv += 1. / m[j];
                double D1_factor = pow(HBAR, 2) * mu1_inv / (24.0 * BOLTZMANN);
                double D2_factor = pow(HBAR, 2) * mu2_inv / (24.0 * BOLTZMANN);

                Q_factors[0][i][j] = Q_factors[0][j][i] = C[0][i][j] * D1_factor * Q1(i, j, lr) * pow(sigma[i][j], -2);
                Q_factors[1][i][j] = Q_factors[1][j][i] = C[1][i][j] * D1_factor * Q1(i, j, la) * pow(sigma[i][j], -2);
                Q_factors[2][i][j] = Q_factors[2][j][i] = C[0][i][j] * pow(D2_factor, 2) * Q2(i, j, lr) * pow(sigma[i][j], -4);
                Q_factors[3][i][j] = Q_factors[3][j][i] = C[1][i][j] * pow(D2_factor, 2) * Q2(i, j, la) * pow(sigma[i][j], -4);
            }
        }
        if (*(std::max_element(FH_order.begin(), FH_order.end())) == 0){
            nterms = 2;
        }
        else if (*(std::max_element(FH_order.begin(), FH_order.end())) == 1){
            nterms = 4;
        }
    }

void QuantumMie::set_temperature(double T){
    if (T == current_temperature) return;
    current_temperature = T;
    for (size_t k = 2; k < nterms; k++){
        for (size_t i = 0; i < Ncomps; i++){
            for (size_t j = i; j < Ncomps; j++){
                C[k][i][j] = C[k][j][i] = Q_factors[k - 2][i][j] / pow(T, k / 2); // Integer division, such that k = {2, 3} => pow(T, 1), k = {4, 5} => pow(T, 2) etc.
            }
        }
    }
    compute_sigma_eps_eff(T);
    compute_vdw_alpha();
}

void QuantumMie::compute_sigma_eps_eff(double T){
    /*
        Search the vector computed_T to find the highest temperature that is lower than T,
        set sigma_min and sigma_eff to the values previously computed at that temperature,
        and set insert_idx to the index one greater than the index of the cached values.

        If the cached temperature is closer than 1 K to the input temperature, set do_insert=false.

        After this block has finished, sigma_min and sigma_eff will have been set to good initial
        guesses for the solvers, and insert_idx, will point to the index at
        which the computed values of sigma_eff and sigma_min should be inserted in the cache to ensure
        that they are sorted by increasing temperature.
    */
    size_t insert_idx = 0;
    bool do_insert = true;
    size_t closest_idx = 0;
    if (computed_T.size() == 0){
        insert_idx = 0;
        sigma_min = sigma;
        sigma_eff = sigma;
    }
    else{
        closest_idx = 0;
        double closest_dist = 1e16;
        if (abs(computed_T.back() - T) < abs(computed_T.front() - T)){ // Search starting from last element
            if (abs(computed_T.back() - T) < 1.){
                closest_idx = computed_T.size() - 1;
                closest_dist = 0.;
            }
            else if (T > computed_T.back()){
                closest_idx = computed_T.size() - 1;
                closest_dist = T - computed_T.back();
                insert_idx = computed_T.size();
            }
            else {
                closest_idx = computed_T.size() - 1;
                for (size_t i = computed_T.size() - 1; i > -1; i--){
                    closest_idx = i;
                    insert_idx = i + 1;
                    if (computed_T[i] < T) {
                        closest_dist = T - computed_T[i];
                        break;
                    }
                }
            }
        }
        else { // Search starting from first element
            if (abs(computed_T[0] - T) < 1.){
                closest_idx = 0;
                closest_dist = 0.;
            }
            else if (T < computed_T[0]){
                sigma_min = sigma;
                sigma_eff = sigma;
                insert_idx = 0;
            }
            else {
                closest_idx = 0;
                for (size_t i = 0; i < computed_T.size(); i++){
                    if (computed_T[i] > T) {
                        closest_dist = computed_T[i - 1] - T;
                        break;
                    }
                    closest_idx = i;
                    insert_idx = i + 1;
                }
            }
        }
        if (closest_dist < 1.) do_insert = false; // Do not insert if we are closer than 1 K to existing saved value.
        sigma_eff = sigma_eff_cache[closest_idx];
        sigma_min = sigma_min_cache[closest_idx];
    }

    Sutherland::compute_sigma_eff();
    Sutherland::compute_epsilon_eff();
    if (do_insert) {
        computed_T.insert(computed_T.begin() + insert_idx, T);
        sigma_eff_cache.insert(sigma_eff_cache.begin() + insert_idx, sigma_eff);
        sigma_min_cache.insert(sigma_min_cache.begin() + insert_idx, sigma_min);
    }
}

std::vector<std::vector<double>> QuantumMie::get_BH_diameters(double T){
    set_temperature(T);
    return Sutherland::get_BH_diameters(T);
}

std::vector<std::vector<double>> QuantumMie::get_collision_diameters(double rho, double T, const std::vector<double>& x){
    set_temperature(T);
    return Sutherland::get_collision_diameters(rho, T, x);
}