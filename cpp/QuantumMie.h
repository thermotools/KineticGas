/*
Author: Vegard Gjeldvik Jervell
Contains: Implementation of the Feynman-Hibbs corrected Mie potential
            Article : Equation of state and force fields for Feynman–Hibbs-corrected Mie fluids.
                        II. Application to mixtures of helium, neon, hydrogen, and deuterium
            doi : https://doi.org/10.1063/1.5136079
Note : This class also overrides the methods of Spherical that call the potential function and derivatives,
        because the potetial is temperature dependent. The only change is that temperature is passed to the
        potential functions.
*/

#pragma once
#include "ExtendedSutherland.h"
#include "global_params.h"
#include <mutex>

class QuantumMie : public ExtSutherland {
    public:
    const std::vector<int> FH_order; // Feynman-Hibbs correction order
    vector3d Q_factors;
    vector1d computed_T;
    vector3d sigma_eff_cache;
    vector3d sigma_min_cache;

    QuantumMie(vector1d mole_weights, vector2d sigma, vector2d eps, vector2d la, vector2d lr, std::vector<int> FH_order,
                bool is_idealgas, bool is_singlecomp);

    double potential(int i, int j, double r, double T){
        set_C_eff(current_rho, T);
        return ExtSutherland::potential(i, j, r);
    }

    double potential_derivative_r(int i, int j, double r, double T){
        set_C_eff(current_rho, T);
        return ExtSutherland::potential_derivative_r(i, j, r);
    }

    double potential_dblderivative_rr(int i, int j, double r, double T){
        set_C_eff(current_rho, T);
        return ExtSutherland::potential_dblderivative_rr(i, j, r);
    }

    vector2d get_sigma_eff(double T){return ExtSutherland::get_sigma_eff(current_rho, T);}
    vector2d get_sigma_min(double T){return ExtSutherland::get_sigma_min(current_rho, T);}
    vector2d get_epsilon_eff(double T){return ExtSutherland::get_epsilon_eff(current_rho, T);}
    
    StatePoint get_transfer_length_point(double rho, double T, const vector1d& x) override {
        return StatePoint(current_rho, T);
    }

private:
    using ExtSutherland::potential;
    using ExtSutherland::potential_derivative_r;
    using ExtSutherland::potential_dblderivative_rr;

    inline double Q1(size_t i, size_t j, const vector2d& lamb){
        return lamb[i][j] * (lamb[i][j] - 1);
    }

    inline double Q2(size_t i, size_t j, const vector2d& lamb){
        return 0.5 * (lamb[i][j] + 2) * (lamb[i][j] + 1) * Q1(i, j, lamb);
    }

};