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
#include "Sutherland.h"
#include "global_params.h"
#include <mutex>

class QuantumMie : public Sutherland{
    public:
    const std::vector<int> FH_order; // Feynman-Hibbs correction order
    vector3d Q_factors;
    vector1d computed_T;
    vector3d sigma_eff_cache;
    vector3d sigma_min_cache;

    QuantumMie(vector1d mole_weights, vector2d sigma, vector2d eps, vector2d la, vector2d lr, std::vector<int> FH_order, bool is_idealgas);

    double omega(int i, int j, int l, int r, double T) override {
        set_temperature(T);
        return Sutherland::omega(i, j, l, r, T);
    }

    double potential(int i, int j, double r, double T){
        set_temperature(T);
        return Sutherland::potential(i, j, r);
    }

    double potential_derivative_r(int i, int j, double r, double T){
        set_temperature(T);
        return Sutherland::potential_derivative_r(i, j, r);
    }

    double potential_dblderivative_rr(int i, int j, double r, double T){
        set_temperature(T);
        return Sutherland::potential_dblderivative_rr(i, j, r);
    }

    inline vector2d get_sigma_eff(double T){set_temperature(T); return sigma_eff;}
    inline vector2d get_sigma_min(double T){set_temperature(T); return sigma_min;}
    inline vector2d get_epsilon_eff(double T){set_temperature(T); return eps_eff;}
    inline vector2d get_vdw_alpha(double T){set_temperature(T); return vdw_alpha;}

    vector2d get_BH_diameters(double T) override;
    vector2d get_collision_diameters(double rho, double T, const std::vector<double>& x) override;
    vector2d model_rdf(double rho, double T, const std::vector<double>& x) override {
        set_temperature(T);
        return Sutherland::model_rdf(rho, T, x);
    }
    inline vector2d saft_rdf(double rho, double T, const std::vector<double>& x, int order=2, bool g2_correction=true) override {
        set_temperature(T);
        return Sutherland::saft_rdf(rho, T, x, order, g2_correction);
    }

    inline vector2d a1_func(double rho, double T, const vector1d& x) override {
        set_temperature(T);
        return Sutherland::a1_func(rho, T, x);
    }

    inline vector2d a_1s_func(double rho, double T, const vector1d& x, const vector2d& lambda) override {
        set_temperature(T);
        return Sutherland::a_1s_func(rho, T, x, lambda);
    }

    inline vector2d B_func(double rho, double T, const vector1d& x, const vector2d& lambda) override {
        set_temperature(T);
        return Sutherland::B_func(rho, T, x, lambda);
    }

    inline vector2d da1_drho_func(double rho, double T, const vector1d& x) override {
        set_temperature(T);
        return Sutherland::da1_drho_func(rho, T, x);
    }

    inline vector2d a2ij_div_chi_func(double rho, double T, const vector1d& x) override {
        set_temperature(T);
        return Sutherland::a2ij_div_chi_func(rho, T, x);
    }

    inline vector2d da2ij_div_chi_drho_func(double rho, double T, const vector1d& x) override {
        set_temperature(T);
        return Sutherland::da2ij_div_chi_drho_func(rho, T, x);
    }

    private:
    double current_temperature = -1.;
    void set_temperature(double T);
    void set_sigma_to_cached(double T, size_t* insert_idx_p);
    void compute_sigma_eps_eff(double T);

    using Sutherland::potential;
    using Sutherland::potential_derivative_r;
    using Sutherland::potential_dblderivative_rr;

    inline double Q1(size_t i, size_t j, const vector2d& lamb){
        return lamb[i][j] * (lamb[i][j] - 1);
    }

    inline double Q2(size_t i, size_t j, const vector2d& lamb){
        return 0.5 * (lamb[i][j] + 2) * (lamb[i][j] + 1) * Q1(i, j, lamb);
    }

    // NOTE: This method is private, because set_temperature MUST called before calling this method.
    // If you need to compute a collision integral, use the public override of omega.
    void precompute_omega(const std::vector<int>& i_vec, const std::vector<int>& j_vec,
                        const std::vector<int>& l_vec, const std::vector<int>& r_vec, double T) override
     {
        for (int idx = 0; idx < i_vec.size(); idx++){
            Sutherland::omega(i_vec[idx], j_vec[idx], l_vec[idx], r_vec[idx], T); // Computed values are stored in omega_map
        }
     }

    void precompute_conductivity_omega(int N, double T) override {
        set_temperature(T);
        KineticGas::precompute_conductivity_omega(N, T);
    }
    void precompute_viscosity_omega(int N, double T) override {
        set_temperature(T);
        KineticGas::precompute_viscosity_omega(N, T);
    }
};