/*
Author: Vegard Gjeldvik Jervell

Contains: The abstract class 'KineticGas', which computes the A_pqrl factors and corresponding square bracket integrals
            as derived by Thompson, Tipton and Lloyalka 
            (Chapman–Enskog solutions to arbitrary order in Sonine polynomials IV: Summational expressions for 
            the diffusion- and thermal conductivity-related bracket integrals, 
            [European Journal of Mechanics - B/Fluids, **28**, 6, pp. 695 - 721, 2009]
            (https://doi.org/10.1016/j.euromechflu.2009.05.002))
            
            Derived classes must implement the collision integral for a given potential model. 
            
            Subclasses must also implement the function 'model_rdf', ie. the radial distribution function
            at contact between particles.

            Finally, the function get_K_factors() must be implemented for all derived models. See the series of papers by 
            M. López de Haro, E. G. D. Cohen, and J. M. Kincaid: The Enskog theory for multicomponent mixtures. I - IV,
            Specifically: I. Linear transport theory
            (https://doi.org/10.1063/1.444985)

            Implementation is found in KineticGas.cpp
*/

#pragma once
#include "Factorial.h"
#include "global_params.h"
#include <vector>
#include <map>
#include <functional>
#include <thread>
#include <math.h>

/*
   To avoid unneccesary evaluations of the collision integrals, this struct is used to represent a point in 
   The five-dimensional (i, j, l, r, T)-space where the collision integral has been evaluated.

   NB: Resolution along the T-axis is 0.1 K, as (by experience) the collision integrals are a weak enough
   function of T to justify using T-values rounded to the nearest .1 K to improve speed, with little cost to precision.
*/
struct OmegaPoint{
    int i, j, l, r, T_dK;
    OmegaPoint(int i, int j, int l, int r, double T) : i{i}, j{j}, l{l}, r{r} {
         T_dK = (int) ((T * 10.0) + 0.5); // Temperature in dK (10^-1 K)
    };

    bool operator<(const OmegaPoint& other) const {
        if (i < other.i) return true;
        else if (i == other.i){
            if (j < other.j) return true;
            else if (j == other.j){
                if (l < other.l) return true;
                else if (l == other.l){
                    if (r < other.r) return true;
                    else if (r == other.r){
                        if (T_dK < other.T_dK) return true;
                    }
                }
            }
        }
        return false;
    }

};

class KineticGas{
    public:
    const unsigned long Ncomps; // must be ulong to be initialised from mole_weights.size()
    const bool is_idealgas;
    std::vector<double> m;
    std::vector<std::vector<double>> M, m0;
    std::map<OmegaPoint, double> omega_map;

    KineticGas(std::vector<double> mole_weights, bool is_idealgas);
    virtual ~KineticGas(){};
    // Collision integrals
    virtual double omega(const int& i, const int& j, const int& l, const int& r, const double& T) = 0;

    // The "distance between particles" at contact.
    virtual std::vector<std::vector<double>> get_contact_diameters(double rho, double T, const std::vector<double>& x) = 0;
    
    /* 
       The radial distribution function "at contact" for the given potential model
       model_rdf is only called if object is initialized with is_idealgas=true
       If a potential model is implemented only for the ideal gas state, its implementation
       of model_rdf should throw an std::invalid_argument error.
    */
    virtual std::vector<std::vector<double>> model_rdf(double T, double rho, const std::vector<double>& mole_fracs) = 0;

    std::vector<double> get_wt_fracs(const std::vector<double> mole_fracs); // Compute weight fractions from mole fractions
    std::vector<std::vector<double>> get_rdf(double rho, double T, const std::vector<double>& mole_fracs) {
        if (is_idealgas) return std::vector<std::vector<double>>(Ncomps, std::vector<double>(Ncomps, 1.));
        return model_rdf(rho, T, mole_fracs);
    }

// ------------------------------------------------------------------------------------------------------------------------ //
// -------------------------------------------------- Square bracket integrals -------------------------------------------- //

    // Linear combination weights by Tompson, Tipton and Lloyalka
    double A(const int& p, const int& q, const int& r, const int& l);
    double A_prime(const int& p, const int& q, const int& r, const int& l, const double& tmp_M1, const double& tmp_M2);
    double A_trippleprime(const int& p, const int& q, const int& r, const int& l);

    // The diffusion and conductivity related square bracket integrals
    double H_ij(const int& p, const int& q, const int& i, const int& j, const double& T); // [S^(p)_{3/2}(U^2_i), S^(q)_{3/2}(U^2_j)]_{ij}
    double H_i(const int& p, const int& q, const int& i, const int& j, const double& T);  // [S^(p)_{3/2}(U^2_i), S^(q)_{3/2}(U^2_i)]_{ij}
    double H_simple(const int& p, const int& q, const int& i, const double& T);           // [S^(p)_{3/2}(U^2_i), S^(q)_{3/2}(U^2_i)]_{i}

    // Linear combination weights by Tompson, Tipton and Lloyalka
    double B_prime(const int& p, const int& q, const int& r, const int& l, const double& M1, const double& M2);
    double B_dblprime(const int& p, const int& q, const int& r, const int& l, const double& M1, const double& M2);
    double B_trippleprime(const int& p, const int& q, const int& r, const int& l);

    // Viscosity related square bracket integrals
    double L_ij(const int& p, const int& q, const int& i, const int& j, const double& T); // [S^(p)_{5/2}(U^2_i), S^(q)_{5/2}(U^2_j)]_{ij}
    double L_i(const int& p, const int& q, const int& i, const int& j, const double& T);  // [S^(p)_{5/2}(U^2_i), S^(q)_{5/2}(U^2_i)]_{ij}
    double L_simple(const int& p, const int& q, const int& i, const double& T);           // [S^(p)_{5/2}(U^2_i), S^(q)_{5/2}(U^2_i)]_{i}

    // Bulk viscosity
    double Lb_ij(int p, int q, int i, int j, double T); // [S^(p)_{1/2}(U^2_i), S^(q)_{1/2}(U^2_j)]_{ij}
    double Lb_i(int p, int q, int i, int j, double T);  // [S^(p)_{1/2}(U^2_i), S^(q)_{1/2}(U^2_i)]_{ij}

// ---------------------------------------------------------------------------------------------------------------------------------------------- //
// ---------------------------------- Matrices and vectors for multicomponent, density corrected solutions -------------------------------------- //

    std::vector<std::vector<double>> get_conductivity_matrix(double rho, double T, const std::vector<double>& x, int N);
    std::vector<double> get_diffusion_vector(double rho, double T, const std::vector<double>& x, int N);
    double get_Lambda_ijpq(int i, int j, int p, int q, double rho, double T, const std::vector<double>& x);
    std::vector<std::vector<double>> get_diffusion_matrix(double rho, double T, const std::vector<double>& x, int N);
    std::vector<double> get_conductivity_vector(double rho, double T, const std::vector<double>& x, int N);
    std::vector<std::vector<double>> get_viscosity_matrix(double rho, double T, const std::vector<double>&x, int N);
    std::vector<double> get_viscosity_vector(double rho, double T, const std::vector<double>& x, int N);
    std::vector<std::vector<double>> get_bulk_viscosity_matrix(double rho, double T, const std::vector<double>&x, int N);
    std::vector<double> get_bulk_viscosity_vector(double rho, double T, double p, const std::vector<double>& x, int N);

    // Eq. (1.2) of 'multicomponent docs'
    std::vector<double> get_K_factors(double rho, double T, const std::vector<double>& mole_fracs);
    // Eq. (5.4) of 'multicomponent docs'
    std::vector<double> get_K_prime_factors(double rho, double T, const std::vector<double>& mole_fracs);
    // Eq. (S.16) of RET for Mie fluids, supporting material.
    std::vector<double> get_K_dblprime_factors(double rho, double T, double p, const std::vector<double>& mole_fracs);

// ----------------------------------------------------------------------------------------------------------------------------------- //
// --------------------------------------- Methods to facilitate multithreading ------------------------------------------------------ //
    /*
        The easiest part of the program to multithread is the computation of collision integrals. Each KineticGas
        object holds its own map of previously computed collision integrals in this->omega_map. Whenever `omega` is called,
        this map is checked to see if the integrals have been previously evaluated. The following methods divide a given
        set of collision integrals between different treads, so that when computing e.g. diffusion coefficients, the
        required collision integrals are first computed and stored by calling `omega` from `precompute_diffusion_omega`,
        then accessed by subsequent calls to `omega` by whatever method needs the integrals. In general, what you want
        to do is:
            1) Determine what collision integrals you need.
            2) Evaluate them by executing `omega` on different threads.
            3) Retrieve them by calling `omega` after all the treads are finished.
    */
    void precompute_omega(const std::vector<int>& i_vec, const std::vector<int>& j_vec,
                        const std::vector<int>& l_vec, const std::vector<int>& r_vec, double T);
    void precompute_conductivity_omega(int N, double T);
    void precompute_diffusion_omega(int N, double T);
    void precompute_th_diffusion_omega(int N, double T);
    void precompute_viscosity_omega(int N, double T);


};

inline int delta(int i, int j) {return (i == j) ? 1 : 0;}