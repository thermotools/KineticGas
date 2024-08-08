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
#include <iostream>

#ifdef NOPYTHON
    #include "eos_interface.h"
    #include <Eigen/Dense>
    #include <nlohmann/json.hpp>
    #include <cppThermopack/thermo.h>
    #include <memory>
    using json = nlohmann::json;
#endif

/*
   To avoid unneccesary evaluations of the collision integrals, this struct is used to represent a point in 
   The five-dimensional (i, j, l, r, T)-space where the collision integral has been evaluated.

   NB: Resolution along the T-axis is 0.1 K, as (by experience) the collision integrals are a weak enough
   function of T to justify using T-values rounded to the nearest .1 K to improve speed, with little cost to precision.
*/
struct StatePoint{
    int T_dK, rho;
    StatePoint(double T) : T_dK{static_cast<int>((T * 100.) + 0.5)} {}
    StatePoint(double T, double rho) : T_dK{static_cast<int>((T * 100.) + 0.5)}, rho{static_cast<int>(rho)}{}

    bool operator<(const StatePoint& other) const {
        if (T_dK < other.T_dK) return true;
        else if (T_dK == other.T_dK){
            if (rho < other.rho) return true;
        }
        return false;
    }
};

struct OmegaPoint{
    int i, j, l, r, T_dK, rho;
    OmegaPoint(int i, int j, int l, int r, double T) : i{i}, j{j}, l{l}, r{r}, rho{0} {
         T_dK = (int) ((T * 10.0) + 0.5); // Temperature in dK (10^-1 K)
    };

    OmegaPoint(int i, int j, int l, int r, double T, double rho_) : OmegaPoint(i, j, l, r, T){
        rho = static_cast<int>(rho_);
    }

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
                        else if (T_dK == other.T_dK){
                            if (rho < other.rho) return true;
                        }
                    }
                }
            }
        }
        return false;
    }
};

enum FrameOfReference{
        CoM,
        CoN,
        CoV,
        solvent,
        zarate,
        zarate_x,
        zarate_w
    };

class KineticGas{
    public:
    KineticGas(std::vector<double> mole_weights, bool is_idealgas, bool is_singlecomp);
    #ifdef NOPYTHON
        KineticGas(std::string comps, bool is_idealgas);
    #endif
    virtual ~KineticGas(){};

    // The transfer lengths related to momentum (MTL) and energy (ETL)
    // Inheriting classes must implement model_[mtl/etl]
    inline vector2d get_mtl(double rho, double T, const vector1d& x){
        set_internals(rho, T, x);
        return model_mtl(rho, T, x);
    }
    inline vector2d get_etl(double rho, double T, const vector1d& x){
        set_internals(rho, T, x);
        return model_etl(rho, T, x);
    }

    // Radial distribution function "at contact". Inheriting classes must implement model_rdf.
    inline vector2d get_rdf(double rho, double T, const vector1d& mole_fracs) {
        if (is_idealgas) return vector2d(Ncomps, vector1d(Ncomps, 1.));
        set_internals(rho, T, mole_fracs);
        return model_rdf(rho, T, mole_fracs);
    }

    #ifdef NOPYTHON
// ---------------------------------------------------------------------------------------------------------------------------------------------- //
// --------------------------------------------- Interfaces to compute transport coefficients --------------------------------------------------- //
// ------------------------ NOTE: Only compiled if compilation flag -DNOPYTHON is used ---------------------------------------- //

        Eigen::MatrixXd interdiffusion(double T, double Vm, const std::vector<double>& x, int N=2, int frame_of_reference=FrameOfReference::CoN, int dependent_idx=-1, int solvent_idx=-1, bool do_compress=true);
        double thermal_conductivity(double T, double Vm, const std::vector<double>& x, int N=2);
        double viscosity(double T, double Vm, const std::vector<double>& x, int N=2);
        Eigen::VectorXd thermal_diffusion_coeff(double T, double Vm, const std::vector<double>& x, int N=2, int frame_of_reference=FrameOfReference::CoN, int dependent_idx=-1, int solvent_idx=-1);
        Eigen::VectorXd thermal_diffusion_ratio(double T, double Vm, const std::vector<double>& x, int N=2);
        Eigen::MatrixXd thermal_diffusion_factor(double T, double Vm, const std::vector<double>& x, int N=2);
        Eigen::MatrixXd interdiffusion_dependent_CoM(double T, double Vm, const std::vector<double>& x, int N=2);

// ------------------------------------------------------------------------------------------------------------------------------------- //
// ----------------- Matrices and vectors for the sets of equations (6-10) in Revised Enskog Theory for Mie fluids  -------------------- //
// ----------------- doi : 10.1063/5.0149865, which are solved to obtain the Sonine polynomial expansion coefficients ------------------ //
// ----------------- for the velocity distribution functions. -------------------------------------------------------------------------- //
// ----------------- The methods compute_* solve the appropriate equations and return the expansion coefficients ----------------------- //
    
        Eigen::VectorXd compute_viscous_expansion_coeff(double rho, double T, const vector1d& x, int N);
        Eigen::VectorXd compute_thermal_expansion_coeff(double rho, double T, const vector1d& x, int N);
        Eigen::VectorXd compute_diffusive_expansion_coeff(double rho, double T, const vector1d& x, int N);
        std::vector<std::vector<std::vector<double>>> reshape_diffusive_expansion_vector(const Eigen::VectorXd& d_ijq);
        Eigen::VectorXd compute_dth_vector(const std::vector<std::vector<std::vector<double>>>& d_ijq, const Eigen::VectorXd& l);
    #endif

    std::vector<std::vector<double>> get_conductivity_matrix(double rho, double T, const std::vector<double>& x, int N);
    std::vector<double> get_diffusion_vector(double rho, double T, const std::vector<double>& x, int N);
    // double get_Lambda_ijpq(int i, int j, int p, int q, double rho, double T, const std::vector<double>& x);
    std::vector<std::vector<double>> get_diffusion_matrix(double rho, double T, const std::vector<double>& x, int N);
    std::vector<double> get_conductivity_vector(double rho, double T, const std::vector<double>& x, int N);
    std::vector<std::vector<double>> get_viscosity_matrix(double rho, double T, const std::vector<double>&x, int N);
    std::vector<double> get_viscosity_vector(double rho, double T, const std::vector<double>& x, int N);

// ----------------------------------------------------------------------------------------------------------------------------------- //
// -------------------------------------------------- Utility methods ---------------------------------------------------------------- //
    /*
        The CoM_to_FoR method is a dispatcher (switchboard) to the other CoM_to_* methods.
        The CoM_to_* methods return the transformation matrix (psi) used to transform diffusion coefficients from the 
        centre of mass (CoM) frame of reference (FoR) to the centre of moles (CoN), centre of volume (CoV) or solvent FoR
    */

    std::vector<double> get_wt_fracs(const std::vector<double> mole_fracs); // Compute weight fractions from mole fractions
    #ifdef NOPYTHON
        Eigen::MatrixXd CoM_to_FoR_matr(double T, double Vm, const std::vector<double>& x, int frame_of_reference, int solvent_idx);
    #endif
    vector1d get_K_factors(double rho, double T, const vector1d& mole_fracs); // Eq. (1.2) of 'multicomponent docs'
    vector1d get_K_prime_factors(double rho, double T, const vector1d& mole_fracs); // Eq. (5.4) of 'multicomponent docs'


// ------------------------------------------------------------------------------------------------------------------------ //
// --------------------------------------- KineticGas internals are below here -------------------------------------------- //
// -------------------------------- End users should not need to care about any of this ----------------------------------- //
protected:

    const size_t Ncomps;
    const bool is_idealgas;
    const bool is_singlecomp;

    std::vector<double> m;
    std::vector<std::vector<double>> M, m0;
    std::map<OmegaPoint, double> omega_map;
    std::map<StatePoint, vector2d> mtl_map;
    std::map<StatePoint, vector2d> etl_map;

    #ifdef NOPYTHON
        const std::vector<json> compdata;
        std::unique_ptr<GenericEoS> eos;
    #endif

    // set_internals is called at the start of all public methods. If a derived class needs to set any internals before running a computation,
    // it should be done by overriding this method.
    inline virtual void set_internals(double rho, double T, const vector1d& x){};
    virtual OmegaPoint get_omega_point(int i, int j, int l, int r, double T){
        return OmegaPoint(i, j, l, r, T);
    }
    virtual StatePoint get_transfer_length_point(double rho, double T, const vector1d& x){
        return StatePoint(T);
    }
    
    // Collision integrals
    virtual double omega(int i, int j, int l, int r, double T) = 0;

    /*
       The radial distribution function "at contact" for the given potential model. model_rdf is only called if object
       is initialized with is_idealgas=true. If a potential model is implemented only for the ideal gas state, its
       implementation of model_rdf should throw an std::invalid_argument error.
    */
    virtual vector2d model_rdf(double rho, double T, const vector1d& mole_fracs) = 0;

    virtual vector2d model_mtl(double rho, double T, const vector1d& x) = 0;
    virtual vector2d model_etl(double rho, double T, const vector1d& x) = 0;

// ----------------------------------------------------------------------------------------------------------------------------------- //
// --------------------------------------- Methods to facilitate multithreading ------------------------------------------------------ //
    /*
        Computing collision integrals and transfer lengths is the most computationally intensive part of computing a transport property,
        given that one does not use correlations for the computation. Therefore, computed collision integrals and transfer lengths are 
        stored in this->omega_map, this->mtl_map, and this->etl_map. The following methods are used to deduce which collision integrals 
        are required to compute a given property, and then split the computation of those integrals among several threads. When 
        computation of the property proceeds, the integrals are then retrieved from the maps.

        To adjust the number of threads used to compute collision integrals, set the compile-time constant Ncores in
        KineticGas_mthr.cpp. In practice, very little is to be gained by increasing this beyond 10, unless you are
        using high Enskog approximation orders (>4), or are working with a large number of components, and have a lot of
        cores available.

        These need to be protected instead of private, because QuantumMie needs to override them to prevent a race condition
        without locking everything with a mutex on every call to omega.
    */
    virtual void precompute_conductivity(int N, double T, bool precompute_etl=true);
    virtual void precompute_viscosity(int N, double T);
    virtual void precompute_omega(const std::vector<int>& i_vec, const std::vector<int>& j_vec,
                        const std::vector<int>& l_vec, const std::vector<int>& r_vec, double T);

private:

    #ifdef NOPYTHON
        Eigen::MatrixXd CoM_to_CoN_matr(double T, double Vm, const std::vector<double>& x);
        Eigen::MatrixXd CoM_to_solvent_matr(double T, double Vm, const std::vector<double>& x, int solvent_idx);
        Eigen::MatrixXd CoM_to_CoV_matr(double T, double Vm, const std::vector<double>& x);
        Eigen::MatrixXd get_zarate_X_matr(const std::vector<double>& x, int dependent_idx);
        Eigen::MatrixXd get_zarate_W_matr(const std::vector<double>& x, int dependent_idx);
        
        std::vector<std::vector<double>> get_chemical_potential_factors(double T, double Vm, const std::vector<double>& x);
        std::vector<double> get_ksi_factors(double T, double Vm, const std::vector<double>& x);
    #endif

    void set_masses(); // Precompute reduced mass of particle pairs (used only on init.)
    void precompute_diffusion(int N, double T); // Forwards call to precompute_conductivity_omega. Override that instead.
    void precompute_th_diffusion(int N, double T); // Forwards call to precompute_conductivity_omega. Override that instead.

// ------------------------------------------------------------------------------------------------------------------------ //
// ---------------------------------------------- Square bracket integrals ------------------------------------------------ //

    // Linear combination weights by Tompson, Tipton and Lloyalka
    double A(int p, int q, int r, int l) const;
    double A_prime(int p, int q, int r, int l, double tmp_M1, double tmp_M2) const;
    double A_trippleprime(int p, int q, int r, int l) const;

    // The diffusion and conductivity related square bracket integrals
    double H_ij(int p, int q, int i, int j, double T); // [S^(p)_{3/2}(U^2_i), S^(q)_{3/2}(U^2_j)]_{ij}
    double H_i(int p, int q, int i, int j, double T);  // [S^(p)_{3/2}(U^2_i), S^(q)_{3/2}(U^2_i)]_{ij}
    double H_simple(int p, int q, int i, double T);    // [S^(p)_{3/2}(U^2_i), S^(q)_{3/2}(U^2_i)]_{i}

    // Linear combination weights by Tompson, Tipton and Lloyalka
    double B_prime(int p, int q, int r, int l, double M1, double M2) const;
    double B_dblprime(int p, int q, int r, int l, double M1, double M2) const;
    double B_trippleprime(int p, int q, int r, int l) const;

    // Viscosity related square bracket integrals
    double L_ij(int p, int q, int i, int j, double T); // [S^(p)_{5/2}(U^2_i), S^(q)_{5/2}(U^2_j)]_{ij}
    double L_i(int p, int q, int i, int j, double T);  // [S^(p)_{5/2}(U^2_i), S^(q)_{5/2}(U^2_i)]_{ij}
    double L_simple(int p, int q, int i, double T);    // [S^(p)_{5/2}(U^2_i), S^(q)_{5/2}(U^2_i)]_{i}
};

inline int delta(int i, int j) {return (i == j) ? 1 : 0;}