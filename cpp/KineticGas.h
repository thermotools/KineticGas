/*
Author: Vegard Gjeldvik Jervell

Contains: 
    The abstract class 'KineticGas', which is the interface to all transport property computations,
    and is responsible for all generic computations that are agnostic to the form of the collision integrals,
    radial distribution function, and transfer lengths.

    The classes constructor is also responsible for fetching parameter sets from the fluid database.

    Derived classes must implement collision integrals, the radial distribution function "at contact",
    and the transfer lengths.
            
References:
    (I) Chapman–Enskog solutions to arbitrary order in Sonine polynomials (I-IV)
            IV : Summational expressions for the diffusion- and thermal conductivity-related bracket integrals, Thompson, Tipton and Lloyalka 
                 European Journal of Mechanics - B/Fluids, 28, 6, pp. 695 - 721 (2009)
                 https://doi.org/10.1016/j.euromechflu.2009.05.002
            
    (II) The Enskog theory for multicomponent mixtures. I. Linear transport theory, 
         M. López de Haro, E. D. G. Cohen, and J. M. Kincaid, J. Chem. Phys. (1983)
         https://doi.org/10.1063/1.444985
    
    (III) Revised Enskog theory for Mie fluids: Prediction of diffusion coefficients, thermal diffusion coefficients,
            viscosities and thermal conductivities, V. G. Jervell and Ø. Wilhelmsen, J. Chem. Phys. (2023)

    (IV) Predicting viscosities and thermal conductivities from dilute gas to dense liquid: 
            Deriving fundamental transfer lengths for momentum and energy exchange in revised Enskog theory,
            V. G. Jervell and Ø. Wilhelmsen, J. Chem. Phys. (2024)
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
#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include <cppThermopack/thermo.h>
#include <memory>
#include "utils.h"
#include "eos_interface.h"

#ifdef PYLIB
    #include <pybind11/pybind11.h>
    namespace py = pybind11;
#endif

using json = nlohmann::json;

class KineticGas{
public:
    KineticGas(vector1d mole_weights, vector2d sigma, vector2d eps, bool is_idealgas, bool is_singlecomp);
    KineticGas(std::string comps, bool is_idealgas);

    virtual ~KineticGas(){};

    /********************************************** TRANSFER LENGTHS ************************************************/
    /*   The transfer lengths related to momentum (MTL) and energy (ETL)                                            */
    /*   Inheriting classes must implement the protected methods model_[mtl/etl]                                    */
    /*   Note: Some generic model_[mtl/etl] methods are implemented in Spherical, so if you are inherriting         */
    /*         Spherical, you can use get_mtl and get_etl once the potential is implemented.                        */
    /****************************************************************************************************************/
    inline vector2d get_mtl(double rho, double T, const vector1d& x){
        set_internals(rho, T, x);
        return model_mtl(rho, T, x);
    }
    inline vector2d get_etl(double rho, double T, const vector1d& x){
        set_internals(rho, T, x);
        return model_etl(rho, T, x);
    }

    // Radial distribution function "at contact". 
    // Inheriting classes must implement the protected method model_rdf.
    inline vector2d get_rdf(double rho, double T, const vector1d& mole_fracs) {
        if (is_idealgas) return vector2d(Ncomps, vector1d(Ncomps, 1.));
        set_internals(rho, T, mole_fracs);
        return model_rdf(rho, T, mole_fracs);
    }

// ---------------------------------------------------------------------------------------------------------------------------------------------- //
// --------------------------------------------- Interfaces to compute transport coefficients --------------------------------------------------- //
// ---------------------------------------------------------------------------------------------------------------------------------------------- // 

    Eigen::MatrixXd interdiffusion(double T, double Vm, const std::vector<double>& x, int N=2, int frame_of_reference=FrameOfReference::CoN, int dependent_idx=-1, int solvent_idx=-1, bool do_compress=true);
    double thermal_conductivity(double T, double Vm, const std::vector<double>& x, int N=2);
    double thermal_diffusivity(double T, double Vm, const vector1d& x, int N=2);
    double viscosity(double T, double Vm, const std::vector<double>& x, int N=2);
    double kinematic_viscosity(double T, double Vm, const vector1d& x, int N=2);
    Eigen::VectorXd thermal_diffusion_coeff(double T, double Vm, const std::vector<double>& x, int N=2, int frame_of_reference=FrameOfReference::CoN, int dependent_idx=-1, int solvent_idx=-1);
    Eigen::VectorXd thermal_diffusion_ratio(double T, double Vm, const std::vector<double>& x, int N=2);
    Eigen::MatrixXd thermal_diffusion_factor(double T, double Vm, const std::vector<double>& x, int N=2);
    Eigen::MatrixXd interdiffusion_dependent_CoM(double T, double Vm, const std::vector<double>& x, int N=2);
    Eigen::VectorXd soret_coefficient(double T, double Vm, const std::vector<double>& x, int N, int dependent_idx=-1);
    std::map<std::string, double> thermal_conductivity_contributions(double T, double Vm, const std::vector<double>& x, int N=2, std::string contribs="tdi");

    // ------------------------------------------------------------------------------------------------------------------- //
    // ----------- TP-interface methods: These just compute molar volume and feed the call to the methods above ---------- //

    inline Eigen::MatrixXd interdiffusion_tp(double T, double p, const vector1d& x, int N=2, int frame_of_reference=FrameOfReference::CoN, int dependent_idx=-1, int solvent_idx=-1, bool do_compress=true){return interdiffusion(T, eos->specific_volume(T, p, sanitize_mole_fracs_eos(x), eos->VAPPH), x, N, frame_of_reference, dependent_idx, solvent_idx, do_compress);}
    inline double thermal_conductivity_tp(double T, double p, const std::vector<double>& x, int N=2){return thermal_conductivity(T, eos->specific_volume(T, p, sanitize_mole_fracs_eos(x), eos->VAPPH), x, N);}
    inline double thermal_diffusivity_tp(double T, double p, const vector1d& x, int N=2){return thermal_diffusivity(T, eos->specific_volume(T, p, sanitize_mole_fracs_eos(x), eos->VAPPH), x, N);}
    inline double viscosity_tp(double T, double p, const std::vector<double>& x, int N=2){return viscosity(T, eos->specific_volume(T, p, sanitize_mole_fracs_eos(x), eos->VAPPH), x, N);}
    inline double kinematic_viscosity_tp(double T, double p, const vector1d& x, int N=2){return kinematic_viscosity(T, eos->specific_volume(T, p, sanitize_mole_fracs_eos(x), eos->VAPPH), x, N);}
    inline Eigen::VectorXd thermal_diffusion_coeff_tp(double T, double p, const std::vector<double>& x, int N=2, int frame_of_reference=FrameOfReference::CoN, int dependent_idx=-1, int solvent_idx=-1){return thermal_diffusion_coeff(T, eos->specific_volume(T, p, sanitize_mole_fracs_eos(x), eos->VAPPH), x, N, frame_of_reference, dependent_idx, solvent_idx);}
    inline Eigen::VectorXd thermal_diffusion_ratio_tp(double T, double p, const std::vector<double>& x, int N=2){return thermal_diffusion_ratio(T, eos->specific_volume(T, p, sanitize_mole_fracs_eos(x), eos->VAPPH), x, N);}
    inline Eigen::MatrixXd thermal_diffusion_factor_tp(double T, double p, const std::vector<double>& x, int N=2){return thermal_diffusion_factor(T, eos->specific_volume(T, p, sanitize_mole_fracs_eos(x), eos->VAPPH), x, N);}
    inline Eigen::VectorXd soret_coefficient_tp(double T, double p, const std::vector<double>& x, int N, int dependent_idx=-1){return soret_coefficient(T, eos->specific_volume(T, p, sanitize_mole_fracs_eos(x), eos->VAPPH), x, N, dependent_idx);}

// ------------------------------------------------------------------------------------------------------------------------------------- //
// ----------------- Matrices and vectors for the sets of equations (6-10) in Ref. (III)  ---------------------------------------------- //
// ----------------- doi : 10.1063/5.0149865, which are solved to obtain the Sonine polynomial expansion coefficients ------------------ //
// ----------------- for the velocity distribution functions. -------------------------------------------------------------------------- //
// ----------------- The methods compute_* solve the appropriate equations and return the expansion coefficients ----------------------- //
    
    Eigen::VectorXd compute_viscous_expansion_coeff(double rho, double T, const vector1d& x, int N);
    Eigen::VectorXd compute_thermal_expansion_coeff(double rho, double T, const vector1d& x, int N);
    Eigen::VectorXd compute_diffusive_expansion_coeff(double rho, double T, const vector1d& x, int N);
    std::vector<std::vector<std::vector<double>>> reshape_diffusive_expansion_vector(const Eigen::VectorXd& d_ijq);
    Eigen::VectorXd compute_dth_vector(const std::vector<std::vector<std::vector<double>>>& d_ijq, const Eigen::VectorXd& l);

    std::vector<std::vector<double>> get_conductivity_matrix(double rho, double T, const std::vector<double>& x, int N);
    std::vector<double> get_conductivity_vector(double rho, double T, const std::vector<double>& x, int N);
    std::vector<double> get_diffusion_vector(double rho, double T, const std::vector<double>& x, int N);
    std::vector<std::vector<double>> get_diffusion_matrix(double rho, double T, const std::vector<double>& x, int N);
    std::vector<std::vector<double>> get_viscosity_matrix(double rho, double T, const std::vector<double>&x, int N);
    std::vector<double> get_viscosity_vector(double rho, double T, const std::vector<double>& x, int N);

    vector1d get_K_factors(double rho, double T, const vector1d& mole_fracs); // Eq. (1.2) of 'multicomponent docs'
    vector1d get_K_prime_factors(double rho, double T, const vector1d& mole_fracs); // Eq. (5.4) of 'multicomponent docs'

// ----------------------------------------------------------------------------------------------------------------------------------- //
// -------------------------------------------------- Utility methods ---------------------------------------------------------------- //
    /*
        The CoM_to_FoR method is a dispatcher (switchboard) to the other CoM_to_* methods.
        The CoM_to_* methods return the transformation matrix (psi) used to transform diffusion coefficients from the 
        centre of mass (CoM) frame of reference (FoR) to the centre of moles (CoN), centre of volume (CoV) or solvent FoR
    */
    virtual Units get_reducing_units(int ci, int cj); // Return a `Units` struct holding the reducing units created from the ci-cj potential parameters.
    std::vector<double> get_wt_fracs(const std::vector<double> mole_fracs); // Compute weight fractions from mole fractions

    Eigen::MatrixXd CoM_to_FoR_matr(double T, double Vm, const std::vector<double>& x, int frame_of_reference, int solvent_idx);
    Eigen::MatrixXd CoM_to_CoN_matr(double T, double Vm, const std::vector<double>& x);
    Eigen::MatrixXd CoM_to_solvent_matr(double T, double Vm, const std::vector<double>& x, int solvent_idx);
    Eigen::MatrixXd CoM_to_CoV_matr(double T, double Vm, const std::vector<double>& x);
    Eigen::MatrixXd get_zarate_X_matr(const std::vector<double>& x, int dependent_idx);
    Eigen::MatrixXd get_zarate_W_matr(const std::vector<double>& x, int dependent_idx);
    
    std::vector<std::vector<double>> get_chemical_potential_factors(double T, double Vm, const std::vector<double>& x);
    std::vector<double> get_ksi_factors(double T, double Vm, const std::vector<double>& x);

    vector1d sanitize_mole_fracs(const vector1d& x);
    vector1d sanitize_mole_fracs_eos(const vector1d& x);
    #ifdef PYLIB
        void set_eos(py::object eos_){
            eos = std::make_unique<GenericEoS>(PyWrapper(eos_));
        }
    #endif
    void set_eos(GenericEoS&& other){
        eos = std::make_unique<GenericEoS>(std::move(other));
    }

    // Different transfer length models, see Ref. (IV)
    void set_transfer_length_model(int model_id);
    std::pair<int, std::string> get_transfer_length_model(); // Return the current transfer length model
    std::map<int, std::string> get_valid_transfer_length_models(); // Get a map of valid models with descriptions

    int frame_of_reference_map(std::string frame_of_ref);

// ------------------------------------------------------------------------------------------------------------------------ //
// --------------------------------------- KineticGas internals are below here -------------------------------------------- //
// -------------------------------- End users should not need to care about any of this ----------------------------------- //
protected:

    const size_t Ncomps;
    const bool is_idealgas;
    const bool is_singlecomp;

    vector1d m; // Particle masses (kg)
    vector2d M, m0, red_mass; // Various combinations of particle masses that show up often
    std::map<OmegaPoint, double> omega_map;
    std::map<StatePoint, vector2d> mtl_map;
    std::map<StatePoint, vector2d> etl_map;

    // In the general case, sigma and eps are scaling parameters for the molecular interaction, 
    // with sigma being the length scale (m) and eps being the energy scale (J).
    // In general, these are just used for convenience to make things non-dimensional. If your potential
    // model does not use them (like HardSphere, which has no energy scale), just set them to dummy-values.
    vector2d sigma, eps;

    std::unique_ptr<GenericEoS> eos;
    const std::vector<json> compdata; // Fluid data for each component, 

    const int default_tl_model_id = TransferLengthModel::EWCA; // Default transfer length model
    int transfer_length_model_id = default_tl_model_id; // Currently active transfer length model

    // set_internals is called at the start of all public methods. 
    // If a derived class needs to set any internals before running a computation,
    // it should be done by overriding this method.
    inline virtual void set_internals(double rho, double T, const vector1d& x){};

    // Implementations of omega, model_mtl, and model_etl, may wish to use the omega_map, mtl_map and etl_map
    // To store values for lazy evaluation (massive speedups)
    // In that case, using these methods to get the state point of the computations ensures that inherriting 
    // classes can override these if they want to handle state points differently 
    // For example, if you are using a density-dependent potential, you will want to include the density in the OmegaPoint.
    // See: Spherical for example.
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
       is initialized with is_idealgas=false. If a potential model is implemented only for the ideal gas state, its
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
        using high Enskog approximation orders (>5), or are working with a large number of components, and have a lot of
        cores available.
    */
    virtual void precompute_conductivity(int N, double T, double rho, bool precompute_etl=true);
    virtual void precompute_viscosity(int N, double T, double rho);
    virtual void precompute_omega(const std::vector<int>& i_vec, const std::vector<int>& j_vec,
                        const std::vector<int>& l_vec, const std::vector<int>& r_vec, double T);

private:

    void set_masses(); // Precompute reduced mass of particle pairs (used only on init.)
    void precompute_diffusion(int N, double T, double rho); // Forwards call to precompute_conductivity_omega. Override that instead.
    void precompute_th_diffusion(int N, double T, double rho); // Forwards call to precompute_conductivity_omega. Override that instead.

// ------------------------------------------------------------------------------------------------------------------------ //
// ---------------------------------------------- Square bracket integrals ------------------------------------------------ //

    // Linear combination weights from Ref. (IV)
    double A(int p, int q, int r, int l) const;
    double A_prime(int p, int q, int r, int l, double tmp_M1, double tmp_M2) const;
    double A_trippleprime(int p, int q, int r, int l) const;

    // The diffusion and conductivity related square bracket integrals, Eq. (4) of Ref. (II)
    double H_ij(int p, int q, int i, int j, double T); // [S^(p)_{3/2}(U^2_i), S^(q)_{3/2}(U^2_j)]_{ij}
    double H_i(int p, int q, int i, int j, double T);  // [S^(p)_{3/2}(U^2_i), S^(q)_{3/2}(U^2_i)]_{ij}
    double H_simple(int p, int q, int i, double T);    // [S^(p)_{3/2}(U^2_i), S^(q)_{3/2}(U^2_i)]_{i}

    // Linear combination weights from Ref. (IV)
    double B_prime(int p, int q, int r, int l, double M1, double M2) const;
    double B_dblprime(int p, int q, int r, int l, double M1, double M2) const;
    double B_trippleprime(int p, int q, int r, int l) const;

    // Viscosity related square bracket integrals, Eq. (4) of Ref. (II)
    double L_ij(int p, int q, int i, int j, double T); // [S^(p)_{5/2}(U^2_i), S^(q)_{5/2}(U^2_j)]_{ij}
    double L_i(int p, int q, int i, int j, double T);  // [S^(p)_{5/2}(U^2_i), S^(q)_{5/2}(U^2_i)]_{ij}
    double L_simple(int p, int q, int i, double T);    // [S^(p)_{5/2}(U^2_i), S^(q)_{5/2}(U^2_i)]_{i}
};

inline int delta(int i, int j) {return (i == j) ? 1 : 0;}