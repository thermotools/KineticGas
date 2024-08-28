/*
Author: Vegard Gjeldvik Jervell
Contains: The abstract class 'Spherical', inheriting from 'KineticGas'. This class implements
        Methods to evaluate the collision integrals for an arbitrary spherically symmetric potential.

        Subclasses must implement the potential, as well as the first and second derivative wrt. distance.

References:
    (I) The Kinetic Gas theory of Mie fluids, V. G. Jervell, Norwegian university of science and technology (2022)
    (II) Revised Enskog theory for Mie fluids: Prediction of diffusion coefficients, thermal diffusion coefficients,
    viscosities and thermal conductivities, V. G. Jervell and Ø. Wilhelmsen, J. Chem. Phys. (2023)
*/

#pragma once
#include "KineticGas.h"
#include "Integration/Integration.h"
#include "global_params.h"
#include <autodiff/forward/dual.hpp>

using dual = autodiff::dual;
using dual2 = autodiff::dual2nd;

class IntegrationParam;

class Spherical : public KineticGas {
    public:
    Spherical(vector1d mole_weights,
                vector2d sigmaij,
                bool is_idealgas, bool is_singlecomp);
    
    Spherical(vector1d mole_weights,
                vector2d sigmaij, vector2d eps,
                bool is_idealgas, bool is_singlecomp);

    Spherical(std::vector<double> mole_weights, 
                std::vector<std::vector<double>> sigmaij,
                bool is_idealgas);
    
    Spherical(std::string comps, bool is_idealgas) : KineticGas(comps, is_idealgas) {}
    
    virtual ~Spherical(){};

    // Potential model, the dual case must be overridden. Overriding the other cases can give improved 
    // efficiency of the code.
    virtual dual2 potential(int i, int j, dual2 r) = 0;

    virtual double potential(int i, int j, double r){
        return potential(i, j, static_cast<dual>(r)).val.val;
    };
    virtual double potential_derivative_r(int i, int j, double r){
        dual2 rd = r;
        const auto func = [&](dual2 r_){return potential(i, j, r_);};
        auto [u0, ur, urr] = autodiff::derivatives(func, autodiff::wrt(rd), autodiff::at(rd));
        return ur;
    }
    virtual double potential_dblderivative_rr(int i, int j, double r){
        dual2 rd = r;
        const auto func = [&](dual2 r_){return potential(i, j, r_);};
        auto [u0, ur, urr] = autodiff::derivatives(func, autodiff::wrt(rd, rd), autodiff::at(rd));
        return urr;
    }

    double omega(int i, int j, int l, int r, double T) override;

    vector2d model_rdf(double rho, double T, const vector1d& mole_fracs) override = 0;
    vector2d get_mtl(double rho, double T, const vector1d& x) override; // Momentum transfer length
    vector2d get_etl(double rho, double T, const vector1d& x) override; // Energy transfer length

    // Different collision diameter models are
    // -1 : Default model
    // 0 : Model presented in Refs. (I) and (II) (see top of file)
    // 1 : Unpublished model for momentum- and energy transfer lengths (MTL and ETL) (2024)
    void set_transfer_length_model(int model_id){
        if (model_id == -1) transfer_length_model_id = default_tl_model_id;
        transfer_length_model_id = model_id;
    }
    const int default_tl_model_id = 1; // Default transfer length model
    int transfer_length_model_id = default_tl_model_id;

    double theta(int i, int j, const double T, const double g, const double b); // Angular coordinate at distance of closest approach.
    double chi(int i, int j, double T, double g, double b); // Deflection angle at given temperature, dimensionless velocity and impact parameter
    double get_R(int i, int j, double T, double g, double b); // Distance of closest approach at given temperature, dimensionless velocity and impact parameter
    // Angular position at given particle separation (r). A plot of r * sin(theta_r) vs. r * cos(theta_r) will show the particle trajectory for a collision with given T, g, b.
    double theta_r(int i, int j, double r, double T, double g, double b);

    double omega_tester(int i, int j, int l, int r, double T, IntegrationParam& param);
    double w_integral_tester(int i, int j, double T, int l, int r, IntegrationParam& param);

    protected:
    // In the general case, sigma is a scaling parameter on the order of the molecular size (m). Its specific physical meaning
    // may be different for different potential models. Likewise, 'eps' is a scaling parameter describing the molecular
    // interaction potential (J).
    vector2d sigma;
    vector2d eps;

    vector2d get_transfer_length(double rho, double T, const vector1d& x, int property);
    vector2d get_collision_diameters(double rho, double T, const vector1d& x); // Obsolete "collision diameter" model

    /*****************************************************************************/
    /**********************         TL MODEL 1              **********************/
    /*****************************************************************************/
    double tl_inner(int i, int j, double T, double g, double I, int property);
    double tl_integrand(int i, int j, double T, double g, double b, double I, double bmax, int property);
    double momentum_transfer(int i, int j, double T, double g, double b);
    double energy_transfer(int i, int j, double T, double g, double b);

    double get_b_max_g(int i, int j, double g, double T); // Find b such that eps < chi(b) < 0, for small eps.
    double get_bmid(int i, int j, double g, double T); // Solve chi = 0

    double tl_weight_integrand(int i, int j, double T, double g, double b, double bmax, int property);
    double tl_weight_inner(int i, int j, double T, double g, int property);
    double get_tl_weight_normalizer(int i, int j, double T, int property);
    double get_tl_weight(int i, int j, double T, double g, double b, double I, double bmax, int property);

    // ------------------------------------------------------------------------------------------------------------------- //
    // -------------------------- Spherical Internals are below here ----------------------------------------------------- //
    // ---------------------- End users should not need to care about anything below -------------------------------------- //
    // ------------------------------------------------------------------------------------------------------------------- //
    double get_R0(int i, int j, double T, double g); // Solve get_R when b = 0

    // This method is responsible for calling get_b_max(T, ierr), and handling eventual failures.
    // Most inheritance should only require overriding the failure handling, not the computation in get_b_max(T, ierr)
    virtual vector2d get_b_max(double T);
    vector2d get_b_max(double T, std::vector<std::vector<int>>& ierr);

    double theta_r(int i, int j, double R, double r, double T, double g, double b); 
    double theta_integral(int i, int j, const double T, const double R, const double g, const double b);

    // private:
    // Helper functions for computing dimentionless collision integrals

    double theta_lim(int i, int j, const double T, const double g);
    double theta_integrand(int i, int j, double T, double r, double g, double b);
    double transformed_theta_integrand(int i, int j, double T, double u, double R, double g, double b);
    double theta_integrand_dblderivative(int i, int j, double T, double r, double g, double b);
    double get_R_rootfunc(int i, int j, double T, double g, double b, double& r);
    double get_R_rootfunc_derivative(int i, int j, double T, double g, double b, double& r);

    double w_integral(int i, int j, double T, int l, int r); // Dimentionless collision integral for spherical potentials
    double w_integrand(int i, int j, double T, double g, double b, int l, int r);
};

class IntegrationParam{
    public:
    Point origin{1e-7, 1e-7};
    Point end{8, 5};
    double dg{0.5}, db{0.03125};
    int refinement_levels_g{4};
    int refinement_levels_b{16};
    double subdomain_dblder_limit{1e-5};

    IntegrationParam(vector1d origin, vector1d end, double dg, double db, int rg, int rb, double dd_lim)
        : origin{origin[0], origin[1]}, end{end[0], end[1]}, dg{dg}, db{db}, refinement_levels_g{rg},
        refinement_levels_b{rb}, subdomain_dblder_limit{dd_lim}
        {}

    void set_end(vector1d end_){end = Point(end_[0], end_[1]);}
    void set_dg(double dg_){dg = dg_;}
    void set_db(double db_){db = db_;}
    void set_rlg(int rlg){refinement_levels_g = rlg;}
    void set_rlb(int rlb){refinement_levels_b = rlb;}
    void set_dd_lim(double lim){subdomain_dblder_limit = lim;}
};