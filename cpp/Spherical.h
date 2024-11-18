/*
Author: Vegard Gjeldvik Jervell
Contains: The abstract class 'Spherical', inheriting from 'KineticGas'. This class implements
        Methods to evaluate the collision integrals for an arbitrary spherically symmetric potential.

        Subclasses must implement the potential, as well as the first and second derivative wrt. distance.

References:
    (I) The Kinetic Gas theory of Mie fluids, V. G. Jervell, Norwegian university of science and technology (2022)
    
    (II) Revised Enskog theory for Mie fluids: Prediction of diffusion coefficients, thermal diffusion coefficients,
            viscosities and thermal conductivities, V. G. Jervell and Ø. Wilhelmsen, J. Chem. Phys. (2023)

    (III) Predicting viscosities and thermal conductivities from dilute gas to dense liquid: 
            Deriving fundamental transfer lengths for momentum and energy exchange in revised Enskog theory,
            V. G. Jervell and Ø. Wilhelmsen, J. Chem. Phys. (2024)
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
    Spherical(vector1d mole_weights, vector2d sigma, vector2d eps, bool is_idealgas, bool is_singlecomp);
    Spherical(std::string comps, bool is_idealgas) : KineticGas(comps, is_idealgas) {}
    
    virtual ~Spherical(){};

    // Potential model, the dual case must be overridden. Overriding the other cases can give improved efficiency.
    virtual dual2 potential(int i, int j, dual2 r) = 0;
    virtual double potential(int i, int j, double r);
    virtual double potential_derivative_r(int i, int j, double r);
    virtual double potential_dblderivative_rr(int i, int j, double r);

    virtual double cross_section(int i, int j, int l, double E);
    double hs_cross_section(int i, int j, int l);
    double reduced_cross_section(int i, int j, int l, double E); // Reduced using corresponding Hard sphere value

    double omega(int i, int j, int l, int r, double T) override;
    vector2d model_rdf(double rho, double T, const vector1d& mole_fracs) override = 0;
    vector2d model_mtl(double rho, double T, const vector1d& x) override; // Momentum transfer length
    vector2d model_etl(double rho, double T, const vector1d& x) override; // Energy transfer length

    // ------------------------------------------------------------------------------------------- //
    // -------------------- GENERIC METHODS DESCRIBING A COLLISION TRAJECTORY -------------------- //
    // ------------------------------------------------------------------------------------------- //
    double theta(int i, int j, const double T, const double g, const double b); // Angular coordinate at distance of closest approach.
    double chi(int i, int j, double T, double g, double b); // Deflection angle at given temperature, dimensionless velocity and impact parameter
    double get_R(int i, int j, double T, double g, double b); // Distance of closest approach at given temperature, dimensionless velocity and impact parameter
    double theta_r(int i, int j, double r, double T, double g, double b); // Angular position at given particle separation (r). A plot of r * sin(theta_r) vs. r * cos(theta_r) will show the particle trajectory for a collision with given T, g, b.
    double theta_r(int i, int j, double R, double r, double T, double g, double b); // Faster, if R is known.

protected:
    // ------------------------------------------------------------------------------------------- //
    // --------------------------------  TRANSFER LENGTH CACHING --------------------------------- //
    /*  
        Because transfer lengths can be expensive to compute, we want to save computed values 
        so that we can re-use them. Importantly, some transfer length models are only functions 
        of temperature, while others are functions of both temperature and density. For optimal 
        efficiency, saved transfer lengths are stored in an std::map<StatePoint, std::vector<std::vector<double>>.

        When a transfer length is computed using get_transfer_length, the StatePoint at which 
        it is being computed is retrieved using get_transfer_length_point, which regulates whether 
        the StatePoint should be a unique value for every (temperature, density) combination,
        or only unique for every temperature (i.e. not a function of density).
    */
    StatePoint get_transfer_length_point(double rho, double T, const vector1d& x) override;
    vector2d get_transfer_length(double rho, double T, const vector1d& x, int property);

    /**************************************************************************************************/
    /**********************         TL MODEL 0 : Collision diameter              **********************/
    /**************************************************************************************************/
    vector2d get_collision_diameters(double rho, double T, const vector1d& x);
    virtual vector2d get_b_max(double T); // Solve Eq. (37) in Ref. (II), and handle failure to solve
    vector2d get_b_max(double T, std::vector<std::vector<int>>& ierr); // Attempt to solve Eq. (37) in Ref (II).

    /*************************************************************************************************************************/
    /**********************         TL MODEL 1 : Exchange weighted closest approach (EWCA)              **********************/
    /*************************************************************************************************************************/
    double tl_ewca(int i, int j, double T, int property); // Evaluate Eq. (13) or (20) in Ref. (III), depending on `property`.
    std::pair<double, double> ewca_inner(int i, int j, double T, double g, int property); // Inner integrals (numerator, denominator) of Eq. (13) or (20) in Ref. (III), depending on `property`.
    double ewca_weight(int i, int j, double T, double g, double b, int property); // Integrand of Eq. (19) or first line of Eq. (11) in Ref. (III), depending on `property`
    double momentum_transfer(int i, int j, double T, double g, double b); // First line of Eq. (11) in Ref. (III).
    double energy_transfer(int i, int j, double T, double g, double b); // Integrand of Eq. (19) in Ref. (III).

    double get_b_max_g(int i, int j, double g, double T, double bmid); // Find b such that eps < chi(b) < 0, for small eps, when (d chi / db) > 0.
    double get_bmid(int i, int j, double g, double T); // Solve chi = 0 for b.

    /********************************************************************************************/
    /**********************         TL MODEL 2 : correlations              **********************/
    /********************************************************************************************/
    vector2d MTL_correlation(double rho, double T); // Eq. (22) in Ref. (III)
    vector2d ETL_correlation(double rho, double T); // Eq. (22) in Ref. (III)

private:
    // ------------------------------------------------------------------------------------------------------------------- //
    // -------------------------------- Spherical Internals are below here ----------------------------------------------- //
    // ---------------------- End users should not need to care about anything below ------------------------------------- //
    // ------------------------------------------------------------------------------------------------------------------- //
    double omega_tester(int i, int j, int l, int r, double T, IntegrationParam& param);
    double w_integral_tester(int i, int j, double T, int l, int r, IntegrationParam& param);

    double get_R0(int i, int j, double T, double g); // Solve get_R when b = 0

    double theta_integral(int i, int j, const double T, const double R, const double g, const double b); // Evaluate Eq. (49) in Ref. (II)

    double theta_lim(int i, int j, const double T, const double g); // Get theta for very large b
    double theta_integrand(int i, int j, double T, double r, double g, double b); // Eq. (48) in Ref. (II)
    double transformed_theta_integrand(int i, int j, double T, double u, double R, double g, double b); // Integrand of Eq. (50) in Ref. (II)
    double theta_integrand_dblderivative(int i, int j, double T, double r, double g, double b);
    double get_R_rootfunc(int i, int j, double T, double g, double b, double& r); // Eq. (45) in Ref. (II)
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