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

using vector1d = std::vector<double>;
using vector2d = std::vector<vector1d>;

class IntegrationParam;

class Spherical : public KineticGas {
    public:

    Spherical(std::vector<double> mole_weights, 
                std::vector<std::vector<double>> sigmaij,
                bool is_idealgas);
    virtual ~Spherical(){};

    // Potential models, these must be overridden in derived classes, in addition to model_rdf and get_contact_diameters.
    virtual double potential(int i, int j, double r) = 0;
    virtual double potential_derivative_r(int i, int j, double r) = 0;
    virtual double potential_dblderivative_rr(int i, int j, double r) = 0;

    double omega(int i, int j, int l, int r, double T) override;
    double omega_tester(int i, int j, int l, int r, double T, IntegrationParam& param);
    double w_integral_tester(int i, int j, double T, int l, int r, IntegrationParam& param);

    std::vector<std::vector<double>> model_rdf(double rho, double T, const std::vector<double>& mole_fracs) override = 0;

    vector2d get_collision_diameters(double rho, double T, const vector1d& x) override;
    vector2d get_collision_diameters_model0(double rho, double T, const vector1d& x);
    vector2d get_collision_diameters_model1(double rho, double T, const vector1d& x);

    // Different collision diameter models are
    // -1 : Default model
    // 0 : Model presented in Refs. (I) and (II) (see top of file)
    // 1 : Unpublished model (2024)
    // 2 : Constant (= sigma)
    // 3 : BH-diameters
    void set_collision_diameter_model(int model_id){
        if (model_id == -1) collision_diameter_model_id = default_cd_model_id;
        collision_diameter_model_id = model_id;
    }
    int collision_diameter_model_id = 0;
    const int default_cd_model_id = 0; // Default collision diameter model

    // ------------------------------------------------------------------------------------------------------------------- //
    // -------------------------- Spherical Internals are below here ----------------------------------------------------- //
    // ---------------------- End users should not need to care about anything below -------------------------------------- //
    // ------------------------------------------------------------------------------------------------------------------- //
    double chi(int i, int j, double T, double g, double b);
    double get_R(int i, int j, double T, double g, double b);
    double get_R0(int i, int j, double T, double g); // Solve get_R when b = 0

    double get_b_max_g(double g, double T);
    double get_bmid(double g, double T);
    double momentum_transfer(double T, double g, double b);

    double cd_weight_inner(double T, double g);
    double get_cd_weight_normalizer(double T);
    double get_cd_weight(double T, double g, double b, double I);

    double cd_integrand(double T, double g, double b, double I);
    double cd_inner(double T, double g, double I);
    double momentum_collision_diameter(double T);

    // In the general case, sigma is a scaling parameter
    // On the order of the molecular size. Its specific physical meaning
    // is different for different potential models. It may also be without physical meaning
    std::vector<std::vector<double>> sigma;

    // This method is responsible for calling get_b_max(T, ierr), and handling eventual failures.
    // Most inheritance should only require overriding the failure handling, not the computation in get_b_max(T, ierr)
    virtual std::vector<std::vector<double>> get_b_max(double T);
    std::vector<std::vector<double>> get_b_max(double T, std::vector<std::vector<int>>& ierr);

    double theta(int i, int j, const double T, const double g, const double b);


    private:
    // Helper functions for computing dimentionless collision integrals

    double theta_lim(int i, int j, const double T, const double g);
    double theta_integral(int i, int j, const double T, const double R, const double g, const double b);
    double theta_integrand(int i, int j, double T, double r, double g, double b);
    double transformed_theta_integrand(int i, int j, double T, double u, double R, double g, double b);
    double theta_integrand_dblderivative(int i, int j, double T, double r, double g, double b);
    double get_R_rootfunc(int i, int j, double T, double g, double b, double& r);
    double get_R_rootfunc_derivative(int i, int j, double T, double g, double b, double& r);

    double w_integral(int i, int j, double T, int l, int r); // Dimentionless collision integral for spherical potentials
    double w_integrand(int i, int j, double T, double g, double b, int l, int r);
    std::function<double(int, int, double, double, double, int, int)> w_integrand_export; // Will bind w_integrand to this function such that it can be passed to the external integration module

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