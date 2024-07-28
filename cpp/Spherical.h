/*
Author: Vegard Gjeldvik Jervell
Contains: The abstract class 'Spherical', inheriting from 'KineticGas'. This class implements
        Methods to evaluate the collision integrals for an arbitrary spherically symmetric potential.

        Subclasses must implement the potential, as well as the first and second derivative wrt. distance.
*/

#pragma once
#include "KineticGas.h"
#include "Integration/Integration.h"
#include "global_params.h"

class IntegrationParam;

class Spherical : public KineticGas {
    public:

    Spherical(std::vector<double> mole_weights, 
                std::vector<std::vector<double>> sigmaij,
                bool is_idealgas);
    
    #ifdef NOPYTHON
        Spherical(std::string comps, bool is_idealgas) : KineticGas(comps, is_idealgas) {}
    #endif
    
    virtual ~Spherical(){};

    // Potential models, these must be overridden in derived classes, in addition to model_rdf and get_contact_diameters.
    virtual double potential(int i, int j, double r) = 0;
    virtual double potential_derivative_r(int i, int j, double r) = 0;
    virtual double potential_dblderivative_rr(int i, int j, double r) = 0;

    double omega(int i, int j, int l, int r, double T) override;
    double omega_tester(int i, int j, int l, int r, double T, IntegrationParam& param);
    double w_integral_tester(int i, int j, double T, int l, int r, IntegrationParam& param);
    std::vector<std::vector<double>> get_collision_diameters(double rho, double T, const std::vector<double>& x) override;
    std::vector<std::vector<double>> model_rdf(double rho, double T, const std::vector<double>& mole_fracs) override = 0;

    // ------------------------------------------------------------------------------------------------------------------- //
    // -------------------------- Spherical Internals are below here ----------------------------------------------------- //
    // ---------------------- End users should not need to care about anything below -------------------------------------- //
    // ------------------------------------------------------------------------------------------------------------------- //
    protected:
    // In the general case, sigma is a scaling parameter
    // On the order of the molecular size. Its specific physical meaning
    // is different for different potential models. It may also be without physical meaning
    std::vector<std::vector<double>> sigma;

    // This method is responsible for calling get_b_max(T, ierr), and handling eventual failures.
    // Most inheritance should only require overriding the failure handling, not the computation in get_b_max(T, ierr)
    virtual std::vector<std::vector<double>> get_b_max(double T);
    std::vector<std::vector<double>> get_b_max(double T, std::vector<std::vector<int>>& ierr);

    double theta(int i, int j, const double T, const double g, const double b);
    double chi(int i, int j, double T, double g, double b);
    double get_R(int i, int j, double T, double g, double b);

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
};

class IntegrationParam{
    public:
    Point origin{1e-7, 1e-7};
    Point end{8, 5};
    double dg{0.5}, db{0.03125};
    int refinement_levels_g{4};
    int refinement_levels_b{16};
    double subdomain_dblder_limit{1e-5};

    void set_dg(double dg_){dg = dg_;}
    void set_db(double db_){db = db_;}
    void set_rlg(int rlg){refinement_levels_g = rlg;}
    void set_rlb(int rlb){refinement_levels_b = rlb;}
    void set_dblder_lim(double lim){subdomain_dblder_limit = lim;}
};