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

class Spherical : public KineticGas {
    public:
    // In the general case, sigma is a scaling parameter
    // On the order of the molecular size. Its specific physical meaning 
    // is different for different potential models. It may also be without physical meaning
    std::vector<std::vector<double>> sigma;

    Spherical(std::vector<double> mole_weights, 
                std::vector<std::vector<double>> sigmaij,
                bool is_idealgas);
    virtual ~Spherical(){};

    // Potential models, these must be explicitly overridden in derived classes, in addition to model_rdf and get_contact_diameters.
    virtual double potential(int i, int j, double r) = 0;
    virtual double potential_derivative_r(int i, int j, double r) = 0;
    virtual double potential_dblderivative_rr(int i, int j, double r) = 0;
    virtual std::vector<std::vector<double>> model_rdf(double rho, double T, const std::vector<double>& mole_fracs) override = 0;
    virtual std::vector<std::vector<double>> get_contact_diameters(double rho, double T, const std::vector<double>& x) override = 0;

    double omega(const int& i, const int& j, const int& l, const int& r, const double& T) override;
    double w_integral(const int& i, const int& j, const double& T, const int& l, const int& r); // Dimentionless collision integral for spherical potentials
    double w_integrand(const int& i, const int& j, const double& T, 
                            const double& g, const double& b, 
                            const int& l, const int& r);
    std::function<double(int, int, double, double, double, int, int)> w_integrand_export; // Will bind w_integrand to this function such that it can be passed to the external integration module

    // Helper functions for computing dimentionless collision integrals
    double theta(int i, int j, const double T, const double g, const double b);
    double theta_lim(int i, int j, const double T, const double g);
    double theta_integral(int i, int j, const double T, const double R, const double g, const double b);
    virtual double theta_integrand(int i, int j, double T, double r, double g, double b);
    virtual double transformed_theta_integrand(int i, int j, double T, double u, double R, double g, double b);
    virtual double theta_integrand_dblderivative(int i, int j, double T, double r, double g, double b);
    double get_R(int i, int j, double T, double g, double b);
    virtual double get_R_rootfunc(int i, int j, double T, double g, double b, double& r);
    virtual double get_R_rootfunc_derivative(int i, int j, double T, double g, double b, double& r);
    double chi(int i, int j, double T, double g, double b);
};