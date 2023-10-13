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
#include "MieKinGas.h"
#include "global_params.h"

class QuantumMie : public MieKinGas{
    public:
    const int FH_order; // Feynman-Hibbs correction order

    QuantumMie(std::vector<double> mole_weights,
        std::vector<std::vector<double>> sigmaij,
        std::vector<std::vector<double>> eps,
        std::vector<std::vector<double>> la,
        std::vector<std::vector<double>> lr,
        int FH_order,
        bool is_idealgas)
        : MieKinGas(mole_weights, sigmaij, eps, la, lr, is_idealgas),
        FH_order{FH_order}
        {}

    inline double Q1(double lambda){
        if (FH_order > 0) return lambda * (lambda - 1);
        return 0;
    }

    inline double Q2(double lambda){
        if (FH_order > 1) return 0.5 * (lambda + 2) * (lambda + 1) * Q1(lambda);
        return 0;
    }

    inline double D(const int& i, const int& j, double T){
        double mu = 1.0 / ( (1.0 / m[i]) + (1.0 / m[j]));
        return pow(HBAR, 2) / (12.0 * mu * BOLTZMANN * T);
    }

    inline void q_mie_error(){throw std::runtime_error("QuantumMie uses a temperature dependent potential, but T was not supplied!");}

    double potential(int i, int j, double r) override {q_mie_error(); return 0.;}
    inline double potential(int i, int j, double r, double T){
        return C[i][j] * eps[i][j]
                * ( pow(sigma[i][j] / r, lr[i][j]) - pow(sigma[i][j] / r, la[i][j])
                    + D(i, j, T) * (Q1(lr[i][j]) * pow(sigma[i][j], lr[i][j]) / pow(r, lr[i][j] + 2)
                                    - Q1(la[i][j]) * pow(sigma[i][j], la[i][j]) / pow(r, la[i][j] + 2))
                    + pow(D(i, j, T), 2) * (Q2(lr[i][j]) * pow(sigma[i][j], lr[i][j]) / pow(r, lr[i][j] + 4)
                                            - Q2(la[i][j]) * pow(sigma[i][j], la[i][j]) / pow(r, la[i][j] + 4))
                    );
    }

    double potential_derivative_r(int i, int j, double r) override {q_mie_error(); return 0.;}
    inline double potential_derivative_r(int i, int j, double r, double T){
        return C[i][j] * eps[i][j] * ((la[i][j] * pow(sigma[i][j], la[i][j]) / pow(r, la[i][j] + 1))
                                        - (lr[i][j] * pow(sigma[i][j], lr[i][j]) / pow(r, lr[i][j] + 1))
                                        + D(i, j, T) * (Q1(lr[i][j]) * (-(lr[i][j] + 2)) * pow(sigma[i][j], lr[i][j]) / pow(r, lr[i][j] + 3)
                                                         - Q1(la[i][j]) * (-(la[i][j] + 2)) * pow(sigma[i][j], la[i][j]) / pow(r, la[i][j] + 3))
                                        + pow(D(i, j, T), 2) * (Q2(lr[i][j]) * (-(lr[i][j] + 4)) * pow(sigma[i][j], lr[i][j]) / pow(r, lr[i][j] + 5)
                                                            - Q2(la[i][j]) * (-(la[i][j] + 4)) * pow(sigma[i][j], la[i][j]) / pow(r, la[i][j] + 5))
                                        );
    }

    double potential_dblderivative_rr(int i, int j, double r) override {q_mie_error(); return 0.;}
    inline double potential_dblderivative_rr(int i, int j, double r, double T){
        return C[i][j] * eps[i][j] * ((lr[i][j] * (lr[i][j] + 1) * pow(sigma[i][j], lr[i][j]) / pow(r, lr[i][j] + 2))
                                    - (la[i][j] * (la[i][j] + 1) * pow(sigma[i][j], la[i][j]) / pow(r, la[i][j] + 2))
                                    + D(i, j, T) * (Q1(lr[i][j]) * (lr[i][j] + 2) * (lr[i][j] + 3) * pow(sigma[i][j], lr[i][j]) / pow(r, lr[i][j] + 4)
                                                    - Q1(la[i][j]) * (la[i][j] + 2) * (lr[i][j] + 3) * pow(sigma[i][j], la[i][j]) / pow(r, la[i][j] + 4))
                                    + pow(D(i, j, T), 2) * (Q2(lr[i][j]) * (lr[i][j] + 4) * (lr[i][j] + 5) * pow(sigma[i][j], lr[i][j]) / pow(r, lr[i][j] + 6)
                                                            - Q2(la[i][j]) * (la[i][j] + 4) * (lr[i][j] + 5) * pow(sigma[i][j], la[i][j]) / pow(r, la[i][j] + 6))
                                    );
    }

    double get_R_rootfunc(int i, int j, double T, double g, double b, double& r) override {
        return (potential(i, j, r, T) / (BOLTZMANN * T * pow(g, 2))) + pow(b / r, 2) - 1;
    }

    double get_R_rootfunc_derivative(int i, int j, double T, double g, double b, double& r) override {
        return (potential_derivative_r(i, j, r, T) / (BOLTZMANN * T * pow(g, 2))) - 2 * pow(b, 2) / pow(r, 3);
    }

    double theta_integrand(int i, int j, double T, double r, double g, double b) override {
        return pow((pow(r, 4) / pow(b, 2)) * (1.0 - potential(i, j, r, T) / (BOLTZMANN * T * pow(g, 2))) - pow(r, 2), -0.5);
    }

    double theta_integrand_dblderivative(int i, int j, double T, double r, double g, double b) override {
        // Expressing the integrand as f = (core)^{-1/2}
        const double a = 1.0 / (pow(b, 2) * BOLTZMANN * T * pow(g, 2));
        const double u = potential(i, j, r, T);
        const double u_prime = potential_derivative_r(i, j, r, T);
        const double u_dblprime = potential_dblderivative_rr(i, j, r, T);
        const double core = pow(r, 4) / pow(b, 2) - a * pow(r, 4) * u - pow(r, 2);
        const double core_prime = 4 * pow(r, 3) / pow(b, 2) - a * (4 * pow(r, 3) * u + pow(r, 4) * u_prime) - 2 * r;
        const double core_dblprime = 12.0 * pow(r, 2) / pow(b, 2) - a * (12 * pow(r, 2) * u + 8 * pow(r, 3) * u_prime + pow(r, 4) * u_dblprime) - 2;

        double val = (3.0 / 4.0) * pow(core, -2.5) * pow(core_prime, 2) - 0.5 * pow(core, - 1.5) * core_dblprime;
        #ifdef DEBUG
            if (val < 0){
                std::printf("\nd3tdr3 at r = %E sigma\n", r / sigma[i][j]);
                std::printf("val = %E\n\n", val);
            }
        #endif
        return val;
    }

    std::vector<std::vector<double>> get_sigma_eff(double T){
        // Newton solver to find sigma_eff by solving potential(sigma_eff) = 0
        std::vector<std::vector<double>> sigma_eff(Ncomps, std::vector<double>(Ncomps, 0));
        for (int i = 0; i < Ncomps; i++){
            for (int j = i; j < Ncomps; j++){
                sigma_eff[i][j] = sigma[i][j];
                while ( abs(potential(i, j, sigma_eff[i][j], T) / (C[i][j] * eps[i][j])) > 1e-10 ){
                    sigma_eff[i][j] -= potential(i, j, sigma_eff[i][j], T) / potential_derivative_r(i, j, sigma_eff[i][j], T);
                }
                sigma_eff[j][i] = sigma_eff[i][j];
            }
        }
        return sigma_eff;
    }

    std::vector<std::vector<double>> get_BH_diameters(double T) override {
        // Gauss-Legendre points taken from SAFT-VR-MIE docs (see: ThermoPack)
        std::vector<std::vector<double>> d_BH(Ncomps, std::vector<double>(Ncomps, 0.0));
        std::vector<std::vector<double>> sigma_eff = get_sigma_eff(T);
        double beta = 1. / (BOLTZMANN * T);
        for (int i = 0; i < Ncomps; i++){
            for (int j = i; j < Ncomps; j++){
                for (int n = 0; n < 10; n++){
                    d_BH[i][j] += mie_rdf_constants::gl_w[n] * (1. - exp(- beta * potential(i, j, sigma_eff[i][j] * (mie_rdf_constants::gl_x[n] + 1) / 2., T)));
                }
                d_BH[i][j] *= sigma_eff[i][j] / 2;
                d_BH[j][i] = d_BH[i][j];
            }
        }
        return d_BH;
    }

};