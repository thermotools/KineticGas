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
    std::vector<std::vector<double>> sigma_0, eps_0; // The potential parameters
    // Position of the minimum, saved between computations of eps_eff for better runtime
    // upon consecutive computations of eps_eff at similar temperatures
    std::vector<std::vector<double>> sigma_min;

    /* ------------------------------- TAKE HEED YE WHO COMETH HERE -----------------------------------*/
    // To make inheriting methods from MieKinGas as painless as possible, this class uses the inheritted
    // Attributes sigma and eps to holde the FH-Mie effective parameters, and uses sigma_0 and eps_0
    // to hold the potential parameters. This allows us to override any method from MieKinGas by simply doing
    //
    //  double my_method(double T){
    //      set_sigma_eff(T); set_epsilon_eff(T);
    //      return MieKinGas::my_method(T);}
    //
    // But carries the caveat of having a "Statefull" class. I'm considering making that less dangerous by making
    // A bunch of helper functions protected to prevent them being called without set_sigma_eff(T) being called first.
    /* ---------------------------------- YE HATH NOW TAKEN HEED --------------------------------------*/

    QuantumMie(std::vector<double> mole_weights,
        std::vector<std::vector<double>> sigmaij,
        std::vector<std::vector<double>> eps,
        std::vector<std::vector<double>> la,
        std::vector<std::vector<double>> lr,
        int FH_order,
        bool is_idealgas)
        : MieKinGas(mole_weights, sigmaij, eps, la, lr, is_idealgas),
        FH_order{FH_order}, sigma_0{sigmaij}, eps_0{eps}, sigma_min{sigmaij}
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
    double potential(int i, int j, double r, double T){
        return C[i][j] * eps_0[i][j] * ( pow(sigma_0[i][j] / r, lr[i][j]) - pow(sigma_0[i][j] / r, la[i][j])
                                        + D(i, j, T) * (Q1(lr[i][j]) * pow(sigma_0[i][j], lr[i][j]) / pow(r, lr[i][j] + 2)
                                                    - Q1(la[i][j]) * pow(sigma_0[i][j], la[i][j]) / pow(r, la[i][j] + 2))
                                        + pow(D(i, j, T), 2) * (Q2(lr[i][j]) * pow(sigma_0[i][j], lr[i][j]) / pow(r, lr[i][j] + 4)
                                                            - Q2(la[i][j]) * pow(sigma_0[i][j], la[i][j]) / pow(r, la[i][j] + 4))
                                        );
    }

    double potential_derivative_r(int i, int j, double r) override {q_mie_error(); return 0.;}
    double potential_derivative_r(int i, int j, double r, double T){
        return C[i][j] * eps_0[i][j] * ((la[i][j] * pow(sigma_0[i][j], la[i][j]) / pow(r, la[i][j] + 1))
                                        - (lr[i][j] * pow(sigma_0[i][j], lr[i][j]) / pow(r, lr[i][j] + 1))
                                        + D(i, j, T) * (Q1(lr[i][j]) * (-(lr[i][j] + 2)) * pow(sigma_0[i][j], lr[i][j]) / pow(r, lr[i][j] + 3)
                                                         - Q1(la[i][j]) * (-(la[i][j] + 2)) * pow(sigma_0[i][j], la[i][j]) / pow(r, la[i][j] + 3))
                                        + pow(D(i, j, T), 2) * (Q2(lr[i][j]) * (-(lr[i][j] + 4)) * pow(sigma_0[i][j], lr[i][j]) / pow(r, lr[i][j] + 5)
                                                            - Q2(la[i][j]) * (-(la[i][j] + 4)) * pow(sigma_0[i][j], la[i][j]) / pow(r, la[i][j] + 5))
                                        );
    }

    double potential_dblderivative_rr(int i, int j, double r) override {q_mie_error(); return 0.;}
    double potential_dblderivative_rr(int i, int j, double r, double T){
        return C[i][j] * eps_0[i][j] * (lr[i][j] * (lr[i][j] + 1) * pow(sigma_0[i][j], lr[i][j]) / pow(r, lr[i][j] + 2)
                                        - (la[i][j] * (la[i][j] + 1) * pow(sigma_0[i][j], la[i][j]) / pow(r, la[i][j] + 2))
                                        + D(i, j, T) * (Q1(lr[i][j]) * (lr[i][j] + 2) * (lr[i][j] + 3) * pow(sigma_0[i][j], lr[i][j]) / pow(r, lr[i][j] + 4)
                                                        - Q1(la[i][j]) * (la[i][j] + 2) * (la[i][j] + 3) * pow(sigma_0[i][j], la[i][j]) / pow(r, la[i][j] + 4))
                                        + pow(D(i, j, T), 2) * (Q2(lr[i][j]) * (lr[i][j] + 4) * (lr[i][j] + 5) * pow(sigma_0[i][j], lr[i][j]) / pow(r, lr[i][j] + 6)
                                            - Q2(la[i][j]) * (la[i][j] + 4) * (la[i][j] + 5) * pow(sigma_0[i][j], la[i][j]) / pow(r, la[i][j] + 6)
                                            )
                                        );
    }

    void set_sigma_eff(double T);
    void set_epsilon_eff(double T);
    inline void set_eff_sigma_eps(double T){set_sigma_eff(T); set_epsilon_eff(T);};
    inline std::vector<std::vector<double>> get_sigma_eff(double T){set_sigma_eff(T); return sigma;};
    inline std::vector<std::vector<double>> get_epsilon_eff(double T){set_epsilon_eff(T); return eps;};
    inline std::vector<std::vector<double>> get_sigma_min(double T){set_epsilon_eff(T); return sigma_min;};

    std::vector<std::vector<double>> get_BH_diameters(double T) override;
    virtual std::vector<std::vector<double>> get_contact_diameters(double rho, double T, const std::vector<double>& x) override;
    std::vector<std::vector<double>> model_rdf(double rho, double T, const std::vector<double>& x) override;

    inline double get_R_rootfunc(int i, int j, double T, double g, double b, double& r) override {
        return (potential(i, j, r, T) / (BOLTZMANN * T * pow(g, 2))) + pow(b / r, 2) - 1;
    }
    inline double get_R_rootfunc_derivative(int i, int j, double T, double g, double b, double& r) override {
        return (potential_derivative_r(i, j, r, T) / (BOLTZMANN * T * pow(g, 2))) - 2 * pow(b, 2) / pow(r, 3);
    }
    inline double theta_integrand(int i, int j, double T, double r, double g, double b) override {
        return pow((pow(r, 4) / pow(b, 2)) * (1.0 - potential(i, j, r, T) / (BOLTZMANN * T * pow(g, 2))) - pow(r, 2), -0.5);
    }
    double theta_integrand_dblderivative(int i, int j, double T, double r, double g, double b) override;

    // The following methods must be overridden to first compute the effective parameters before calling
    // the parent method, until then they are private (should not be used).
    private:
    using MieKinGas::get_b_max; // (double T);
    using MieKinGas::rdf_HS; // (double rho, double T, const std::vector<double>& x);
    using MieKinGas::rdf_g1_func; // (double rho, double T, const std::vector<double>& x);
    using MieKinGas::rdf_g2_func; // (double rho, double T, const std::vector<double>& x,
                                  //              const std::vector<std::vector<double>>& d_BH,
                                  //              const std::vector<std::vector<double>>& x0);
    // using MieKinGas::rdf_g2_func; // (double rho, double T, const std::vector<double>& x);
    using MieKinGas::a_1s_func; // (double rho, double T, const std::vector<double>& x,
                                //                const std::vector<std::vector<double>>& lambda);
    using MieKinGas::da1s_drho_func; // (double rho, double T, const std::vector<double>& x,
                                     //           const std::vector<std::vector<double>>& lambda);

    using MieKinGas::a1ij_func; // (double rho, double T, const std::vector<double>& x);
    using MieKinGas::da1ij_drho_func; // (double rho, double T, const std::vector<double>& x);
    using MieKinGas::a2ij_func; // (double rho, double T, const std::vector<double>& x);
    using MieKinGas::da2ij_drho_func; // (double rho, double T, const std::vector<double>& x);

    using MieKinGas::gamma_corr; // (double zeta_x, double T);


};