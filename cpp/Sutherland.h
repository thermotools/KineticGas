#pragma once
#include "MieKinGas.h"

class Sutherland : public MieKinGas{

    std::vector<std::vector<std::vector<double>>> coeff, lambda
    std::vector<std::vector<double>> sigma_p, eps_p
    int nterms;

    Sutherland(std::vector<double> mole_weights, std::vector<std::vector<double>> sigmaij,
               std::vector<std::vector<double>> epsilonij, std::vector<std::vector<std::vector<double>>> lambda,
               std::vector<std::vector<std::vector<double>>> coeff, bool is_idealgas)
               : coeff{coeff}, lambda{lambda}, sigma_p{sigmaij}, eps_p{epsilonij}, mole_weights{mole_weights},
                    is_idealgas{is_idealgas}
               {
                    nterms = coeff[0][0].size();
                    set_sigma_eff();
                    set_epsilon_eff();
                    set_alpha();
               }

    void set_sigma_eff();
    std::vector<std::vector<double>> get_sigma_eff();
    void set_epsilon_eff();
    std::vector<std::vector<double>> get_epsilon_eff();
    void set_alpha();

    double potential(int i, int j, double r) override;
    double potential_derivative_r(int i, int j, double r) override;
    double potential_dblderivative_rr(int i, int j, double r) override;
}