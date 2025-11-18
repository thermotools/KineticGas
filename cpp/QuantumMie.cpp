#include "QuantumMie.h"
#include <cppThermopack/saftvrqmie.h>
#include <iostream>

QuantumMie::QuantumMie(vector1d mole_weights, vector2d sigma, vector2d eps, vector2d la, vector2d lr, std::vector<int> FH_order,
                        bool is_idealgas, bool is_singlecomp)
        : ExtSutherland(mole_weights, sigma, eps, 6, is_idealgas, is_singlecomp), FH_order{FH_order}
{
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = i; j < Ncomps; j++){
            C[0][i][j] = C[0][j][i] = (lr[i][j] / (lr[i][j] - la[i][j])) * pow(lr[i][j] / la[i][j], (la[i][j] / (lr[i][j] - la[i][j])));
            C[1][i][j] = C[1][j][i] = - C[0][i][j];

            lambda[0][i][j] = lambda[0][j][i] = lr[i][j];
            lambda[1][i][j] = lambda[1][j][i] = la[i][j];
            lambda[2][i][j] = lambda[2][j][i] = lr[i][j] + 2.;
            lambda[3][i][j] = lambda[3][j][i] = la[i][j] + 2.;
            lambda[4][i][j] = lambda[4][j][i] = lr[i][j] + 4.;
            lambda[5][i][j] = lambda[5][j][i] = la[i][j] + 4.;

            beta_exp[2][i][j] = beta_exp[2][j][i] = beta_exp[3][i][j] = beta_exp[3][j][i] = 1;
            beta_exp[4][i][j] = beta_exp[4][j][i] = beta_exp[5][i][j] = beta_exp[5][j][i] = 2;

            double mu1_inv{0.}, mu2_inv{0.};
            if (FH_order[i] > 0) mu1_inv += 1. / m[i];
            if (FH_order[j] > 0) mu1_inv += 1. / m[j];
            if (FH_order[i] > 1) mu2_inv += 1. / m[i];
            if (FH_order[j] > 1) mu2_inv += 1. / m[j];
            double D1_factor = pow(HBAR, 2) * mu1_inv / 24.0;
            double D2_factor = pow(HBAR, 2) * mu2_inv / 24.0;

            C[2][i][j] = C[2][j][i] = C[0][i][j] * D1_factor * Q1(i, j, lr) * pow(sigma[i][j], -2) * pow(eps[i][j], -1);
            C[3][i][j] = C[3][j][i] = C[1][i][j] * D1_factor * Q1(i, j, la) * pow(sigma[i][j], -2) * pow(eps[i][j], -1);
            C[4][i][j] = C[4][j][i] = C[0][i][j] * pow(D2_factor, 2) * Q2(i, j, lr) * pow(sigma[i][j], -4) * pow(eps[i][j], -2);
            C[5][i][j] = C[5][j][i] = C[1][i][j] * pow(D2_factor, 2) * Q2(i, j, la) * pow(sigma[i][j], -4) * pow(eps[i][j], -2);
        }
    }
    int max_FH_order = *(std::max_element(FH_order.begin(), FH_order.end()));
    for (double fh : FH_order){
        if (fh != max_FH_order) throw std::runtime_error("All FH orders must be equal!");
    }

    if (max_FH_order == 0){
        nterms = 2;
    }
    else if (max_FH_order == 1){
        nterms = 4;
    }
}

QuantumMie::QuantumMie(std::string comps, int FH_order_, bool is_idealgas)
    : ExtSutherland(comps, 2 * (FH_order_ + 1), is_idealgas), FH_order(Ncomps, FH_order_)
{   
    std::string potential_key = "Mie-FH" + std::to_string(FH_order_);
    for (size_t i = 0; i < Ncomps; i++){
        const auto pdata = compdata[i][potential_key]["default"];
        sigma[i][i] = pdata["sigma"];
        eps[i][i] = static_cast<double>(pdata["eps_div_k"]) * BOLTZMANN;
        lambda[0][i][i] = pdata["lambda_r"];
        lambda[1][i][i] = pdata["lambda_a"];
    }
    mix_sigma();
    mix_epsilon();
    mix_exponents(lambda[0]); mix_exponents(lambda[1]);
    init_FH_terms();
    Saftvrqmie svrqm(comps, FH_order_);
    for (size_t i = 0; i < ((is_singlecomp) ? 1 : Ncomps); i++){
        svrqm.set_pure_fluid_param(i + 1, 1., sigma[i][i], eps[i][i] / BOLTZMANN, lambda[1][i][i], lambda[0][i][i]);
    }
    eos = std::make_unique<GenericEoS>(ThermoWrapper(std::move(svrqm)));
}

void QuantumMie::init_FH_terms(){
    vector2d lr = lambda[0]; vector2d la = lambda[1];
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = i; j < Ncomps; j++){
            C[0][i][j] = C[0][j][i] = (lr[i][j] / (lr[i][j] - la[i][j])) * pow(lr[i][j] / la[i][j], (la[i][j] / (lr[i][j] - la[i][j])));
            C[1][i][j] = C[1][j][i] = - C[0][i][j];

            if (nterms > 2){
                lambda[2][i][j] = lambda[2][j][i] = lr[i][j] + 2.;
                lambda[3][i][j] = lambda[3][j][i] = la[i][j] + 2.;
            }
            if (nterms > 4){
                lambda[4][i][j] = lambda[4][j][i] = lr[i][j] + 4.;
                lambda[5][i][j] = lambda[5][j][i] = la[i][j] + 4.;
            }

            if (nterms > 2) beta_exp[2][i][j] = beta_exp[2][j][i] = beta_exp[3][i][j] = beta_exp[3][j][i] = 1;
            if (nterms > 4) beta_exp[4][i][j] = beta_exp[4][j][i] = beta_exp[5][i][j] = beta_exp[5][j][i] = 2;

            double mu1_inv{0.}, mu2_inv{0.};
            if (FH_order[i] > 0) mu1_inv += 1. / m[i];
            if (FH_order[j] > 0) mu1_inv += 1. / m[j];
            if (FH_order[i] > 1) mu2_inv += 1. / m[i];
            if (FH_order[j] > 1) mu2_inv += 1. / m[j];
            double D1_factor = pow(HBAR, 2) * mu1_inv / 24.0;
            double D2_factor = pow(HBAR, 2) * mu2_inv / 24.0;

            if (nterms > 2){
                C[2][i][j] = C[2][j][i] = C[0][i][j] * D1_factor * Q1(i, j, lr) * pow(sigma[i][j], -2) * pow(eps[i][j], -1);
                C[3][i][j] = C[3][j][i] = C[1][i][j] * D1_factor * Q1(i, j, la) * pow(sigma[i][j], -2) * pow(eps[i][j], -1);
            }
            if (nterms > 4){
                C[4][i][j] = C[4][j][i] = C[0][i][j] * pow(D2_factor, 2) * Q2(i, j, lr) * pow(sigma[i][j], -4) * pow(eps[i][j], -2);
                C[5][i][j] = C[5][j][i] = C[1][i][j] * pow(D2_factor, 2) * Q2(i, j, la) * pow(sigma[i][j], -4) * pow(eps[i][j], -2);
            }
        }
    }
}