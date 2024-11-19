#include "QuantumMie.h"
#include <iostream>
#include <mutex>
#include <thread>

QuantumMie::QuantumMie(vector1d mole_weights, vector2d sigma, vector2d eps, vector2d la, vector2d lr, std::vector<int> FH_order,
                        bool is_idealgas, bool is_singlecomp)
        : ExtSutherland(mole_weights, sigma, eps, 6, is_idealgas, is_singlecomp), FH_order{FH_order}, Q_factors(4, vector2d(Ncomps, vector1d(Ncomps, 0.)))
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
    if (*(std::max_element(FH_order.begin(), FH_order.end())) == 0){
        nterms = 2;
    }
    else if (*(std::max_element(FH_order.begin(), FH_order.end())) == 1){
        nterms = 4;
    }
}