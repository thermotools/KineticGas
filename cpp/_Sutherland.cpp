#include "Sutherland.h"

void Sutherland::compute_sigma_eff(){
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            sigma[i][j] = sigma_p[i][j];
            while (abs(potential(i, j, sigma[i][j])) / eps_p[i][j] > 1e-6){
                sigma[i][j] -= potential(i, j, sigma[i][j]) / potential_derivative_r(i, j, sigma[i][j]);
            }
            sigma[j][i] = sigma[i][j];
        }
    }
}


void Sutherland::compute_epsilon_eff(){
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            sigma_min = sigma_eff[i][j];
            while (abs(potential_derivative_r(i, j, sigma[i][j])) / eps_p[i][j] > 1e-6){
                sigma_min -= potential_derivative_r(i, j, sigma_min) / potential_dblderivative_rr(i, j, sigma_min);
            }
            eps[i][j] = eps[j][i] = - potential(i, j, sigma_min);
        }
    }
}


void Sutherland::set_alpha(){
    alpha = std::vector<std::vector<double>>(Ncomps, std::vector<double>(Ncomps, 0.));
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            for (int n = 0; n < nterms; n++){
                alpha[i][j] -= (1 / eps[i][j]) * (1. / (lambda[i][j][n] - 3.)) * pow(1 / sigma[i][j], lambda[i][j][n]);
            }
            alpha[j][i] = alpha[i][j];
        }
    }
}

double Sutherland::potential(int i, int j, double r) override {
    double p{0.0};
    for (int n = 0; n < nterms; n++){
        p += coeff[i][j][n] * pow(1 / sigma_p[i][j], lambda[i][j][n]);
    }
    return eps_p[i][j] * p;
}

double Sutherland::potential_derivative_r(int i, int j, double r) override {
    double p{0.0};
    for (int n = 0; n < nterms; n++){
        p -= coeff[i][j][n] * lambda[i][j][n] * pow(1. / sigma_p[i][j], lambda[i][j][n] + 1.);
    }
    return eps_p[i][j] * p;
}

double Sutherland::potential_dblderivative_rr(int i, int j, double r) override {
    double p{0.0};
    for (int n = 0; n < nterms; n++){
        p -= coeff[i][j][n] * lambda[i][j][n] * (lambda[i][j][n] + 1.) * pow(1. / sigma_p[i][j], lambda[i][j][n] + 2.);
    }
    return eps_p[i][j] * p;
}