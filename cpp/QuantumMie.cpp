#include "QuantumMie.h"
#include <iostream>

void QuantumMie::set_sigma_eff(double T){
    // Newton solver to find sigma_eff by solving potential(sigma_eff) = 0
    // Uses previously computed sigma_eff as initial guess, to improve runtime upon consecutive calls at the same temp.
    // At init, sigma_eff = sigma, so the first call uses the potential parameter as the initial guess.
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            while (potential_derivative_r(i, j, sigma[i][j], T) > 0){
                sigma[i][j] *= 0.9;
            }
            while (abs(potential(i, j, sigma[i][j], T)) / eps_0[i][j] > 1e-10){
                sigma[i][j] -= potential(i, j, sigma[i][j], T) / potential_derivative_r(i, j, sigma[i][j], T);
            }
            sigma[j][i] = sigma[i][j];

            // After converging the first root, make new initial guesses for the rest.
            if ((i == 0) && (j == 0)){
                for (int ii = 0; ii < Ncomps; ii++){
                    for (int jj = ii; jj < Ncomps; jj++){
                        sigma[ii][jj] = sigma_0[ii][jj] * sigma[0][0] / sigma_0[0][0];
                    }
                }
            }
        }
    }
}

void QuantumMie::set_epsilon_eff(double T){
    double step;
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            while (potential_dblderivative_rr(i, j, sigma_min[i][j], T) < 0){
                sigma_min[i][j] *= 0.9;
            }

            step = potential_derivative_r(i, j, sigma_min[i][j], T) / potential_dblderivative_rr(i, j, sigma_min[i][j], T);
            while ((abs(potential_derivative_r(i, j, sigma_min[i][j], T)) / eps_0[i][j] > 1e-5) && (abs(step) / sigma_0[i][j] > 1e-5)){
                sigma_min[i][j] -= step;
                step = potential_derivative_r(i, j, sigma_min[i][j], T) / potential_dblderivative_rr(i, j, sigma_min[i][j], T);
            }
            sigma_min[j][i] = sigma_min[i][j];
            eps[i][j] = eps[j][i] = - potential(i, j, sigma_min[i][j], T);

            // After converging the first position, make new initial guesses for the rest.
            if ((i == 0) && (j == 0)){
                for (int ii = 0; ii < Ncomps; ii++){
                    for (int jj = ii; jj < Ncomps; jj++){
                        sigma_min[ii][jj] = sigma_0[ii][jj] * sigma_min[0][0] / sigma_0[0][0];
                    }
                }
            }
        }
    }
}

std::vector<std::vector<double>> QuantumMie::model_rdf(double rho, double T, const std::vector<double>& x){
    set_eff_sigma_eps(T);
    return MieKinGas::model_rdf(rho, T, x);
}

std::vector<std::vector<double>> QuantumMie::get_BH_diameters(double T){
    std::cout << "Calling barker henderson!" << std::endl;
    set_eff_sigma_eps(T);
    return MieKinGas::get_BH_diameters(T);
}

std::vector<std::vector<double>> QuantumMie::get_contact_diameters(double rho, double T, const std::vector<double>& x){
    set_eff_sigma_eps(T);
    return MieKinGas::get_contact_diameters(rho, T, x);
}

double QuantumMie::theta_integrand_dblderivative(int i, int j, double T, double r, double g, double b){
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