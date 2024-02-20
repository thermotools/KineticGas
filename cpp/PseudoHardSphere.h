/*
Author: Vegard Gjeldvik Jervell
Contains: The PseudoHardSphere class. This class implements a potential with continious first and second derivatives
            That rises steeply at r = sigma, and goes to zero at r > sigma.
*/
#pragma once
#include "Spherical.h"

class PseudoHardSphere : public Spherical {
    public:

    PseudoHardSphere(std::vector<double> mole_weights,
        std::vector<std::vector<double>> sigmaij,
        bool is_idealgas)
        : Spherical(mole_weights, sigmaij, is_idealgas) {}
    
    double potential(int i, int j, double r) override {
        // To get this, start with a potential that has a second derivative f''(r) = (sigma / r)^22 + A
        // Then integrate the function and require that f''(sigma) = f'(sigma) = f(sigma) = 0
        if (r > sigma[i][j]){
            return 0.0;
        }
        return (pow(sigma[i][j] / (r), 20) - (20.0 * 21.0 / 2) * pow(r / sigma[i][j], 2) + 20.0 * 22.0 * (r / sigma[i][j]) + 20.0 * ((21.0 / 2.0) - 22.0) - 1.0) / BOLTZMANN; // Force continiuous function
    }
    
    double potential_derivative_r(int i, int j, double r) override {
        if (r > sigma[i][j]){
            return 0.0;
        }
        return (- 20.0 * pow(sigma[i][j], 20) / pow((r), 21) - 20.0 * 21.0 * (r) / pow(sigma[i][j], 2) + 20.0 * 22.0 / sigma[i][j]) / BOLTZMANN; // Force continiuous first derivative
    }
    
    double potential_dblderivative_rr(int i, int j, double r) override {
        if (r > sigma[i][j]){
            return 0.0;
        }
        return (20.0 * 21.0 * pow(sigma[i][j], 20) / pow((r), 22) - 20.0 * 21.0 / pow(sigma[i][j], 2) ) / BOLTZMANN; // Force continiuous second derivative
    }

    std::vector<std::vector<double>> model_rdf(double rho, double T, const std::vector<double>& xi) override {
        std::vector<double> Zi(3);
        for (int i = 1; i < 4; i++){
            for (int j = 0; j < Ncomps; j++){
                Zi[i - 1] += rho * xi[j] * pow(sigma[j][j], i);
            }
            Zi[i - 1] *= (PI / 6.);
        }
        double Z = 1 - Zi[2];
        std::vector<std::vector<double>> rdf(Ncomps, std::vector<double>(Ncomps));
        for (int i = 0; i < Ncomps; i++){
            for (int j = i; j < Ncomps; j++){
                rdf[i][j] = (pow(Z, 2) 
                            + (3. * sigma[i][i] * sigma[j][j] / (sigma[i][i] + sigma[j][j])) * Z * Zi[1] 
                            + 2. * pow(sigma[i][i] * sigma[j][j] / (sigma[i][i] + sigma[j][j]), 2) * pow(Zi[1], 2) 
                            ) / pow(Z, 3);
                rdf[j][i] = rdf[i][j];
            }
        }
        return rdf;
    }

    std::vector<std::vector<double>> get_collision_diameters(double rho, double T, const std::vector<double>& x) override {return sigma;}
};