/*
Author: Vegard Gjeldvik Jervell
Contains: The analytical collision integrals and deflection angle for a hard sphere potential.
*/
#pragma once
#include "KineticGas.h"
#include "Factorial.h"
#include "global_params.h"

class HardSphere : public KineticGas {
    public: 

    HardSphere(std::vector<double> mole_weights,
        std::vector<std::vector<double>> sigma,
        bool is_idealgas, bool is_singlecomp)
        : KineticGas(mole_weights, sigma, 
            vector2d(mole_weights.size(), vector1d(mole_weights.size(), 1)),  // Use dummy value for energy parameter
            is_idealgas, is_singlecomp) 
    {}
    
    double omega(int i, int j, int l, int r, double T) override {
        double w = w_integral(i, j, T, l, r); 
        if (i == j) return pow(sigma.at(i).at(j), 2) * sqrt((PI * BOLTZMANN * T) / m.at(i)) * w;
        return 0.5 * pow(sigma.at(i).at(j), 2) * sqrt(2 * PI * BOLTZMANN * T / (m0[i][j] * M[i][j] * M[j][i])) * w;
    }

    double w_integral(int i, int j, double T, int l, int r){
        long long f = Fac(r + 1).eval();
        if (l % 2 == 0){
            return 0.25 * (2 - ((1.0 / (l + 1)) * 2)) * f;
        }
        return 0.5 * f;
    }

    double chi(int i, int j, double T, double g, double b){
        if (b >= sigma[i][j]) return 0;
        return acos(1 - 2 * (1 - pow(b / sigma[i][j], 2)));
    }
    
    std::vector<std::vector<double>> model_rdf(double rho, double T, const std::vector<double>& xi) override {
        std::vector<double> Zi(3, 0.);
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

    std::vector<std::vector<double>> model_mtl(double rho, double T, const std::vector<double>& x) override {return sigma;}
    std::vector<std::vector<double>> model_etl(double rho, double T, const std::vector<double>& x) override {return sigma;}

};