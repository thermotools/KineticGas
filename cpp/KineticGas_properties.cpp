#include "KineticGas.h"
#include <Eigen/Dense>
#include <json/json.hpp>
#include <iostream>

using vector1d = std::vector<double>;
using vector2d = std::vector<vector1d>;
using vector3d = std::vector<vector2d>;

Eigen::MatrixXd stdvector_to_EigenMatrix(const std::vector<std::vector<double>>& in_matr){
    Eigen::MatrixXd out_matr(in_matr.size(), in_matr[0].size());
    for (int i = 0; i < in_matr.size(); i++){
        for (int j = 0; j < in_matr[i].size(); j++){
            out_matr(i, j) = in_matr[i][j];
        }
    }
    return out_matr;
}

Eigen::VectorXd stdvector_to_EigenVector(const std::vector<double>& in_vec){
    Eigen::VectorXd out_vec(in_vec.size());
    for (int i = 0; i < in_vec.size(); i++){
        out_vec(i) = in_vec[i];
    }
    return out_vec;
}

double KineticGas::viscosity(double T, double Vm, const std::vector<double>& x, int N){
    double rho{AVOGADRO / Vm};

    Eigen::MatrixXd visc_matr{stdvector_to_EigenMatrix(get_viscosity_matrix(rho, T, x, N))};
    Eigen::VectorXd visc_vec{stdvector_to_EigenVector(get_viscosity_vector(rho, T, x, N))};

    Eigen::VectorXd expansion_coeff{visc_matr.partialPivLu().solve(visc_vec)};

    std::vector<std::vector<double>> cd{get_collision_diameters(rho, T, x)};
    std::vector<std::vector<double>> rdf{get_rdf(rho, T, x)};
    std::vector<double> K_prime{get_K_prime_factors(rho, T, x)};

    double eta_prime{0.0};
    for (int i = 0; i < Ncomps; i++){
        eta_prime += K_prime[i] * x[i] * expansion_coeff(i);
    }
    eta_prime *= BOLTZMANN * T / 2.0;

    if (is_idealgas) return eta_prime; // The second term is only included if (is_idealgas == True)

    double eta_dblprime{0.0};
    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j < Ncomps; j++){
            eta_dblprime += sqrt(m[i] * m[j] / (m[i] + m[j])) * x[i] * x[j] * pow(cd[i][j], 4) * rdf[i][j];
        }
    }
    eta_dblprime *= 4.0 * pow(rho, 2) * sqrt(2.0 * PI * BOLTZMANN * T) / 15.0;

    return eta_prime + eta_dblprime;

}

double KineticGas::thermal_conductivity(double T, double Vm, const vector1d& x, int N){
        double rho = AVOGADRO / Vm;
        
        Eigen::MatrixXd cond_matr{stdvector_to_EigenMatrix(get_conductivity_matrix(rho, T, x, N))};
        Eigen::VectorXd cond_vec{stdvector_to_EigenVector(get_conductivity_vector(rho, T, x, N))};
        Eigen::VectorXd th_expansion_coeff{cond_matr.partialPivLu().solve(cond_vec)};

        vector2d rdf = get_rdf(rho, T, x);
        vector1d K = get_K_factors(rho, T, x);
        vector2d cd = get_collision_diameters(rho, T, x);

        vector3d diff_expansion_coeff = reshape_diffusive_expansion_vector(compute_diffusive_expansion_coeff(rho, T, x, N));

        double lambda_dblprime = 0.;
        if (is_idealgas){ // lambda_dblprime is only nonzero when density corrections are present, and vanishes at infinite dilution
            for (size_t i = 0; i < Ncomps; i++){
                for (size_t j = 0; j < Ncomps; j++){
                    lambda_dblprime += pow(rho, 2) * sqrt(2 * PI * m[i] * m[j] * BOLTZMANN * T / (m[i] + m[j])) 
                                        * (x[i] * x[j]) / (m[i] + m[j]) * pow(cd[i][j], 4) * rdf[i][j];
                }
            }
            lambda_dblprime *= (4. * BOLTZMANN / 3.);
        }

        double lambda_prime = 0.;
        Eigen::VectorXd dth{compute_dth_vector(diff_expansion_coeff, th_expansion_coeff)};
        for (size_t i = 0; i < Ncomps; i++){
            double tmp = 0.;
            for (size_t k = 0; k < Ncomps; k++){
                tmp += diff_expansion_coeff[1][i][k] * dth(k);
            }
            lambda_prime += x[i] * K[i] * (th_expansion_coeff(Ncomps + i) - tmp);
        }
        lambda_prime *= (5. * BOLTZMANN / 4.);

        double lamba_internal = 0.;
        double lamb_int_f = 1.32e3;
        double eta0 = viscosity(T, 1e6, x, N);
        double avg_mol_weight = 0.;
        for (size_t i = 0; i < Ncomps; i++) avg_mol_weight += x[i] * m[i];
        avg_mol_weight *= AVOGADRO * 1e3;
        double Cp_id = (is_singlecomp) ? eos->idealenthalpysingle(T, 1, true).dt() : 0.;
        if (!is_singlecomp){
            for (size_t i = 0; i < Ncomps; i++){
                Cp_id += x[i] * eos->idealenthalpysingle(T, i + 1, true).dt();
                std::cout << "Cp(" << i << ") : " <<  eos->idealenthalpysingle(T, i + 1, true) << std::endl;
            }
        }
        double Cp_factor = (Cp_id - 5. * GAS_CONSTANT / 2.) / avg_mol_weight;
        std::cout << "Cp_factor : " << Cp_factor << std::endl;
        lamba_internal = lamb_int_f * eta0 * Cp_factor;

        const double cond = lambda_prime + lambda_dblprime + lamba_internal;
        return cond;
}

Eigen::MatrixXd KineticGas::CoM_to_FoR_matr(double T, double Vm, const std::vector<double>& x, int frame_of_reference, int solvent_idx=-1){
    switch (frame_of_reference) {
    case FrameOfReference::CoM:
        return Eigen::MatrixXd::Identity(Ncomps, Ncomps);;
    case FrameOfReference::CoN:
        return CoM_to_CoN_matr(T, Vm, x);
    case FrameOfReference::CoV:
        return CoM_to_CoV_matr(T, Vm, x);
    case FrameOfReference::solvent:
        return CoM_to_solvent_matr(T, Vm, x, solvent_idx);
    default:
        throw std::runtime_error("Unknown (not-implemented) Frame of Reference.");
    }
}

Eigen::MatrixXd KineticGas::CoM_to_CoN_matr(double T, double Vm, const std::vector<double>& x){
    Eigen::MatrixXd matr = Eigen::MatrixXd::Identity(Ncomps, Ncomps);
    std::vector<double> wt_fracs = get_wt_fracs(x);
    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j < Ncomps; j++){
            matr(i, j) += x[i] * ((wt_fracs[j] / x[j]) - 1);
        }
    }
    return matr;
}

Eigen::MatrixXd KineticGas::CoM_to_solvent_matr(double T, double Vm, const std::vector<double>& x, int solvent_idx){
    Eigen::MatrixXd matr = Eigen::MatrixXd::Identity(Ncomps, Ncomps);
    std::vector<double> wt_fracs = get_wt_fracs(x);
    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j < Ncomps; j++){
            matr(i, j) += x[i] * ((wt_fracs[j] / x[j]) - (delta(j, solvent_idx) / x[solvent_idx]));
        }
    }
    return matr;
}
Eigen::MatrixXd KineticGas::CoM_to_CoV_matr(double T, double Vm, const std::vector<double>& x){
    double p = eos->pressure_tv(T, Vm, x);
    std::vector<double> dvdn = eos->specific_volume(T, p, x, eos->VAPPH, false, false, true).dn();
    Eigen::MatrixXd psi = Eigen::MatrixXd::Identity(Ncomps, Ncomps);
    std::vector<double> wt_frac = get_wt_fracs(x);
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = 0; j < Ncomps; j++){
            psi(i, j) += x[i] * ((wt_frac[j] / x[j]) - (dvdn[j] / Vm));
        }
    }
    return psi;
}

Eigen::VectorXd KineticGas::compute_diffusive_expansion_coeff(double rho, double T, const vector1d& x, int N){
    Eigen::MatrixXd diff_matr{stdvector_to_EigenMatrix(get_diffusion_matrix(rho, T, x, N))};
    Eigen::VectorXd diff_vec{stdvector_to_EigenVector(get_diffusion_vector(rho, T, x, N))};
    return diff_matr.partialPivLu().solve(diff_vec);
}

vector3d KineticGas::reshape_diffusive_expansion_vector(const Eigen::VectorXd& d_ijq){
    unsigned long N{d_ijq.size() / (Ncomps * Ncomps)};
    vector3d d_qij_matr(N, vector2d(Ncomps, vector1d(Ncomps, 0.)));
    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j < Ncomps; j++){
            for (int q = 0; q < N; q++){
                d_qij_matr[q][i][j] = d_ijq(N * Ncomps * j + Ncomps * q + i);
            }
        }
    }
    return d_qij_matr;
}

Eigen::VectorXd KineticGas::compute_dth_vector(const vector3d& d_qij, const Eigen::VectorXd& l){
    Eigen::MatrixXd d_ij(Ncomps, Ncomps);
    Eigen::VectorXd l_i(Ncomps);

    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j < Ncomps; j++){
            d_ij(i, j) = d_qij[0][i][j];
        }
        l_i(i) = l(i);
    }

    return d_ij.partialPivLu().solve(l_i);
}

vector2d KineticGas::get_chemical_potential_factors(double T, double Vm, const std::vector<double>& x){
    vector2d dmudrho(Ncomps, vector1d(Ncomps, 0));
    if (is_singlecomp){
        const vector2d dmudn_pure = eos->chemical_potential_tv(T, Vm, {1.}, false, false, true).dn();
        const double RT = GAS_CONSTANT * T;
        const double dmudrho_pure = Vm * dmudn_pure[0][0];
        const double rho = 1 / Vm;
        dmudrho[0][0] = (dmudrho_pure + RT * x[1] / (rho * x[0])) / AVOGADRO;
        dmudrho[0][1] = dmudrho[1][0] = (dmudrho_pure - RT / rho) / AVOGADRO;
        dmudrho[1][1] = (dmudrho_pure + RT * x[0] / (rho * x[1])) / AVOGADRO;
    }
    else{
        const vector2d dmudn = eos->chemical_potential_tv(T, Vm, x, false, false, true).dn();
        for (size_t i = 0; i < Ncomps; i++){
            for (size_t j = 0; j < Ncomps; j++){
                dmudrho[i][j] = Vm * dmudn[i][j] / AVOGADRO;
            }
        }
    }

    vector2d Eij(Ncomps, vector1d(Ncomps, 0));
    for (size_t i = 0; i < Ncomps; i++){
        double ni = x[i] / Vm;
        for (size_t j = 0; j < Ncomps; j++){
            Eij[i][j] = (ni / (BOLTZMANN * T)) * dmudrho[i][j];
        }
    }
    return Eij;
}