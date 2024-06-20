#include "KineticGas.h"
#include <Eigen/Dense>
#include <json/json.hpp>
#include <iostream>

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

double KineticGas::viscosity(double T, double Vm, std::vector<double>& x, int N){
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

std::vector<std::vector<std::vector<double>>> KineticGas::reshape_diffusive_expansion_vector(const Eigen::VectorXd& d_ijq){
    unsigned long N{d_ijq.size() / (Ncomps * Ncomps)};
    std::vector<std::vector<std::vector<double>>> d_ijq_matr(Ncomps, std::vector<std::vector<double>>(Ncomps, std::vector<double>(N, 0.)));
    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j < Ncomps; j++){
            for (int q = 0; q < Ncomps; q++){
                d_ijq_matr[i][j][q] = d_ijq(N * Ncomps * j + Ncomps * q + i);
            }
        }
    }
    return d_ijq_matr;
}

Eigen::VectorXd KineticGas::compute_dth_vector(const std::vector<std::vector<std::vector<double>>>& d_ijq, const Eigen::VectorXd& l){
    Eigen::MatrixXd d_ij(Ncomps, Ncomps);
    Eigen::VectorXd l_i(Ncomps);

    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j < Ncomps; j++){
            d_ij(i, j) = d_ijq[i][j][0];
        }
        l_i(i) = l(i);
    }

    return d_ij.partialPivLu().solve(l_i);
}