#include "KineticGas.h"
#include <Eigen/Dense>
#include <nlohmann/json.hpp>
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

    std::vector<std::vector<double>> mtl{get_mtl(rho, T, x)};
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
            eta_dblprime += sqrt(m[i] * m[j] / (m[i] + m[j])) * x[i] * x[j] * pow(mtl[i][j], 4) * rdf[i][j];
        }
    }
    eta_dblprime *= 4.0 * pow(rho, 2) * sqrt(2.0 * PI * BOLTZMANN * T) / 15.0;

    return eta_prime + eta_dblprime;

}

double KineticGas::thermal_conductivity(double T, double Vm, const vector1d& x, int N){
    #ifdef DEBUG
        if (!eos.get()) throw std::runtime_error("EoS is not set (in get_chemical_potential_factors)");
    #endif
    double rho = AVOGADRO / Vm;
    
    Eigen::VectorXd th_expansion_coeff{compute_thermal_expansion_coeff(rho, T, x, N)};

    vector2d rdf = get_rdf(rho, T, x);
    vector1d K = get_K_factors(rho, T, x);
    vector2d etl = get_etl(rho, T, x);

    double lambda_dblprime = 0.;
    if (!is_idealgas){ // lambda_dblprime is only nonzero when density corrections are present, and vanishes at infinite dilution
        for (size_t i = 0; i < Ncomps; i++){
            for (size_t j = 0; j < Ncomps; j++){
                lambda_dblprime += pow(rho, 2) * sqrt(2 * PI * m[i] * m[j] * BOLTZMANN * T / (m[i] + m[j])) 
                                    * (x[i] * x[j]) / (m[i] + m[j]) * pow(etl[i][j], 4) * rdf[i][j];
            }
        }
        lambda_dblprime *= (4. * BOLTZMANN / 3.);
    }

    double lambda_prime = 0.;
    Eigen::VectorXd dth = Eigen::VectorXd::Zero(Ncomps);
    vector3d diff_expansion_coeff(N, vector2d(Ncomps, vector1d(Ncomps, 0.)));
    if (!is_singlecomp){
        diff_expansion_coeff = reshape_diffusive_expansion_vector(compute_diffusive_expansion_coeff(rho, T, x, N));
        dth = compute_dth_vector(diff_expansion_coeff, th_expansion_coeff);
    }
    for (size_t i = 0; i < Ncomps; i++){
        double tmp = 0.;
        if (!is_singlecomp){ // tmp is a Thermal diffusion related term, which is zero for a pure component
            for (size_t k = 0; k < Ncomps; k++){
                tmp += diff_expansion_coeff[1][i][k] * dth(k);
            }
        }
        lambda_prime += x[i] * K[i] * (th_expansion_coeff(Ncomps + i) - tmp);
    }
    lambda_prime *= (5. * BOLTZMANN / 4.);

    double lamba_internal = 0.;
    double lamb_int_f = 1.32e3;
    double eta0 = viscosity(T, 1e6, x, N);
    double avg_mol_weight = 0.;
    for (size_t i = 0; i < Ncomps; i++) {avg_mol_weight += x[i] * m[i];}
    avg_mol_weight *= AVOGADRO * 1e3;
    double Cp_id = (is_singlecomp) ? eos->Cp_ideal(T, 1) : 0.;
    if (!is_singlecomp){
        for (size_t i = 0; i < Ncomps; i++){
            Cp_id += x[i] * eos->Cp_ideal(T, i + 1);
        }
    }
    double Cp_factor = (Cp_id - 5. * GAS_CONSTANT / 2.) / avg_mol_weight;
    lamba_internal = lamb_int_f * eta0 * Cp_factor;

    const double cond = lambda_prime + lambda_dblprime + lamba_internal;
    return cond;
}

Eigen::MatrixXd compress_diffusion_matrix(const Eigen::MatrixXd& D_in, int dependent_idx){
    size_t N = D_in.cols();
    Eigen::MatrixXd D_out(N - 1, N - 1);
    size_t Mi = 0;
    for (size_t i = 0; i < N; i++){
        if (i == dependent_idx) continue;
        size_t Mj = 0;
        for (size_t j = 0; j < N; j++){
            if (j == dependent_idx) continue;
            D_out(Mi, Mj) = D_in(i, j);
            Mj++;
        }
        Mi++;
    }
    return D_out;
}

Eigen::MatrixXd KineticGas::interdiffusion(double T, double Vm, const vector1d& x, int N, int frame_of_reference, 
                                    int dependent_idx, int solvent_idx, bool do_compress){
    if (dependent_idx < 0 && frame_of_reference == FrameOfReference::solvent) dependent_idx = solvent_idx;
    while (dependent_idx < 0) dependent_idx += Ncomps;
    
    if (frame_of_reference == FrameOfReference::zarate_x){
        Eigen::MatrixXd D = interdiffusion(T, Vm, x, N, FrameOfReference::CoN, dependent_idx);
        return compress_diffusion_matrix(D, dependent_idx);
    }
    else if (frame_of_reference == FrameOfReference::zarate){
        Eigen::MatrixXd X = get_zarate_X_matr(x, dependent_idx);
        Eigen::MatrixXd Dx = interdiffusion(T, Vm, x, N, FrameOfReference::zarate_x, dependent_idx);
        return X * Dx * X.inverse();
    }
    else if (frame_of_reference == FrameOfReference::zarate_w){
        Eigen::MatrixXd W = get_zarate_W_matr(x, dependent_idx);
        Eigen::MatrixXd Dz = interdiffusion(T, Vm, x, N, FrameOfReference::zarate, dependent_idx);
        return W * Dz * W.inverse();
    }

    Eigen::MatrixXd psi = CoM_to_FoR_matr(T, Vm, x, frame_of_reference, solvent_idx);
    Eigen::MatrixXd D_dep = interdiffusion_dependent_CoM(T, Vm, x, N);
    D_dep = psi * D_dep;
    
    Eigen::MatrixXd D{D_dep};
    vector1d ksi = get_ksi_factors(T, Vm, x);
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = 0; j < Ncomps; j++){
            D(i, j) -= (ksi[j] / ksi[dependent_idx]) * D_dep(i, dependent_idx);
        }
    }
    return (do_compress) ? compress_diffusion_matrix(D, dependent_idx) : D;
}

Eigen::VectorXd KineticGas::thermal_diffusion_coeff(double T, double Vm, const vector1d& x, int N, 
                                                int frame_of_reference, int dependent_idx, int solvent_idx){
    if (dependent_idx < 0 && frame_of_reference == FrameOfReference::solvent) dependent_idx = solvent_idx;
    while (dependent_idx < 0) dependent_idx += Ncomps;
    
    double rho = AVOGADRO / Vm;
    vector3d d = reshape_diffusive_expansion_vector(compute_diffusive_expansion_coeff(rho, T, x, N));
    Eigen::VectorXd l = compute_thermal_expansion_coeff(rho, T, x, N);
    vector1d ksi = get_ksi_factors(T, Vm, x);
    vector2d rdf = get_rdf(rho, T, x);
    vector2d etl = get_etl(rho, T, x);

    Eigen::MatrixXd b = Eigen::MatrixXd::Identity(Ncomps, Ncomps);
    if (!is_idealgas){
        for (size_t j = 0; j < Ncomps; j++){
            for (size_t k = 0; k < Ncomps; k++){
                b(j, k) += (4. * PI / 3) * rho * x[k] * pow(etl[j][k], 3) * M[j][k] * rdf[j][k]; 
            }
        }
    }
    Eigen::VectorXd DT(Ncomps);
    for (size_t i = 0; i < Ncomps; i++){
        double outer_sum = 0.;
        for (size_t j = 0; j < Ncomps; j++){
            double inner_sum = 0.;
            for (size_t k = 0; k < Ncomps; k++){
                inner_sum += b(j, k);
            }
            outer_sum += d[0][i][j] * x[j] * inner_sum;
        }
        DT(i) = (x[i] / 2.) * (l[i] - outer_sum);
    }

    bool use_zarate{false};
    if (frame_of_reference == FrameOfReference::zarate){
        frame_of_reference = FrameOfReference::CoN;
        use_zarate = true;
    }

    Eigen::MatrixXd psi = CoM_to_FoR_matr(T, Vm, x, frame_of_reference, solvent_idx);
    DT = psi * DT;
    Eigen::MatrixXd Dij = psi * interdiffusion_dependent_CoM(T, Vm, x, N);

    for (size_t i = 0; i < Ncomps; i++){
        double tmp = 0.;
        for (size_t k = 0; k < Ncomps; k++){
            for (size_t m = 0; m < Ncomps; m++){
                tmp += rho * x[m] * b(m, k);
            }
        }
        DT(i) += (Dij(i, dependent_idx) / ksi[dependent_idx]) * tmp;
    }
    DT /= AVOGADRO;

    if (use_zarate){
        Eigen::VectorXd DT_indep(Ncomps - 1);
        Eigen::MatrixXd D_indep = interdiffusion(T, Vm, x, N, frame_of_reference, dependent_idx, solvent_idx);

        size_t DT_idx = 0;
        Eigen::VectorXd x_factor(Ncomps - 1);
        for (size_t i = 0; i < Ncomps; i++){
            if (i == dependent_idx) continue;
            x_factor[DT_idx] = x[i];
            DT_indep[DT_idx] = DT[i];
            DT_idx += 1;
        }
        Eigen::MatrixXd X = get_zarate_X_matr(x, dependent_idx);
        double c = 1. / Vm;
        DT = (- c * X).partialPivLu().solve((DT_indep + c * D_indep * x_factor) / T);
    }
    return DT;
}

Eigen::VectorXd KineticGas::thermal_diffusion_ratio(double T, double Vm, const vector1d& x, int N){
    double rho = AVOGADRO / Vm;
    Eigen::VectorXd DT = thermal_diffusion_coeff(T, Vm, x, N, FrameOfReference::CoM, -1);
    Eigen::MatrixXd Dij = interdiffusion(T, Vm, x, N, FrameOfReference::CoM, -1, -1, false);
    vector2d rdf = get_rdf(rho, T, x);
    vector2d etl = get_etl(rho, T, x);
    vector1d ksi = get_ksi_factors(T, Vm, x);
    Eigen::MatrixXd A(Ncomps, Ncomps);

    for (size_t i = 0; i < Ncomps - 1; i++){
        for (size_t j = 0; j < Ncomps; j++){
            A(i, j) = - Dij(i, j) * x[j] * (1 / Vm);
        }
        A(Ncomps - 1, i) = x[i] * ksi[i];
    }
    A(Ncomps - 1, Ncomps - 1) = x[Ncomps - 1] * ksi[Ncomps - 1];

    if (is_idealgas){
        DT(Ncomps - 1) = 1;
    }
    else{
        DT(Ncomps - 1) = 0.;
        for (size_t i = 0; i < Ncomps; i++){
            for (size_t j = 0; j < Ncomps; j++){
                int k_delta = (i == j) ? 1 : 0;
                DT(Ncomps - 1) += x[i] * (k_delta + (4. * PI / 3.) * rho * x[j] * pow(etl[i][j], 3) * M[i][j] * rdf[i][j]);
            }
        }
    }
    return A.partialPivLu().solve(DT);
}

Eigen::MatrixXd KineticGas::thermal_diffusion_factor(double T, double Vm, const vector1d& x, int N){
    Eigen::VectorXd kT = thermal_diffusion_ratio(T, Vm, x, N);
    Eigen::MatrixXd alpha(Ncomps, Ncomps);
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = 0; j < Ncomps; j++){
            alpha(i, j) = kT(i) - kT(j);
        }
    }
    return alpha;
}

Eigen::MatrixXd KineticGas::interdiffusion_dependent_CoM(double T, double Vm, const std::vector<double>& x, int N){
    double rho = AVOGADRO / Vm;
    vector3d d = reshape_diffusive_expansion_vector(compute_diffusive_expansion_coeff(rho, T, x, N));
    vector2d E = get_chemical_potential_factors(T, Vm, x);
    Eigen::MatrixXd Dij(Ncomps, Ncomps);
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = 0; j < Ncomps; j++){
            for (size_t k = 0; k < Ncomps; k++){
                Dij(i, j) += d[0][i][k] * E[k][j];
            }
            Dij(i, j) *= x[i] / (2 * rho);
        }
    }
    return Dij;
}

Eigen::MatrixXd KineticGas::CoM_to_FoR_matr(double T, double Vm, const vector1d& x, int frame_of_reference, int solvent_idx=-1){
    switch (frame_of_reference) {
    case FrameOfReference::CoM:
        return Eigen::MatrixXd::Identity(Ncomps, Ncomps);
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
Eigen::MatrixXd KineticGas::CoM_to_CoN_matr(double T, double Vm, const vector1d& x){
    Eigen::MatrixXd matr = Eigen::MatrixXd::Identity(Ncomps, Ncomps);
    std::vector<double> wt_fracs = get_wt_fracs(x);
    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j < Ncomps; j++){
            matr(i, j) += x[i] * ((wt_fracs[j] / x[j]) - 1);
        }
    }
    return matr;
}
Eigen::MatrixXd KineticGas::CoM_to_solvent_matr(double T, double Vm, const vector1d& x, int solvent_idx){
    Eigen::MatrixXd matr = Eigen::MatrixXd::Identity(Ncomps, Ncomps);
    std::vector<double> wt_fracs = get_wt_fracs(x);
    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j < Ncomps; j++){
            matr(i, j) += x[i] * ((wt_fracs[j] / x[j]) - (delta(j, solvent_idx) / x[solvent_idx]));
        }
    }
    return matr;
}
Eigen::MatrixXd KineticGas::CoM_to_CoV_matr(double T, double Vm, const vector1d& x){
    #ifdef DEBUG
        if (!eos.get()) throw std::runtime_error("EoS is not set (in CoM_to_CoV_matr)");
    #endif
    double p = eos->pressure_tv(T, Vm, x);
    std::vector<double> dvdn = eos->dvdn(T, p, x, eos->VAPPH);
    Eigen::MatrixXd psi = Eigen::MatrixXd::Identity(Ncomps, Ncomps);
    std::vector<double> wt_frac = get_wt_fracs(x);
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = 0; j < Ncomps; j++){
            psi(i, j) += x[i] * ((wt_frac[j] / x[j]) - (dvdn[j] / Vm));
        }
    }
    return psi;
}
Eigen::MatrixXd KineticGas::get_zarate_X_matr(const vector1d& x, int dependent_idx){
    while (dependent_idx < 0) dependent_idx += Ncomps;
    Eigen::MatrixXd X(Ncomps - 1, Ncomps - 1);
    size_t Mi = 0;
    for (size_t i = 0; i < Ncomps; i++){
        if (i == dependent_idx) continue;
        X(Mi, Mi) = x[i];
        size_t Mj = 0;
        for (size_t j = 0; j < Ncomps; j++){
            if (j == dependent_idx) continue;
            X(Mi, Mj) -= x[i] * x[j];
            Mj++;
        }
        Mi++;
    }
    return X;
}
Eigen::MatrixXd KineticGas::get_zarate_W_matr(const vector1d& x, int dependent_idx){
    vector1d wt_fracs = get_wt_fracs(x);
    return get_zarate_X_matr(wt_fracs, dependent_idx);
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
    if (is_singlecomp){ // dth_vector is related to thermal diffusion, which vanishes for a pure component.
        return Eigen::VectorXd::Zero(Ncomps);
    }
    Eigen::MatrixXd d_ij(Ncomps, Ncomps);
    Eigen::VectorXd l_i(Ncomps);

    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j < Ncomps; j++){
            d_ij(i, j) = d_qij[0][i][j] / d_qij[0][0][0];
        }
        l_i(i) = l(i);
    }

    return d_ij.partialPivLu().solve(l_i / d_qij[0][0][0]);
}
Eigen::VectorXd KineticGas::compute_thermal_expansion_coeff(double rho, double T, const vector1d& x, int N){
    Eigen::MatrixXd cond_matr{stdvector_to_EigenMatrix(get_conductivity_matrix(rho, T, x, N))};
    Eigen::VectorXd cond_vec{stdvector_to_EigenVector(get_conductivity_vector(rho, T, x, N))};
    return cond_matr.partialPivLu().solve(cond_vec);
}
vector2d KineticGas::get_chemical_potential_factors(double T, double Vm, const std::vector<double>& x){
    #ifdef DEBUG
        if (!eos.get()) throw std::runtime_error("EoS is not set (in get_chemical_potential_factors)");
    #endif
    vector2d dmudrho(Ncomps, vector1d(Ncomps, 0));
    if (is_singlecomp){
        // const vector2d dmudn_pure = eos->dmudn(T, Vm, {1.});
        // const double RT = GAS_CONSTANT * T;
        // const double dmudrho_pure = Vm * dmudn_pure[0][0];
        // const double rho = 1 / Vm;
        // dmudrho[0][0] = (dmudrho_pure + RT * x[1] / (rho * x[0])) / AVOGADRO;
        // dmudrho[0][1] = dmudrho[1][0] = (dmudrho_pure - RT / rho) / AVOGADRO;
        // dmudrho[1][1] = (dmudrho_pure + RT * x[0] / (rho * x[1])) / AVOGADRO;
        const double RT = GAS_CONSTANT * T;
        dmudrho[0][0] = (Vm*RT+Vm*RT*x[1]/x[0]) / AVOGADRO;
        dmudrho[0][1] = (Vm*RT-Vm*RT) / AVOGADRO;
        dmudrho[1][0] = (Vm*RT-Vm*RT) / AVOGADRO;
        dmudrho[1][1] = (Vm*RT+Vm*RT*x[0]/x[1]) / AVOGADRO;
        std::cout << dmudrho[0][0] << std::endl;
        std::cout << dmudrho[1][0] << std::endl;
        std::cout << dmudrho[1][1] << std::endl;
    }
    else{
        const vector2d dmudn = eos->dmudn(T, Vm, x);
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
vector1d KineticGas::get_ksi_factors(double T, double Vm, const vector1d& x){
    vector2d E = get_chemical_potential_factors(T, Vm, x);
    vector1d ksi(Ncomps, 0.);
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = 0; j < Ncomps; j++){
            ksi[i] += E[j][i];
        }
    }
    return ksi;
}