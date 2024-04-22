#include "Sutherland.h"

double Sutherland::potential(int i, int j, double r){
    double p{0.0};
    for (size_t k = 0; k < nterms; k++){
        p += C[k][i][j] * pow(sigma[i][j] / r, lambda[k][i][j]);
    }
    return eps[i][j] * p;
}

double Sutherland::potential_derivative_r(int i, int j, double r){
    double p{0.0};
    for (size_t k = 0; k < nterms; k++){
        p -= C[k][i][j] * lambda[k][i][j] * pow(sigma[i][j], lambda[k][i][j]) / pow(r, lambda[k][i][j] + 1);
    }
    return eps[i][j] * p;
}

double Sutherland::potential_dblderivative_rr(int i, int j, double r){
    double p{0.0};
    for (size_t k = 0; k < nterms; k++){
        p += C[k][i][j] * lambda[k][i][j] * (lambda[k][i][j] + 1.) * pow(sigma[i][j], lambda[k][i][j]) / pow(r, lambda[k][i][j] + 2);
    }
    return eps[i][j] * p;
}

vector2d Sutherland::get_b_max(double T){
    std::vector<std::vector<int>> ierr(Ncomps, std::vector<int>(Ncomps, 1));
    vector2d b_max = Spherical::get_b_max(T, ierr);
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = i; j < Ncomps; j++){
            if (ierr[i][j]) {
                std::cout << "Could not compute upper integration limit for collision diameter (" << i << ", " << j
                            << ") using Barker-Henderson diameter as fallback value." << std::endl;
                vector2d d_BH = get_BH_diameters(T);
                b_max[j][i] = b_max[i][j] = d_BH[i][j];
            }
        }
    }
    return b_max;
}

void Sutherland::compute_sigma_eff(){
    // Assumes that sigma_eff is already set to a reasonable initial guess (i.e. a value smaller than the true value of sigma_min)
    // If sigma_eff > sigma_min (where sigma_min denotes the true value) when this function is called, the solver will diverge.
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            while (abs(potential(i, j, sigma_eff[i][j])) / eps[i][j] > 1e-6){
                sigma_eff[i][j] -= potential(i, j, sigma_eff[i][j]) / potential_derivative_r(i, j, sigma_eff[i][j]);
            }
            sigma_eff[j][i] = sigma_eff[i][j];
        }
    }
}

void Sutherland::compute_epsilon_eff(){
    // Assumes that sigma_min is already set to a reasonable initial guess (i.e. a value smaller than the maximum of the first derivative of the potential)
    // If sigma_min is larger than this value when the function is called, the solver will diverge.
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            while (abs(potential_derivative_r(i, j, sigma_min[i][j])) * sigma_eff[i][j] / eps[i][j] > 1e-10){
                sigma_min[i][j] -= potential_derivative_r(i, j, sigma_min[i][j]) / potential_dblderivative_rr(i, j, sigma_min[i][j]);
            }
            eps_eff[i][j] = eps_eff[j][i] = - potential(i, j, sigma_min[i][j]);
        }
    }
}

void Sutherland::compute_vdw_alpha(){
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = 0; j < Ncomps; j++){
            vdw_alpha[i][j] = 0.;
            for (size_t k = 0; k < nterms; k++){
                vdw_alpha[i][j] -= C[k][i][j] * eps[i][j] * pow(sigma[i][j] / sigma_eff[i][j], lambda[k][i][j]) / (lambda[k][i][j] - 3.);
            }
            vdw_alpha[i][j] /= eps_eff[i][j];
        }
    }
}

vector2d Sutherland::saft_rdf(double rho, double T, const std::vector<double>& x, int order, bool g2_correction){
    double beta = (1. / (BOLTZMANN * T));
    vector2d d_BH = get_BH_diameters(T);
    vector2d x_eff = get_xeff(d_BH);
    vector2d gHS = rdf_g0_func(rho, x, d_BH);
    vector2d g1 = (order > 0) ? rdf_g1_func(rho, x, d_BH) : vector2d(Ncomps, vector1d(Ncomps, 0.0));
    vector2d g2 = (order > 1) ? rdf_g2_func(rho, T, x, d_BH, x_eff, g2_correction) : vector2d(Ncomps, vector1d(Ncomps, 0.0));
    vector2d g(Ncomps, vector1d(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            g[i][j] = gHS[i][j] * exp(beta * eps[i][j] * (g1[i][j] / gHS[i][j]) + pow(beta * eps[i][j], 2) * (g2[i][j] / gHS[i][j]));
            g[j][i] = g[i][j];
        }
    }
    return g;
}

vector2d Sutherland::rdf_g0_func(double rho, const vector1d& x, const vector2d& d_BH){
    double zeta_x = zeta_x_func(rho, x, d_BH);
    double k0 = - log(1 - zeta_x) + (42. * zeta_x - 39. * pow(zeta_x, 2) + 9. * pow(zeta_x, 3) - 2 * pow(zeta_x, 4)) / (6 * pow(1 - zeta_x, 3));
    double k1 = (pow(zeta_x, 4) + 6. * pow(zeta_x, 2) - 12. * zeta_x) / (2. * pow(1 - zeta_x, 3));
    double k2 = - 3. * pow(zeta_x, 2) / (8. * pow(1 - zeta_x, 2));
    double k3 = (- pow(zeta_x, 4.) + 3. * pow(zeta_x, 2) + 3. * zeta_x) / (6. * pow(1 - zeta_x, 3));
    
    double x_eff;
    vector2d rdf(Ncomps, std::vector<double>(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j <= i; j++){
            x_eff = sigma_eff[i][j] / d_BH[i][j];
            rdf[i][j] = exp(k0 + k1 * x_eff + k2 * pow(x_eff, 2) + k3 * pow(x_eff, 3));
            rdf[j][i] = rdf[i][j];
        }
    }
    return rdf;
}

vector2d Sutherland::rdf_g0_func(double rho, double T, const vector1d& x){
    vector2d d_BH = get_BH_diameters(T);
    return rdf_g0_func(rho, x, d_BH);
}

vector2d Sutherland::rdf_g1_func(double rho, const vector1d& x, const vector2d& d_BH){
    vector2d g1(Ncomps, std::vector<double>(Ncomps, 0.));
    vector2d dadr = da1_drho_func(rho, x, d_BH);
    vector2d x_eff = get_xeff(d_BH);
    vector2d x0 = get_x0(d_BH);
    double zeta_x = zeta_x_func(rho, x, d_BH);

    for (size_t k = 0; k < nterms; k++){
        vector2d a1s_k = a_1s_func(rho, x, zeta_x, d_BH, lambda[k]);
        vector2d B_k = B_func(rho, x, zeta_x, x_eff, d_BH, lambda[k]);
        for (int i = 0; i < Ncomps; i++){
            for (int j = i; j < Ncomps; j++){
                g1[i][j] += lambda[k][i][j] * C[k][i][j] * pow(x0[i][j], lambda[k][i][j]) * (a1s_k[i][j] + B_k[i][j]) / rho;
            }
        }
    }

    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = i; j < Ncomps; j++){
            g1[i][j] += 3. * dadr[i][j];
            g1[i][j] *= (1.0 / (2. * PI * eps[i][j] * pow(d_BH[i][j], 3)));
            g1[j][i] = g1[i][j];
        }
    }
    return g1;
}

vector2d Sutherland::rdf_g1_func(double rho, double T, const std::vector<double>& x){
    std::vector<std::vector<double>> d_BH = get_BH_diameters(T);
    return rdf_g1_func(rho, x, d_BH);
}

vector2d Sutherland::rdf_g2_func(double rho, double T, const vector1d& x, const vector2d& d_BH, const vector2d& x_eff, bool g2_correction){
    vector2d g2(Ncomps, std::vector<double>(Ncomps));

    const vector2d x0 = get_x0(d_BH);
    double zeta_x = zeta_x_func(rho, x, d_BH);
    double K_HS = K_HS_func(zeta_x);
    for (size_t k = 0; k < nterms; k++){
        for (size_t l = 0; l < nterms; l++){
            vector2d lambda_kl = get_lambda_kl(k, l);
            vector2d a1s_kl = a_1s_func(rho, x, zeta_x, d_BH, lambda_kl);
            vector2d B_kl = B_func(rho, x, zeta_x, x_eff, d_BH, lambda_kl);
            for (size_t i = 0; i < Ncomps; i++){
                for (size_t j = i; j < Ncomps; j++){
                    g2[i][j] -= K_HS * eps[i][j] * C[k][i][j] * C[l][i][j] * lambda[k][i][j] * pow(x0[i][j], lambda_kl[i][j]) * (a1s_kl[i][j] + B_kl[i][j]) / rho;
                }
            }
        }
    }
    vector2d da2ij_div_chi_drho = da2ij_div_chi_drho_func(rho, x, K_HS, d_BH, x_eff);
    double zeta_x_HS = zeta_x_func(rho, x, sigma_eff);
    vector2d gamma_c = g2_correction ? gamma_corr(zeta_x_HS, T) : vector2d(Ncomps, vector1d(Ncomps, 0.0));
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = i; j < Ncomps; j++){
            g2[i][j] += 3. * da2ij_div_chi_drho[i][j];
            g2[i][j] *= (1 + gamma_c[i][j]) * (1. / (2. * PI * pow(eps[i][j], 2) * pow(d_BH[i][j], 3))) ;
            g2[j][i] = g2[i][j];
        }
    }
    return g2;
}

vector2d Sutherland::rdf_g2_func(double rho, double T, const vector1d& x, bool g2_correction){
    vector2d d_BH = get_BH_diameters(T);
    vector2d x_eff = get_xeff(d_BH);
    return rdf_g2_func(rho, T, x, d_BH, x_eff, g2_correction);
}

vector2d Sutherland::get_lambda_kl(size_t k, size_t l){
    vector2d lambda_kl(Ncomps, vector1d(Ncomps, 0.));
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = 0; j < Ncomps; j++){
            lambda_kl[j][i] = lambda_kl[i][j] = lambda[k][i][j] + lambda[l][i][j];
        }
    }
    return lambda_kl;
}

vector2d Sutherland::gamma_corr(double zeta_x, double T){
    vector2d gamma(Ncomps, std::vector<double>(Ncomps));
    constexpr double phi[5] = {10., 10., 0.57, -6.7, -8.};

    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            double theta = exp(eps_eff[i][j] / (BOLTZMANN * T)) - 1.;
            gamma[i][j] = phi[0] * (1 - tanh(phi[1] * (phi[2] - vdw_alpha[i][j]))) 
                        * zeta_x * theta * exp(phi[3] * zeta_x + phi[4] * pow(zeta_x, 2));
            gamma[j][i] = gamma[i][j];
        }
    }
    return gamma;
}

vector2d Sutherland::get_BH_diameters(double T){
    // 20-point Gauss-Legendre from 0.5 sigma to 1 sigma. Using constant value of 1 for integrand within (0, 0.5) sigma.
    std::vector<std::vector<double>> d_BH(Ncomps, std::vector<double>(Ncomps, 0.0));
    double beta = 1. / (BOLTZMANN * T);
    for (int i = 0; i < Ncomps; i++){
        d_BH[i][i] = 2.;
        for (int n = 0; n < 20; n++){
            d_BH[i][i] += rdf_constants.gl_w[n] * (1. - exp(- beta * potential(i, i, sigma_eff[i][i] * (rdf_constants.gl_x[n] / 4. + 3. / 4.))));
        }
        d_BH[i][i] *= sigma_eff[i][i] / 4.;
    }
    for (int i = 0; i < Ncomps - 1; i++){
        for (int j = i + 1; j < Ncomps; j++){
            d_BH[i][j] = (d_BH[i][i] + d_BH[j][j]) / 2.0;
            d_BH[j][i] = d_BH[i][j];
        }
    }
    return d_BH;
}

vector2d Sutherland::get_x0(const vector2d& d_BH){
    vector2d x0(Ncomps, std::vector<double>(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j <= i; j++){
            x0[i][j] = sigma[i][j] / d_BH[i][j];
            x0[j][i] = x0[i][j];
        }
    }
    return x0;
}

vector2d Sutherland::get_xeff(const vector2d& d_BH){
    vector2d x_eff(Ncomps, std::vector<double>(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j <= i; j++){
            x_eff[i][j] = sigma_eff[i][j] / d_BH[i][j];
            x_eff[j][i] = x_eff[i][j];
        }
    }
    return x_eff;
}

vector2d Sutherland::da2ij_div_chi_drho_func(double rho, const vector1d& x, double K_HS, const vector2d& d_BH, const vector2d& x_eff){
    double zeta_x = zeta_x_func(rho, x, d_BH);
    double dzxdr = dzetax_drho_func(x, d_BH);
    double dKHS_drho = dKHS_drho_func(zeta_x, dzxdr);
    vector2d rdf_chi_HS = rdf_chi_func(rho, x);
    vector2d a2ij = a2ij_func(rho, x, K_HS, rdf_chi_HS, d_BH, x_eff);
    vector2d dchi_HS_drho = drdf_chi_drho_func(rho, x);
    vector2d da2ij_drho(Ncomps, vector1d(Ncomps, 0.0));
    const vector2d x0 = get_x0(d_BH);

    for (size_t k = 0; k < nterms; k++){
        for (size_t l = 0; l < nterms; l++){
            vector2d lambda_kl = get_lambda_kl(k, l);
            vector2d da1s_drho = da1s_drho_func(rho, x, d_BH, lambda_kl);
            vector2d dB_drho = dBdrho_func(rho, x, zeta_x, x_eff, d_BH, lambda_kl);
            for (size_t i = 0; i < Ncomps; i++){
                for (size_t j = i; j < Ncomps; j++){
                    da2ij_drho[i][j] += .5 * K_HS * eps[i][j] * C[k][i][j] * C[l][i][j]
                                    * pow(x0[i][j], lambda_kl[i][j]) * (da1s_drho[i][j] + dB_drho[i][j]);

                }
            }
        }
    }
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = i; j < Ncomps; j++){
            da2ij_drho[i][j] += dKHS_drho * (a2ij[i][j] / (K_HS * (1 + rdf_chi_HS[i][j])));
            da2ij_drho[j][i] = da2ij_drho[i][j];
        }
    }
    return da2ij_drho;
}

vector2d Sutherland::a2ij_func(double rho, const vector1d& x, double K_HS, const vector2d& rdf_chi_HS, const vector2d& d_BH, const vector2d& x_eff){
    vector2d a2ij(Ncomps, std::vector<double>(Ncomps, 0.0));
    double zeta_x = zeta_x_func(rho, x, d_BH);
    const vector2d x0 = get_x0(d_BH);
    for (size_t k = 0; k < nterms; k++){
        for (size_t l = 0; l < nterms; l++){
            vector2d lambda_kl = get_lambda_kl(k, l);
            vector2d a1s = a_1s_func(rho, x, zeta_x, d_BH, lambda_kl);
            vector2d B = B_func(rho, x, zeta_x, x_eff, d_BH, lambda_kl);
            for (size_t i = 0; i < Ncomps; i++){
                for (size_t j = i; j < Ncomps; j++){
                    a2ij[i][j] += C[k][i][j] * C[l][i][j] * pow(x0[i][j], lambda_kl[i][j]) * (a1s[i][j] + B[i][j]);
                }
            }
        }
    }
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = i; j < Ncomps; j++){
            a2ij[i][j] *= 0.5 * K_HS * (1 + rdf_chi_HS[i][j]) * eps[i][j];
            a2ij[j][i] = a2ij[i][j];
        }
    }
    return a2ij;
}

vector2d Sutherland::rdf_chi_func(double rho, const vector1d& x){
    vector2d rdf_chi(Ncomps, vector1d(Ncomps));
    double zeta_x = zeta_x_func(rho, x, sigma_eff);
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            vector1d f = f_corr(vdw_alpha[i][j]);
            rdf_chi[i][j] = f[0] * zeta_x + f[1] * pow(zeta_x, 5) + f[2] * pow(zeta_x, 8);
            rdf_chi[j][i] = rdf_chi[i][j];
        }
    }
    return rdf_chi;
}

vector1d Sutherland::f_corr(double alpha){
    std::vector<double> f(3);
    for (int i = 0; i < 3; i++){
        double num{0.0}, denom{1.0};
        for (int n = 0; n <= 3; n++){
            num += rdf_constants.phi[n][i] * pow(alpha, n);
        }
        for (int n = 4; n <= 6; n++){
            denom += rdf_constants.phi[n][i] * pow(alpha, n - 3);
        }
        f[i] = num / denom;
    }
    return f;
}

vector2d Sutherland::drdf_chi_drho_func(double rho, const vector1d& x){
    vector2d drdf_chi_drho(Ncomps, vector1d(Ncomps, 0.0));
    double zeta_x = zeta_x_func(rho, x, sigma_eff);
    double dzdr = dzetax_drho_func(x, sigma_eff);
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            std::vector<double> f = f_corr(vdw_alpha[i][j]);
            drdf_chi_drho[i][j] = dzdr * (f[0] + 5. * f[1] * pow(zeta_x, 4) + 8. * f[2] * pow(zeta_x, 7));
            drdf_chi_drho[j][i] = drdf_chi_drho[i][j];
        }
    }
    return drdf_chi_drho;
}

double Sutherland::zeta_x_func(double rho, const vector1d& x, const vector2d& d_BH){
    double zeta{0.0};
    for (int i = 0; i < x.size(); i++){
        for (int j = 0; j < x.size(); j++){
            zeta += x[i] * x[j] * pow(d_BH[i][j], 3);
        }
    }
    zeta *= PI * rho / 6.0;
    return zeta;
}

double Sutherland::zeta_eff_func(double rho,  const vector1d& x, double zeta_x, double lambdaijk){
    std::vector<double> c_coeffs(4, 0.0);
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            c_coeffs[i] += rdf_constants.C_coeff_matr[i][j] * pow(lambdaijk, - j);
        }
    }
    double zeta_eff{0.0};
    for (int i = 0; i < 4; i++){
        zeta_eff += c_coeffs[i] * pow(zeta_x, i + 1);
    }
    return zeta_eff;
}

double Sutherland::dzeta_eff_drho_func(double rho, const std::vector<double>& x, const vector2d& d_BH, double lambdaijk){

    std::vector<double> c_coeffs(4, 0.0);
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            c_coeffs[i] += rdf_constants.C_coeff_matr[i][j] * pow(lambdaijk, - j);
        }
    }
    double zeta_x = zeta_x_func(rho, x, d_BH);
    double dzdrho{0.0};
    for (int i = 0; i < 4; i++){
        dzdrho += c_coeffs[i] * (i + 1) * pow(zeta_x, i) * dzetax_drho_func(x, d_BH);
    }
    return dzdrho;
}

double Sutherland::dzetax_drho_func(const vector1d& x, const vector2d& d_BH){
    double dzdrho{0.0};
    for (int i = 0; i < x.size(); i++){
        for (int j = 0; j < x.size(); j++){
            dzdrho += x[i] * x[j] * pow(d_BH[i][j], 3);
        }
    }
    dzdrho *= PI / 6.0;
    return dzdrho;
}

vector2d Sutherland::a_1s_func(double rho, const vector1d& x, double zeta_x, const vector2d& d_BH, const vector2d& lambda_k){
    double zeta_eff;
    vector2d a1s(Ncomps, std::vector<double>(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            zeta_eff = zeta_eff_func(rho, x, zeta_x, lambda_k[i][j]);
            a1s[i][j] = (-2. * rho * PI * eps[i][j] * pow(d_BH[i][j], 3) / (lambda_k[i][j] - 3)) * (1. - zeta_eff/2.) / pow(1 - zeta_eff, 3);
            a1s[j][i] = a1s[i][j];
        }
    }
    return a1s;
}

vector2d Sutherland::a_1s_func(double rho, double T, const vector1d& x, const vector2d& lambda_k){
    vector2d d_BH = get_BH_diameters(T);
    double zeta_x = zeta_x_func(rho, x, d_BH);
    return a_1s_func(rho, x, zeta_x, d_BH, lambda_k);
}

vector2d Sutherland::B_func(double rho, const std::vector<double>& x, double zeta_x, const vector2d& x_eff, const vector2d& d_BH, const vector2d& lambda_k){
    vector2d I{I_func(x_eff, lambda_k)}, J{J_func(x_eff, lambda_k)};
    vector2d B(Ncomps, vector1d(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            B[i][j] = 2. * PI * rho * pow(d_BH[i][j], 3) * eps[i][j] 
                        * (I[i][j] * (1. - zeta_x / 2.) / pow(1. - zeta_x, 3)
                            - J[i][j] * 9. * zeta_x * (1. + zeta_x) / (2. * pow(1. - zeta_x, 3)));
            B[j][i] = B[i][j];
        }
    }
    return B;
}

vector2d Sutherland::B_func(double rho, const vector1d& x, const vector2d& d_BH, const vector2d& lambda_k){
    double zeta_x = zeta_x_func(rho, x, d_BH);
    vector2d x_eff = get_xeff(d_BH);
    return B_func(rho, x, zeta_x, x_eff, d_BH, lambda_k);
}

vector2d Sutherland::B_func(double rho, double T, const vector1d& x, const vector2d& lambda_k){
    std::vector<std::vector<double>> d_BH = get_BH_diameters(T);
    return B_func(rho, x, d_BH, lambda_k);
}

vector2d Sutherland::I_func(const vector2d& xeff, const vector2d& lambda_k){
    vector2d I(Ncomps, std::vector<double>(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            I[i][j] = - (pow(xeff[i][j], 3. - lambda_k[i][j]) - 1.) / (lambda_k[i][j] - 3.);
            I[j][i] = I[i][j];
        }
    }
    return I;
}

vector2d Sutherland::J_func(const vector2d& xeff, const vector2d& lambda_k){
    vector2d J(Ncomps, std::vector<double>(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            J[i][j] = - (pow(xeff[i][j], 4. - lambda_k[i][j]) * (lambda_k[i][j] - 3.) 
                            - pow(xeff[i][j], 3. - lambda_k[i][j]) * (lambda_k[i][j] - 4.) - 1.)
                        / ((lambda_k[i][j] - 3.) * (lambda_k[i][j] - 4.));
            J[j][i] = J[i][j];
        }
    }
    return J;
}

vector2d Sutherland::a1_func(double rho, double T, const vector1d& x){
    vector2d d_BH = get_BH_diameters(T);
    return a1_func(rho, x, d_BH);
}

vector2d Sutherland::a1_func(double rho, const vector1d& x, const vector2d& d_BH){
    vector2d a1(Ncomps, vector1d(Ncomps, 0.0));
    vector2d x_eff = get_xeff(d_BH);
    vector2d x0 = get_x0(d_BH);
    double zeta_x = zeta_x_func(rho, x, d_BH);
    for (size_t k = 0; k < nterms; k++){
        vector2d a1s_k = a_1s_func(rho, x, zeta_x, d_BH, lambda[k]);
        vector2d B_k = B_func(rho, x, zeta_x, x_eff, d_BH, lambda[k]);
        for (size_t i = 0; i < Ncomps; i++){
            for (size_t j = i; j < Ncomps; j++){
                a1[i][j] -= C[k][i][j] * pow(x0[i][j], lambda[k][i][j]) * (a1s_k[i][j] + B_k[i][j]);
                a1[j][i] = a1[i][j];
            }
        }
    }
    return a1;
}

vector2d Sutherland::da1_drho_func(double rho, const vector1d& x, const vector2d& d_BH){
    vector2d da1drho(Ncomps, std::vector<double>(Ncomps));
    const vector2d x_eff = get_xeff(d_BH);
    const vector2d x0 = get_x0(d_BH);
    double zeta_x = zeta_x_func(rho, x, d_BH);
    for (size_t k = 0; k < nterms; k++){
        vector2d da1s_drho = da1s_drho_func(rho, x, d_BH, lambda[k]);
        vector2d dB_drho = dBdrho_func(rho, x, zeta_x, x_eff, d_BH, lambda[k]);
        for (int i = 0; i < Ncomps; i++){
            for (int j = i; j < Ncomps; j++){
                da1drho[i][j] -= C[k][i][j] * pow(x0[i][j], lambda[k][i][j]) * (da1s_drho[i][j] + dB_drho[i][j]);
                da1drho[j][i] = da1drho[i][j];
            }
        }
    }
    return da1drho;
}

vector2d Sutherland::da1s_drho_func(double rho, const vector1d& x, const vector2d& d_BH, const vector2d& lambda_k){
    vector2d  da1sdrho(Ncomps, std::vector<double>(Ncomps));
    double zeta_x = zeta_x_func(rho, x, d_BH);
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            double ze = zeta_eff_func(rho, x, zeta_x, lambda_k[i][j]);
            double dzedrho = dzeta_eff_drho_func(rho, x, d_BH, lambda_k[i][j]);
            da1sdrho[i][j] = -2. * (PI * eps[i][j] * pow(d_BH[i][j], 3) / (lambda_k[i][j] - 3)) 
                                    * (((1. - ze / 2.) / pow(1 - ze, 3)) 
                                        + rho * dzedrho * (2.5 - ze) / pow(1 - ze, 4));
            da1sdrho[j][i] = da1sdrho[i][j];
        }
    }
    return da1sdrho;
}

vector2d Sutherland::dBdrho_func(double rho, const vector1d& x, double zeta_x, const vector2d& x_eff, const vector2d& d_BH, const vector2d& lambda_k){
    vector2d B = B_func(rho, x, zeta_x, x_eff, d_BH, lambda_k);
    vector2d I = I_func(x_eff, lambda_k);
    vector2d J = J_func(x_eff, lambda_k);
    double dzxdrho = dzetax_drho_func(x, d_BH);
    vector2d dBdrho(Ncomps, std::vector<double>(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            dBdrho[i][j] = B[i][j] / rho
                            + 2. * PI * rho * eps[i][j] * pow(d_BH[i][j], 3)
                            * (I[i][j] * (- 0.5 * (1 - zeta_x) + 3. * (1. - 0.5 * zeta_x)) / pow(1 - zeta_x, 4)
                                - J[i][j] * (9. * (1. + 2. * zeta_x) * (1. - zeta_x) + 27. * zeta_x * (1. + zeta_x)) / (2. * pow(1 - zeta_x, 4))
                                ) * dzxdrho;
            dBdrho[j][i] = dBdrho[i][j];
        }
    }
    return dBdrho;
}

vector2d Sutherland::dBdrho_func(double rho, double T, const vector1d& x, const vector2d& lambda_k)
{
    vector2d d_BH = get_BH_diameters(T);
    double zeta_x = zeta_x_func(rho, x, d_BH);
    vector2d x_eff = get_xeff(d_BH);

    return dBdrho_func(rho, x, zeta_x, x_eff, d_BH, lambda_k);
}