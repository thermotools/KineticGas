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

vector2d Sutherland::get_contact_diameters(double rho, double T, const vector1d& x){
    throw std::runtime_error("Collision diameters not implemented for Sutherland!");
}

void Sutherland::compute_sigma_eff(){
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            sigma_eff[i][j] = sigma[i][j];
            while (abs(potential(i, j, sigma_eff[i][j])) / eps[i][j] > 1e-6){
                sigma_eff[i][j] -= potential(i, j, sigma_eff[i][j]) / potential_derivative_r(i, j, sigma_eff[i][j]);
            }
            sigma_eff[j][i] = sigma_eff[i][j];
        }
    }
}

void Sutherland::compute_epsilon_eff(){
    sigma_min = sigma_eff;
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            while (abs(potential_derivative_r(i, j, sigma_min[i][j])) * sigma_eff[i][j] / eps[i][j] > 1e-10){
                sigma_min[i][j] -= potential_derivative_r(i, j, sigma_min[i][j]) / potential_dblderivative_rr(i, j, sigma_min[i][j]);
            }
            eps_eff[i][j] = eps_eff[j][i] = - potential(i, j, sigma_min[i][j]);
        }
    }
}

vector2d Sutherland::model_rdf(double rho, double T, const vector1d& x){
    double beta = (1. / (BOLTZMANN * T));
    vector2d d_BH = get_BH_diameters(T);
    vector2d x0 = get_x0(d_BH);
    vector2d gHS = rdf_g0_func(rho, x, d_BH);
    vector2d g1 = rdf_g1_func(rho, x, d_BH);
    // vector2d g2 = rdf_g2_func(rho, T, x, d_BH, x0);
    vector2d g(Ncomps, vector1d(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            g[i][j] = gHS[i][j] * exp(beta * eps[i][j] * (g1[i][j] / gHS[i][j])); //  + pow(beta * eps[i][j], 2) * (g2[i][j] / gHS[i][j]));
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
    vector2d x0 = get_x0(d_BH);
    double zeta_x = zeta_x_func(rho, x, d_BH);

    for (size_t k = 0; k < nterms; k++){
        vector2d a1s_k = a_1s_func(rho, x, zeta_x, d_BH, lambda[k]);
        vector2d B_k = B_func(rho, x, zeta_x, x0, d_BH, lambda[k]);
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

/*
vector2d rdf_g2_func(double rho, double T, const vector1d& x, const vector2d& d_BH, const vector2d& x0){
    vector2d g2(Ncomps, std::vector<double>(Ncomps));
    vector2d _la(Ncomps, std::vector<double>(Ncomps));
    vector2d _lr(Ncomps, std::vector<double>(Ncomps));
    vector2d la_lr(Ncomps, std::vector<double>(Ncomps));
    
    double zeta_x = zeta_x_func(rho, x, d_BH);
    double K_HS = K_HS_func(zeta_x);
    vector2d rdf_chi = rdf_chi_func(rho, x, d_BH);
    vector2d rdf_chi_HS = rdf_chi_func(rho, x, sigma);
    vector2d drdf_chi_drho = drdf_chi_drho_func(rho, x, d_BH);
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            _la[i][j] = _la[j][i] = 2 * la[i][j];
            _lr[i][j] = _lr[j][i] = 2 * lr[i][j];
            la_lr[i][j] = la_lr[j][i] = la[i][j] + lr[i][j];
        }
    }
    vector2d a1s_la = a_1s_func(rho, x, d_BH, _la);
    vector2d a1s_lr = a_1s_func(rho, x, d_BH, _lr);
    vector2d a1s_la_lr = a_1s_func(rho, x, d_BH, la_lr);
    vector2d B_la = B_func(rho, x, d_BH, _la);
    vector2d B_lr = B_func(rho, x, d_BH, _lr);
    vector2d B_la_lr = B_func(rho, x, d_BH, la_lr);
    vector2d da2ij_div_chi_drho = da2ij_div_chi_drho_func(rho, x, K_HS, d_BH, x0);
    vector2d a2ij = a2ij_func(rho, x, K_HS, rdf_chi_HS, d_BH, x0);
    double zeta_x_HS = zeta_x_func(rho, x, sigma);
    vector2d gamma_c = gamma_corr(zeta_x_HS, T);
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            g2[i][j] = (1 + gamma_c[i][j])
                        * (1. / (2 * PI * pow(eps[i][j], 2) * pow(d_BH[i][j], 3))) 
                        * (3 * da2ij_div_chi_drho[i][j] 
                            + eps[i][j] * K_HS * pow(C[i][j], 2) 
                                * (- lr[i][j] * pow(x0[i][j], 2 * lr[i][j]) * (a1s_lr[i][j] + B_lr[i][j]) / rho
                                    + la_lr[i][j] * pow(x0[i][j], la_lr[i][j]) * (a1s_la_lr[i][j] + B_la_lr[i][j]) / rho
                                    -  la[i][j] * pow(x0[i][j], 2 * la[i][j]) * (a1s_la[i][j] + B_la[i][j]) / rho
                                    )
                            );
            g2[j][i] = g2[i][j];
        }
    }
    return g2;
}


vector2d rdf_g2_func(double rho, double T, const vector1d& x){
    vector2d d_BH = get_BH_diameters(T);
    vector2d x0 = get_x0(d_BH);
    return rdf_g2_func(rho, T, x, d_BH, x0);
}
*/

vector2d Sutherland::get_BH_diameters(double T){
    // Gauss-Legendre points taken from SAFT-VR-MIE docs (see: ThermoPack)
    vector2d d_BH(Ncomps, std::vector<double>(Ncomps, 0.0));
    double beta = 1. / (BOLTZMANN * T);
    for (int i = 0; i < Ncomps; i++){
        for (int n = 0; n < 10; n++){
            d_BH[i][i] += sutherland_rdf::gl_w[n] * (1. - exp(- beta * potential(i, i, sigma[i][i] * (sutherland_rdf::gl_x[n] + 1) / 2.)));
        }
        d_BH[i][i] *= sigma[i][i] / 2;
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
            c_coeffs[i] += sutherland_rdf::C_coeff_matr[i][j] * pow(lambdaijk, - j);
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
            c_coeffs[i] += sutherland_rdf::C_coeff_matr[i][j] * pow(lambdaijk, - j);
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

vector2d Sutherland::da1_drho_func(double rho, const vector1d& x, const vector2d& d_BH){
    vector2d da1drho(Ncomps, std::vector<double>(Ncomps));
    const vector2d x_eff = get_xeff(d_BH);
    double zeta_x = zeta_x_func(rho, x, d_BH);
    for (size_t k = 0; k < nterms; k++){
        vector2d da1s_drho = da1s_drho_func(rho, x, d_BH, lambda[k]);
        vector2d dB_drho = dBdrho_func(rho, x, zeta_x, x_eff, d_BH, lambda[k]);
        for (int i = 0; i < Ncomps; i++){
            for (int j = 0; j < Ncomps; j++){
                da1drho[i][j] -= C[k][i][j] * pow(x_eff[i][j], lambda[k][i][j]) * (da1s_drho[i][j] + dB_drho[i][j]);
            }
        }
    }
    return da1drho;
}

vector2d Sutherland::da1s_drho_func(double rho, const vector1d& x, const vector2d& d_BH, const vector2d& lambda_k){
    vector2d  da1sdrho(Ncomps, std::vector<double>(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            double ze = zeta_eff_func(rho, x, d_BH, lambda_k[i][j]);
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
            dBdrho[i][j] = B[i][j] / rho + 2. * PI * eps[i][j] * d_BH[i][j] * dzxdrho 
                            * (I[i][j] * (2.5 - zeta_x) / pow(1 - zeta_x, 4) 
                                - 4.5 * J[i][j] * (1 + 3 * zeta_x - pow(zeta_x, 2)) / pow(1 - zeta_x, 2));
        }
    }
    return dBdrho;
}