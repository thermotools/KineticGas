#include "ExtendedSutherland.h"

vector2d dual_to_double(const vector2d2& vin){
    vector2d vout;
    for (auto it_outer = vin.begin(); it_outer != vin.end(); ++it_outer){
        vout.push_back(vector1d());
        for (auto it = it_outer->begin(); it != it_outer->end(); ++it){
            vout.back().push_back(static_cast<double>(*it));
        }
    }
    return vout;
}

// ------------------------------------------------------------------------------ //
// -------------------------- Constructors and helpers -------------------------- // 
// ------------------------------------------------------------------------------ //

ExtSutherland::ExtSutherland(std::string comps, size_t nterms, bool is_idealgas)
    : Spherical(comps, is_idealgas),
    C(nterms, vector2d(Ncomps, vector1d(Ncomps))), 
    lambda(nterms, vector2d(Ncomps, vector1d(Ncomps, 0.))), 
    beta_exp(nterms, vector2d(Ncomps, vector1d(Ncomps, 0.))), 
    rho_exp(nterms, vector2d(Ncomps, vector1d(Ncomps, 0.))), 
    nterms{nterms},
    sigma_eff(Ncomps, vector1d2(Ncomps)), 
    r_min(Ncomps, vector1d2(Ncomps)), 
    eps_eff(Ncomps, vector1d2(Ncomps)), 
    vdw_alpha(Ncomps, vector1d2(Ncomps)),
    C_eff(nterms, vector2d2(Ncomps, vector1d2(Ncomps)))
{}

ExtSutherland::ExtSutherland(vector1d mole_weights, vector2d sigma, vector2d eps, 
            vector3d C, vector3d lambda, vector3d beta_exp, vector3d rho_exp, 
            bool is_idealgas, bool is_singlecomp)
    : Spherical(mole_weights, sigma, eps, is_idealgas, is_singlecomp), 
    C{C}, 
    lambda{lambda}, 
    beta_exp(beta_exp), 
    rho_exp{rho_exp}, nterms{C.size()},
    sigma_eff(Ncomps, vector1d2(Ncomps)), 
    r_min(Ncomps, vector1d2(Ncomps)), 
    eps_eff(Ncomps, vector1d2(Ncomps)), 
    vdw_alpha(Ncomps, vector1d2(Ncomps)),
    C_eff(C.size(), vector2d2(Ncomps, vector1d2(Ncomps)))
{init_effective_params();}
    
ExtSutherland::ExtSutherland(vector1d mole_weights, vector2d sigma, vector2d eps, 
                        size_t nterms, bool is_idealgas, bool is_singlecomp)
    : Spherical(mole_weights, sigma, eps, is_idealgas, is_singlecomp), 
    C(nterms, vector2d(Ncomps, vector1d(Ncomps))), 
    lambda(nterms, vector2d(Ncomps, vector1d(Ncomps, 0.))), 
    beta_exp(nterms, vector2d(Ncomps, vector1d(Ncomps, 0.))), 
    rho_exp(nterms, vector2d(Ncomps, vector1d(Ncomps, 0.))), 
    nterms{nterms},
    sigma_eff(Ncomps, vector1d2(Ncomps)), 
    r_min(Ncomps, vector1d2(Ncomps)), 
    eps_eff(Ncomps, vector1d2(Ncomps)), 
    vdw_alpha(Ncomps, vector1d2(Ncomps)),
    C_eff(nterms, vector2d2(Ncomps, vector1d2(Ncomps)))
{init_effective_params();}

void ExtSutherland::init_effective_params(){
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = 0; j < Ncomps; j++){
            sigma_eff[i][j] = sigma_eff[j][i] = r_min[i][j] = r_min[j][i] = sigma[i][j];
            eps_eff[i][j] = eps_eff[j][i] = eps[i][j];
            for (size_t k = 0; k < nterms; k++){
                C_eff[k][i][j] = C[k][i][j];
            }
        }
    }
}

void ExtSutherland::mix_sigma(){
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = i; j < Ncomps; j++){
            sigma[i][j] = sigma[j][i] = 0.5 * (sigma[i][i] + sigma[j][j]);
        }
    }
}

void ExtSutherland::mix_epsilon(){
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = i; j < Ncomps; j++){
            eps[i][j] = eps[j][i] = sqrt(eps[i][i] * eps[j][j]);
        }
    }
}

void ExtSutherland::mix_exponents(std::vector<std::vector<double>>& expo){
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = i; j < Ncomps; j++){
            expo[i][j] = expo[j][i] = 3. + sqrt((expo[i][i] - 3.) * (expo[j][j] - 3.));
        }
    }
}

// ----------------------------------------------------------------------- //
// -------------------------- Potential Methods -------------------------- // 
// ----------------------------------------------------------------------- //

dual2 ExtSutherland::potential(int i, int j, dual2 r, dual2 T, dual2 rho) const {
    if (T > 1e4) throw std::runtime_error("Using density instead of T?");
    dual2 p{0.0};
    dual2 beta = eps[i][j] / (BOLTZMANN * T);
    dual2 rho_r = rho * pow(sigma[i][j], 3);
    for (size_t k = 0; k < nterms; k++){
        p += C[k][i][j] * pow(rho_r, rho_exp[k][i][j]) * pow(beta, beta_exp[k][i][j]) * pow(sigma[i][j] / r, lambda[k][i][j]);
    }
    return eps[i][j] * p;
}

double ExtSutherland::potential(int i, int j, double r) const {
    double p{0.0};
    for (size_t k = 0; k < nterms; k++){
        p += static_cast<double>(C_eff[k][i][j]) * pow(sigma[i][j] / r, lambda[k][i][j]);
    }
    return eps[i][j] * p;
}

double ExtSutherland::potential(int i, int j, double r, double T, double rho) const {
    if (T > 1e4) throw std::runtime_error("Using density instead of T?");
    double p{0.0};
    double beta = eps[i][j] / (BOLTZMANN * T);
    double rho_r = rho * pow(sigma[i][j], 3);
    for (size_t k = 0; k < nterms; k++){
        p += C[k][i][j] * pow(rho_r, rho_exp[k][i][j]) * pow(beta, beta_exp[k][i][j])* pow(sigma[i][j] / r, lambda[k][i][j]);
    }
    return eps[i][j] * p;
}

dual2 ExtSutherland::potential_r(int i, int j, dual2 r, dual2 T, dual2 rho) const {
    dual2 p{0.0};
    dual2 beta = eps[i][j] / (BOLTZMANN * T);
    dual2 rho_r = rho * pow(sigma[i][j], 3);
    for (size_t k = 0; k < nterms; k++){
        p -= C[k][i][j] * pow(rho_r, rho_exp[k][i][j]) * pow(beta, beta_exp[k][i][j]) * lambda[k][i][j] * pow(sigma[i][j] / r, lambda[k][i][j] + 1.) / sigma[i][j];
    }
    return eps[i][j] * p;
}

dual2 ExtSutherland::potential_rr(int i, int j, dual2 r, dual2 T, dual2 rho) const {
    dual2 p{0.0};
    dual2 beta = eps[i][j] / (BOLTZMANN * T);
    dual2 rho_r = rho * pow(sigma[i][j], 3);
    for (size_t k = 0; k < nterms; k++){
        p += C[k][i][j] * pow(rho_r, rho_exp[k][i][j]) * pow(beta, beta_exp[k][i][j]) * lambda[k][i][j] * (lambda[k][i][j] + 1.) * pow(sigma[i][j] / r, lambda[k][i][j] + 2) / pow(sigma[i][j], 2);
    }
    return eps[i][j] * p;
}

dual2 ExtSutherland::potential(int i, int j, dual2 r) const {
    return potential(i, j, r, current_T, current_rho);
}

dual2 ExtSutherland::potential_r(int i, int j, dual2 r) const {
    return potential_r(i, j, r, current_T, current_rho);
}

dual2 ExtSutherland::potential_rr(int i, int j, dual2 r) const {
    return potential_rr(i, j, r, current_T, current_rho);
}

// -------------------------------------------------------------------------------------- //
// -------------------------- Handling of effective parameters -------------------------- // 
// -------------------------------------------------------------------------------------- //

size_t ExtSutherland::set_effective_params(dual2 rho, dual2 T){
    if ((rho == current_rho) && (T == current_T)) return 0;
    size_t r;
    if ((rho != current_rho) && (T != current_T)) {r = 1;}
    else if (rho != current_rho) {r = 2;}
    else if (T != current_T) {r = 3;}
    set_C_eff(rho, T);
    set_sigma_eff(rho, T);
    set_epsilon_eff(rho, T);
    set_vdw_alpha(rho, T);
    current_rho = static_cast<double>(rho); current_T = static_cast<double>(T);
    return r;
}

void ExtSutherland::set_C_eff(dual2 rho, dual2 T){
    for (size_t k = 0; k < nterms; k++){
        for (size_t i = 0; i < Ncomps; i++){
            for (size_t j = i; j < Ncomps; j++){
                dual2 rho_r = rho * pow(sigma[i][j], 3);
                dual2 beta = eps[i][j] / (BOLTZMANN * T);
                dual2 C_eff_i = C[k][i][j];
                if (rho_exp[k][i][j] != 0){
                    C_eff_i *= pow(rho_r, rho_exp[k][i][j]);
                }
                if (beta_exp[k][i][j] != 0){
                    C_eff_i *= pow(beta, beta_exp[k][i][j]);
                }
                C_eff[k][j][i] = C_eff[k][i][j] = C_eff_i;
            }
        }
    }
    C_set = true;
}

void ExtSutherland::set_sigma_eff(dual2 rho, dual2 T){
    dual2 u, ur;
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            if ((rho != current_rho) || (T != current_T)){
                sigma_eff[i][j] = sigma[i][j];
            }
            do {
                u = potential(i, j, sigma_eff[i][j], T, rho);
                ur = potential_r(i, j, sigma_eff[i][j], T, rho);
                sigma_eff[i][j] -= u / ur;
            } while (abs(u) / eps[i][j] > 1e-6);
            sigma_eff[j][i] = sigma_eff[i][j];
        }
    }
}

void ExtSutherland::set_epsilon_eff(dual2 rho, dual2 T){
    dual2 ur, urr;
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            if ((rho != current_rho) || (T != current_T)){
                r_min[i][j] = sigma[i][j];
            }
            do {
                ur = potential_r(i, j, r_min[i][j], T, rho);
                urr = potential_rr(i, j, r_min[i][j], T, rho);
                r_min[i][j] -= ur / urr;
            } while (abs(ur) * r_min[i][j] / eps[i][j] > 1e-10);
            eps_eff[i][j] = eps_eff[j][i] = - potential(i, j, r_min[i][j], T, rho);
        }
    }
}

void ExtSutherland::set_vdw_alpha(dual2 rho, dual2 T){
    // Assumes that C_eff, sigma_eff and eps_eff are already set. Use the public method if you want to compute vdw_alpha.
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = 0; j < Ncomps; j++){
            vdw_alpha[i][j] = 0.;
            for (size_t k = 0; k < nterms; k++){
                vdw_alpha[i][j] -= C_eff[k][i][j] * eps[i][j] * pow(sigma[i][j] / sigma_eff[i][j], lambda[k][i][j]) / (lambda[k][i][j] - 3.);
            }
            vdw_alpha[i][j] /= eps_eff[i][j];
        }
    }
}

vector2d ExtSutherland::get_b_max(double T){
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

// ---------------------------------------------------------------------------------- //
// -------------------------- Radial distribution function -------------------------- // 
// ---------------------------------------------------------------------------------- //

vector2d2 ExtSutherland::get_BH_diameters(dual2 rho, dual2 T){
    // 20-point Gauss-Legendre from 0.5 sigma to 1 sigma. Using constant value of 1 for integrand within (0, 0.5) sigma.
    vector2d2 d_BH(Ncomps, vector1d2(Ncomps, 0.0));
    dual2 beta = 1. / (BOLTZMANN * T);
    for (size_t i = 0; i < Ncomps; i++){
        d_BH[i][i] = 2.;
        for (size_t n = 0; n < 20; n++){
            d_BH[i][i] += rdf_constants.gl_w[n] * (1. - exp(- beta * potential(i, i, sigma_eff[i][i] * (rdf_constants.gl_x[n] / 4. + 3. / 4.), T, rho)));
        }
        d_BH[i][i] *= sigma_eff[i][i] / 4.;
    }
    for (size_t i = 0; i < Ncomps - 1; i++){
        for (size_t j = i + 1; j < Ncomps; j++){
            d_BH[i][j] = (d_BH[i][i] + d_BH[j][j]) / 2.0;
            d_BH[j][i] = d_BH[i][j];
        }
    }
    return d_BH;
}

vector2d ExtSutherland::get_BH_diameters(double rho, double T){
    set_effective_params(rho, T);
    vector2d2 dBH_dual = get_BH_diameters(static_cast<dual2>(rho), static_cast<dual2>(T));
    return dual_to_double(dBH_dual);
}

vector2d ExtSutherland::get_BH_diameters(double T){
    // 20-point Gauss-Legendre from 0.5 sigma to 1 sigma. Using constant value of 1 for integrand within (0, 0.5) sigma.
    vector2d d_BH(Ncomps, vector1d(Ncomps, 0.0));
    double beta = 1. / (BOLTZMANN * T);
    for (int i = 0; i < Ncomps; i++){
        d_BH[i][i] = 2.;
        for (int n = 0; n < 20; n++){
            d_BH[i][i] += rdf_constants.gl_w[n] * (1. - exp(- beta * potential(i, i, static_cast<double>(sigma_eff[i][i]) * (rdf_constants.gl_x[n] / 4. + 3. / 4.))));
        }
        d_BH[i][i] *= static_cast<double>(sigma_eff[i][i]) / 4.;
    }
    for (int i = 0; i < Ncomps - 1; i++){
        for (int j = i + 1; j < Ncomps; j++){
            d_BH[i][j] = (d_BH[i][i] + d_BH[j][j]) / 2.0;
            d_BH[j][i] = d_BH[i][j];
        }
    }
    return d_BH;
}

vector2d ExtSutherland::saft_rdf(double rho_d, double T_d, const vector1d& x, int order, bool g2_correction){
    dual2 T = T_d;
    dual2 rho = rho_d;
    set_effective_params(rho, T);
    dual2 beta = (1. / (BOLTZMANN * T));
    vector2d2 d_BH = get_BH_diameters(rho, T);
    vector2d2 x_eff = get_xeff(d_BH);
    vector2d2 g0 = rdf_g0_func(rho, x, d_BH, x_eff);
    vector2d1 g1 = (order > 0) ? rdf_g1_func(rho, T, x, d_BH, x_eff) : vector2d1(Ncomps, vector1d1(Ncomps, 0.0));
    vector2d1 g2 = (order > 1) ? rdf_g2_func(rho, T, x, d_BH, x_eff, g2_correction) : vector2d1(Ncomps, vector1d1(Ncomps, 0.0));
    vector2d g(Ncomps, vector1d(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            g[i][j] = static_cast<double>(static_cast<dual2>(g0[i][j] * exp(beta * eps[i][j] * (g1[i][j] / g0[i][j]) + pow(beta * eps[i][j], 2) * (g2[i][j] / g0[i][j])) ));
            g[j][i] = g[i][j];
        }
    }
    return g;
}

vector3d ExtSutherland::get_rdf_terms(double rho, double T, const vector1d& x){
    vector3d terms;
    terms.push_back(saft_rdf(rho, T, x, 0));
    terms.push_back(saft_rdf(rho, T, x, 1));
    terms.push_back(saft_rdf(rho, T, x, 2, false));
    terms.push_back(saft_rdf(rho, T, x, 2, true));
    return terms;
}

vector2d2 ExtSutherland::rdf_g0_func(dual2 rho, const vector1d& x, const vector2d2& d_BH, const vector2d2& x_eff){
    dual2 zeta_x = zeta_x_func(rho, x, d_BH);
    dual2 k0 = - log(1 - zeta_x) + (42. * zeta_x - 39. * pow(zeta_x, 2) + 9. * pow(zeta_x, 3) - 2 * pow(zeta_x, 4)) / (6 * pow(1 - zeta_x, 3));
    dual2 k1 = (pow(zeta_x, 4) + 6. * pow(zeta_x, 2) - 12. * zeta_x) / (2. * pow(1 - zeta_x, 3));
    dual2 k2 = - 3. * pow(zeta_x, 2) / (8. * pow(1 - zeta_x, 2));
    dual2 k3 = (- pow(zeta_x, 4.) + 3. * pow(zeta_x, 2) + 3. * zeta_x) / (6. * pow(1 - zeta_x, 3));
    
    vector2d2 rdf(Ncomps, vector1d2(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            rdf[i][j] = exp(k0 + k1 * x_eff[i][j] + k2 * pow(x_eff[i][j], 2) + k3 * pow(x_eff[i][j], 3));
            rdf[j][i] = rdf[i][j];
        }
    }
    return rdf;
}

dual2 ExtSutherland::zeta_x_func(dual2 rho, const vector1d& x, const vector2d2& d_BH){
    dual2 zeta{0.0};
    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j < Ncomps; j++){
            zeta += x[i] * x[j] * pow(d_BH[i][j], 3);
        }
    }
    zeta *= PI * rho / 6.0;
    return zeta;
}

vector2d1 ExtSutherland::rdf_g1_func(dual2 rho, dual2 T, const vector1d& x, const vector2d2& d_BH, const vector2d2& x_eff){
    vector2d1 g1(Ncomps, vector1d1(Ncomps, 0.));
    vector2d2 x0 = get_x0(d_BH);
    dual2 zeta_x = zeta_x_func(rho, x, d_BH);

    for (size_t k = 0; k < nterms; k++){
        for (int i = 0; i < Ncomps; i++){
            for (int j = i; j < Ncomps; j++){
                dual2 a1s_k = a_1s_func(i, j, rho, x, zeta_x, d_BH, lambda[k]);
                dual2 B_k = B_func(i, j, rho, x, zeta_x, x_eff, d_BH, lambda[k]);
                dual2 term = lambda[k][i][j] * C_eff[k][i][j] * pow(x0[i][j], lambda[k][i][j]) * (a1s_k + B_k) / rho;
                g1[i][j] += term.val;
            }
        }
    }

    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = i; j < Ncomps; j++){
            auto a1 = [&](dual2 rho_){return a1_func(i, j, rho_, T, x);};
            g1[i][j] += 3. * derivative(a1, wrt(rho), at(rho));
            g1[i][j] *= (1.0 / (2. * PI * eps[i][j] * pow(d_BH[i][j].val, 3)));
            g1[j][i] = g1[i][j];
        }
    }
    return g1;
}

dual2 ExtSutherland::a_1s_func(int i, int j, dual2 rho, const vector1d& x, dual2 zeta_x, const vector2d2& d_BH, const vector2d& lambda_k){
    dual2 zeta_eff;
    dual2 a1s;
    zeta_eff = zeta_eff_func(rho, x, zeta_x, lambda_k[i][j]);
    a1s = (-2. * rho * PI * eps[i][j] * pow(d_BH[i][j], 3) / (lambda_k[i][j] - 3)) * (1. - zeta_eff/2.) / pow(1 - zeta_eff, 3);
    return a1s;
}

dual2 ExtSutherland::zeta_eff_func(dual2 rho,  const vector1d& x, dual2 zeta_x, double lambdaijk){
    vector1d c_coeffs(4, 0.0);
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            c_coeffs[i] += rdf_constants.C_coeff_matr[i][j] * pow(lambdaijk, - j);
        }
    }
    dual2 zeta_eff{0.0};
    for (int i = 0; i < 4; i++){
        zeta_eff += c_coeffs[i] * pow(zeta_x, i + 1);
    }
    return zeta_eff;
}

dual2 ExtSutherland::B_func(int i, int j, dual2 rho, const vector1d& x, dual2 zeta_x, const vector2d2& x_eff, const vector2d2& d_BH, const vector2d& lambda_k){
    dual2 I = I_func(i, j, x_eff, lambda_k);
    dual2 J = J_func(i, j, x_eff, lambda_k);
    dual2 B = 2. * PI * rho * pow(d_BH[i][j], 3) * eps[i][j] 
                * (I * (1. - zeta_x / 2.) / pow(1. - zeta_x, 3)
                    - J * 9. * zeta_x * (1. + zeta_x) / (2. * pow(1. - zeta_x, 3)));
    return B;
}

dual2 ExtSutherland::I_func(int i, int j, const vector2d2& xeff, const vector2d& lambda_k){
    dual2 I = - (pow(xeff[i][j], 3. - lambda_k[i][j]) - 1.) / (lambda_k[i][j] - 3.);
    return I;
}

dual2 ExtSutherland::J_func(int i, int j, const vector2d2& xeff, const vector2d& lambda_k){
    dual2 J = - (pow(xeff[i][j], 4. - lambda_k[i][j]) * (lambda_k[i][j] - 3.) 
                    - pow(xeff[i][j], 3. - lambda_k[i][j]) * (lambda_k[i][j] - 4.) - 1.)
                / ((lambda_k[i][j] - 3.) * (lambda_k[i][j] - 4.));
    return J;
}

dual2 ExtSutherland::a1_func(int i, int j, dual2 rho, dual2 T, const vector1d& x){
    vector2d2 d_BH = get_BH_diameters(rho, T);
    return a1_func(i, j, rho, x, d_BH);
}

dual2 ExtSutherland::a1_func(int i, int j, dual2 rho, const vector1d& x, const vector2d2& d_BH){
    dual2 a1 = 0.;
    vector2d2 x_eff = get_xeff(d_BH);
    vector2d2 x0 = get_x0(d_BH);
    dual2 zeta_x = zeta_x_func(rho, x, d_BH);
    for (size_t k = 0; k < nterms; k++){
        dual2 a1s_k = a_1s_func(i, j, rho, x, zeta_x, d_BH, lambda[k]);
        dual2 B_k = B_func(i, j, rho, x, zeta_x, x_eff, d_BH, lambda[k]);
        a1 -= C_eff[k][i][j] * pow(x0[i][j], lambda[k][i][j]) * (a1s_k + B_k);
    }
    return a1;
}

vector2d1 ExtSutherland::rdf_g2_func(dual2 rho, dual2 T, const vector1d& x, const vector2d2& d_BH, 
                                    const vector2d2& x_eff, bool g2_correction){
    vector2d1 g2(Ncomps, vector1d1(Ncomps));

    const vector2d2 x0 = get_x0(d_BH);
    dual2 zeta_x = zeta_x_func(rho, x, d_BH);
    dual2 K_HS = K_HS_func(zeta_x);
    for (size_t k = 0; k < nterms; k++){
        for (size_t l = 0; l < nterms; l++){
            vector2d lambda_kl = get_lambda_kl(k, l);
            for (size_t i = 0; i < Ncomps; i++){
                for (size_t j = i; j < Ncomps; j++){
                    dual2 a1s_kl = a_1s_func(i, j, rho, x, zeta_x, d_BH, lambda_kl);
                    dual2 B_kl = B_func(i, j, rho, x, zeta_x, x_eff, d_BH, lambda_kl);
                    dual2 term = K_HS * eps[i][j] * C_eff[k][i][j] * C_eff[l][i][j] * lambda[k][i][j] * pow(x0[i][j], lambda_kl[i][j]) * (a1s_kl + B_kl) / rho;
                    g2[i][j] -= term.val;
                }
            }
        }
    }
    // vector2d da2ij_div_chi_drho = da2ij_div_chi_drho_func(rho, x, K_HS, d_BH, x_eff);
    dual2 zeta_x_HS = zeta_x_func(rho, x, sigma_eff);
    vector2d2 gamma_c = g2_correction ? gamma_corr(zeta_x_HS, T) : vector2d2(Ncomps, vector1d2(Ncomps, 0.0));
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = i; j < Ncomps; j++){
            auto a2_div_chi = [&](dual2 rho_){return a2_div_chi_func(i, j, rho_, T, x);};
            g2[i][j] += 3. * derivative(a2_div_chi, wrt(rho), at(rho));
            dual2 mul_factor = (1 + gamma_c[i][j]) * (1. / (2. * PI * pow(eps[i][j], 2) * pow(d_BH[i][j], 3)));
            g2[i][j] *= mul_factor.val;
            g2[j][i] = g2[i][j];
        }
    }
    return g2;
}

vector2d2 ExtSutherland::gamma_corr(dual2 zeta_x, dual2 T){
    vector2d2 gamma(Ncomps, vector1d2(Ncomps));
    constexpr double phi[5] = {10., 10., 0.57, -6.7, -8.};
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            dual2 theta = exp(eps_eff[i][j] / (BOLTZMANN * T)) - 1.;
            gamma[i][j] = phi[0] * (1 - tanh(phi[1] * (phi[2] - vdw_alpha[i][j]))) 
                        * zeta_x * theta * exp(phi[3] * zeta_x + phi[4] * pow(zeta_x, 2));
            gamma[j][i] = gamma[i][j];
        }
    }
    return gamma;
}

vector2d ExtSutherland::get_lambda_kl(size_t k, size_t l){
    vector2d lambda_kl(Ncomps, vector1d(Ncomps, 0.));
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = 0; j < Ncomps; j++){
            lambda_kl[j][i] = lambda_kl[i][j] = lambda[k][i][j] + lambda[l][i][j];
        }
    }
    return lambda_kl;
}

vector2d2 ExtSutherland::get_x0(const vector2d2& d_BH){
    vector2d2 x0(Ncomps, vector1d2(Ncomps, 0.));
    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j <= i; j++){
            x0[i][j] = sigma[i][j] / d_BH[i][j];
            x0[j][i] = x0[i][j];
        }
    }
    return x0;
}

vector2d2 ExtSutherland::get_xeff(const vector2d2& d_BH){
    vector2d2 x_eff(Ncomps, vector1d2(Ncomps, 0.));
    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j <= i; j++){
            x_eff[i][j] = sigma_eff[i][j] / d_BH[i][j];
            x_eff[j][i] = x_eff[i][j];
        }
    }
    return x_eff;
}

dual2 ExtSutherland::a2_div_chi_func(int i, int j, dual2 rho, dual2 T, const vector1d& x){
    vector2d2 d_BH = get_BH_diameters(rho, T);
    dual2 zeta_x = zeta_x_func(rho, x, d_BH);
    dual2 K_HS = K_HS_func(zeta_x);
    vector2d2 x_eff = get_xeff(d_BH);
    return a2_div_chi_func(i, j, rho, x, K_HS, d_BH, x_eff);
}

dual2 ExtSutherland::a2_div_chi_func(int i, int j, dual2 rho, const vector1d& x, dual2 K_HS, const vector2d2& d_BH, const vector2d2& x_eff){
    dual2 a2ij = 0.;
    dual2 zeta_x = zeta_x_func(rho, x, d_BH);
    const vector2d2 x0 = get_x0(d_BH);
    for (size_t k = 0; k < nterms; k++){
        for (size_t l = 0; l < nterms; l++){
            vector2d lambda_kl = get_lambda_kl(k, l);
            dual2 a1s = a_1s_func(i, j, rho, x, zeta_x, d_BH, lambda_kl);
            dual2 B = B_func(i, j, rho, x, zeta_x, x_eff, d_BH, lambda_kl);
            a2ij += C_eff[k][i][j] * C_eff[l][i][j] * pow(x0[i][j], lambda_kl[i][j]) * (a1s + B);
        }
    }
    a2ij *= 0.5 * K_HS * eps[i][j];
    return a2ij;
}
vector2d2 ExtSutherland::rdf_chi_func(dual2 rho, dual2 T, const vector1d& x){
    vector2d2 rdf_chi(Ncomps, vector1d2(Ncomps));
    dual2 zeta_x = zeta_x_func(rho, x, sigma_eff);
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            vector1d2 f = f_corr(vdw_alpha[i][j]);
            rdf_chi[i][j] = f[0] * zeta_x + f[1] * pow(zeta_x, 5) + f[2] * pow(zeta_x, 8);
            rdf_chi[j][i] = rdf_chi[i][j];
        }
    }
    return rdf_chi;
}

vector1d2 ExtSutherland::f_corr(dual2 alpha){
    vector1d2 f(3);
    for (int i = 0; i < 3; i++){
        dual2 num{0.0}, denom{1.0};
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