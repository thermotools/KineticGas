/*
Author : Vegard Gjeldvik Jervell
Contains : Implementation of functions used to compute the radial distribution function at contact for a Mie fluid
See : J. Chem. Phys. 139, 154504 (2013); https://doi.org/10.1063/1.4819786
*/

#include "MieKinGas.h"

using namespace mie_rdf_constants;

double MieKinGas::omega(int i, int j, int l, int r, double T){
    if ((l <= 2) && (r >= l) && (r <= 3) && (abs(la[i][j] - 6.) < 1e-10)){ // Use Correlation by Fokin, Popov and Kalashnikov, High Temperature, Vol. 37, No. 1 (1999)
        // NOTE: There is a typo in Eq. (4b) in the article, ln(1/m) should be ln(m).
        // See: I. H. Bell et. al. J. Chem. Eng. Data (2020), https://doi.org/10.1021/acs.jced.9b00455
        // The correlation implemented here has been verified vs. tabulated data.

        // The correlation gives the logarithm of the reduced collision integral.
        // The collision integral is reduced using the hard-sphere value, i.e. lnomega_star = log(omega / omega_hs)
        double T_star = T * BOLTZMANN / eps[i][j];
        if (T_star > 0.4) return omega_correlation(i, j, l, r, T_star) * omega_hs(i, j, l, r, T);
    }
    return Spherical::omega(i, j, l, r, T);
}

double MieKinGas::omega_correlation(int i, int j, int l, int r, double T_star){
    // Correlation by Fokin, Popov and Kalashnikov, High Temperature, Vol. 37, No. 1 (1999)
    // NOTE: There is a typo in Eq. (4b) in the article, ln(1/m) should be ln(m).
    // See: I. H. Bell et. al. J. Chem. Eng. Data (2020), https://doi.org/10.1021/acs.jced.9b00455
    // The correlation implemented here has been verified vs. tabulated data.

    // The correlation gives the logarithm of the reduced collision integral.
    // The collision integral is reduced using the hard-sphere value, i.e. lnomega_star = log(omega / omega_hs)
    if (l == r){
        double lnomega_star = - (2. / lr[i][j]) * log(T_star);
        double a_m;
        for (int n = 1; n < 7; n++){
            a_m = omega_correlation_factors[l - 1][n - 1][0]
                    + (omega_correlation_factors[l - 1][n - 1][1] / lr[i][j])
                    + pow(lr[i][j], -2.) * (omega_correlation_factors[l - 1][n - 1][2]
                                             + omega_correlation_factors[l - 1][n - 1][3] * log(lr[i][j]));

            lnomega_star += a_m * pow(T_star, (1. - n) * 0.5);
        }
        if (l == 2) { lnomega_star += log(1. - (2. / (3. * lr[i][j]))); }
        return exp(lnomega_star);
    }
    return omega_recursive_factor(i, j, l, r - 1, T_star) * omega_correlation(i, j, l, r - 1, T_star);
}

double MieKinGas::omega_recursive_factor(int i, int j, int l, int r, double T_star){
    // Higher order collision integrals can be computed from the derivative of lower order ones using the recursion
    // Given in Fokin, Popov and Kalashnikov, High Temperature, Vol. 37, No. 1 (1999)
    // See also: Hirchfelder, Curtiss & Bird, Molecular Theory of Gases and Liquids.
    // For reduced integrals : omega(l, r + 1) / omega(l, r) = 1 + (d omega(l, r) / d lnT^*) / (r + 2)
    if ((l > 2) || (r > 3)) {throw std::runtime_error("No recursive factor available!");}
    double dlnomega_dlnT = - (2 / lr[i][j]);
    double a_m;
    for (int n = 2; n < 7; n++){
        a_m = omega_correlation_factors[l - 1][n - 1][0] + (omega_correlation_factors[l - 1][n - 1][1] / lr[i][j])
                + pow(lr[i][j], -2.) * (omega_correlation_factors[l - 1][n - 1][2]
                                        + omega_correlation_factors[l - 1][n - 1][3] * log(lr[i][j]));
        dlnomega_dlnT += a_m * pow(T_star, - (n - 1.) / 2.) * (1. - n) / 2.;
    }
    return 1. + dlnomega_dlnT / (r + 2);
}

std::vector<std::vector<double>> MieKinGas::saft_rdf(double rho, double T, const std::vector<double>& x, int order, bool g2_correction){
        double beta = (1. / (BOLTZMANN * T));
        std::vector<std::vector<double>> d_BH = get_BH_diameters(T);
        std::vector<std::vector<double>> x0 = get_x0(d_BH);
        std::vector<std::vector<double>> gHS = rdf_HS(rho, x, d_BH);
        std::vector<std::vector<double>> g1 = (order > 0) ? rdf_g1_func(rho, x, d_BH) : std::vector<std::vector<double>>(Ncomps, std::vector<double>(Ncomps, 0.));
        std::vector<std::vector<double>> g2 = (order > 1) ? rdf_g2_func(rho, T, x, d_BH, x0, g2_correction) : std::vector<std::vector<double>>(Ncomps, std::vector<double>(Ncomps, 0.));
        std::vector<std::vector<double>> g(Ncomps, std::vector<double>(Ncomps, 0.0));
        for (int i = 0; i < Ncomps; i++){
            for (int j = i; j < Ncomps; j++){
                g[i][j] = gHS[i][j] * exp(beta * eps[i][j] * (g1[i][j] / gHS[i][j]) + pow(beta * eps[i][j], 2) * (g2[i][j] / gHS[i][j]));
                g[j][i] = g[i][j];
            }
        }
        return g;
    }
std::vector<std::vector<double>> MieKinGas::rdf_HS(double rho, const std::vector<double>& x, const std::vector<std::vector<double>>& d_BH){
    double zeta = zeta_x_func(rho, x, d_BH);
    double k0 = - log(1 - zeta) + (42. * zeta - 39. * pow(zeta, 2) + 9. * pow(zeta, 3) - 2 * pow(zeta, 4)) / (6 * pow(1 - zeta, 3));
    double k1 = (pow(zeta, 4) + 6. * pow(zeta, 2) - 12. * zeta) / (2. * pow(1 - zeta, 3));
    double k2 = - 3. * pow(zeta, 2) / (8. * pow(1 - zeta, 2));
    double k3 = (- pow(zeta, 4.) + 3. * pow(zeta, 2) + 3. * zeta) / (6. * pow(1 - zeta, 3));
    
    double x0;
    std::vector<std::vector<double>> rdf(Ncomps, std::vector<double>(Ncomps));

    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j <= i; j++){
            x0 = sigma[i][j] / d_BH[i][j];
            rdf[i][j] = exp(k0 + k1 * x0 + k2 * pow(x0, 2) + k3 * pow(x0, 3));
            rdf[j][i] = rdf[i][j];
        }
    }
    return rdf;
}

std::vector<std::vector<double>> MieKinGas::rdf_HS(double rho, double T, const std::vector<double>& x){
    std::vector<std::vector<double>> d_BH = get_BH_diameters(T);
    return rdf_HS(rho, x, d_BH);
}

std::vector<std::vector<double>> MieKinGas::rdf_g1_func(double rho, const std::vector<double>& x,
                                                        const std::vector<std::vector<double>>& d_BH)
{
    std::vector<std::vector<double>> g1(Ncomps, std::vector<double>(Ncomps));
    std::vector<std::vector<double>> dadr = da1ij_drho_func(rho, x, d_BH);
    std::vector<std::vector<double>> x0 = get_x0(d_BH);
    std::vector<std::vector<double>> a1s_la = a_1s_func(rho, x, d_BH, la);
    std::vector<std::vector<double>> a1s_lr = a_1s_func(rho, x, d_BH, lr);
    std::vector<std::vector<double>> B_la = B_func(rho, x, d_BH, la);
    std::vector<std::vector<double>> B_lr = B_func(rho, x, d_BH, lr);
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            g1[i][j] = (1.0 / (2. * PI * eps[i][j] * pow(d_BH[i][j], 3))) 
                        * (3. * dadr[i][j] - C[i][j] * la[i][j] * pow(x0[i][j], la[i][j]) * (a1s_la[i][j] + B_la[i][j]) / rho
                                            + C[i][j] * lr[i][j] * pow(x0[i][j], lr[i][j]) * (a1s_lr[i][j] + B_lr[i][j]) / rho);
            g1[j][i] = g1[i][j];
        }
    }
    return g1;
}

std::vector<std::vector<double>> MieKinGas::rdf_g1_func(double rho, double T, const std::vector<double>& x){
    std::vector<std::vector<double>> d_BH = get_BH_diameters(T);
    return rdf_g1_func(rho, x, d_BH);
}

std::vector<std::vector<double>> MieKinGas::rdf_g2_func(double rho, double T, const std::vector<double>& x,
                                                        const std::vector<std::vector<double>>& d_BH,
                                                        const std::vector<std::vector<double>>& x0,
                                                        bool g2_correction)
{
    std::vector<std::vector<double>> g2(Ncomps, std::vector<double>(Ncomps));
    std::vector<std::vector<double>> _la(Ncomps, std::vector<double>(Ncomps));
    std::vector<std::vector<double>> _lr(Ncomps, std::vector<double>(Ncomps));
    std::vector<std::vector<double>> la_lr(Ncomps, std::vector<double>(Ncomps));
    
    double zeta_x = zeta_x_func(rho, x, d_BH);
    double K_HS = K_HS_func(zeta_x);
    std::vector<std::vector<double>> rdf_chi_HS = rdf_chi_func(rho, x);
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            _la[i][j] = _la[j][i] = 2 * la[i][j];
            _lr[i][j] = _lr[j][i] = 2 * lr[i][j];
            la_lr[i][j] = la_lr[j][i] = la[i][j] + lr[i][j];
        }
    }
    std::vector<std::vector<double>> a1s_la = a_1s_func(rho, x, d_BH, _la);
    std::vector<std::vector<double>> a1s_lr = a_1s_func(rho, x, d_BH, _lr);
    std::vector<std::vector<double>> a1s_la_lr = a_1s_func(rho, x, d_BH, la_lr);
    std::vector<std::vector<double>> B_la = B_func(rho, x, d_BH, _la);
    std::vector<std::vector<double>> B_lr = B_func(rho, x, d_BH, _lr);
    std::vector<std::vector<double>> B_la_lr = B_func(rho, x, d_BH, la_lr);
    std::vector<std::vector<double>> da2ij_div_chi_drho = da2ij_div_chi_drho_func(rho, x, K_HS, d_BH, x0);
    std::vector<std::vector<double>> a2ij = a2ij_func(rho, x, K_HS, rdf_chi_HS, d_BH, x0);
    double zeta_x_HS = zeta_x_func(rho, x, sigma);
    std::vector<std::vector<double>> gamma_c = g2_correction ? gamma_corr(zeta_x_HS, T) : std::vector<std::vector<double>>(Ncomps, std::vector<double>(Ncomps, 0.));
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

std::vector<std::vector<double>> MieKinGas::rdf_g2_func(double rho, double T, const std::vector<double>& x, bool g2_correction){
    std::vector<std::vector<double>> d_BH = get_BH_diameters(T);
    std::vector<std::vector<double>> x0 = get_x0(d_BH);
    return rdf_g2_func(rho, T, x, d_BH, x0, g2_correction);
}

std::vector<std::vector<double>> MieKinGas::get_b_max(double T){
        // bmax[i][j] in units of sigma[i][j]
        // The maximum value of the impact parameter at which deflection angle (chi) is positive
        const double g_avg = sqrt(4. / PI); // Average dimensionless relative velocity
        std::vector<std::vector<double>> bmax(Ncomps, std::vector<double>(Ncomps, 0));
        double b, db, chi_val;
        for (int i = 0; i < Ncomps; i++){
            for (int j = i; j < Ncomps; j++){
                b = 0.5;
                db = 0.1;
                chi_val = chi(i, j, T, g_avg, b * sigma[i][j]);
                while (abs(db) > 1e-5){
                    if ((chi_val < 0 && db > 0) || (chi_val > 0 && db < 0)){
                        db *= -0.5;
                    }
                    b += db;
                    chi_val = chi(i, j, T, g_avg, b * sigma[i][j]);
                    if (b > 10) {
                        std::printf("Could not compute contact diameter from deflection angle\n");
                        std::printf("Possibly because of insufficient integral resolution in theta integral\n");
                        std::printf("Returning BH diameters as fallback value for upper integration limit for contact diameter.\n");
                        bmax = get_BH_diameters(T);
                        for (int bi = 0; bi < Ncomps; bi++){
                            for (int bj = bi; bj < Ncomps; bj++){
                                bmax[bi][bj] /= sigma[bi][bj];
                                bmax[bj][bi] = bmax[bi][bj];
                            }
                        }
                        return bmax;
                    }
                }
                bmax[j][i] = bmax[i][j] = b;
            }
        }
        return bmax;
    }

std::vector<std::vector<double>> MieKinGas::get_BH_diameters(double T){
    // 20-point Gauss-Legendre from 0.5 sigma to 1 sigma. Using constant value of 1 for integrand within (0, 0.5) sigma.
    std::vector<std::vector<double>> d_BH(Ncomps, std::vector<double>(Ncomps, 0.0));
    double beta = 1. / (BOLTZMANN * T);
    for (int i = 0; i < Ncomps; i++){
        d_BH[i][i] = 2.;
        for (int n = 0; n < 20; n++){
            d_BH[i][i] += gl_w[n] * (1. - exp(- beta * potential(i, i, sigma[i][i] * (gl_x[n] / 4. + 3. / 4.))));
        }
        d_BH[i][i] *= sigma[i][i] / 4.;
    }
    for (int i = 0; i < Ncomps - 1; i++){
        for (int j = i + 1; j < Ncomps; j++){
            d_BH[i][j] = (d_BH[i][i] + d_BH[j][j]) / 2.0;
            d_BH[j][i] = d_BH[i][j];
        }
    }
    return d_BH;
}

std::vector<std::vector<double>> MieKinGas::get_x0(const std::vector<std::vector<double>>& d_BH){
    std::vector<std::vector<double>> x0(Ncomps, std::vector<double>(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j <= i; j++){
            x0[i][j] = sigma[i][j] / d_BH[i][j];
            x0[j][i] = x0[i][j];
        }
    }
    return x0;
}

std::vector<std::vector<double>> MieKinGas::a_1s_func(double rho,
                                                        const std::vector<double>& x,
                                                        const std::vector<std::vector<double>>& d_BH,
                                                        const std::vector<std::vector<double>>& lambda)
{
    double zeta_eff;
    std::vector<std::vector<double>> a1s(Ncomps, std::vector<double>(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            zeta_eff = zeta_eff_func(rho, x, d_BH, lambda[i][j]);
            a1s[i][j] = (-2. * rho * PI * eps[i][j] * pow(d_BH[i][j], 3) / (lambda[i][j] - 3))
                        * (1. - zeta_eff/2.) / pow(1. - zeta_eff, 3);
            a1s[j][i] = a1s[i][j];
        }
    }
    return a1s;
}

std::vector<std::vector<double>> MieKinGas::a_1s_func(double rho, double T, const std::vector<double>& x,
                                                const std::vector<std::vector<double>>& lambda)
{
    std::vector<std::vector<double>> d_BH = get_BH_diameters(T);
    return a_1s_func(rho, x, d_BH, lambda);
}

std::vector<std::vector<double>> MieKinGas::da1s_drho_func(double rho,
                                                            const std::vector<double>& x,
                                                            const std::vector<std::vector<double>>& d_BH,
                                                            const std::vector<std::vector<double>>& lambda)
{
    std::vector<std::vector<double>>  da1sdrho(Ncomps, std::vector<double>(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            double ze = zeta_eff_func(rho, x, d_BH, lambda[i][j]);
            double dzedrho = dzeta_eff_drho_func(rho, x, d_BH, lambda[i][j]);
            da1sdrho[i][j] = -2. * (PI * eps[i][j] * pow(d_BH[i][j], 3) / (lambda[i][j] - 3)) 
                                    * (((1. - ze / 2.) / pow(1 - ze, 3)) 
                                        + rho * dzedrho * (2.5 - ze) / pow(1 - ze, 4));
            da1sdrho[j][i] = da1sdrho[i][j];
        }
    }
    return da1sdrho;
}

std::vector<std::vector<double>> MieKinGas::da1s_drho_func(double rho, double T, const std::vector<double>& x,
                                                const std::vector<std::vector<double>>& lambda)
{
    std::vector<std::vector<double>> d_BH = get_BH_diameters(T);
    return da1s_drho_func(rho, x, d_BH, lambda);
}

std::vector<std::vector<double>> MieKinGas::I_func(const std::vector<std::vector<double>>& x0,
                                        const std::vector<std::vector<double>>& lambda)
{
    std::vector<std::vector<double>> I(Ncomps, std::vector<double>(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            I[i][j] = - (pow(x0[i][j], 3. - lambda[i][j]) - 1.) / (lambda[i][j] - 3.);
            I[j][i] = I[i][j];
        }
    }
    return I;
}

std::vector<std::vector<double>> MieKinGas::J_func(const std::vector<std::vector<double>>& x0,
                                        const std::vector<std::vector<double>>& lambda)
{
    std::vector<std::vector<double>> J(Ncomps, std::vector<double>(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            J[i][j] = - (pow(x0[i][j], 4. - lambda[i][j]) * (lambda[i][j] - 3.) 
                            - pow(x0[i][j], 3. - lambda[i][j]) * (lambda[i][j] - 4.) - 1.)
                        / ((lambda[i][j] - 3.) * (lambda[i][j] - 4.));
            J[j][i] = J[i][j];
        }
    }
    return J;
}

std::vector<std::vector<double>> MieKinGas::B_func(double rho, const std::vector<double>& x,
                                        const std::vector<std::vector<double>>& d_BH,
                                        const std::vector<std::vector<double>>& lambda)
{
    double zeta_x = zeta_x_func(rho, x, d_BH);
    std::vector<std::vector<double>> x0 = get_x0(d_BH);
    std::vector<std::vector<double>> I{I_func(x0, lambda)}, J{J_func(x0, lambda)};
    std::vector<std::vector<double>> B(Ncomps, std::vector<double>(Ncomps));
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

std::vector<std::vector<double>> MieKinGas::B_func(double rho, double T, const std::vector<double>& x,
                                                    const std::vector<std::vector<double>>& lambda)
{
    std::vector<std::vector<double>> d_BH = get_BH_diameters(T);
    return B_func(rho, x, d_BH, lambda);
}

std::vector<std::vector<double>> MieKinGas::dBdrho_func(double rho, const std::vector<double>& x,
                                        const std::vector<std::vector<double>>& d_BH,
                                        const std::vector<std::vector<double>>& lambda)
{
    std::vector<std::vector<double>> B = B_func(rho, x, d_BH, lambda);
    std::vector<std::vector<double>> x0 = get_x0(d_BH);
    std::vector<std::vector<double>> I = I_func(x0, lambda);
    std::vector<std::vector<double>> J = J_func(x0, lambda);
    std::vector<std::vector<double>> dBdrho(Ncomps, std::vector<double>(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            double dzxdrho = dzetax_drho_func(x, d_BH);
            double zx = zeta_x_func(rho, x, d_BH);
            dBdrho[i][j] = B[i][j] / rho
                            + 2. * PI * rho * eps[i][j] * pow(d_BH[i][j], 3)
                            * (I[i][j] * (- 0.5 * (1 - zx) + 3. * (1. - 0.5 * zx)) / pow(1 - zx, 4)
                                - J[i][j] * (9. * (1. + 2. * zx) * (1. - zx) + 27. * zx * (1. + zx)) / (2. * pow(1 - zx, 4))
                                ) * dzxdrho;
        }
    }
    return dBdrho;
}

std::vector<std::vector<double>> MieKinGas::dBdrho_func(double rho, double T, const std::vector<double>& x,
                                                    const std::vector<std::vector<double>>& lambda)
{
    std::vector<std::vector<double>> d_BH = get_BH_diameters(T);
    return dBdrho_func(rho, x, d_BH, lambda);
}

std::vector<std::vector<double>> MieKinGas::a1ij_func(double rho, double T, const std::vector<double>& x)
{
    std::vector<std::vector<double>> d_BH = get_BH_diameters(T);
    std::vector<std::vector<double>> x0 = get_x0(d_BH);
    std::vector<std::vector<double>> a1s_la = a_1s_func(rho, x, d_BH, la);
    std::vector<std::vector<double>> a1s_lr = a_1s_func(rho, x, d_BH, lr);
    std::vector<std::vector<double>> B_la = B_func(rho, x, d_BH, la);
    std::vector<std::vector<double>> B_lr = B_func(rho, x, d_BH, lr);

    std::vector<std::vector<double>> a1ij(Ncomps, std::vector<double>(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            a1ij[i][j] = C[i][j] * (pow(x0[i][j], la[i][j]) * (a1s_la[i][j] + B_la[i][j])
                                    - pow(x0[i][j], lr[i][j]) * (a1s_lr[i][j] + B_lr[i][j]));
            a1ij[j][i] = a1ij[i][j];
        }
    }
    return a1ij;
}

std::vector<std::vector<double>> MieKinGas::da1ij_drho_func(double rho, const std::vector<double>& x,
                                const std::vector<std::vector<double>>& d_BH)
{
    std::vector<std::vector<double>> da1drho(Ncomps, std::vector<double>(Ncomps));
    std::vector<std::vector<double>> x0 = get_x0(d_BH);
    std::vector<std::vector<double>> da1sdr_la = da1s_drho_func(rho, x, d_BH, la);
    std::vector<std::vector<double>> da1sdr_lr = da1s_drho_func(rho, x, d_BH, lr);
    std::vector<std::vector<double>> dBdr_la = dBdrho_func(rho, x, d_BH, la);
    std::vector<std::vector<double>> dBdr_lr = dBdrho_func(rho, x, d_BH, lr);
    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j < Ncomps; j++){
            da1drho[i][j] += C[i][j] * (pow(x0[i][j], la[i][j]) * (da1sdr_la[i][j] + dBdr_la[i][j])
                                            - pow(x0[i][j], lr[i][j]) * (da1sdr_lr[i][j] + dBdr_lr[i][j]));
        }
    }
    return da1drho;
}

std::vector<std::vector<double>> MieKinGas::da1ij_drho_func(double rho, double T, const std::vector<double>& x)
{
    std::vector<std::vector<double>> d_BH = get_BH_diameters(T);
    return da1ij_drho_func(rho, x, d_BH);
}

std::vector<std::vector<double>> MieKinGas::rdf_chi_func(double rho, const std::vector<double>& x)
{
    std::vector<std::vector<double>> rdf_chi(Ncomps, std::vector<double>(Ncomps));
    double zeta_x = zeta_x_func(rho, x, sigma);
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            std::vector<double> f = f_corr(alpha[i][j]);
            rdf_chi[i][j] = f[0] * zeta_x + f[1] * pow(zeta_x, 5) + f[2] * pow(zeta_x, 8);
            rdf_chi[j][i] = rdf_chi[i][j];
        }
    }
    return rdf_chi;
}

std::vector<std::vector<double>> MieKinGas::drdf_chi_drho_func(double rho, const std::vector<double>& x)
{
    std::vector<std::vector<double>> drdf_chi_drho(Ncomps, std::vector<double>(Ncomps));
    double zeta_x = zeta_x_func(rho, x, sigma);
    double dzdr = dzetax_drho_func(x, sigma);
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            std::vector<double> f = f_corr(alpha[i][j]);
            drdf_chi_drho[i][j] = dzdr * (f[0] + 5. * f[1] * pow(zeta_x, 4) + 8. * f[2] * pow(zeta_x, 7));
            drdf_chi_drho[j][i] = drdf_chi_drho[i][j]; 
        }
    }
    return drdf_chi_drho;
}

std::vector<std::vector<double>> MieKinGas::a2ij_func(double rho, const std::vector<double>& x, double K_HS, 
                                                const std::vector<std::vector<double>>& rdf_chi_HS,
                                                const std::vector<std::vector<double>>& d_BH,
                                                const std::vector<std::vector<double>>& x0)
{
    std::vector<std::vector<double>> _la(Ncomps, std::vector<double>(Ncomps)), 
                                    _lr(Ncomps, std::vector<double>(Ncomps)), 
                                    la_lr(Ncomps, std::vector<double>(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            _la[i][j] = _la[j][i] = 2 * la[i][j];
            _lr[i][j] = _lr[j][i] = 2 * lr[i][j];
            la_lr[i][j] = la_lr[j][i] = la[i][j] + lr[i][j];
        }
    }
    std::vector<std::vector<double>> a1s_la = a_1s_func(rho, x, d_BH, _la);
    std::vector<std::vector<double>> a1s_lr = a_1s_func(rho, x, d_BH, _lr);
    std::vector<std::vector<double>> a1s_la_lr = a_1s_func(rho, x, d_BH, la_lr);
    std::vector<std::vector<double>> B_la = B_func(rho, x, d_BH, _la);
    std::vector<std::vector<double>> B_lr = B_func(rho, x, d_BH, _lr);
    std::vector<std::vector<double>> B_la_lr = B_func(rho, x, d_BH, la_lr);
    std::vector<std::vector<double>> a2ij(Ncomps, std::vector<double>(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            a2ij[i][j] = .5 * K_HS * (1 + rdf_chi_HS[i][j]) * eps[i][j] * pow(C[i][j], 2)
                        * (pow(x0[i][j], _la[i][j]) * (a1s_la[i][j] + B_la[i][j]) 
                            - 2 * pow(x0[i][j], la_lr[i][j]) * (a1s_la_lr[i][j] + B_la_lr[i][j])
                            + pow(x0[i][j], _lr[i][j]) * (a1s_lr[i][j] + B_lr[i][j]));
            a2ij[j][i] = a2ij[i][j];
        }
    }
    return a2ij;
}
std::vector<std::vector<double>> MieKinGas::a2ij_func(double rho, double T, const std::vector<double>& x)
{
    std::vector<std::vector<double>> d_BH = get_BH_diameters(T);
    std::vector<std::vector<double>> x0 = get_x0(d_BH);
    std::vector<std::vector<double>> rdf_chi_HS = rdf_chi_func(rho, x);
    double zeta_x = zeta_x_func(rho, x, d_BH);
    double K_HS = K_HS_func(zeta_x);
    return a2ij_func(rho, x, K_HS, rdf_chi_HS, d_BH, x0);
}


std::vector<std::vector<double>> MieKinGas::da2ij_drho_func(double rho, const std::vector<double>& x, double K_HS, 
                                                const std::vector<std::vector<double>>& d_BH,
                                                const std::vector<std::vector<double>>& x0)
{
    std::vector<std::vector<double>> _la(Ncomps, std::vector<double>(Ncomps)), 
                                    _lr(Ncomps, std::vector<double>(Ncomps)), 
                                    la_lr(Ncomps, std::vector<double>(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            _la[i][j] = _la[j][i] = 2 * la[i][j];
            _lr[i][j] = _lr[j][i] = 2 * lr[i][j];
            la_lr[i][j] = la_lr[j][i] = la[i][j] + lr[i][j];
        }
    }
    double zeta_x = zeta_x_func(rho, x, d_BH);
    double dzxdr = dzetax_drho_func(x, d_BH);
    double dKHS_drho = dKHS_drho_func(zeta_x, dzxdr);
    std::vector<std::vector<double>> rdf_chi_HS = rdf_chi_func(rho, x);
    std::vector<std::vector<double>> a2ij = a2ij_func(rho, x, K_HS, rdf_chi_HS, d_BH, x0);
    std::vector<std::vector<double>> dchi_HS_drho = drdf_chi_drho_func(rho, x);
    std::vector<std::vector<double>> da1s_la = da1s_drho_func(rho, x, d_BH, _la);
    std::vector<std::vector<double>> da1s_lr = da1s_drho_func(rho, x, d_BH, _lr);
    std::vector<std::vector<double>> da1s_la_lr = da1s_drho_func(rho, x, d_BH, la_lr);
    std::vector<std::vector<double>> dB_la = dBdrho_func(rho, x, d_BH, _la);
    std::vector<std::vector<double>> dB_lr = dBdrho_func(rho, x, d_BH, _lr);
    std::vector<std::vector<double>> dB_la_lr = dBdrho_func(rho, x, d_BH, la_lr);
    std::vector<std::vector<double>> da2ij_drho(Ncomps, std::vector<double>(Ncomps));

    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            da2ij_drho[i][j] = dKHS_drho * (a2ij[i][j] / K_HS) + dchi_HS_drho[i][j] * (a2ij[i][j] / (1 + rdf_chi_HS[i][j]))
                        + .5 * K_HS * (1 + rdf_chi_HS[i][j]) * eps[i][j] * pow(C[i][j], 2)
                        * (pow(x0[i][j], _la[i][j]) * (da1s_la[i][j] + dB_la[i][j]) 
                            - 2 * pow(x0[i][j], la_lr[i][j]) * (da1s_la_lr[i][j] + dB_la_lr[i][j])
                            + pow(x0[i][j], _lr[i][j]) * (da1s_lr[i][j] + dB_lr[i][j]));
            da2ij_drho[j][i] = da2ij_drho[i][j];
        }
    }
    return da2ij_drho;

}

std::vector<std::vector<double>> MieKinGas::da2ij_div_chi_drho_func(double rho, const std::vector<double>& x, double K_HS, 
                                                const std::vector<std::vector<double>>& d_BH,
                                                const std::vector<std::vector<double>>& x0)
{
    std::vector<std::vector<double>> _la(Ncomps, std::vector<double>(Ncomps)), 
                                    _lr(Ncomps, std::vector<double>(Ncomps)), 
                                    la_lr(Ncomps, std::vector<double>(Ncomps));
    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            _la[i][j] = _la[j][i] = 2 * la[i][j];
            _lr[i][j] = _lr[j][i] = 2 * lr[i][j];
            la_lr[i][j] = la_lr[j][i] = la[i][j] + lr[i][j];
        }
    }
    double zeta_x = zeta_x_func(rho, x, d_BH);
    double dzxdr = dzetax_drho_func(x, d_BH);
    double dKHS_drho = dKHS_drho_func(zeta_x, dzxdr);
    std::vector<std::vector<double>> rdf_chi_HS = rdf_chi_func(rho, x);
    std::vector<std::vector<double>> a2ij = a2ij_func(rho, x, K_HS, rdf_chi_HS, d_BH, x0);
    std::vector<std::vector<double>> dchi_HS_drho = drdf_chi_drho_func(rho, x);
    std::vector<std::vector<double>> da1s_la = da1s_drho_func(rho, x, d_BH, _la);
    std::vector<std::vector<double>> da1s_lr = da1s_drho_func(rho, x, d_BH, _lr);
    std::vector<std::vector<double>> da1s_la_lr = da1s_drho_func(rho, x, d_BH, la_lr);
    std::vector<std::vector<double>> dB_la = dBdrho_func(rho, x, d_BH, _la);
    std::vector<std::vector<double>> dB_lr = dBdrho_func(rho, x, d_BH, _lr);
    std::vector<std::vector<double>> dB_la_lr = dBdrho_func(rho, x, d_BH, la_lr);
    std::vector<std::vector<double>> da2ij_drho(Ncomps, std::vector<double>(Ncomps));

    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            da2ij_drho[i][j] = dKHS_drho * (a2ij[i][j] / (K_HS * (1 + rdf_chi_HS[i][j])))
                            + .5 * K_HS * eps[i][j] * pow(C[i][j], 2)
                            * (pow(x0[i][j], _la[i][j]) * (da1s_la[i][j] + dB_la[i][j]) 
                                - 2 * pow(x0[i][j], la_lr[i][j]) * (da1s_la_lr[i][j] + dB_la_lr[i][j])
                                + pow(x0[i][j], _lr[i][j]) * (da1s_lr[i][j] + dB_lr[i][j]));
            da2ij_drho[j][i] = da2ij_drho[i][j];
        }
    }
    return da2ij_drho;
}

std::vector<std::vector<double>> MieKinGas::da2ij_drho_func(double rho, double T, const std::vector<double>& x)
{
    std::vector<std::vector<double>> d_BH = get_BH_diameters(T);
    std::vector<std::vector<double>> x0 = get_x0(d_BH);
    double zeta_x = zeta_x_func(rho, x, d_BH);
    double K_HS = K_HS_func(zeta_x);
    return da2ij_drho_func(rho, x, K_HS, d_BH, x0);
}

std::vector<std::vector<double>> MieKinGas::gamma_corr(double zeta_x, double T){
    std::vector<std::vector<double>> gamma(Ncomps, std::vector<double>(Ncomps));
    constexpr double phi[5] = {10., 10., 0.57, -6.7, -8.};

    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            double theta = exp(eps[i][j] / (BOLTZMANN * T)) - 1.;
            gamma[i][j] = phi[0] * (1 - tanh(phi[1] * (phi[2] - alpha[i][j]))) 
                        * zeta_x * theta * exp(phi[3] * zeta_x + phi[4] * pow(zeta_x, 2));
            gamma[j][i] = gamma[i][j];

        }
    }
    return gamma;
}

double MieKinGas::zeta_x_func(double rho,
                const std::vector<double>& x,
                const std::vector<std::vector<double>>& d_BH)
{
    double zeta{0.0};
    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j < Ncomps; j++){
            zeta += x[i] * x[j] * pow(d_BH[i][j], 3);
        }
    }
    zeta *= PI * rho / 6.0;
    return zeta;
}

double MieKinGas::dzetax_drho_func(const std::vector<double>& x,
                const std::vector<std::vector<double>>& d_BH)
{
    double dzdrho{0.0};
    for (int i = 0; i < x.size(); i++){
        for (int j = 0; j < x.size(); j++){
            dzdrho += x[i] * x[j] * pow(d_BH[i][j], 3);
        }
    }
    dzdrho *= PI / 6.0;
    return dzdrho;
}

double MieKinGas::zeta_eff_func(double rho,
                const std::vector<double>& x,
                const std::vector<std::vector<double>>& d_BH,
                double lambdaij)
{

    std::vector<double> c_coeffs(4, 0.0);
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            c_coeffs[i] += C_coeff_matr[i][j] * pow(lambdaij, - j);
        }
    }
    double zeta_x = zeta_x_func(rho, x, d_BH);
    double zeta_eff{0.0};
    for (int i = 0; i < 4; i++){
        zeta_eff += c_coeffs[i] * pow(zeta_x, i + 1);
    }
    return zeta_eff;
}

double MieKinGas::dzeta_eff_drho_func(double rho,
                const std::vector<double>& x,
                const std::vector<std::vector<double>>& d_BH,
                double lambdaij)
{

    std::vector<double> c_coeffs(4, 0.0);
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            c_coeffs[i] += C_coeff_matr[i][j] * pow(lambdaij, - j);
        }
    }
    double zeta_x = zeta_x_func(rho, x, d_BH);
    double dzdrho{0.0};
    for (int i = 0; i < 4; i++){
        dzdrho += c_coeffs[i] * (i + 1) * pow(zeta_x, i) * dzetax_drho_func(x, d_BH);
    }
    return dzdrho;
}

std::vector<double> MieKinGas::f_corr(double alpha){

    std::vector<double> f(3);
    for (int i = 0; i < 3; i++){
        double num{0.0}, denom{1.0};
        for (int n = 0; n <= 3; n++){
            num += phi[n][i] * pow(alpha, n);
        }
        for (int n = 4; n <= 6; n++){
            denom += phi[n][i] * pow(alpha, n - 3);
        }
        f[i] = num / denom;
    }
    return f;
}