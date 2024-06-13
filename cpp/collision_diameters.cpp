#include "Spherical.h"
#include "Integration/Integration.h"

vector2d Spherical::get_collision_diameters(double rho, double T, const std::vector<double>& x){
    int T_idx = static_cast<int>(T * 100 + 0.5); // Using T in centi kelvin for lookup.
    const std::map<int, vector2d>::iterator pos = collision_diameter_map.find(T_idx);
    vector2d cd;
    if (pos == collision_diameter_map.end()){
        switch (collision_diameter_model_id) {
        case 0:
            cd = get_collision_diameters_model0(rho, T, x); break;
        case 1:
        case 2:
        case 3:
            cd = get_collision_diameters_model1(rho, T, x); break;
        case 4: {
            double T_coeffs[4] = {0., 0., 0., 0.};
            constexpr double rho_coeffs[4][5] = {{48.1090, -68.2761, 18.4683, 9.1140, 1.2057},
                                                 {29.1735, -116.0489, 137.6770, -58.2916, 5.6222},
                                                 {-6.6618, 13.0363, -8.0000, 1.4533, 0.5971},
                                                 {-1.1381, 2.3210, -1.5264, 0.3194, 0.8324}
                                                };
            double rho_r;
            double Tr;
            cd = vector2d(Ncomps, vector1d(Ncomps, 0.));
            for (int i = 0; i < Ncomps; i++){
                for (int j = i; j < Ncomps; j++){
                    rho_r = rho * pow(sigma[i][j], 3);
                    Tr = T * BOLTZMANN / eps[i][j];
                    for (int ti = 0; ti < 4; ti++){
                        T_coeffs[ti] = 0.;
                        for (int ri = 0; ri < 5; ri++){
                            T_coeffs[ti] += rho_coeffs[ti][ri] * pow(rho_r, 4 - ri);
                        }
                    }
                    cd[i][j] = ((0.1 * T_coeffs[0] * Tr / pow(Tr - 0.1 * T_coeffs[1], 2)) - 1e-2 * T_coeffs[2] * Tr + T_coeffs[3]) * sigma[i][j];
                    cd[j][i] = cd[i][j];
                }
            }
            return cd;
            break;
            }
        default:
            throw std::runtime_error("Invalid collision diameter model id!"); // collision_diameter_model_id is private, so it shouldn't be possible to end up here.
        }
        collision_diameter_map[T_idx] = cd;
    }
    else{
        cd = pos->second;
    }
    return cd;
}

/**********************************************************************************************/
/***********************                      MODEL 0                  ************************/
/**********************************************************************************************/

vector2d Spherical::get_collision_diameters_model0(double rho, double T, const std::vector<double>& x){
    // Evaluate the integral of Eq. (40) in RET for Mie fluids (doi: 10.1063/5.0149865)
    // Note: For models with is_idealgas==true, zeros are returned, as there is no excluded volume at infinite dilution.
    // Weights and nodes for 6-point Gauss-Legendre quadrature
    constexpr int n_gl_points = 6;
    constexpr double gl_points[n_gl_points] = {-0.93246951, -0.66120939, -0.23861919, 0.23861919, 0.66120939, 0.93246951};
    constexpr double gl_weights[n_gl_points] = {0.17132449, 0.36076157, 0.46791393, 0.46791393, 0.36076157, 0.17132449};

    const double g_avg = sqrt(4. / PI); // Average dimensionless relative speed of collision
    std::vector<std::vector<double>> avg_R(Ncomps, std::vector<double>(Ncomps, 0.));
    if (is_idealgas) {
        return avg_R;
    }
    std::vector<std::vector<double>> bmax = get_b_max(T);
    double point, weight;

    for (int i = 0; i < Ncomps; i++){
        for (int j = i; j < Ncomps; j++){
            for (int p = 0; p < n_gl_points; p++){
                point = gl_points[p] * (bmax[i][j] / 2.) + bmax[i][j] / 2.;
                weight = gl_weights[p] * bmax[i][j] / 2.;
                avg_R[i][j] += get_R(i, j, T, g_avg, point * sigma[i][j]) * weight / bmax[i][j];
            }
            avg_R[j][i] = avg_R[i][j];
        }
    }
    return avg_R;
}

vector2d Spherical::get_b_max(double T, std::vector<std::vector<int>>& ierr){
        // bmax[i][j] in units of sigma[i][j]
        // The maximum value of the impact parameter at which deflection angle (chi) is positive
        const double g_avg = sqrt(4. / PI); // Average dimensionless relative velocity
        std::vector<std::vector<double>> bmax(Ncomps, std::vector<double>(Ncomps, 0));
        double b, db, chi_val;
        for (int i = 0; i < Ncomps; i++){
            for (int j = i; j < Ncomps; j++){
                ierr[i][j] = ierr[j][i] = 0;
                b = 0.5;
                db = 0.1;
                chi_val = chi(i, j, T, g_avg, b * sigma[i][j]);
                while (abs(db) > 1e-5){
                    if ((chi_val < 0 && db > 0) || (chi_val > 0 && db < 0)){
                        db *= -0.5;
                    }
                    b += db;
                    chi_val = chi(i, j, T, g_avg, b * sigma[i][j]);
                    if (b > 10) { // Unable to converge, this failure should be handled in get_b_max(double T)
                        ierr[i][j] = 1;
                        break;
                    }
                }
                bmax[j][i] = bmax[i][j] = b;
            }
        }
        return bmax;
    }

vector2d Spherical::get_b_max(double T){
    std::vector<std::vector<int>> ierr(Ncomps, std::vector<int>(Ncomps, 1));
    std::vector<std::vector<double>> b_max = get_b_max(T, ierr);
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = i; j < Ncomps; j++){
            if (ierr[i][j]) {
                std::cout << "Could not compute upper integration limit for collision diameter (" << i << ", " << j
                            << ") using sigma as fallback value." << std::endl;
                b_max[j][i] = b_max[i][j] = sigma[i][j];
            }
        }
    }
    return b_max;
}

/**********************************************************************************************/
/***********************                      MODEL 1                  ************************/
/**********************************************************************************************/

inline double dimless_relative_vdf(double g){
    return sqrt(2. / PI) * pow(g, 2) * exp(- 0.5 * pow(g, 2));
}

vector2d Spherical::get_collision_diameters_model1(double rho, double T, const std::vector<double>& x){
    if (is_singlecomp){
        const double cd = momentum_collision_diameter(0, 0, T);
        return vector2d(Ncomps, vector1d(Ncomps, cd));
    }
    vector2d cd(Ncomps, vector1d(Ncomps, 0.));
    for (int i = 0; i < Ncomps; i++){
        for (int j = 0; j < Ncomps; j++){
            cd[i][j] = momentum_collision_diameter(i, j, T);
        }
    }
    return cd;
}

double Spherical::momentum_collision_diameter(int i, int j, double T){
    const double I = get_cd_weight_normalizer(i, j, T);
    const double g0 = 0.;
    const double g1 = 3.5;
    const double gmax = 5.;
    const auto integrand = [&](double g){return cd_inner(i, j, T, g, I);};
    double cd = simpson(integrand, g0, g1, 15);
    cd += simpson(integrand, g1, gmax, 10);
    return cd;
}

double Spherical::cd_inner(int i, int j, double T, double g, double I){
    const double bmax = get_b_max_g(i, j, g, T);
    const double bmid = get_bmid(i, j, g, T);
    const double b0 = 0.;
    const auto integrand = [&](double b){return cd_integrand(i, j, T, g, b, I, bmax);};
    double val = simpson(integrand, b0, bmid, 15);
    val += simpson(integrand, bmid, bmax, 20);
    return val;
}

double Spherical::cd_weight_integrand(int i, int j, double T, double g, double b, double bmax){
    switch (collision_diameter_model_id) {
        case 1:
            return momentum_transfer(i, j, T, g, b) * dimless_relative_vdf(g) * g * 2 * PI * b;
        case 2:
            return momentum_transfer(i, j, T, g, b) * dimless_relative_vdf(g) * g * 2 * PI * b;
        case 3:
            return momentum_transfer(i, j, T, g, b) * dimless_relative_vdf(g) * g * 2 * PI * b;
        default:
            throw std::runtime_error("Invalid collision diameter model id!");
    }
}

double Spherical::get_cd_weight(int i, int j, double T, double g, double b, double I, double bmax){
    return cd_weight_integrand(i, j, T, g, b, bmax) / I;
}

double Spherical::get_cd_weight_normalizer(int i, int j, double T){
    double g0 = 0;
    double g1 = 3.5;
    double g_max = 5.;
    double I = 0.;
    switch (collision_diameter_model_id){
    case 1:
        // {
        // const auto integrand = [&](double g, double b){return cd_weight_integrand(i, j, T, g, b * sigma[i][j], 5.5 * sigma[i][j]);};
        // Point origin{1e-7, 1e-7};
        // Point end{5.5, 5.5};
        // double dg{0.5}, db{0.03125};
        // int refinement_levels_g{4};
        // int refinement_levels_b{16};
        // double subdomain_dblder_limit{1e-5};
//
        // I = integrate2d(origin, end,
        //                     dg, db,
        //                     refinement_levels_g, refinement_levels_b,
        //                     subdomain_dblder_limit,
        //                     integrand);
        // return I * sigma[i][j];
        // }
    case 2:
    case 3:
        {
        const auto integrand = [&](double g){return cd_weight_inner(i, j, T, g);};
        I = simpson(integrand, g0, g1, 15);
        I += simpson(integrand, g1, g_max, 10);
        break;
        }
    default:
        throw std::runtime_error("Invalid CD model!");
    }
    return I;
}

double Spherical::cd_weight_inner(int i, int j, double T, double g){
    double bmax = get_b_max_g(i, j, g, T);
    double bmid = get_bmid(i, j, g, T);
    const auto integrand = [&](double b){return cd_weight_integrand(i, j, T, g, b, bmax);};
    double val = simpson(integrand, 0., bmid, 15);
    val += simpson(integrand, bmid, bmax, 20);
    return val;
}

double Spherical::cd_integrand(int i, int j, double T, double g, double b, double I, double bmax){
    const double wt = get_cd_weight(i, j, T, g, b, I, bmax);
    double r;
    switch (collision_diameter_model_id){
        case 1:
            r = momentum_transfer_length(i, j, T, g, b * sigma[i][j]);
            if (isnan(r)){
                std::cout << "Model 1 gave NAN at " << T << ", " << g << ", " << b << std::endl;
                r = get_R(i, j, T, g, b * sigma[i][j]);
            }
            break;
        case 0:
        case 2:
        case 3:
            r = get_R(i, j, T, g, b * sigma[i][j]); break;
        default:
            throw std::runtime_error("Invalid collision diameter model!");
    }
    return wt * r;
}

double Spherical::momentum_transfer(int i, int j, double T, double g, double b){
    double chi_val = chi(i, j, T, g, b * sigma[i][j]);
    double red_mass = m[i] * m[j] / (m[i] + m[j]);
    double U = sqrt(2 * BOLTZMANN * T / red_mass) * g;
    double dp;
    switch (collision_diameter_model_id){
        case 1:
        case 2:
            dp = red_mass * U * sqrt(2 * (1 - cos(chi_val))) * abs(sin(chi_val / 2.));
            break;
        case 3:
        case 4:
            dp = red_mass * U * abs(cos(chi_val) - sin(chi_val) - 1); // Energy transfer
            break;
       default:
            throw std::runtime_error("Invalid collision diameter model!");
    }
    return dp;
}

double Spherical::get_b_max_g(int i, int j, double g, double T){
    double b, db, chi_val;
    if (g == 0) return 10.;
    b = 5.;
    db = -0.5;
    chi_val = chi(i, j, T, g, b * sigma[i][j]);
    while (abs(db) > 1e-3){
        double scale = g * sqrt(2 * (1 - cos(chi_val))) * abs(sin(chi_val / 2.));
        if (scale > 1e-3){
            b -= db;
            db *= 0.5;
        }
        b += db;
        chi_val = chi(i, j, T, g, b * sigma[i][j]);
    }
    return b;
}

double Spherical::get_bmid(int i, int j, double g, double T){
    if (g == 0) return 5.;
    double b = 0.5;
    double db = 0.1;
    double chi_val = chi(i, j, T, g, b * sigma[i][j]);
    while (abs(db) > 1e-5){
        if ((chi_val < 0 && db > 0) || (chi_val > 0 && db < 0)){
            db *= -0.5;
        }
        b += db;
        chi_val = chi(i, j, T, g, b * sigma[i][j]);
        if (b > 10) {
            std::cout << "bmid did not converge!" << std::endl;
            throw std::runtime_error("Unable to converge bmid!");
            break;
        }
    }
    return b;
}

double Spherical::momentum_transfer_length_weight(int i, int j, double r, double chi_val, double T, double g, double b){
    const double R = get_R(i, j, T, g, b);
    if (b == 0.){
        return 0.;
    }
    const double mu = m[i] * m[j] / (m[i] + m[j]);
    const double u0 = sqrt(2 * BOLTZMANN * T / mu) * g;
    const double ur = sqrt(pow(u0, 2) * (1. - pow(b / r, 2)) - 2. * potential(0, 0, r) / mu);

    const double theta = theta_r(i, j, R, r, T, g, b);
    // const double theta_p = PI - chi_val - theta_n;
    const double F = - potential_derivative_r(i, j, r);
    // const double dpdt = - 2 * mu * u0 * sin(2 * theta); //- (2. * u0 + (1. / u0)) * sin(theta) * abs(cos(theta));
    // const double dtdr = theta_integrand(i, j, T, r, g, b);
    // std::cout << "In weight(" << r / R << ", " << theta_n / PI << " / " << chi_val / PI << " / " << theta_p / PI << ") : "
    //     << R / r << ", "
    //     << dpdt << ", "
    //     << dtdr * sigma[i][j] << ", "
    //     << std::endl;
    return abs(F / ur * (cos(theta) + cos(PI - chi_val - theta))); // dpdt * dtdr; //  dpdt * dtdr; //
}

double Spherical::dpdt(int i, int j, double theta_n, double chi_val, double T, double g){
    const double mu = m[i] * m[j] / (m[i] + m[j]);
    const double u0 = sqrt(2 * BOLTZMANN * T / mu) * g;
    // const double theta_p = PI - chi_val - theta_n;
    return - 2 * mu * u0 * sin(2 * theta_n);
}

double Spherical::momentum_transfer_length(int i, int j, double T, double g, double b){
    const double R = get_R(i, j, T, g, b);
    if ((b == 0.) || (g == 0.)){
        return R;
    }
    const double dh{2.5e-2}, tol{1e-6}; // dh{5e-4}, tol{1e-10};
    const double chi_val = chi(i, j, T, g, b);
    const auto w_integrand = [&](double h){return (R / pow(h, 2)) * momentum_transfer_length_weight(i, j, R / h, chi_val, T, g, b);};
    const auto integrand = [&](double h){return (R / h) * w_integrand(h);};
    const double I = tanh_sinh(w_integrand, dh, tol);
    return tanh_sinh(integrand, dh, tol) / I;
    // return tanh_sinh(integrand, 7.5e-3) / I;
}

