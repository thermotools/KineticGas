#include "Spherical.h"
#include "Integration/Integration.h"

vector2d Spherical::get_collision_diameters(double rho, double T, const std::vector<double>& x){
    int T_idx = static_cast<int>(T * 100 + 0.5); // Using T in centi kelvin for lookup.
    const std::map<int, double>::iterator pos = collision_diameter_map.find(T_idx);
    double cd{0.};
    if (pos == collision_diameter_map.end()){
        switch (collision_diameter_model_id) {
        case 0:
            cd = get_collision_diameters_model0(rho, T, x)[0][0]; break;
        case 1:
        case 2:
        case 3:
            cd = get_collision_diameters_model1(rho, T, x)[0][0]; break;
        default:
            throw std::runtime_error("Invalid collision diameter model id!"); // collision_diameter_model_id is private, so it shouldn't be possible to end up here.
        }
        collision_diameter_map[T_idx] = cd;
    }
    else{
        cd = pos->second;
    }
    return vector2d(Ncomps, vector1d(Ncomps, cd));
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
    if (m[0] != m[1]){
        throw std::runtime_error("Collision diameter model 1 only implemented for equal-mass particles!.");
    }
    const double cd = momentum_collision_diameter(T);
    return vector2d(Ncomps, vector1d(Ncomps, cd));
}

double Spherical::momentum_collision_diameter(double T){
    const double I = get_cd_weight_normalizer(T);
    const double g0 = 0.;
    const double g1 = 3.5;
    const double gmax = 5.;
    const auto integrand = [&](double g){return cd_inner(T, g, I);};
    double cd = simpson(integrand, g0, g1, 15);
    cd += simpson(integrand, g1, gmax, 10);
    return cd;
}

double Spherical::cd_inner(double T, double g, double I){
    const double bmax = get_b_max_g(g, T);
    const double bmid = get_bmid(g, T);
    const double b0 = 0.;
    const auto integrand = [&](double b){return cd_integrand(T, g, b, I, bmax);};
    double val = simpson(integrand, b0, bmid, 15);
    val += simpson(integrand, bmid, bmax, 20);
    return val;
}

double Spherical::cd_weight_integrand(double T, double g, double b, double bmax){
    switch (collision_diameter_model_id) {
        case 1:
            return momentum_transfer(T, g, b) * dimless_relative_vdf(g) * g * 2 * PI * b;
        case 2:
            return momentum_transfer(T, g, b) * dimless_relative_vdf(g) * g * 2 * PI * b;
        case 3:
            return momentum_transfer(T, g, b) * dimless_relative_vdf(g) * g * 2 * PI * b;
        default:
            throw std::runtime_error("Invalid collision diameter model id!");
    }
}

double Spherical::get_cd_weight(double T, double g, double b, double I, double bmax){
    return cd_weight_integrand(T, g, b, bmax) / I;
}

double Spherical::get_cd_weight_normalizer(double T){
    double g0 = 0;
    double g1 = 3.5;
    double g_max = 5.;
    const auto integrand = [&](double g){return cd_weight_inner(T, g);};
    double I = simpson(integrand, g0, g1, 15);
    I += simpson(integrand, g1, g_max, 10);
    return I;
}

double Spherical::cd_weight_inner(double T, double g){
    double bmax = get_b_max_g(g, T);
    double bmid = get_bmid(g, T);
    const auto integrand = [&](double b){return cd_weight_integrand(T, g, b, bmax);};
    double val = simpson(integrand, 0., bmid, 15);
    val += simpson(integrand, bmid, bmax, 20);
    return val;
}

double Spherical::cd_integrand(double T, double g, double b, double I, double bmax){
    const double wt = get_cd_weight(T, g, b, I, bmax);
    const double R = get_R(0, 0, T, g, b * sigma[0][0]);
    return wt * R;
}

double Spherical::momentum_transfer(double T, double g, double b){
    double chi_val = chi(0, 0, T, g, b * sigma[0][0]);
    double U = sqrt(4 * BOLTZMANN * T / m[0]) * g;
    double dp;
    switch (collision_diameter_model_id){
        case 1:
            dp = U * sqrt(2 * (1 - cos(chi_val)));
            break;
        case 2:
            dp = U * sqrt(2 * (1 - cos(chi_val))) * abs(sin(chi_val / 2.));
            break;
        case 3:
            dp = U * sqrt(2 * (1 - cos(chi_val))) * sin(chi_val / 2.);
            break;
        default:
            throw std::runtime_error("Invalid collision diameter model!");
    }
    return dp;
}

double Spherical::get_b_max_g(double g, double T){
    double b, db, chi_val;
    if (g == 0) return 10.;
    b = 5.;
    db = -0.5;
    chi_val = chi(0, 0, T, g, b * sigma[0][0]);
    while (abs(db) > 1e-3){
        double scale = g * sqrt(2 * (1 - cos(chi_val))) * abs(sin(chi_val / 2.));
        if (scale > 1e-3){
            b -= db;
            db *= 0.5;
        }
        b += db;
        chi_val = chi(0, 0, T, g, b * sigma[0][0]);
    }
    return b;
}

double Spherical::get_bmid(double g, double T){
    if (g == 0) return 5.;
    double b = 0.5;
    double db = 0.1;
    double chi_val = chi(0, 0, T, g, b * sigma[0][0]);
    while (abs(db) > 1e-5){
        if ((chi_val < 0 && db > 0) || (chi_val > 0 && db < 0)){
            db *= -0.5;
        }
        b += db;
        chi_val = chi(0, 0, T, g, b * sigma[0][0]);
        if (b > 10) {
            std::cout << "bmid did not converge!" << std::endl;
            throw std::runtime_error("Unable to converge bmid!");
            break;
        }
    }
    return b;
}

/**********************************************************************************************/
/***********************                      MODEL 2                  ************************/
/**********************************************************************************************/

vector2d Spherical::get_collision_diameters_model2(double T){
    const double I = get_cd_weight_normalizer_2(T);
    const auto integrand = [&](double g){return cd_inner_2(T, g, I);};
    double cd = simpson(integrand, 1e-6, 5., 30);
    return vector2d(Ncomps, vector1d(Ncomps, cd));
}

double Spherical::cd_inner_2(double T, double g, double I){
    const auto integrand = [&](double R){return get_cd_weight_2(T, g, R, I) * R;};
    const double Rmin = get_R_min(T, g);
    return simpson(integrand, Rmin, 5., 30);
}

double Spherical::get_cd_weight_2(double T, double g, double R, double I){
    return momentum_transfer_R(T, g, R) * ideal_rdf(T, R * sigma[0][0]) * dimless_relative_vdf(g) / I;
}

double Spherical::get_R_min(double T, double g){
    double r = sigma[0][0];
    double dr = 0.01 * sigma[0][0];
    while (potential(0, 0, r) < BOLTZMANN * T * pow(g, 2)){
        r -= dr;
    }
    dr *= 0.5;
    r += dr;
    while (potential(0, 0, r) > BOLTZMANN * T * pow(g, 2)){
        r += dr;
    }
    return r;
}

double Spherical::cd_weight_inner_2(double T, double g){
    const double R_min = get_R_min(T, g);
    const auto integrand = [&](double R){return momentum_transfer_R(T, g, R) * ideal_rdf(T, R * sigma[0][0]) * dimless_relative_vdf(g);};
    return simpson(integrand, R_min, 5., 30);
}

double Spherical::get_cd_weight_normalizer_2(double T){
    const auto integrand = [&](double g){return cd_weight_inner_2(T, g);};
    return simpson(integrand, 1e-6, 5., 30);
}

double Spherical::ideal_rdf(double T, double r){
    return exp(- potential(0, 0, r) / (BOLTZMANN * T));
}

double Spherical::momentum_transfer_R(double T, double g, double R){
    double chi_val = chi_R(T, g, R * sigma[0][0]);
    double U = sqrt(4 * BOLTZMANN * T / m[0]) * g;
    return U * sqrt(2 * (1 - cos(chi_val)));
}

double Spherical::chi_R(double T, double g, double R){
    return PI - 2 * theta_R(T, g, R);
}

double Spherical::theta_R(double T, double g, double R){
    double b = R * sqrt(1 - (potential(0, 0, R) / (BOLTZMANN * T * pow(g, 2))));
    return theta_integral(0, 0, T, R, g, b) - theta_lim(0, 0, T, g) + PI / 2;
}