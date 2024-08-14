/*
Author: Vegard Gjeldvik Jervell
Contains: Transfer length models (see: <unpublished> 2024). 
        Primary control flow is as follows: 
            - Use get_mtl and get_etl to compute the transfer lengths you need.
            
            - These forward the call to get_transfer_length, which handles dispatch to the appropriate 
              transfer length model: collision diameters, (old) or transfer lengths (new)
                - Use set_transfer_length_model to select the model beforehand.

            - In the case of the (new) models, all methods used in the computation are the same, except for the
              function computing the weight, so the parameter "property", which should be one of the values in the 
              enum "transfer_lengths", is passed through the call chain until the weight is computed, and used in
              tl_weight_integrand to select the correct weight function.
*/
#include "Spherical.h"
#include "Integration/Integration.h"

enum transfer_lengths{
    MTL = 0,
    ETL
};

vector2d Spherical::model_mtl(double rho, double T, const vector1d& x){
    StatePoint point = get_transfer_length_point(rho, T, x);
    const auto pos = mtl_map.find(point);

    if (pos != mtl_map.end()) return pos->second;

    vector2d mtl = get_transfer_length(rho, T, x, transfer_lengths::MTL);
    mtl_map[point] = mtl;
    return mtl;
}

vector2d Spherical::model_etl(double rho, double T, const vector1d& x){
    StatePoint point = get_transfer_length_point(rho, T, x);
    const auto pos = etl_map.find(point);

    if (pos != etl_map.end()) return pos->second;

    vector2d etl = get_transfer_length(rho, T, x, transfer_lengths::ETL);
    etl_map[point] = etl;
    return etl;
}

vector2d Spherical::get_transfer_length(double rho, double T, const vector1d& x, int property){
    switch (transfer_length_model_id) {
        case 0:
            return get_collision_diameters(rho, T, x);
        case 1: {
            double g0 = 0.;
            double g1 = 3.5;
            double gmax = 5.;
            if (is_singlecomp){
                double I = get_tl_weight_normalizer(0, 0, T, property);
                const auto integrand = [&](double g){return tl_inner(0, 0, T, g, I, property);};
                const double tl = simpson(integrand, g0, g1, 15) + simpson(integrand, g1, gmax, 10);
                return vector2d(Ncomps, vector1d(Ncomps, tl));
            }
            vector2d tl(Ncomps, vector1d(Ncomps, 0.));
            std::vector<std::thread> threads;
            for (int i = 0; i < Ncomps; i++){
                for (int j = i; j < Ncomps; j++){
                    threads.push_back(std::thread(
                        [&](double& tl_ij, const int ci, const int cj){
                            double I = get_tl_weight_normalizer(ci, cj, T, property);
                            const auto integrand = [&](double g){return tl_inner(ci, cj, T, g, I, property);};
                            tl_ij = simpson(integrand, g0, g1, 15) + simpson(integrand, g1, gmax, 10);
                        }, std::ref(tl[i][j]), i, j));
                }
            }
            for (auto it = threads.begin(); it != threads.end(); ++it){
                it->join();
            }
            for (int i = 0; i < Ncomps; i++){
                for (int j = i; j < Ncomps; j++){
                    tl[j][i] = tl[i][j];
                }
            }
            return tl;
        }
        case 2: {
            switch (property){
            case transfer_lengths::MTL:
                return MTL_correlation(rho, T);
            case transfer_lengths::ETL:
                return ETL_correlation(rho, T);
            default:
                throw std::runtime_error("Invalid transfer length type!");
            }
        }
        default:
            throw std::runtime_error("Invalid transfer length model!");
    }
}


/**************************************************************************************************/
/**********************         TL MODEL 0 : Collision diameter              **********************/
/**************************************************************************************************/

vector2d Spherical::get_collision_diameters(double rho, double T, const std::vector<double>& x){
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

/*************************************************************************************************************************/
/**********************         TL MODEL 1 : Exchange weighted closest approach (EWCA)              **********************/
/*************************************************************************************************************************/

inline double dimless_relative_vdf(double g){
    return sqrt(2. / PI) * pow(g, 2) * exp(- 0.5 * pow(g, 2));
}

double Spherical::tl_inner(int i, int j, double T, double g, double I, int property){
    const double bmax = get_b_max_g(i, j, g, T);
    const double bmid = get_bmid(i, j, g, T);
    const double b0 = 0.;
    const auto integrand = [&](double b){return get_tl_weight(i, j, T, g, b, I, bmax, property) * get_R(i, j, T, g, b * sigma[i][j]);};
    double val = simpson(integrand, b0, bmid, 15);
    val += simpson(integrand, bmid, bmax, 20);
    return val;
}

double Spherical::tl_weight_integrand(int i, int j, double T, double g, double b, double bmax, int property){
    switch (property) {
        case transfer_lengths::MTL:
            return momentum_transfer(i, j, T, g, b) * dimless_relative_vdf(g) * g * 2 * PI * b;
        case transfer_lengths::ETL:
            return energy_transfer(i, j, T, g, b) * dimless_relative_vdf(g) * g * 2 * PI * b;
        default:
            throw std::runtime_error("Invalid collision diameter model id!");
    }
}

double Spherical::get_tl_weight(int i, int j, double T, double g, double b, double I, double bmax, int property){
    return tl_weight_integrand(i, j, T, g, b, bmax, property) / I;
}

double Spherical::get_tl_weight_normalizer(int i, int j, double T, int property){
    double g0 = 0;
    double g1 = 3.5;
    double g_max = 5.;
    const auto integrand = [&](double g){return tl_weight_inner(i, j, T, g, property);};
    double I = simpson(integrand, g0, g1, 15);
    I += simpson(integrand, g1, g_max, 10);
    return I;
}

double Spherical::tl_weight_inner(int i, int j, double T, double g, int property){
    double bmax = get_b_max_g(i, j, g, T);
    double bmid = get_bmid(i, j, g, T);
    const auto integrand = [&](double b){return tl_weight_integrand(i, j, T, g, b, bmax, property);};
    double val = simpson(integrand, 0., bmid, 15);
    val += simpson(integrand, bmid, bmax, 20);
    return val;
}

double Spherical::momentum_transfer(int i, int j, double T, double g, double b){
    double chi_val = chi(i, j, T, g, b * sigma[i][j]);
    double red_mass = m[i] * m[j] / (m[i] + m[j]);
    double U = sqrt(2 * BOLTZMANN * T / red_mass) * g;
    return red_mass * U * sqrt(2 * (1 - cos(chi_val))) * abs(sin(chi_val / 2.));
}

double Spherical::energy_transfer(int i, int j, double T, double g, double b){
    double chi_val = chi(i, j, T, g, b * sigma[i][j]);
    double red_mass = m[i] * m[j] / (m[i] + m[j]);
    double U = sqrt(2 * BOLTZMANN * T / red_mass) * g;
    return red_mass * U * abs(cos(chi_val) - sin(chi_val) - 1);
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

/********************************************************************************************/
/**********************         TL MODEL 2 : correlations              **********************/
/********************************************************************************************/

vector2d Spherical::MTL_correlation(double rho, double T){
    double rho_r = rho * pow(sigma[0][0], 3);
    double Tr = T * BOLTZMANN / eps[0][0];

    static constexpr double ai[2] = {1.033, 0.058};
    static constexpr double bi[2] = {- 0.244, 0.325};
    static constexpr double ci[2] = {0.270, 2.413};

    double a = ai[0] + ai[1] * rho_r;
    double b = bi[0] + bi[1] * a;
    double c = ci[0] * sinh(ci[1] * rho_r);
    double mtl = a - b * log(Tr) + c * exp(- Tr) / Tr;
    return vector2d(Ncomps, vector1d(Ncomps, mtl * sigma[0][0]));
}

vector2d Spherical::ETL_correlation(double rho, double T){
    throw std::runtime_error("ETL correlation not implemented!");
    return MTL_correlation(rho, T);
}