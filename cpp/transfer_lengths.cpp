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

// vector2d Spherical::get_transfer_length(double rho, double T, const vector1d& x, int property){
//     switch (transfer_length_model_id) {
//         case 0:
//             return get_collision_diameters(rho, T, x);
//         case 1: {
//             double g0 = 0.;
//             double g1 = 3.5;
//             double gmax = 5.;
//             if (is_singlecomp){
//                 double I = get_tl_weight_normalizer(0, 0, T, property);
//                 const auto integrand = [&](double g){return tl_inner(0, 0, T, g, I, property);};
//                 const auto wt = [&](double g){return tl_weight_inner(0, 0, T, g, property);};
//                 std::pair<double, double> I1, I2;
//                 std::thread t1([&](std::pair<double, double>& I_){I_ = weighted_simpson(integrand, wt, g0, g1, 15);}, std::ref(I1));
//                 std::thread t2([&](std::pair<double, double>& I_){I_ = weighted_simpson(integrand, wt, g1, g_max, 10);}, std::ref(I2));
//                 t1.join(); t2.join();
//                 double F = I1.first + I2.first;
//                 double W = I1.second + I2.second;
//                 const double tl = F / W;
//                 // std::vector<std::thread> threads;
//                 // constexpr size_t ncores = 8;
//                 // double dg01 = (g1 - g0) / ncores;
//                 // double dg11 = (gmax - g1) / ncores;
//                 // vector1d tl_parts(ncores);
//                 // for (size_t ci = 0; ci < ncores; ci++){
//                 //     std::cout << "Started thread " << ci << std::endl;
//                 //     threads.push_back(std::thread(
//                 //         [&](double& tl_part){
//                 //             // double I = get_tl_weight_normalizer(0, 0, T, property);
//                 //             // const auto integrand = [&](double g){return tl_inner(0, 0, T, g, I, property);};
//                 //             tl_part = simpson(integrand, g0 + dg01 * ci, g0 + dg01 * (ci + 1), 2) + simpson(integrand, g1 + dg11 * ci, g1 + dg11 * (ci + 1), 2);
//                 //         }, std::ref(tl_parts[ci])));
//                 // }
//                 // int ci = 0;
//                 // for (auto it = threads.begin(); it != threads.end(); ++it){
//                 //     it->join();
//                 //     std::cout << "Joined thread " << ci++ << std::endl;
//                 // }
//                 // double tl = 0;
//                 // for (size_t ci = 0; ci < ncores; ci++){
//                 //     tl += tl_parts[ci];
//                 // }
//                 return vector2d(Ncomps, vector1d(Ncomps, tl));
//             }
//             vector2d tl(Ncomps, vector1d(Ncomps, 0.));
//             std::vector<std::thread> threads;
//             for (int i = 0; i < Ncomps; i++){
//                 for (int j = i; j < Ncomps; j++){
//                     threads.push_back(std::thread(
//                         [&](double& tl_ij, const int ci, const int cj){
//                             double I = get_tl_weight_normalizer(ci, cj, T, property);
//                             const auto integrand = [&](double g){return tl_inner(ci, cj, T, g, I, property);};
//                             tl_ij = simpson(integrand, g0, g1, 15) + simpson(integrand, g1, gmax, 10);
//                         }, std::ref(tl[i][j]), i, j));
//                 }
//             }
//             for (auto it = threads.begin(); it != threads.end(); ++it){
//                 it->join();
//             }
//             for (int i = 0; i < Ncomps; i++){
//                 for (int j = i; j < Ncomps; j++){
//                     tl[j][i] = tl[i][j];
//                 }
//             }
//             return tl;
//         }
//         case 2: {
//             switch (property){
//             case transfer_lengths::MTL:
//                 return MTL_correlation(rho, T);
//             case transfer_lengths::ETL:
//                 return ETL_correlation(rho, T);
//             default:
//                 throw std::runtime_error("Invalid transfer length type!");
//             }
//         }
//         default:
//             throw std::runtime_error("Invalid transfer length model!");
//     }
// }

vector2d Spherical::get_transfer_length(double rho, double T, const vector1d& x, int property){
    switch (transfer_length_model_id) {
        case 0:
            return get_collision_diameters(rho, T, x);
        case 1: {
            if (is_singlecomp){
                double tl = tl_ewca(0, 0, T, property);
                return vector2d(Ncomps, vector1d(Ncomps, tl));
            }
            vector2d tl(Ncomps, vector1d(Ncomps, 0.));
            std::vector<std::thread> threads;
            for (int i = 0; i < Ncomps; i++){
                for (int j = i; j < Ncomps; j++){
                    threads.push_back(std::thread(
                        [&](double& tl_ij, const int ci, const int cj){
                            tl_ij = tl_ewca(ci, cj, T, property);
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

double Spherical::tl_ewca(int i, int j, double T, int property){
    double g0 = 0.;
    double g1 = 3.5;
    double gmax = 5.;
    const auto integrand = [&](double g) -> std::pair<double, double> {return ewca_inner(i, j, T, g, property);};
    std::pair<double, double> I1, I2;
    std::thread t1([&](std::pair<double, double>& I_){I_ = weighted_simpson(integrand, g0, g1, 15);}, std::ref(I1));
    I2 = weighted_simpson(integrand, g1, gmax, 10);
    t1.join();
    double F = I1.first + I2.first;
    double W = I1.second + I2.second;
    return F / W;
}

std::pair<double, double> Spherical::ewca_inner(int i, int j, double T, double g, int property){
    const double bmid = get_bmid(i, j, g, T);
    const double bmax = get_b_max_g(i, j, g, T, bmid);
    
    const double b0 = 0.;
    const auto func = [&](double b){return get_R(i, j, T, g, b * sigma[i][j]);};
    const auto wt = [&](double b){return ewca_weight(i, j, T, g, b, property);};
    std::pair<double, double> I1, I2;
    I1 = weighted_simpson(func, wt, b0, bmid, 15);
    I2 = weighted_simpson(func, wt, bmid, bmax, 20);
    double F = I1.first + I2.first;
    double W = I1.second + I2.second;
    return std::pair<double, double>(F, W);
}

double Spherical::ewca_weight(int i, int j, double T, double g, double b, int property){
    switch (property) {
        case transfer_lengths::MTL:
            return momentum_transfer(i, j, T, g, b) * dimless_relative_vdf(g) * g * b;
        case transfer_lengths::ETL:
            return energy_transfer(i, j, T, g, b) * dimless_relative_vdf(g) * g * b;
        default:
            throw std::runtime_error("Invalid collision diameter model id!");
    }
}

double Spherical::momentum_transfer(int i, int j, double T, double g, double b){
    double chi_val = chi(i, j, T, g, b * sigma[i][j]);
    double red_mass = m[i] * m[j] / (m[i] + m[j]);
    double U = sqrt(2 * BOLTZMANN * T / red_mass) * g;
    return g * sqrt(2 * (1 - cos(chi_val))) * abs(sin(chi_val / 2.));
}

double Spherical::energy_transfer(int i, int j, double T, double g, double b){
    double chi_val = chi(i, j, T, g, b * sigma[i][j]);
    double red_mass = m[i] * m[j] / (m[i] + m[j]);
    double U = sqrt(2 * BOLTZMANN * T / red_mass) * g;
    return g * abs(cos(chi_val) - sin(chi_val) - 1);
}

double Spherical::get_b_max_g(int i, int j, double g, double T, double bmid){
    double b, db, chi_val;
    if (g == 0) return 10.;
    b = 2 * bmid;
    db = 0.5;
    chi_val = chi(i, j, T, g, b * sigma[i][j]);
    while (abs(chi_val + 1e-5) > 1e-3){
        b += db;
        chi_val = chi(i, j, T, g, b * sigma[i][j]);
    }
    return b;
}

double Spherical::get_bmid(int i, int j, double g, double T){
    if (g == 0) return 5.;
    double b = 0.9;
    double db = 0.1;
    double chi_val = chi(i, j, T, g, b * sigma[i][j]);
    double next_chi_val;
    double dchidb = 10; // Init with placeholder to start loop
    double step_tol = 1e-6;
    double chi_tol = 1e-6;
    while ((abs(chi_val) > chi_tol) && (abs(db) > step_tol)){
        if ((chi_val < 0 && db > 0) || (chi_val > 0 && db < 0)){
            db *= -0.9;
        }
        next_chi_val = chi(i, j, T, g, (b + db) * sigma[i][j]);
        while (abs(next_chi_val) >= abs(chi_val)){
            db *= 0.1;
            next_chi_val = chi(i, j, T, g, (b + db) * sigma[i][j]);
            if (abs(db) < step_tol) break;
        }
        if (abs(db) > step_tol){
            b += db;
            dchidb = (chi_val - next_chi_val) / db;
            chi_val = next_chi_val;
            db = - chi_val / dchidb;
        }
        if ((b > 10) || (b < 0)) {
            std::cout << "bmid did not converge!" << std::endl;
            throw std::runtime_error("Unable to converge bmid!");
            break;
        }
    }
    return b;
}

/********************************************************************************************/
/**********************         TL MODEL 2 : correlations              **********************/
/***      NOTE: These correlations have been fitted to reproduce the viscosity and       ****/ 
/***      thermal conductivity of argon, and should not be used for anything else.       ****/
/********************************************************************************************/

vector2d Spherical::MTL_correlation(double rho, double T){
    double rho_r = rho * pow(sigma[0][0], 3);
    double Tr = T * BOLTZMANN / eps[0][0];

    static constexpr double s[5] = {1.029, 0.091, 0.615, 1.074, - 0.603};

    double mtl = s[0] - s[1] * log(Tr) + s[2] * (rho_r / Tr) * exp(- (s[3] + s[4] * rho_r) * sqrt(Tr) );
    return vector2d(Ncomps, vector1d(Ncomps, mtl * sigma[0][0]));
}

vector2d Spherical::ETL_correlation(double rho, double T){
    double rho_r = rho * pow(sigma[0][0], 3);
    double Tr = T * BOLTZMANN / eps[0][0];

    const double a = 1.054;
    const double bi[2] = {0.100, - 0.023};
    const double ci[2] = {- 1.166, 2.389};
    const double c0 = 3. * bi[0];
    const double c1 = ci[0];
    const double c2 = ci[1];
    const double c3 = (c1 - c0);
    const double b = bi[0] + bi[1] * rho_r;
    const double c = c0 + c1 * rho_r + c2 * pow(rho_r, 2) + c3 * pow(rho_r, 3);
    const double etl = a - b * log(Tr) + c / pow(Tr, 1.5);
    return vector2d(Ncomps, vector1d(Ncomps, etl * sigma[0][0]));
}