#include "Spherical.h"

Spherical::Spherical(std::vector<double> mole_weights,
                    std::vector<std::vector<double>> sigmaij,
                    bool is_idealgas) 
                    : KineticGas(mole_weights, is_idealgas), sigma{sigmaij}{


    w_integrand_export = std::bind(&Spherical::w_integrand, this,
                                        std::placeholders::_1, std::placeholders::_2,
                                        std::placeholders::_3, std::placeholders::_4,
                                        std::placeholders::_5, std::placeholders::_6,
                                        std::placeholders::_7);

}

double Spherical::omega(int i, int j, int l, int r, double T){
    OmegaPoint point{i, j, l, r, T}, sympoint{j, i, l, r, T};
    const std::map<OmegaPoint, double>::iterator pos = omega_map.find(point);
    if (pos == omega_map.end()){
        double w = w_integral(i, j, T, l, r);
        double val;
        if (i == j) val = pow(sigma[i][j], 2) * sqrt((PI * BOLTZMANN * T) / m[i]) * w;
        else val = 0.5 * pow(sigma[i][j], 2) * sqrt(2 * PI * BOLTZMANN * T / (m0[i][j] * M[i][j] * M[j][i])) * w;
        omega_map[point] = val;
        omega_map[sympoint] = val; // Collision integrals are symmetric wrt. particle indices.
        return val;
    }
    return pos->second;
}

double Spherical::omega_tester(int i, int j, int l, int r, double T, IntegrationParam& param){
    double w = w_integral_tester(i, j, T, l, r, param);
    if (i == j) return pow(sigma[i][j], 2) * sqrt((PI * BOLTZMANN * T) / m[i]) * w;
    return 0.5 * pow(sigma[i][j], 2) * sqrt(2 * PI * BOLTZMANN * T / (m0[i][j] * M[i][j] * M[j][i])) * w;
}

double Spherical::w_integral_tester(int i, int j, double T, int l, int r, IntegrationParam& param){
    double I = integrate2d(param.origin, param.end,
                        param.dg, param.db,
                        param.refinement_levels_g, param.refinement_levels_b,
                        param.subdomain_dblder_limit,
                        i, j, T, l, r,
                        w_integrand_export);

    return I;
}

double Spherical::w_integral(int i, int j, double T, int l, int r){
    /*
    Evaulate the dimensionless collision integral

    See: The Kinetic Gas Theory of Mie fluids (V. G. Jervell, Norwegian University of Science and Technology, 2022)
    https://ntnuopen.ntnu.no/ntnu-xmlui/handle/11250/3029213
    For details on the integration routine implemented in integrate2d.
    */
    Point origin{1e-7, 1e-7};
    Point end{8, 5};
    double dg{0.5}, db{0.03125};
    int refinement_levels_g{4};
    int refinement_levels_b{16};
    double subdomain_dblder_limit{1e-5};

    double I = integrate2d(origin, end,
                        dg, db,
                        refinement_levels_g, refinement_levels_b,
                        subdomain_dblder_limit,
                        i, j, T, l, r,
                        w_integrand_export);
    
    return I;
}

double Spherical::w_integrand(int i, int j, double T, 
                                        double g, double b,
                                        int l, int r){ // Using b = b / sigma to better scale the axes. Multiply the final integral by sigma.
    const double chi_val = chi(i, j, T, g, b * sigma[i][j]);
    return 2 * exp(- pow(g, 2)) * pow(g, 2.0 * r + 3.0) * (1 - pow(cos(chi_val), l)) * b;
};

std::vector<std::vector<double>> Spherical::get_b_max(double T, std::vector<std::vector<int>>& ierr){
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

std::vector<std::vector<double>> Spherical::get_b_max(double T){
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

std::vector<std::vector<double>> Spherical::get_collision_diameters(double rho, double T, const std::vector<double>& x){
    // Evaluate the integral of Eq. (40) in RET for Mie fluids (doi: 10.1063/5.0149865)
    // Note: For models with is_idealgas==true, zeros are returned, as there is no excluded volume at infinite dilution.
    // Weights and nodes for 6-point Gauss-Legendre quadrature
    constexpr int n_gl_points = 6;
    constexpr double gl_points[n_gl_points] = {-0.93246951, -0.66120939, -0.23861919, 0.23861919, 0.66120939, 0.93246951};
    constexpr double gl_weights[n_gl_points] = {0.17132449, 0.36076157, 0.46791393, 0.46791393, 0.36076157, 0.17132449};

    const double g_avg = sqrt(1. / PI); // Average dimensionless relative speed of collision
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

/*
Contains functions for computing the collision integrals
Common variables are:
    i, j : Species indices
    T : Temperature [K]
    l, r : Collision integral indices (in omega and w_<potential_type> functions)
    g : Dimentionless relative velocity of molecules
    b : "Impact paramter", closest distance of approach of the centers of mass if molecules were non-interacting point particles
    theta : Angular polar coordinate. Theta = 0 when the particles are infinitely far away from each other, before the collision,
        theta(t1) is the angle between the line between the particles before interaction begins (at infinite distance), and the line between the particles at t = t1
    r : Radial polar coordinate. Distance from the center of mass of one particle to the other.
    R : Actual distance of closest approach. Corresponds to dr/dtheta = 0
    chi : Deflection angle, corresponds to pi - 2 * theta(R)
*/

#include "KineticGas.h"
#include "Integration/Integration.h"

// Helper funcions for computing dimentionless collision integrals

double Spherical::theta(int i, int j, const double T, const double g, const double b){
    // Compute deflection angle for a collision
    if (b / sigma[i][j] > 10) return PI / 2;
    if (b / sigma[i][j] < 1e-3) return 0;
    double R = get_R(i, j, T, g, b);
    return theta_integral(i, j, T, R, g, b); // - theta_lim(i, j, T, g) + PI / 2;
}

double Spherical::theta_lim(int i, int j, const double T, const double g){
    double b = 10 * sigma[i][j];
    double R = get_R(i, j, T, g, b);
    return theta_integral(i, j, T, R, g, b);
}

double Spherical::theta_integral(int i, int j, double T, double R, double g, double b){
    constexpr double h{7.5e-3};
    double I{0.0};
    int k{1};
    double u = tanh(PI * sinh(k * h) / 2.);
    double w = (PI / 2.) * h * cosh(k * h) / pow(cosh(PI * sinh(k * h) / 2.), 2);
    double f = transformed_theta_integrand(i, j, T, u, R, g, b);
    while (abs(f * w) > 1e-8){
        I += w * f;
        k+=1;
        u = tanh(PI * sinh(k * h) / 2.);
        w = (PI / 2) * h * cosh(k * h) / pow(cosh(PI * sinh(k * h) / 2.), 2);
        f = transformed_theta_integrand(i, j, T, u, R, g, b);
        if (isnan(f) || isnan(w) || isinf(f) || isinf(w)) break;
    }
    return I;
}

double Spherical::theta_integrand(int i, int j, double T, double r, double g, double b){
    return pow((pow(r, 4) / pow(b, 2)) * (1.0 - potential(i, j, r) / (BOLTZMANN * T * pow(g, 2))) - pow(r, 2), -0.5);
}

double Spherical::transformed_theta_integrand(int i, int j, double T, double u, double R, double g, double b){
    // Transformed by the substitution u = R / r
    return (R / pow(u, 2)) * theta_integrand(i, j, T, R / u, g, b);
}

double Spherical::theta_integrand_dblderivative(int i, int j, double T, double r, double g, double b){
    // Expressing the integrand as f = (core)^{-1/2}
    const double a = 1.0 / (pow(b, 2) * BOLTZMANN * T * pow(g, 2));
    const double u = potential(i, j, r);
    const double u_prime = potential_derivative_r(i, j, r);
    const double u_dblprime = potential_dblderivative_rr(i, j, r);
    const double core = pow(r, 4) / pow(b, 2) - a * pow(r, 4) * u - pow(r, 2);
    const double core_prime = 4 * pow(r, 3) / pow(b, 2) - a * (4 * pow(r, 3) * u + pow(r, 4) * u_prime) - 2 * r;
    const double core_dblprime = 12.0 * pow(r, 2) / pow(b, 2) - a * (12 * pow(r, 2) * u + 8 * pow(r, 3) * u_prime + pow(r, 4) * u_dblprime) - 2;

    double val = (3.0 / 4.0) * pow(core, -2.5) * pow(core_prime, 2) - 0.5 * pow(core, - 1.5) * core_dblprime;
    return val;
}

double Spherical::get_R_rootfunc(int i, int j, double T, double g, double b, double& r){
    return (potential(i, j, r) / (BOLTZMANN * T * pow(g, 2))) + pow(b / r, 2) - 1;
}

double Spherical::get_R_rootfunc_derivative(int i, int j, double T, double g, double b, double& r){
    return (potential_derivative_r(i, j, r) / (BOLTZMANN * T * pow(g, 2))) - 2 * pow(b, 2) / pow(r, 3);
}

double Spherical::get_R(int i, int j, double T, double g, double b){
    // Newtons method
    double tol = 1e-5; // Relative to sigma[i][j]
    double init_guess_factor = 1.0;
    double r = init_guess_factor * b;
    double f = get_R_rootfunc(i, j, T, g, b, r);
    double dfdr = get_R_rootfunc_derivative(i, j, T, g, b, r);
    double next_r = r - f / dfdr;
    while (abs((r - next_r) / sigma[i][j]) > tol){
        if (next_r < 0){
            init_guess_factor *= 0.95;
            r = init_guess_factor * b;
        }
        else if (f < 0 && f / dfdr < 0){
            init_guess_factor *= 0.95;
            r = init_guess_factor * b;
        }
        else{
            r = next_r;
        }
        f = get_R_rootfunc(i, j, T, g, b, r);
        dfdr = get_R_rootfunc_derivative(i, j, T, g, b, r);
        next_r = r - f / dfdr;
    }
    return next_r;
}

double Spherical::chi(int i, int j, double T, double g, double b){
    if (b / sigma[i][j] > 10) return 0;
    double t = theta(i, j, T, g, b);
    return PI - 2.0 * t;
}


