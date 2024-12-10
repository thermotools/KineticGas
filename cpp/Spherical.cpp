#include "Spherical.h"
#include "Integration/Integration.h"

Spherical::Spherical(vector1d mole_weights,
                    vector2d sigma,
                    vector2d eps,
                    bool is_idealgas, bool is_singlecomp) 
                    : KineticGas(mole_weights, sigma, eps, is_idealgas, is_singlecomp)
                    {}

double Spherical::potential(int i, int j, double r){
    return potential(i, j, static_cast<dual2>(r)).val.val;
}

double Spherical::potential_derivative_r(int i, int j, double r){
    dual2 rd = r;
    const auto func = [&](dual2 r_){return potential(i, j, r_);};
    auto [u0, ur, urr] = autodiff::derivatives(func, autodiff::wrt(rd), autodiff::at(rd));
    return ur;
}

double Spherical::potential_dblderivative_rr(int i, int j, double r){
    dual2 rd = r;
    const auto func = [&](dual2 r_){return potential(i, j, r_);};
    auto [u0, ur, urr] = autodiff::derivatives(func, autodiff::wrt(rd, rd), autodiff::at(rd));
    return urr;
}

StatePoint Spherical::get_transfer_length_point(double rho, double T, const vector1d& x){
    if (transfer_length_model_id == TransferLengthModel::correlation){
        return StatePoint(T, rho);
    }
    return StatePoint(T);
}

double Spherical::omega(int i, int j, int l, int r, double T){
    OmegaPoint point = get_omega_point(i, j, l, r, T);
    OmegaPoint sympoint = get_omega_point(j, i, l, r, T);
    const std::map<OmegaPoint, double>::iterator pos = omega_map.find(point);
    
    if (pos != omega_map.end()) return pos->second;

    double w = w_integral(i, j, T, l, r);
    double val;
    if (i == j) val = pow(sigma[i][j], 2) * sqrt((PI * BOLTZMANN * T) / m[i]) * w;
    else val = 0.5 * pow(sigma[i][j], 2) * sqrt(2 * PI * BOLTZMANN * T / (m0[i][j] * M[i][j] * M[j][i])) * w;
    omega_map[point] = val;
    omega_map[sympoint] = val; // Collision integrals are symmetric wrt. particle indices.
    if (is_singlecomp){
        for (int ci = 0; ci < Ncomps; ci++){
            for (int cj = 0; cj < Ncomps; cj++){
                if (((ci == i) && (cj == j)) || ((ci == j) && (cj == i))) continue;
                OmegaPoint purepoint = get_omega_point(ci, cj, l, r, T);
                omega_map[purepoint] = val;
            }
        }
    }
    return val;
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
    if (T < 70){
        dg /= 4.; 
        end = Point(6, 3.5);
    }
    if (T < 50){
        refinement_levels_b /= 2;
    }
    if (T < 35){
        dg /= 2;
        refinement_levels_g *= 2;
        db /= 2;
        end = Point(3.5, 3.5);
    }
    if (T < 20){
        dg /= 2;
        refinement_levels_g *= 2;
        db /= 2;
        refinement_levels_b *= 2;
    }

    const auto w_integrand_export = [&](double g, double b){return w_integrand(i, j, T, g, b, l, r);};
    double I = integrate2d(origin, end,
                        dg, db,
                        refinement_levels_g, refinement_levels_b,
                        subdomain_dblder_limit,
                        w_integrand_export);
    
    return I;
}

double Spherical::w_integrand(int i, int j, double T, 
                                        double g, double b,
                                        int l, int r){ // Using b = b / sigma to better scale the axes. Multiply the final integral by sigma.
    if ((b > 2.5) && (g > 1.5)) return 0;
    const double chi_val = chi(i, j, T, g, b * sigma[i][j]);
    return 2 * exp(- pow(g, 2)) * pow(g, 2.0 * r + 3.0) * (1 - pow(cos(chi_val), l)) * b;
};

double Spherical::cross_section(int i, int j, int l, double E){
    /*
        Collision cross-sections, \int_0^{\infty} (1 - \cos^{l}(\chi)) b \d b
        - i, j : Component indices
        - l : Cross section moment
        - E : Dinemsionless collision energy (E / (k_B T) = g_{ij}^2)

        Note: Because other internal functions work with dimensionless velocity, rather than collision energy,
              we must set a dummy temperature to convert the energy to a collision velocity.
    */
    if (l == 0) return 0;
    E *= eps[i][j];
    const double T = 300;
    const double g = sqrt(E / (BOLTZMANN * T));
    const auto integrand = [&](double b){
        double chi_val = chi(i, j, T, g, b * sigma[i][j]);
        return (1. - pow(cos(chi_val), l)) * b;
    };
    return pow(sigma[i][j], 2) * simpson_inf(integrand, 1e-6, 0.5); // Integration carried out in reduced units.
}

double Spherical::hs_cross_section(int i, int j, int l){
    if (l % 2 == 1) return PI * pow(sigma[i][j], 2);
    return (1 - 1. / (1 + l)) * PI * pow(sigma[i][j], 2);
}

double Spherical::reduced_cross_section(int i, int j, int l, double E){
    double Q = cross_section(i, j, l, E); 
    double Q_hs = hs_cross_section(i, j, l);
    return Q / Q_hs; 
}

double Spherical::second_virial(int i, int j, double T){
    set_internals(0, T, {1.});
    const auto integrand = [&](double r){return pow(r, 2) * (exp(- potential(i, j, r * sigma[i][j]) / (BOLTZMANN * T)) - 1);}; // Integration variable is r / sigma[i][j];
    const double r0 = 0.1;
    const double r1 = 2;
    double I = (pow(r0, 3) / 3) + simpson(integrand, r0, r1, 100);
    I += simpson_inf(integrand, r1, 3);
    return - 0.5 * 4 * PI * pow(sigma[i][j], 3) * AVOGADRO * I;
}

/*
Contains functions for computing the collision integrals
Common variables are:
    i, j : Species indices
    T : Temperature [K]
    l, r : Collision integral indices
    g : Dimentionless relative velocity of molecules, reduced using : 0.5 * red_mass[i][j] * v^2 = g^2 * k_B T (See: Eq. (3) of Ref. (II) in header file.) 
    b : "Impact paramter", closest distance of approach of the centers of mass if molecules were non-interacting point particles
    theta : Angular polar coordinate. Theta = 0 when the particles are infinitely far away from each other, before the collision,
        theta(t1) is the angle between the line between the particles before interaction begins (at infinite distance), and the line between the particles at t = t1
    r : Radial polar coordinate. Distance from the center of mass of one particle to the other.
    R : Actual distance of closest approach. Corresponds to dr/dtheta = 0
    chi : Deflection angle, corresponds to pi - 2 * theta(R)
*/

// Helper funcions for computing dimentionless collision integrals

double Spherical::theta(int i, int j, const double T, const double g, const double b){
    /* Compute angular collision coordinate at distance of closest approach
        - i, j : Component indices
        - T : Temperature (K)
        - g : Dimensionless relative velocity
        - b : Impact parameter (m)
    */
    if (b / sigma[i][j] > 100) return PI / 2;
    if (b / sigma[i][j] < 1e-3) return 0;
    double R = get_R(i, j, T, g, b);
    return theta_integral(i, j, T, R, g, b) - theta_lim(i, j, T, g) + PI / 2;
}

double Spherical::theta_r(int i, int j, double r, double T, double g, double b){
    /* Compute Angular collision coordinate at specified distance
        - i, j : Component indices
        - r : Separation distance (m)
        - T : Temperature (K)
        - g : Dimensionless relative velocity
        - b : Impact parameter (m)
    */
    double R = get_R(i, j, T, g, b);
    return theta_r(i, j, R, r, T, g, b);
}

double Spherical::theta_r(int i, int j, double R, double r, double T, double g, double b){
    /* Compute Angular collision coordinate at specified distance, when distance of closest approach is known
        - i, j : Component indices
        - R : Distance of closest approach (m)
        - r : Separation distance (m)
        - T : Temperature (K)
        - g : Dimensionless relative velocity
        - b : Impact parameter (m)
    */
    if (b / sigma[i][j] > 100) return PI / 2;
    if (b / sigma[i][j] < 1e-3) return 0;
    const double dh{7.5e-3};
    const double tol{1e-6};

    const auto integrand = [&](double r_prime){return theta_integrand(i, j, T, r_prime, g, b);};
    const auto integrand1 = [&](double u){return (r / pow(u, 2)) * integrand(r / u);}; // Substitution u = r / r'
    const double r_lim = (1 + 5e-2) * R + (1 / (5 * g + 1)) * R;

    if (r < r_lim) {
        const auto lim_integrand = [&](double u){return (r_lim / pow(u, 2)) * integrand(r_lim / u);};
        const double correction = simpson(lim_integrand, 1e-12, 0.8, 100) + simpson(lim_integrand, 0.8, 1., 100);
        const double cutoff_error = tanh_sinh(lim_integrand, dh, tol) - correction;
        const double i1 = tanh_sinh(integrand1, dh, tol) - cutoff_error;
        return i1;
    }

    const double i1 = simpson(integrand1, 1e-12, 0.8, 100);
    const double i2 = simpson(integrand1, 0.8, 1., 100);
    return i1 + i2;
}

double Spherical::theta_lim(int i, int j, const double T, const double g){
    double b = 100 * sigma[i][j];
    double R = get_R(i, j, T, g, b);
    return theta_integral(i, j, T, R, g, b);
}

double Spherical::theta_integral(int i, int j, double T, double R, double g, double b){
    constexpr double dh{1.0e-3};
    const auto integrand = [&](double u){return (R / pow(u, 2)) * theta_integrand(i, j, T, R / u, g, b);};
    return tanh_sinh(integrand, dh);
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
    /* Compute distance of closest approach (m)
        - i, j : Component indices
        - T : Temperature (K)
        - g : Dimensionless relative velocity
        - b : Impact parameter (m)
    */
    if ((b / sigma[i][j] < 1e-5) || (g < 1e-6)) return get_R0(i, j, T, g); // Treat separately (set b = 0, and handle g = 0 case)
    // Newtons method
    double tol = 1e-5; // Relative to sigma[i][j]
    double init_guess_factor = 1.0;
    double r = init_guess_factor * b;
    double f = get_R_rootfunc(i, j, T, g, b, r);
    double dfdr = get_R_rootfunc_derivative(i, j, T, g, b, r);
    double next_r = r - f / dfdr;
    #ifdef DEBUG
        int niter = 0;
    #endif
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
        #ifdef DEBUG
            if (niter++ > 10000) throw std::runtime_error("get_R exceeded 10000 iterations! (T, g, b) : " 
                                                        + std::to_string(T) + ", " + std::to_string(g) + ", " + std::to_string(b / sigma[i][j]) + " : " + std::to_string(next_r / sigma[i][j]));
            if (isnan(next_r)) throw std::runtime_error("Encountered NAN in get_R! (T, g, b) : " 
                                                        + std::to_string(T) + ", " + std::to_string(g) + ", " + std::to_string(b));
        #endif
    }
    return next_r;
}

double Spherical::get_R0(int i, int j, double T, double g){
    double tol = 1e-6;
    double Ek = BOLTZMANN * T * pow(g, 2);
    double r = 0.9 * sigma[i][j];
    double f = potential(i, j, r) - Ek;
    while (f < 0){
        r *= 0.8;
        f = potential(i, j, r) - Ek;
    }
    double dfdr = potential_derivative_r(i, j, r);
    while (abs(f / (dfdr * sigma[i][j])) > tol){
        r -= f / dfdr;
        f = potential(i, j, r) - Ek;
        dfdr = potential_derivative_r(i, j, r);
    }
    return r;
}

double Spherical::chi(int i, int j, double T, double g, double b){
    /* Compute deflection angle
        - i, j : Component indices
        - T : Temperature (K)
        - g : Dimensionless relative velocity
        - b : Impact parameter (m)
    */
    if (b / sigma[i][j] > 100) return 0;
    double t = theta(i, j, T, g, b);
    return PI - 2.0 * t;
}


// Testing functions

double Spherical::omega_tester(int i, int j, int l, int r, double T, IntegrationParam& param){
    double w = w_integral_tester(i, j, T, l, r, param);
    if (i == j) return pow(sigma[i][j], 2) * sqrt((PI * BOLTZMANN * T) / m[i]) * w;
    return 0.5 * pow(sigma[i][j], 2) * sqrt(2 * PI * BOLTZMANN * T / (m0[i][j] * M[i][j] * M[j][i])) * w;
}

double Spherical::w_integral_tester(int i, int j, double T, int l, int r, IntegrationParam& param){
    const auto w_integrand_export = [&](double g, double b){return w_integrand(i, j, T, g, b, l, r);};
    double I = integrate2d(param.origin, param.end,
                        param.dg, param.db,
                        param.refinement_levels_g, param.refinement_levels_b,
                        param.subdomain_dblder_limit,
                        w_integrand_export);

    return I;
}