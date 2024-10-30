#include "Quantum.h"
#include "Integration/Integration.h"
#include <autodiff/forward/dual.hpp>

using dual = autodiff::dual;

double nested_sum(const std::function<double(size_t)>& term, const size_t upper, const size_t depth, const size_t maxdepth){
    double tot = 0.;
    if (depth == maxdepth) {
        for (size_t a = 0; a < upper; a++){
            tot += term(a);
        }
        return tot;
    }

    for (size_t a = 0; a < upper; a++){
        tot += term(a - maxdepth + depth) * nested_sum(term, a, depth + 1, maxdepth);
    }
    return tot;
}

double bessel_deriv(int n, double r, size_t kind){
    switch (kind) {
    case 1:
        if (n == 0) return - std::sph_bessel(1, r);
        return std::sph_bessel(n - 1, r) - (n + 1) * std::sph_bessel(n, r) / r;
    case 2:
        if (n == 0) return - std::sph_neumann(1, r);
        return std::sph_neumann(n - 1, r) - (n + 1) * std::sph_neumann(n, r) / r;
    default:
        throw std::runtime_error("Invalid spherical_bessel derivative kind!\n");
    }
}

dual spherical_bessel(dual r, int n, size_t kind){
    const auto prefactor = [&](int k) -> dual {return (Fac(n + k) / (Fac(k) * Fac(n - k))).eval() / pow(2 * r, k);};
    dual jcos{0.}, jsin{0.};
    switch (kind) {
        case 1:
            for (int k = abs(((n - 1) % 4) - 0); k <= n; k += 4){ jcos += prefactor(k);}
            for (int k = abs(((n - 1) % 4) - 1); k <= n; k += 4){ jsin -= prefactor(k);}
            for (int k = abs(((n - 1) % 4) - 2); k <= n; k += 4){ jcos -= prefactor(k);}
            for (int k = abs(((n - 1) % 4) - 3); k <= n; k += 4){ jsin += prefactor(k);}
            break;
        case 2:
            for (int k = abs((n % 4) - 0); k <= n; k += 4) { jcos += prefactor(k);}
            for (int k = abs((n % 4) - 1); k <= n; k += 4) { jsin -= prefactor(k);}
            for (int k = abs((n % 4) - 2); k <= n; k += 4) { jcos -= prefactor(k);}
            for (int k = abs((n % 4) - 2); k <= n; k += 4) { jsin += prefactor(k);}
            if ((n + 1) % 2 == 1) {jcos *= -1; jsin *= -1;}
            break;
        default: throw std::runtime_error("Invalid spherical_bessel kind : " + std::to_string(kind));
    }
    return (jcos * cos(r) + jsin * sin(r)) / r;
}

Quantum::Quantum(std::string comps) : Spherical(comps, true) {
    half_spin = std::vector<size_t>(Ncomps, 0);
    sigma = vector2d(Ncomps, vector1d(Ncomps, 0.));
    eps = vector2d(Ncomps, vector1d(Ncomps, 0.));
    for (size_t i = 0; i < Ncomps; i++){
        half_spin[i] = static_cast<size_t>(static_cast<double>(compdata[i]["spin"]) * 2 + 0.5);
        sigma[i][i] = 1;
        eps[i][i] = 1;
    }
}

dual2 Quantum::potential(int i, int j, dual2 r){
    return 4 * (pow(1. / r, 12) - pow(1. / r, 6));
}

double Quantum::potential(int i, int j, double r){
    return 4 * (pow(1. / r, 12) - pow(1. / r, 6));
}
double Quantum::potential_derivative_r(int i, int j, double r){
    return 4 * (6 * pow(1. / r, 7) - 12 * pow(1. / r, 13));
}
// double Quantum::potential_rr(int i, int j, double r){
//     return 0;
// }

double Quantum::JKWB_phase_shift(int i, int j, int l, double E){
    double red_mass = 1;
    double k = sqrt(2. * red_mass * E); //  / HBAR;
    double T = 50;
    double g = sqrt(E / (BOLTZMANN * T));
    double b = (l + 0.5) / k;
    
    double R = get_R(i, j, T, g, b);
    if (abs(R - b) < 1e-6) return 0;
    std::cout << "JKWB : " << b << ", " << R;
    while (1 - pow(b / R, 2) - potential(i, j, R) / E < 0) {
        R += 1e-6 * sigma[i][j];
        std::cout << " => " << R;
    }
    std::cout << std::endl;
    const auto integrand1 = [&](double r){return sqrt(1 - pow(b / r, 2) - potential(i, j, r) / E) - sqrt(1 - pow(b / r, 2));};
    double I1, I2;
    if (b < R){
        const auto integrand2 = [&](double r){return sqrt(1 - pow(b / r, 2));};
        I1 = simpson_inf(integrand1, R, 1.5 * R);
        I2 = - simpson(integrand2, b, R, 10);
    }
    else {
        const auto integrand2 = [&](double r){return sqrt(1 - pow(b / r, 2) - potential(i, j, r) / E);};
        I1 = simpson_inf(integrand1, b, 1.5 * b);
        I2 = simpson(integrand2, R, b, 10);
    }
    std::cout << " I1, I2 : " << I1 << ", " << I2 << std::endl;
    return k * (I1 + I2);
}

vector2d Quantum::wave_function(int i, int j, int l, const double E, const double r_end, const double step_size){
    const double k2 = 1.; // 2. * (m[i] * m[j] / (m[i] + m[j]))/ (pow(HBAR, 2));
    const double k = sqrt(k2 * E);
    const double s2 = pow(step_size, 2) / 12.;

    const auto g_fun = [&](double r_i){return (l * (l + 1)) / pow(r_i, 2) + k2 * (potential(i, j, r_i) - E);};
    const auto numerov_step = [&](vector1d& r_n, vector1d& psi_n, vector1d& g_n) -> void {
        size_t n_step = r_n.size();
        r_n.push_back(r_n.back() + step_size);
        g_n.push_back(g_fun(r_n.back()));

        double next_psi = (psi_n[n_step - 1] * (2. + 10. * s2 * g_n[n_step - 1]) - psi_n[n_step - 2] * (1. - s2 * g_n[n_step - 2])) / (1 - s2 * g_n[n_step]);
        psi_n.push_back(next_psi);
    };

    const auto local_phase_shift = [&](const vector1d& r, const vector1d& psi, const vector1d& g_vals) -> double {
        const size_t n_step = r.size() - 3;
        double A1 = 0.5 * (psi[n_step + 1] - psi[n_step - 1]);
        double A2 = 0.5 * (psi[n_step + 2] - psi[n_step - 2]);
        double B1 = s2 * (g_vals[n_step + 1] * psi[n_step + 1] - g_vals[n_step - 1] * psi[n_step - 1]);
        double B2 = s2 * (g_vals[n_step + 2] * psi[n_step + 2] - g_vals[n_step - 2] * psi[n_step - 2]);
        double dpsi = (16. / (21. * step_size)) * (- A1 + (37. * A2 / 32.) - (37. * B1 / 5.) - (17. * B2 / 40.)); // First derivative to order h^9

        double gamma = (dpsi / psi[n_step]) - (1. / r[n_step]);

        double jl = std::sph_bessel(l, k * r[n_step]);
        double djl = bessel_deriv(l, k * r[n_step], 1);
        double yl = std::sph_neumann(l, k * r[n_step]);
        double dyl = bessel_deriv(l, k * r[n_step], 2);

        return atan((k * djl - gamma * jl) / (k * dyl - gamma * yl)); // local phase shift
    }; 

    double r0 = 5;
    vector1d r = {r0, r0 + step_size};
    vector1d g_vals = {g_fun(r[0]), g_fun(r[1])};
    vector1d psi = {0., step_size};
    vector1d phase_shifts;
    for (size_t i = 0; i < 2; i++) numerov_step(r, psi, g_vals);
    while (r.back() < r_end){
        numerov_step(r, psi, g_vals);
        phase_shifts.push_back(local_phase_shift(r, psi, g_vals));
    }
    r = vector1d(r.begin() + 2, r.end() - 2);
    psi = vector1d(psi.begin() + 2, psi.end() - 2);
    vector2d out = {r, psi, phase_shifts};
    return out;
}

double Quantum::phase_shift(int i, int j, int l, double E){
    const double step_size{1e-6 * pow(E, 1. / 3.)};
    const double s2 = pow(step_size, 2) / 12.;
    // E *= BOLTZMANN;
    double k2 = 1.; // 2. * (m[i] * m[j] / (m[i] + m[j]))/ (pow(HBAR, 2));
    double k = sqrt(k2 * E);

    const auto g_fun = [&](double r_i){return (l * (l + 1)) / pow(r_i, 2) + k2 * (potential(i, j, r_i) - E);};
    const auto numerov_step = [&](std::array<double, 5>& r_n, std::array<double, 5>& psi_n, std::array<double, 5>& g_n) -> void {
        size_t n_step = 4;
        for (size_t i = 0; i < n_step; i++) {r_n[i] = r_n[i + 1]; psi_n[i] = psi_n[i + 1]; g_n[i] = g_n[i + 1];}
        r_n[n_step] = r_n[n_step - 1] + step_size;
        g_n[n_step] = g_fun(r_n[n_step]);

        double next_psi = (psi_n[n_step - 1] * (2. + 10. * s2 * g_n[n_step - 1]) - psi_n[n_step - 2] * (1. - s2 * g_n[n_step - 2])) / (1 - s2 * g_n[n_step]);
        psi_n[n_step] = next_psi;
    };

    double r0 = 0.5;
    std::array<double, 5> r = {0, 0, 0, r0, r0 + step_size};
    std::array<double, 5> g_vals = {0, 0, 0, g_fun(r[3]), g_fun(r[4])};
    std::array<double, 5> psi = {0, 0, 0, 0, step_size};

    const auto local_phase_shift = [&](const std::array<double, 5>& r, const std::array<double, 5>& psi, const std::array<double, 5>& g_vals) -> double {
        const size_t n_step = 2;
        double A1 = 0.5 * (psi[n_step + 1] - psi[n_step - 1]);
        double A2 = 0.5 * (psi[n_step + 2] - psi[n_step - 2]);
        double B1 = s2 * (g_vals[n_step + 1] * psi[n_step + 1] - g_vals[n_step - 1] * psi[n_step - 1]);
        double B2 = s2 * (g_vals[n_step + 2] * psi[n_step + 2] - g_vals[n_step - 2] * psi[n_step - 2]);
        double dpsi = (16. / (21. * step_size)) * (- A1 + (37. * A2 / 32.) - (37. * B1 / 5.) - (17. * B2 / 40.)); // First derivative to order h^9

        double gamma = (dpsi / psi[n_step]) - (1. / r[n_step]);

        double jl = std::sph_bessel(l, k * r[n_step]);
        double djl = bessel_deriv(l, k * r[n_step], 1);
        double yl = std::sph_neumann(l, k * r[n_step]);
        double dyl = bessel_deriv(l, k * r[n_step], 2);

        return atan((k * djl - gamma * jl) / (k * dyl - gamma * yl)); // local phase shift
    }; 

    for (size_t i = 0; i < 5; i++) numerov_step(r, psi, g_vals); // Need to fill the arrays before we can start computing phase shifts.
    double prev_delta;
    double new_delta = 2. * PI;
    do {
        prev_delta = new_delta;
        numerov_step(r, psi, g_vals);
        new_delta = local_phase_shift(r, psi, g_vals);
    } while (abs(new_delta - prev_delta) > 1e-8);
    return new_delta;
}

double Quantum::cross_section_A(int n, int l, size_t k){
    switch (k) {
    case 0:
        return 1;
    case 1:
        return static_cast<double>(n * (2 * pow(l, 2) + 2 * (n * l - l) - n + 1)) / static_cast<double>(2 * (2 * l - 1) * (2 * (l + n) - 1));
    case 2:
        return (static_cast<double>(n * (n - 1)) / 8.) 
                * static_cast<double>(4. * pow(l, 4) + (n - 3) * (8. * pow(l, 3) + (3 - 8 * l) * (n - 2)) + 4 * pow(l, 2) * (pow(n, 2) - 8 * n + 13))
                / static_cast<double>((2 * l - 3) * (2 * l - 1) * (2 * l + 2 * n - 5) * (2 * l + 2 * n - 3));
    default:
        const auto A_term = [&](size_t a) -> double {return static_cast<double>(pow(l + a, 2)) 
                                                                / static_cast<double>((2 * (l + a) - 1) * (2 * (l + a) + 1));
                                                            };
        return nested_sum(A_term, n - k, 1, k);
    }
}

double Quantum::cross_section_kernel(int i, int j, double n, double l, double E){
    if (n == 0){
        return (2 * l + 1) * pow(sin(phase_shift(i, j, l, E)), 2);
    }
    double A, B;
    double q = 0;;
    for (size_t k = 0; (k < 3) && (2 * k + 1 <= n); k++){
        B = (2 * l + 1) * pow(sin(phase_shift(i, j, l, E) - phase_shift(i, j, l + n - 2 * k, E)), 2);
        double B_tmp = 0;
        for (size_t p = 1; p <= n - 2 * k; p++){
            B_tmp += (l + p) / (2 * (l + p) - 1);
        }
        B *= B_tmp;
        A = cross_section_A(n, l, k);
        q += A * B;
    }
    return q;
}

double Quantum::cross_section(int i, int j, int n, double E){
    if (is_singlecomp) i = j;
    
    double even_prefactor, odd_prefactor;
    if (i != j) {
        even_prefactor = odd_prefactor = 1.;
    }
    else {
        double S = half_spin[i] * 2;
        even_prefactor = (S + 1) / (2 * S + 1);
        odd_prefactor = S / (2 * S + 1);
        if (half_spin[i] % 2 != 0) { // Odd Half-Integer spin, Fermions: Swap prefactors
            double tmp = even_prefactor;
            even_prefactor = odd_prefactor;
            odd_prefactor = tmp;
        }
    }
    double Q = 0;
    if (even_prefactor > 0.){
        double Q_even = 0;
        int l_even = 0;
        double q_even;
        do {
            q_even = cross_section_kernel(i, j, n, l_even, E);
            Q_even += q_even;
            l_even += 2;
        } while (abs(q_even) > 1e-6 * Q_even);
        Q += even_prefactor * Q_even;
    } 
    if (odd_prefactor > 0.){
        double Q_odd = 0;
        int l_odd = 1;
        double q_odd;
        do {
            q_odd = cross_section_kernel(i, j, n, l_odd, E);
            Q_odd += q_odd;
            l_odd += 2;
        } while (abs(q_odd) > 1e-6 * Q_odd);
        Q += odd_prefactor * Q_odd;
    } 
    double kappa_mul_E_squared = 2. * m[i] * m[j] / ((m[i] + m[j]) * HBAR); // The factor E has been included in omega instead
    Q *= 4. * PI / kappa_mul_E_squared;
    return Q;
}

double Quantum::quantum_omega(int i, int j, int n, int s, double T){
    double beta = 1. / (T * BOLTZMANN);
    double sfac = 1;
    for (int si = 2; si <= s + 1; si++){
        sfac *= si;
    }
    const auto kernel = [&](double E){return cross_section(i, j, n, E) * exp(- beta * E) * pow(beta * E, s + 2) / sfac;};
    return simpson(kernel, 0, 10, 50);
}
