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

double sph_bessel(int n, double r){
    if (r < 2000) return std::sph_bessel(n, r);
    return sin(r - n * PI / 2) / r;
}

double sph_neumann(int n, double r){
    if (r < 2000) return std::sph_neumann(n, r);
    return - cos(r - n * PI / 2) / r;
}

double bessel_deriv(int n, double r, size_t kind){
    switch (kind) {
    case 1:
        if (n == 0) return - sph_bessel(1, r);
        return sph_bessel(n - 1, r) - (n + 1) * sph_bessel(n, r) / r;
    case 2:
        if (n == 0) return - sph_neumann(1, r);
        return sph_neumann(n - 1, r) - (n + 1) * sph_neumann(n, r) / r;
    default:
        throw std::runtime_error("Invalid spherical_bessel derivative kind!\n");
    }
}

Quantum::Quantum(std::string comps) : Spherical(comps, true) {
    half_spin = std::vector<unsigned int>(Ncomps, 0);
    for (size_t i = 0; i < Ncomps; i++){
        half_spin[i] = static_cast<size_t>(static_cast<double>(compdata[i]["spin"]) * 2 + 0.5);
    }
}

vector2d Quantum::get_de_boer(){
    vector2d de_boer(Ncomps, vector1d(Ncomps));
    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = 0; j < Ncomps; j++){
            de_boer[i][j] = get_de_boer(i, j);
        }
    }
    return de_boer;
}
double Quantum::get_de_boer(int i, int j){
    return PLANCK / (sigma[i][j] * sqrt(2. * red_mass[i][j] * eps[i][j]));
}

void Quantum::set_de_boer_mass(int i, double de_boer){
    m[i] = pow(PLANCK / (sigma[i][i] * de_boer), 2) / eps[i][i];
    set_masses(); // Set- methods are responsible for clearing caches
}

double Quantum::de_broglie_wavelength(int i, double T){
    return PLANCK / sqrt(2. * PI * m[i] * BOLTZMANN * T);
}

double Quantum::JKWB_upper_E_limit(int i, int j){
    double deBoer = get_de_boer(i, j);
    double margin = 5; // Some suitably large number ...
    double E_red = pow(margin * deBoer, 2) + 1;
    return E_red;
}

double Quantum::JKWB_phase_shift(int i, int j, int l, double E){
    E *= eps[i][j];
    double k = sqrt(2. * red_mass[i][j] * E) / HBAR;
    double T = 50;
    double g = sqrt(E / (BOLTZMANN * T));
    double b = ((l + 0.5) / k) / sigma[i][j];
    
    double R = get_R(i, j, T, g, b * sigma[i][j]) / sigma[i][j];
    while (1 - pow(b / R, 2) - potential(i, j, R * sigma[i][j]) / E < 0) {
        R += 1e-3;
    }
    
    const auto integrand1 = [&](double r){return sqrt(1 - pow(b / r, 2) - potential(i, j, r * sigma[i][j]) / E) - sqrt(1 - pow(b / r, 2));};
    if (abs(R - b) < 1e-6){
        double lower_lim = (R > b) ? R : b;
        double I = k * simpson_inf(integrand1, lower_lim, 1.5 * lower_lim) * sigma[i][j];
        return I;
    }

    double I1, I2;
    if (b < R){
        const auto integrand2 = [&](double r){return sqrt(1 - pow(b / r, 2));};
        I1 = simpson_inf(integrand1, R, 1.5 * R);
        I2 = - simpson(integrand2, b, R, 10);
    }
    else {
        const auto integrand2 = [&](double r){
            double I = sqrt(1 - pow(b / r, 2) - potential(i, j, r * sigma[i][j]) / E);
            return I;
        };
        I1 = simpson_inf(integrand1, b, 1.5 * b);
        I2 = simpson(integrand2, R, b, 10);
    }
    return k * (I1 + I2) * sigma[i][j];
}

vector2d Quantum::wave_function(int i, int j, int l, double E, double r_end, double step_size){
    E *= eps[i][j]; step_size *= sigma[i][j]; r_end *= sigma[i][j];
    const double k2 = 2. * red_mass[i][j] / (pow(HBAR, 2));
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

        double jl = sph_bessel(l, k * r[n_step]);
        double djl = bessel_deriv(l, k * r[n_step], 1);
        double yl = sph_neumann(l, k * r[n_step]);
        double dyl = bessel_deriv(l, k * r[n_step], 2);

        return atan((k * djl - gamma * jl) / (k * dyl - gamma * yl)); // local phase shift
    }; 

    double r0 = 0.5 * sigma[i][j];
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

double Quantum::quantum_phase_shift(int i, int j, int l, double E){
    // Input is reduced energy (E = E(Joule) / eps[i][j])
    const double step_size = ((E < 1e-3) ? 1e-4 : 1e-3 * pow(E, 1. / 3.)) * sigma[i][j];
    E *= eps[i][j];
    const double s2 = pow(step_size, 2) / 12.;
    double k2 = (2. * red_mass[i][j] / pow(HBAR, 2));
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

    double r0 = 0.5 * sigma[i][j];
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

        double jl = sph_bessel(l, k * r[n_step]);
        double djl = bessel_deriv(l, k * r[n_step], 1);
        double yl = sph_neumann(l, k * r[n_step]);
        double dyl = bessel_deriv(l, k * r[n_step], 2);

        return atan((k * djl - gamma * jl) / (k * dyl - gamma * yl)); // local phase shift
    }; 

    size_t nsteps_init = static_cast<size_t>((5 * sigma[i][j] - r0) / step_size);
    for (size_t i = 0; i < nsteps_init; i++) numerov_step(r, psi, g_vals); // Do some steps computing phase shifts.
    double prev_delta = 2 * PI;
    double new_delta = - prev_delta;
    double current_err;
    // ---------------------------------- ITERATION TO SOLVE FOR PHASE SHIFT --------------------------------- //
    // --- We integrate the wave function outwards, and at regular intervals compute the local phase shift (LPS). 
    // --- When the change between two checks is below a tolerance, we terminate.
    // --- Because the LPS oscillates, the change in the LPS between steps is periodic, with the same period as 
    // --- the wave function. Therefore, we do as follows:
    // --- (1) Iterate out to the first minimum in change, computing the LPS at each step. 
    // --- (2) Iterate out to the first maximum in change, computing the LPS at each step.
    // --- (3) Iterate one half period of the wave function before re-computing the LPS. Terminate when the change between two periods 
    // ---      is below the tolerance.
    // --- This ensures that we are checking the LPS at the point in the period wher it changes most rapidly, to ensure that we
    // --- do not terminate the iteration prematurely because we are at a point in the period where it happens to change slowly. 
    
    // Find minimum in change (check at each step)
    do { 
        current_err = abs(new_delta - prev_delta); 
        prev_delta = new_delta;
        numerov_step(r, psi, g_vals);
        new_delta = local_phase_shift(r, psi, g_vals);
    } while (abs(new_delta - prev_delta) < current_err); 

    // Find maximum in change (check at each step)
    do {
        current_err = abs(new_delta - prev_delta); 
        prev_delta = new_delta;
        numerov_step(r, psi, g_vals);
        new_delta = local_phase_shift(r, psi, g_vals);
    } while (abs(new_delta - prev_delta) > current_err);

    // Find distance where change is below tolerance, only check once per period.
    do {
        size_t nsteps_period = static_cast<size_t>(PI / (k * step_size));
        prev_delta = new_delta;
        for (size_t i = 0; i < nsteps_period; i++) numerov_step(r, psi, g_vals);
        new_delta = local_phase_shift(r, psi, g_vals);
    } while (abs(new_delta - prev_delta) > 1e-5);
    return new_delta;
}

double Quantum::phase_shift(int i, int j, int l, double E){
    if (E > 50) return JKWB_phase_shift(i, j, l, E);
    return quantum_phase_shift(i, j, l, E);
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
    // Input is reduced energy (E = E(Joule) / eps[i][j])
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

double Quantum::cross_section(int i, int j, const int n, const double E){
    // Input is reduced energy (E = E(Joule) / eps[i][j])
    if (is_singlecomp) i = j;
    
    double even_prefactor, odd_prefactor;
    if (i != j) {
        even_prefactor = odd_prefactor = 1.;
    }
    else {
        double S = half_spin[i] / 2.;
        even_prefactor = 1; // (S + 1) / (2 * S + 1);
        odd_prefactor = 0; // S / (2 * S + 1);
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

    double k2 = 2. * red_mass[i][j] * E * eps[i][j] / pow(HBAR, 2);
    Q *= 4. * PI / k2;
    return Q;
}

double Quantum::classical_cross_section(int i, int j, int l, double E){
    return Spherical::cross_section(i, j, l, E);
}

double Quantum::quantum_omega(int i, int j, int n, int s, double T){
    double beta = 1. / (T * BOLTZMANN);
    // double sfac = 1;
    // for (int si = 2; si <= s + 1; si++){
    //     sfac *= si;
    // }
    std::cout << "Computing : " << n << ", " << s << ", " << T << std::endl;
    const auto kernel = [&](double Eb){
        return exp(- Eb) * pow(Eb, s + 1) * cross_section(i, j, n, Eb / (beta * eps[i][j]));
    };
    double maxpoint = s + 1.;
    double I = simpson_inf(kernel, 1e-3, maxpoint / 2);
    std::cout << "Final (" << n << ", " << s << ", " << T << ") : " << I << std::endl;
    return sqrt(PI / (2. * beta * red_mass[i][j])) * I;
}

double Quantum::classical_omega(int i, int j, int l, int r, double T){
    std::cout << "Classical omega ..." << std::endl;
    double val = Spherical::omega(i, j, l, r, T);
    std::cout << "Finished " << l << ", " << r << ", " << T << std::endl;
    return val;
}

double Quantum::omega(int i, int j, int l, int r, double T) {
    OmegaPoint point = get_omega_point(i, j, l, r, T);
    OmegaPoint sympoint = get_omega_point(j, i, l, r, T);
    const std::map<OmegaPoint, double>::iterator pos = omega_map.find(point);
    if (pos == omega_map.end()){
        double val = (quantum_is_active) ? quantum_omega(i, j, l, r, T) : classical_omega(i, j, l, r, T);
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
    return pos->second;
}

void Quantum::set_quantum_active(bool active){
    if (active != quantum_is_active){
        clear_all_caches();
    }
    quantum_is_active = active;
}