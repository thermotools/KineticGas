/* 
A bunch of quantum mechanics pls help
    
References:
Second virial : 
    Kilpatrick et al., Second virial coefficients of He3 and He4, Phys. Rev. 94 (5), 1954
        DOI: https://doi.org/10.1103/PhysRev.94.1103 

Quantal cross sections, Phase shifts, JKWB approximation : 
    Meeks et al., On the quantum cross sections in dilute gases, J. Chem. Phys. 100 (5), 1994
        DOI: https://doi.org/10.1063/1.466370

    Lang et al., Thermophysical properties of argon gas from improved two-body interaction potential, Phys. Rev. A 109 (2024)
        DOI: 10.1103/PhysRevA.109.052803

    Blatt, Practical points concerning the solution of the Schrödinger equation, J. Comp. Phys. 1 (3), 1967
        DOI: https://doi.org/10.1016/0021-9991(67)90046-0

    Munn et al., Some Aspects of the Quantal and Semiclassical Calculation of Phase shifts and 
        cross sections for molecular scattering and transport, J. Chem. Phys. 41 (12), 1964
        DOI: https://doi.org/10.1063/1.1725845

    Mehl et al. Ab initio transport coefficients of gaseous hydrogen, Int. J. Thermophys. 31, 2010
        DOI: https://doi.org/10.1007/s10765-009-0697-9

    Hurly & Mehl, He4 Thermophysical properties: New ab initio calculations, J. Res. Natl. Inst. Stand. Technol. 112, 2007
        DOI: 10.6028/jres.112.006

See Also:
    Aziz & Slaman, An Analysis of the ITS-90 Relations for the Non-Ideality of 3He and 4He: Recommended Relations 
        Based on a New Interatomic Potential for Helium, Metrologia 27 (211), 1990
        DOI: 10.1088/0026-1394/27/4/005

    Sidharth, Phase shifts in the collision of massive particles, J. Math. Phys. 30 (3), 1989
        DOI: https://doi.org/10.1063/1.528381
*/ 

#include "Quantum.h"
#include "Integration/Integration.h"
#include <autodiff/forward/dual.hpp>
#include <fstream>
#include <sstream>
#include <shared_mutex>

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

// Thse wrappers existed because I was handling asymptotic cases (r >> n) using an 
// asymptotic expansion. It turns out that isn't nessecary, and introduces numerical error
// which causes instabilities in the phase shifts.
// If we need to go to very large r at some point, that will need to be handled.
inline double sph_bessel(int n, double r){
    return std::sph_bessel(n, r);
    // return sin(r - n * PI / 2) / r;
}

inline double sph_neumann(int n, double r){
    return std::sph_neumann(n, r);
    // return - cos(r - n * PI / 2) / r;
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

Quantum::Quantum(std::string comps_) 
    : Spherical(comps_, true), 
    E_bound(Ncomps, std::vector<vector2d>(Ncomps)),
    half_spin(Ncomps, 0),
    spin(Ncomps, 0),
    rot_ground_state(Ncomps, 0)
{   
    for (size_t i = 0; i < Ncomps; i++){
        const auto qdata = compdata[i]["Quantum"];
        bool has_quantum_params = static_cast<bool>(qdata["has_quantum_params"]);
        if (!has_quantum_params) {
            std::cout << "WARNING : Component " << comps[i] << " does not support Quantum!\n";
            std::cout << "\tDeactivating Quantum as default behaviour. You can use `set_quantum_active` if you really want to ...\n";
            quantum_is_active = false;
            quantum_supported = false;
            continue;
        }
        half_spin[i] = static_cast<unsigned int>(qdata["half_spin"]);
        spin[i] = static_cast<double>(half_spin[i]) / 2.;
        rot_ground_state[i] = static_cast<unsigned int>(qdata["rot_ground_state"]);
        bool has_E_bound_file = static_cast<bool>(qdata["has_E_bound_file"]);
        if (has_E_bound_file){
            E_bound[i][i] = get_E_bound_from_file(comps[i]);
        }
        else {
            E_bound[i][i] = qdata["E_bound_div_k"];
        }

        if (E_bound[i][i].size() == 0){
            std::cout << "WARNING : No bound state energies supplied for component " << i << std::endl; 
        }

        for (size_t v = 0; v < E_bound[i][i].size(); v++){
            for (size_t l = 0; l < E_bound[i][i][v].size(); l++){
                E_bound[i][i][v][l] *= BOLTZMANN;
            }
        }
    }
    std::string phase_shift_filepath = get_fluid_dir() + "/E_bound/" + comps[0] + ".phase_shifts";
    std::ifstream phase_shift_file(phase_shift_filepath);
    if (phase_shift_file.is_open()){
        json phase_shift_data = json::parse(phase_shift_file);
        absolute_phase_shift_map = phase_shift_data.get<std::map<int, vector2d>>();
    } 
    else {
        std::cout << "No phase shift data found for component : " << comps[0] << std::endl;
    }
}

vector2d Quantum::get_E_bound_from_file(const std::string& comp){
    std::string filepath = get_fluid_dir() + "/E_bound/" + comp + ".dat";
    std::ifstream file(filepath);
    if (!file.is_open()) throw std::runtime_error("Failed to open E_bound file: " + filepath);
    vector2d data;
    std::string line;
    size_t ri{0}, ci{0};
    while (std::getline(file, line)){
        vector1d row;
        std::stringstream ss(line);
        std::string val;
        while (std::getline(ss, val, ',')){
            try {
                row.push_back(std::stod(val)); // Convert to double
            } catch (const std::invalid_argument& e) {
                std::cout << "In file : " << filepath << "\n\tInvalid value on (row, col) : (" << ri << ", " << ci << ")\n\tValue : " << val << std::endl;
                throw e;
            }
            ci++;
        }
        data.push_back(row);
        ri++; ci = 0;
    }
    return data;
}

void Quantum::clear_all_caches(){
    Spherical::clear_all_caches();
    phase_shift_map.clear();
}

int Quantum::get_interaction_statistics(int i, int j){
    if ( (is_singlecomp) ) i = j;
    if ( i != j ) return StatisticType::Boltzmann;
    if ( (half_spin[i] % 2) == 0 ) return StatisticType::BoseEinstein;
    return StatisticType::FermiDirac;
}

std::array<double, 2> Quantum::get_symmetry_weights(int i, int j){
    
    std::array<double, 2> wts;
    if (is_singlecomp) i = j;
    if (i != j){
        wts[0] = wts[1] = 0.5;
    }
    else {
        wts[0] = ((spin[i] + 1) * (rot_ground_state[i] + 1) + spin[i] * rot_ground_state[i]) / ((2 * spin[i] + 1) * (2 * rot_ground_state[i] + 1));
        wts[1] = ((spin[i] + 1) * rot_ground_state[i] + spin[i] * (rot_ground_state[i] + 1)) / ((2 * spin[i] + 1) * (2 * rot_ground_state[i] + 1));
    }
    if (get_interaction_statistics(i, j) == StatisticType::FermiDirac){
        std::swap(wts[0], wts[1]);
    }
    return wts;
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

double Quantum::JKWB_upper_E_limit(int i, int j){
    double deBoer = get_de_boer(i, j);
    double margin = 5; // Some suitably large number ...
    double E_red = pow(margin * deBoer, 2) + 1;
    return E_red;
}

double Quantum::JKWB_phase_shift(int i, int j, int l, double E){
    if (E < 1e-10) return 0.;
    E *= eps[i][j];
    double k = sqrt(2. * red_mass[i][j] * E) / HBAR;
    double T = 50;
    double g = sqrt(E / (BOLTZMANN * T));
    double b = ((l + 0.5) / k) / sigma[i][j];
    
    constexpr bool verbose = false; // (k * sigma[i][j] > 65) && (k * sigma[i][j] < 80);
    if (verbose) std::cout << "JKWB : " << l << ", " << k * sigma[i][j] << std::endl;

    double R = get_R(i, j, T, g, b * sigma[i][j]) / sigma[i][j]; // r_classical_forbidden(i, j, l, E / eps[i][j])
    
    // const auto integrand = [&](double r){
    //     double v1 = 1 - pow(b / r, 2);
    //     double v2 = v1 - potential(i, j, r * sigma[i][j]) / E;
    //     double I = (v2 > 0) ? sqrt(v2) : 0;
    //     if (v1 > 0) I -= sqrt(v1);
    //     return I;
    //     };
    // 
    // double r0 = std::min(R, b);
    // double r1 = std::max(std::max(R, b), 1.5 * r0);
    // double I = simpson(integrand, r0, r1, 50);
    // double r2 = 1.1 * r1;
    // double Ival = integrand(r1);
    // while (abs(integrand(r2) / Ival) < 0.5) {
    //     r2 = 0.5 * (r1 + r2);
    //     std::cout << "Reduced r2 : " << r2 / r1 - 1. << std::endl;
    // }
    // double I_tail = simpson_inf(integrand, r1, r2, 1e-8, I);
    // if (verbose) std::cout << "Integrated : " << r0 << " => " << r1 << " => " << r2 << " : " << I << ", " << I_tail - I << std::endl;
    // return k * sigma[i][j] * (I_tail);

    const auto R_fun = [&](double r){return 1 - pow(b / r, 2) - potential(i, j, r * sigma[i][j]) / E;};
    // if ( (b > 100) && (potential(i, j, R * sigma[i][j]) / E > - 1e-12) ){
    //     R = b;
    // }
    double s0 = get_potential_root(i, j) / sigma[i][j];
    if (b > s0){
        double R2 = R;
        if (verbose) std::cout << "Increase R ? " << R_fun(R2) << std::endl;
        // while (R_fun(R2) < 0) {
        //     if (verbose) std::cout << "Increasing R2 : " << R2 << ", " << b << ", " << R2 / b - 1 << ", " << R_fun(R2) << std::endl;
        //     R2 = 0.5 * (b + R2);
        //     if (abs(R2 / b - 1) < 1e-10) {
        //         if (R_fun(b) < 0) {
        //             throw std::runtime_error("JKWB phase shift: Could not find integration limit...");
        //         }
        //         else {
        //             R2 = b; break;
        //         }
        //     }
        // }
        // if (verbose) std::cout << "Increased R2 : " << R2 << ", " << b << ", " << R2 / b - 1 << ", " << R_fun(R2) << std::endl;
        // if (R2 != R) R = bracket_positive(R_fun, R, R2);
        if (R_fun(R) < 0) R = bracket_positive(R_fun, R, b);
        if (verbose) std::cout << "Increasing R : " << R2 << ", " << b << ", " << R2 / b - 1 << ", " << R_fun(R) << std::endl;
    }
    else if (R_fun(R) < 0) {
        const auto dR_fun = [&](double r){return (2. / b) * pow(b / r, 3) - potential_derivative_r(i, j, r * sigma[i][j]) * (sigma[i][j] / E);};
        if (verbose) std::cout << "R_fun < 0 : " << R_fun(R) << std::endl;
        if (verbose) std::cout << "Increasing R (> b) : " << R << ", " << b << ", " << R / b - 1 << std::endl;
        double R2 = R;
        while (R_fun(R2) < 0){R2 += 1e-2;}
        if (verbose) std::cout << "Increasing R (> b) : " << R << ", " << b << ", " << R / b - 1 << std::endl;
        R = bracket_positive(R_fun, R, R2);
        if (verbose) std::cout << "Increasing R (> b) : " << R << ", " << b << ", " << R / b - 1 << std::endl;
        if (verbose) std::cout << "R_fun > 0 : " << R_fun(R) << std::endl;
    }
    
    const auto integrand1 = [&](double r){return sqrt(1 - pow(b / r, 2) - potential(i, j, r * sigma[i][j]) / E) - sqrt(1 - pow(b / r, 2));};
    // if (abs(R - b) < 1e-12){
    //     double lower_lim = (R > b) ? R : b;
    //     double I = k * simpson_inf(integrand1, lower_lim, 1.5 * lower_lim) * sigma[i][j];
    //     return I;
    // }
    if (verbose) std::cout << "JKWB : " << l << ", " << k * sigma[i][j] << ", " << b << ", " << R << ", " << R / b - 1 << std::endl;
    double I1, I2;
    if (b <= R){
        const auto integrand2 = [&](double r){return sqrt(1 - pow(b / r, 2));};
        double R2 = 1.5 * R;
        double I0 = integrand1(R);
        while (integrand1(R2) < 0.5 * I0 && (R2 / R - 1) > 1e-6){
            R2 = 0.5 * (R + R2);
            if (verbose) std::cout << "Reduced upper R limit : " << R2 / R - 1 << ", " << integrand1(R2) << std::endl;
        }
        if (verbose) std::cout << "Upper R limit : " << R2 / R - 1 << ", " << integrand1(R2) << std::endl;
        I1 = simpson_inf(integrand1, R, R2); // (R2 / R - 1 > 1e-6) ?  : simpson(integrand1, R, 3 * R2, 100);
        I2 = (b == R) ? 0 : - simpson(integrand2, b, R, 10);
    }
    else {
        const auto integrand2 = [&](double r){
            double I = sqrt(1 - pow(b / r, 2) - potential(i, j, r * sigma[i][j]) / E);
            return I;
        };
        double b2 = 1.5 * b;
        double I0 = integrand1(b);
        while (integrand1(b2) < 0.5 * I0 && (b2 / b - 1) > 1e-3){
            b2 = 0.5 * (b + b2);
            if (verbose) std::cout << "Reduced upper b limit : " << b2 / b - 1 << ", " << integrand1(b2) << std::endl;
        }
        if (verbose) std::cout << "Upper b limit : " << b2 / b - 1 << ", " << integrand1(b2) << ", " << I0 << std::endl;
        I1 = simpson_inf(integrand1, b, b2); // (b2 / b - 1 > 1e-3) ?  : simpson(integrand1, b, 3 * b2, 100);
        I2 = simpson(integrand2, R, b, 10);
    }
    if (verbose) std::cout << "\t =>" << I1 << ", " << I2 << " : " << k * (I1 + I2) * sigma[i][j] << std::endl;
    return k * (I1 + I2) * sigma[i][j];
}

vector2d Quantum::wave_function(int i, int j, int l, double E, double r_end, double step_size){
    E *= eps[i][j]; step_size *= sigma[i][j]; r_end *= sigma[i][j];
    const double k2 = 2. * red_mass[i][j] / (pow(HBAR, 2));
    const double k = sqrt(k2 * E);
    if (step_size <= 0){
        step_size = 0.025 * std::min(1 / k, sigma[i][j]);
    }
    const double s2 = pow(step_size, 2) / 12.;

    const auto g_fun = [&](double r_i){return (l * (l + 1)) / pow(r_i, 2) + k2 * (potential(i, j, r_i) - E);};
    const auto numerov_step = [&](vector1d& r_n, vector1d& psi_n, vector1d& g_n) -> void {
        size_t n_step = r_n.size();
        r_n.push_back(r_n.back() + step_size);
        g_n.push_back(g_fun(r_n.back()));

        double next_psi = (psi_n[n_step - 1] * (2. + 10. * s2 * g_n[n_step - 1]) - psi_n[n_step - 2] * (1. - s2 * g_n[n_step - 2])) / (1 - s2 * g_n[n_step]);
        psi_n.push_back(next_psi);
    };

    const auto local_phase_shift = [&](const vector1d& r, const vector1d& psi, const vector1d& g_vals, bool cout) -> double {
        const size_t n_step = r.size() - 3;
        double A1 = 0.5 * (psi[n_step + 1] - psi[n_step - 1]);
        double A2 = 0.5 * (psi[n_step + 2] - psi[n_step - 2]);
        double B1 = s2 * (g_vals[n_step + 1] * psi[n_step + 1] - g_vals[n_step - 1] * psi[n_step - 1]);
        double B2 = s2 * (g_vals[n_step + 2] * psi[n_step + 2] - g_vals[n_step - 2] * psi[n_step - 2]);
        double dpsi = (16. / (21. * step_size)) * (- A1 + (37. * A2 / 32.) - (37. * B1 / 5.) - (17. * B2 / 40.)); // First derivative to order h^9
        // double d2psi = (psi[n_step + 1] - 2 * psi[n_step] + psi[n_step - 1]) / pow(step_size, 2);
        // double d2psi = (psi[n_step + 1] + psi[n_step - 1] - (15. / 8.) * psi[n_step] - (1. / 16.) * (psi[n_step + 2] + psi[n_step - 2])) * (4. / (3. * pow(step_size, 2)));
        
        double gamma = (dpsi / psi[n_step]) - (1. / r[n_step]);

        double jl = sph_bessel(l, k * r[n_step]);
        double djl = (l == 0) ? - sph_bessel(1, k * r[n_step]) : sph_bessel(l - 1, k * r[n_step]) - (l + 1) * jl / (k * r[n_step]); // bessel_deriv(l, k * r[n_step], 1);
        double yl = sph_neumann(l, k * r[n_step]);
        double dyl = (l == 0) ? - sph_neumann(1, k * r[n_step]) : sph_neumann(l - 1, k * r[n_step]) - (l + 1) * yl / (k * r[n_step]); // bessel_deriv(l, k * r[n_step], 2);

        if (cout) std::cout << "Calc phase shift : " << r.back() << ", " << k << ", " << gamma << ", " << jl << ", " << yl << ", " << djl << ", " << dyl << std::endl;
        return atan((k * djl - gamma * jl) / (k * dyl - gamma * yl)); // local phase shift // pow(HBAR, 2) * (d2psi / psi[n_step]) / (2. * red_mass[i][j]); // 
    }; 
    
    const double L_unt = sigma[i][j];
    const double E_unt = eps[i][j];
    const auto r0_fun = [&](double r_val){return (potential(i, j, r_val * L_unt) + l * (l + 1) / (k2 * pow(r_val * sigma[i][j], 2)) - E) / E_unt;};
    const auto r0_fun_d = [&](double r_val){return (L_unt * potential_derivative_r(i, j, r_val * L_unt) 
                                                    - 2 * L_unt * l * (l + 1) / (k2 * pow(r_val * sigma[i][j], 3))) / E_unt;};

    double r0_guess = 0.3;
    while (r0_fun(r0_guess) > 0){r0_guess *= 1.1;}
    while (r0_fun(r0_guess) < 0){r0_guess *= 0.9;}
    std::cout << "r0 guess : " << r0_guess << std::endl;
    double r_cf = newton(r0_fun, r0_fun_d, r0_guess) * sigma[i][j];
    double r0 = std::max(r_cf - 3 * (2 * PI / k), r_cf / 2);

    const auto particles_are_free = [&](double ri){return (abs(potential(i, j, ri) / E) < 1e-7)
                                                       && (abs(potential_derivative_r(i, j, ri) / (k * E)) < 1e-7);};

    if (particles_are_free(r_cf)) r_end = r_cf + 30 * (2 * PI / k);
    // r_end = r_cf + 3 * sigma[i][j];

    std::cout << "Checks : " << abs(potential(i, j, r_cf) / E) << ", " << abs(potential_derivative_r(i, j, r_cf) / (k * E)) << std::endl;
    std::cout << "Starting : " << r0 / sigma[i][j] << ", " << r_cf / sigma[i][j] << ", " << step_size / sigma[i][j] << ", " << (2 * PI / k) / sigma[i][j] << ", " << r_end / sigma[i][j] << std::endl;
    vector1d r = {r0, r0 + step_size};
    vector1d g_vals = {g_fun(r[0]), g_fun(r[1])};
    vector1d psi = {0, step_size};
    vector1d phase_shifts;

    int n = (l <= rot_ground_state[i]) ? 1 : 0;
    double delta, prev_delta;
    for (size_t i = 0; i < 2; i++) numerov_step(r, psi, g_vals);
    delta = local_phase_shift(r, psi, g_vals, false);
    while ( (!particles_are_free(r.back()) ) || (r.back() < r_end) ){
        prev_delta = delta;
        numerov_step(r, psi, g_vals);
        delta = local_phase_shift(r, psi, g_vals, false);
        phase_shifts.push_back(delta);
    }
    std::cout << "Terminated at : " << abs(potential(i, j, r.back()) / E) 
              << ", " << abs(potential_derivative_r(i, j, r.back()) / (k * E)) 
              << ", " << r.back() / sigma[i][j] << " / " << r_end / sigma[i][j] << std::endl;
    delta = local_phase_shift(r, psi, g_vals, true);
    std::cout << "Result : " << delta << std::endl;
    r = vector1d(r.begin() + 2, r.end() - 2);
    psi = vector1d(psi.begin() + 2, psi.end() - 2);
    vector2d out = {r, psi, phase_shifts};
    return out;
}

double Quantum::quantum_phase_shift(int i, int j, int l, double E, double& r_lev){
    // Input is reduced energy (E = E(Joule) / eps[i][j])
    E *= eps[i][j];
    double k2 = (2. * red_mass[i][j] / pow(HBAR, 2));
    double k = sqrt(k2 * E);
    const double step_size = 0.025 * std::min(1 / k, sigma[i][j]);
    const double s2 = pow(step_size, 2) / 12.; 

    // std::cout << "Start phase shift : " << E / eps[i][j] << ", " << step_size / sigma[i][j] << ", " << sigma[i][j] * k << std::endl;

    const auto g_fun = [&](double r_i){return (l * (l + 1)) / pow(r_i, 2) + k2 * (potential(i, j, r_i) - E);};
    const auto numerov_step = [&](std::array<double, 5>& r_n, std::array<double, 5>& psi_n, std::array<double, 5>& g_n) -> void {
        size_t n_step = 4;
        for (size_t i = 0; i < n_step; i++) {r_n[i] = r_n[i + 1]; psi_n[i] = psi_n[i + 1]; g_n[i] = g_n[i + 1];}
        r_n[n_step] = r_n[n_step - 1] + step_size;
        g_n[n_step] = g_fun(r_n[n_step]);

        psi_n[n_step] = (psi_n[n_step - 1] * (2. + 10. * s2 * g_n[n_step - 1]) - psi_n[n_step - 2] * (1. - s2 * g_n[n_step - 2])) / (1 - s2 * g_n[n_step]);
    };

    /*
        Determine the boundary of the classically forbidden region (solve r0_fun = 0).
        The wavefunction decays exponentially inside the classically forbidden region, so we set the 
        integration starting point a couple wavelengths inside this region.

        Note: If the classical boundary is at a very large separation, we can safely set the phase shift to zero.
        These collisions correspond to situations with very large angular momentum (large  l (quantum) or b (classical)),
        and very low energy. Thus, in a classical perspective, these are essentially collisions where the particles trajectories are 
        completely dominated by the centrifugal potential, meaning that they essentially behave as free particles, and we have zero 
        phase shift.

        In other words: When the classical boundary is at a very large distance (where the interaction potential is negligible), 
        and the wavefunction decays exponentially inside this boundary, there is a near-zero probability of finding the particles 
        at a separation where the interaction potential is non-negligible. Thus, we can safely treat them as free particles, 
        and assign them a vanishing phase-shift.
    */
    double r_cf = r_classical_forbidden(i, j, l, E / eps[i][j]);
    // True if the effect of the potential is negligible
    const auto particles_are_free = [&](double ri){return (abs(potential(i, j, ri) / E) < 1e-7)
                                                       && (abs(potential_derivative_r(i, j, ri) / (k * E)) < 1e-7);};

    // std::cout << "Checks : " << abs(potential(i, j, r_cf) / E) << ", " << abs(potential_derivative_r(i, j, r_cf) / (k * E)) << std::endl;

    if (particles_are_free(r_cf)) {r_lev = (l + 0.5) / k; return 0;}
    double r0 = std::max(r_cf - 3 * (2 * PI / k), r_cf / 2);

    // std::cout << "Starting :\n\tStart: " << r0 / sigma[i][j] 
    //                    << "\n\tR_cf : " << r_cf / sigma[i][j] 
    //                    << "\n\tStep : " << step_size / sigma[i][j] 
    //                    << "\n\tWave : " << (2 * PI / k) / sigma[i][j] 
    //                    << std::endl;

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
        // double d2psi = (psi[n_step + 1] - 2 * psi[n_step] + psi[n_step - 1]) / pow(step_size, 2);
        // double d2psi = (psi[n_step + 1] + psi[n_step - 1] - (15. / 8.) * psi[n_step] - (1. / 16.) * (psi[n_step + 2] + psi[n_step - 2])) * (4. / (3. * pow(step_size, 2)));
        
        double gamma = (dpsi / psi[n_step]) - (1. / r[n_step]);

        double jl = sph_bessel(l, k * r[n_step]);
        double djl = (l == 0) ? - sph_bessel(1, k * r[n_step]) : sph_bessel(l - 1, k * r[n_step]) - (l + 1) * jl / (k * r[n_step]); // bessel_deriv(l, k * r[n_step], 1);
        double yl = sph_neumann(l, k * r[n_step]);
        double dyl = (l == 0) ? - sph_neumann(1, k * r[n_step]) : sph_neumann(l - 1, k * r[n_step]) - (l + 1) * yl / (k * r[n_step]); // bessel_deriv(l, k * r[n_step], 2);
        // std::cout << "Calc phase shift : " << r.back() << ", " << k << ", " << gamma << ", " << jl << ", " << yl << ", " << djl << ", " << dyl << std::endl;
        return atan((k * djl - gamma * jl) / (k * dyl - gamma * yl)); // local phase shift // pow(HBAR, 2) * (d2psi / psi[n_step]) / (2. * red_mass[i][j]); // 
    }; 

    size_t nsteps_init = 2;
    double delta, prev_delta;
    int nl = get_levinson_multiple(i, j, l);
    r_lev = (l + 0.5) / k;
    int node_count = 0;
    bool levinson_node_found = false;

    for (size_t i = 0; i < nsteps_init; i++) {
        numerov_step(r, psi, g_vals);
        if (psi[3] * psi[4] < 0 && r.back() > 0.9 * sigma[i][j]){
            node_count++;
            if (node_count == nl + 2){
                double a = (psi[4] - psi[3]) / (r[4] - r[3]);
                double b = psi[4] - a * r[4];
                r_lev = - b / a;
                levinson_node_found = true;
            }
        }
    }
    delta = local_phase_shift(r, psi, g_vals);
    do {
        numerov_step(r, psi, g_vals);
        if (!levinson_node_found && psi[3] * psi[4] < 0 && r.back() > 0.9 * sigma[i][j]){
            node_count++;
            if (node_count == nl + 2){
                double a = (psi[4] - psi[3]) / (r[4] - r[3]);
                double b = psi[4] - a * r[4];
                r_lev = - b / a;
                levinson_node_found = true;
            }
        }
        // if (r.back() <= rmax_node_count && psi[3] * psi[4] < 0) r_lev += 1;
    } while ( ! particles_are_free(r.back())); // || (potential_derivative_r(i, j, r.back()) < 0) ); // while (l * (l + 1) / ( E * k2 * pow(r[4], 2)) > 1e-3); // (abs(new_delta - prev_delta) > 1e-5); // while (error > 1e-5); // 
    // std::cout << "Converged at r = " << r[4] / sigma[i][j] << std::endl;
    // std::cout << "Terminate phase shift : " << abs(potential(i, j, r.back()) / E) 
    //           << ", " << abs(potential_derivative_r(i, j, r.back()) / (k * E))  << ", " << r.back() / sigma[i][j] << std::endl;
    delta = local_phase_shift(r, psi, g_vals);
    // std::cout << "Result : (" << delta / PI << " + " << n << ") pi" << std::endl;
    return delta;
}

double Quantum::phase_shift(int i, int j, int l, double E){    
    if (l >= JKWB_l_limit || E > JKWB_E_limit) return JKWB_phase_shift(i, j, l, E);

    double r_lev;
    double delta = quantum_phase_shift(i, j, l, E, r_lev);
    return delta;
}

double Quantum::r_classical_forbidden(int i, int j, int l, double E){
    E *= eps[i][j];
    const double k2 = (2. * red_mass[i][j] / pow(HBAR, 2));
    const double L_unt = sigma[i][j];
    const double E_unt = eps[i][j];
    const auto r0_fun = [&](double r_val){return (potential(i, j, r_val * L_unt) + l * (l + 1) / (k2 * pow(r_val * L_unt, 2)) - E) / E_unt;};
    const auto r0_fun_d = [&](double r_val){return (L_unt * potential_derivative_r(i, j, r_val * L_unt) 
                                                    - 2 * L_unt * l * (l + 1) / (k2 * pow(r_val * L_unt, 3))) / E_unt;};

    double r0_guess = 0.3;
    while (r0_fun(r0_guess) > 0){r0_guess *= 1.1;}
    while (r0_fun(r0_guess) < 0){r0_guess *= 0.9;}
    return newton(r0_fun, r0_fun_d, r0_guess) * sigma[i][j];
}

double Quantum::absolute_phase_shift(int i, int j, int l, double E, double prev_delta){
    double delta = phase_shift(i, j, l, E);
    if ((E > JKWB_E_limit) || (l > JKWB_l_limit)) return delta;

    if (l <= rot_ground_state[i]) delta += PI;
    
    while (abs(delta - PI - prev_delta) < abs(delta - prev_delta)) delta -= PI;
    while (abs(delta + PI - prev_delta) < abs(delta - prev_delta)) delta += PI;
    return delta;
}

double Quantum::absolute_phase_shift(int i, int j, int l, double E){
    if ((E > JKWB_E_limit) || (l > JKWB_l_limit)) return JKWB_phase_shift(i, j, l, E);
    double delta = 0;
    // return phase_shift(i, j, l, E);
    double prev_delta;
    double Ei = 1e-6;
    double dE = (l == 0) ? 1e-2 : 5e-2;
    int n = 0; // (l <= rot_ground_state[i]) ? 1 : 0;
    std::cout << "E, l : " << E << " / " << Ei << ", " << l << std::endl;
    while (Ei <= E){
        prev_delta = delta;
        delta = absolute_phase_shift(i, j, l, Ei, delta);
        if ((delta < - PI / 4) and (prev_delta > PI / 4)) {
            n += 1;
        }
        else if ((delta > PI / 4) && (prev_delta < - PI / 4)){
            n -= 1;
        }
        // std::cout << "Ei, delta : " << E << ", " << Ei << ", " << delta << std::endl;
        Ei += dE;
    }
    Ei = E;
    prev_delta = delta;
    delta = absolute_phase_shift(i, j, l, Ei, delta);
    if ((delta < - PI / 4) and (prev_delta > PI / 4)) {
        n += 1;
    }
    else if ((delta > PI / 4) && (prev_delta < - PI / 4)){
        n -= 1;
    }
    // delta = absolute_phase_shift(i, j, l, E, delta);
    return delta + n * PI;
}

vector2d Quantum::absolute_phase_shifts(int i, int j, int l, double k_max){
    trace_absolute_phase_shifts(i, j, l, k_max);
    return absolute_phase_shift_map[l];
}

void Quantum::dump_phase_shift_map(){
    std::string filepath = get_fluid_dir() + "/E_bound/" + comps[0] + ".phase_shifts";
    std::ofstream file(filepath);
    if (!file.is_open()) throw std::runtime_error("Failed to open dumpfile: " + filepath);
    json dump_data(absolute_phase_shift_map);
    // dump_data["Total"] = stored_total_phase_shifts;
    file << std::setw(4) << dump_data << std::endl;
    std::cout << "Dumped absolute phase shifts to " << filepath << std::endl;
}

void Quantum::clear_phase_shift_maps(){
    absolute_phase_shift_map.clear();
    stored_total_phase_shifts.clear();
}

vector1d Quantum::get_levinson_r(int i, int j, int l, const vector1d k_vals){
    vector1d r_levinson;
    for (double k : k_vals) {
        r_levinson.push_back(get_levinson_multiple(i, j, l) + 1);
        if (k == 0) {
            r_levinson[0] = (l + 0.5) * sigma[i][j] / k;
            continue;
        }
        double E = pow(k * HBAR / sigma[i][j], 2) / (2. * red_mass[i][j] * eps[i][j]);
        quantum_phase_shift(i, j, l, E, r_levinson.back());
    }
    return r_levinson;
}

void Quantum::fill_absolute_phase_shifts(int i, int j, int l, double next_k, int& n, vector1d& k_vals, vector1d& phase_shifts, vector1d& r_levinson){
    const auto E_from_k = [&](double k_val){return pow(k_val * HBAR, 2) / (2. * red_mass[i][j] * eps[i][j]);};
    double r_lev;
    double next_delta = quantum_phase_shift(i, j, l, E_from_k(next_k / sigma[i][j]), r_lev);
    double delta, extrapol_delta;
    int local_n = n;
    do {
        local_n = n;
        const int n_step = k_vals.size() - 1;
        double k_step = next_k - k_vals.back();
        double dk1 = k_vals[n_step] - k_vals[n_step - 1];
        double dk2 = k_vals[n_step] - k_vals[n_step - 2];
        double dpdk = (phase_shifts[n_step] - phase_shifts[n_step - 1]) / dk1;
        double d2pdk2 = (phase_shifts[n_step - 2] - (1 - dk2 / dk1) * phase_shifts[n_step] - (dk2 / dk1) * phase_shifts[n_step - 1]) / (0.5 * dk2 * (dk2 - dk1));
        extrapol_delta = phase_shifts.back() + k_step * dpdk + pow(k_step, 2) * d2pdk2 / 2;
        double trial_delta = next_delta + local_n * PI;
        // std::cout << "Computing next(" << next_k << ", " << n << ") : " << phase_shifts.back() / PI << ", " << trial_delta / PI << ", " << extrapol_delta / PI
        //           << " ( " << dpdk << " / " << d2pdk2 << " )" << std::endl; 

        if ( (extrapol_delta > phase_shifts.back()) && (trial_delta + 1e-6 < phase_shifts.back()) ){
            local_n += 1;
            // std::cout << "Increased n to " << local_n << std::endl;
        }
        else if ((extrapol_delta < phase_shifts.back()) && (trial_delta - 1e-6 > phase_shifts.back())){
            local_n -= 1;
            // std::cout << "Decreased n to " << local_n << std::endl;
        }
        if (local_n != n){
            // std::cout << "\tNext(" << next_k << ", " << n << ") : " << phase_shifts.back() / PI << std::endl;
        }
        delta = next_delta + local_n * PI;
        if (k_step < 1e-6 || n > 15){
            // std::cout << "Reached maximum recursion ... " << k_vals.back() << ", " << k_step << std::endl;
            break; throw std::runtime_error("Unable to resolve absolute phase shift!");
        }
        if (abs(extrapol_delta - delta) > 0.01){
            // std::cout << "Backstepping to : " << k_vals.back() + k_step / 2 << std::endl;
            fill_absolute_phase_shifts(i, j, l, k_vals.back() + k_step / 2, n, k_vals, phase_shifts, r_levinson);
        }
    } while (abs(extrapol_delta - delta) > 0.01);
    // std::cout << "Completed recursion : " << next_k << std::endl;
    n = local_n;
    k_vals.push_back(next_k);
    phase_shifts.push_back(delta);
    r_levinson.push_back(r_lev);
}

void Quantum::fill_absolute_phase_shifts_tail(int i, int j, int l, double next_k, int& n, vector1d& k_vals, vector1d& phase_shifts){
    const auto E_from_k = [&](double k_val){return pow(k_val * HBAR, 2) / (2. * red_mass[i][j] * eps[i][j]);};
    double next_delta = phase_shift(i, j, l, E_from_k(next_k / sigma[i][j]));
    double delta, extrapol_delta;
    int local_n = n;
    constexpr bool verbose = false;
    if (verbose) std::cout << "Filling tail (" << next_k << ", " << n << ") : " << phase_shifts.back() / PI << std::endl;
    do {
        local_n = n;
        const int n_step = k_vals.size() - 1;
        double k_step = next_k - k_vals.back();
        double dk1 = k_vals[n_step] - k_vals[n_step - 1];
        double dk2 = k_vals[n_step] - k_vals[n_step - 2];
        double dpdk = (phase_shifts[n_step] - phase_shifts[n_step - 1]) / dk1;
        double d2pdk2 = (phase_shifts[n_step - 2] - (1 - dk2 / dk1) * phase_shifts[n_step] - (dk2 / dk1) * phase_shifts[n_step - 1]) / (0.5 * dk2 * (dk2 - dk1));
        extrapol_delta = phase_shifts.back() + k_step * dpdk + pow(k_step, 2) * d2pdk2 / 2;
        
        while (extrapol_delta < PI * (local_n - 0.5)) local_n -= 1;
        while (extrapol_delta > PI * (local_n + 0.5)) local_n += 1;

        if (local_n != n){
            if (verbose) std::cout << "Updated n (" << n << " => " << local_n << " )" << std::endl;
        }

        delta = next_delta + local_n * PI;
        if (verbose) {
            std::cout << "Computing next(" << next_k << ", " << local_n << ") : " << phase_shifts.back() / PI << ", " << delta / PI << ", " << extrapol_delta / PI
                   << ", " << abs(delta - extrapol_delta) / PI << " ( " << dpdk << " / " << d2pdk2 << " )" << std::endl; 
        }
        
        if (k_step < 1e-3 || n > 150){
            break; throw std::runtime_error("Unable to resolve absolute phase shift!");
        }
        if (abs(extrapol_delta - delta) > 0.2 * PI){
            if (verbose) std::cout << "Backstepping to : " << k_vals.back() + k_step / 2 << std::endl;
            fill_absolute_phase_shifts_tail(i, j, l, k_vals.back() + k_step / 2, n, k_vals, phase_shifts);
        }
    } while (abs(extrapol_delta - delta) > 0.2 * PI);
    if (verbose) std::cout << "Completed recursion : " << next_k << std::endl;
    n = local_n;
    k_vals.push_back(next_k);
    phase_shifts.push_back(delta);
}

int Quantum::get_levinson_multiple(int i, int j, int l){
    vector2d& Eb = E_bound[i][j];
    int nl = 0;
    for (size_t vib_i = 0; vib_i < Eb.size(); vib_i++){
        if (Eb[vib_i].size() > l) nl++;
    }
    return nl;
}

void Quantum::trace_absolute_phase_shifts(int i, int j, int l, double k_max){
    constexpr bool verbose = false;
    const auto E_from_k = [&](double k_val){return pow(k_val * HBAR, 2) / (2. * red_mass[i][j] * eps[i][j]);};
    int n = get_levinson_multiple(i, j, l);
    {
        std::lock_guard<std::mutex> loc(abs_phase_shift_map_mutex);
        const auto pos = absolute_phase_shift_map.find(l);
        if (pos == absolute_phase_shift_map.end()){
            vector1d k_vals = {0};
            vector1d phase_shifts = {PI * n};
            absolute_phase_shift_map[l] = vector2d({k_vals, phase_shifts});
            std::cout << "Compute for l, k_max = " << l << ", " << k_max << std::endl; 
        }
        else {
            if (absolute_phase_shift_map[l][0].back() >= k_max){
                std::cout << "Retrieve for l, k_max = " << l << ", " << k_max << std::endl;
            }
            else {
                std::cout << "Continue l = " << l << " : " << absolute_phase_shift_map[l][0].back() << " =>" << k_max << std::endl;
            }
        }
    }

    vector1d& k_vals = absolute_phase_shift_map[l][0];
    vector1d& phase_shifts = absolute_phase_shift_map[l][1];
    double r_lev;
    vector1d r_levinson(k_vals.size());

    if (k_vals.back() >= k_max) return;

    if (l >= JKWB_l_limit){
        std::cout << "Use JKWB (" << l << " > " << JKWB_l_limit << ")" << std::endl;
        double dk = k_max / 100.0;
        while ((k_vals.back() + dk) / k_max < 1 - 1e-12){
            k_vals.push_back(k_vals.back() + dk);
            phase_shifts.push_back(JKWB_phase_shift(i, j, l, E_from_k(k_vals.back() / sigma[i][j])));
        }
        k_vals.push_back(k_max);
        phase_shifts.push_back(JKWB_phase_shift(i, j, l, E_from_k(k_vals.back() / sigma[i][j])));
        return;
    }

    double k = (k_vals.size() == 1) ? 0.01 : k_vals.back();
    double dk = 0.001;

    if (k_vals.size() > 1){
        while (phase_shifts.back() > PI * (n + 0.5)) n += 1;
        while (phase_shifts.back() < PI * (n - 0.5)) n -= 1;
    }
    double delta{0};

    
    if (verbose) std::cout << "Starting trace : " << l << ", " << k_max << std::endl;
    if (verbose) std::cout << "n0, d0 : " << n << ", " << k_vals[0] << ", " << phase_shifts[0] / PI << std::endl;
    if (k_vals.size() == 1) r_levinson[0] = NAN;
    double prev_delta = quantum_phase_shift(i, j, l, E_from_k(k / sigma[i][j]), r_levinson.back());
    if (verbose) std::cout << "n, d : " << n << ", " << prev_delta / PI << std::endl;
    for (; k_vals.size() < 5;){
        k += dk;
        r_lev = n;
        delta = quantum_phase_shift(i, j, l, E_from_k(k / sigma[i][j]), r_lev);
        if (prev_delta < - PI / 4 && delta > 0){
            n -= 1;
        }
        else if (prev_delta > PI / 4 && delta < 0){
            n += 1;
        }
        // if (verbose) std::cout << "n, d : " << n << ", " << delta / PI << ", " << n + delta / PI << std::endl;
        prev_delta = delta;
        k_vals.push_back(k);
        phase_shifts.push_back(delta + n * PI);
        r_levinson.push_back(r_lev);
    }
    int n_step = k_vals.size() - 1;
    
    const Eigen::Vector3d b1 = {1, 0, 0};
    const Eigen::Vector3d b2 = {0, 1, 0};
    const Eigen::Vector3d b3 = {0, 0, 1};

    std::array<double, 3> dk_vals;
    dk_vals[0] = k_vals[2] - k_vals[1];
    dk_vals[1] = k_vals[3] - k_vals[1];
    dk_vals[2] = k_vals[4] - k_vals[1];
    Eigen::Matrix3d A(3, 3);
    for (unsigned int dki = 0; dki < 3; dki++){
        for (unsigned int fi = 0; fi < 3; fi++){
            A(fi, dki) = pow(dk_vals[dki], fi + 1);
        }
    }
    Eigen::PartialPivLU<Eigen::Matrix3d> A_piv(A);
    Eigen::Vector3d c = A_piv.solve(b1);
    double c0 = - (c(0) + c(1) + c(2));
    double dpdk = (c0 * phase_shifts[1] + c(0) * phase_shifts[2] + c(1) * phase_shifts[3] + c(2) * phase_shifts[4]);
    double d2pdk2 = 4 * (c0 * phase_shifts[1] + c(0) * phase_shifts[2] + c(1) * phase_shifts[3] + c(2) * phase_shifts[4]);
    double d3pdk3 = 6 * (c0 * phase_shifts[1] + c(0) * phase_shifts[2] + c(1) * phase_shifts[3] + c(2) * phase_shifts[4]);

    dk_vals[0] = k_vals[n_step] - k_vals[n_step - 1];
    dk_vals[1] = k_vals[n_step] - k_vals[n_step - 2];
    dk_vals[2] = k_vals[n_step] - k_vals[n_step - 3];
    for (unsigned int dki = 0; dki < 3; dki++){
        for (unsigned int fi = 0; fi < 3; fi++){
            A(fi, dki) = pow(dk_vals[dki], fi + 1);
        }
    }
    c = A_piv.compute(A).solve(b1);
    c0 = - (c(0) + c(1) + c(2));
    dpdk = - (c0 * phase_shifts[n_step] + c(0) * phase_shifts[n_step - 1] + c(1) * phase_shifts[n_step - 2] + c(2) * phase_shifts[n_step - 3]);
    
    // c = A_piv.solve(b2);
    // c0 = - (c(0) + c(1) + c(2));
    // double d2pdk2 = 2 * (c0 * phase_shifts[n_step] + c(0) * phase_shifts[n_step - 1] + c(1) * phase_shifts[n_step - 2] + c(2) * phase_shifts[n_step - 3]);

    c = A_piv.solve(b3);
    c0 = - (c(0) + c(1) + c(2));
    d3pdk3 = - 6 * (c0 * phase_shifts[n_step] + c(0) * phase_shifts[n_step - 1] + c(1) * phase_shifts[n_step - 2] + c(2) * phase_shifts[n_step - 3]);

    bool is_linear = (   (pow(dk, 3) / 6) * abs(d3pdk3) < 1e-3 
                      && (dpdk < 0) 
                      && (phase_shifts.back() < - PI / 8) 
                     );

    double extrapol_tol = 1e-2;
    double max_dk = 5e-1;
    double min_dk = 1e-2;
    if (l == 14){
        max_dk /= 2; min_dk /= 2; extrapol_tol /= 2;
    }
    bool resonance_found{false}, resonance_passed{false};
    bool prepare_jump{false}, done_jump{false};
    size_t n_jump = 0;
    bool worst_is_over{false};
    double dr0_k, dr0_k2, dr0_k3;
    // extrapol_tol = 1e-5;
    // max_dk = 1e-2;
    // min_dk = 1e-6;
    while ((k + max_dk < k_max) && (!is_linear)){
        // std::array<double, 3> extrapol_coeff = quadric_extrapolate_coeff(k_vals, phase_shifts);
        n_step = k_vals.size() - 1;

        dk_vals[0] = k_vals[n_step] - k_vals[n_step - 1];
        dk_vals[1] = k_vals[n_step] - k_vals[n_step - 2];
        dk_vals[2] = k_vals[n_step] - k_vals[n_step - 3];
        for (unsigned int dki = 0; dki < 3; dki++){
            for (unsigned int fi = 0; fi < 3; fi++){
                A(fi, dki) = pow(dk_vals[dki], fi + 1);
            }
        }
        A_piv.compute(A);
        c = A_piv.solve(b1);
        c0 = - (c(0) + c(1) + c(2));
        dpdk = - (c0 * phase_shifts[n_step] + c(0) * phase_shifts[n_step - 1] + c(1) * phase_shifts[n_step - 2] + c(2) * phase_shifts[n_step - 3]);
        
        dr0_k = - (c0 * r_levinson[n_step] + c(0) * r_levinson[n_step - 1] + c(1) * r_levinson[n_step - 2] + c(2) * r_levinson[n_step - 3]);
        c = A_piv.solve(b2);
        c0 = - (c(0) + c(1) + c(2));
        dr0_k2 = 2 * (c0 * r_levinson[n_step] + c(0) * r_levinson[n_step - 1] + c(1) * r_levinson[n_step - 2] + c(2) * r_levinson[n_step - 3]);
        d2pdk2 = 2 * (c0 * phase_shifts[n_step] + c(0) * phase_shifts[n_step - 1] + c(1) * phase_shifts[n_step - 2] + c(2) * phase_shifts[n_step - 3]);

        c = A_piv.solve(b3);
        c0 = - (c(0) + c(1) + c(2));
        d3pdk3 = - 6 * (c0 * phase_shifts[n_step] + c(0) * phase_shifts[n_step - 1] + c(1) * phase_shifts[n_step - 2] + c(2) * phase_shifts[n_step - 3]);
        dr0_k3 = - 6 * (c0 * r_levinson[n_step] + c(0) * r_levinson[n_step - 1] + c(1) * r_levinson[n_step - 2] + c(2) * r_levinson[n_step - 3]);

        double tol_dk = pow(6 * extrapol_tol / abs(d3pdk3), 1. / 3.);
        dk = std::max(min_dk, std::min(tol_dk, max_dk));
        // std::cout << "Using dk = " << dk << ", " << k << ", " << phase_shifts.back() << std::endl;
        k += dk;
        if (dk != max_dk && !resonance_found){
            resonance_found = true;
            resonance_passed = false;
            if (verbose) std::cout << "Resonance found at (l, k) = " << l << ", " << k - dk << ", " << phase_shifts.back() / PI << std::endl;
        }
        if (dk == max_dk && resonance_found && !resonance_passed){
            resonance_passed = true;
            if (verbose) std::cout << "Resonance passed at (l, k) = " << l << ", " << k - dk << ", " << phase_shifts.back() / PI << std::endl;
        }
        // if (resonance_passed){
        //     min_dk = 0.1; max_dk = 0.5; extrapol_tol = 1e-3;
        // }
        // if (resonance_found && !resonance_passed){
        //     std::cout << "In resonance : " << k << ", " << dk << ", " << phase_shifts.back() / PI << std::endl;
        // }
        // if (r_levinson.back() == n && !prepare_jump){
        //     prepare_jump = true;
        //     if (verbose) std::cout << "Prepare jump at " << l << ", " << k << std::endl;
        // }

        double prev_r0 = r_levinson.back();
        int prev_n = n;
        double prev_delta = phase_shifts.back();
        fill_absolute_phase_shifts(i, j, l, k, n, k_vals, phase_shifts, r_levinson);
        dk_vals[0] = k_vals[n_step] - k_vals[n_step - 1];
        dk_vals[1] = k_vals[n_step] - k_vals[n_step - 2];
        dk_vals[2] = k_vals[n_step] - k_vals[n_step - 3];
        for (unsigned int dki = 0; dki < 3; dki++){
            for (unsigned int fi = 0; fi < 3; fi++){
                A(fi, dki) = pow(dk_vals[dki], fi + 1);
            }
        }
        A_piv.compute(A);
        c = A_piv.solve(b1);
        c0 = - (c(0) + c(1) + c(2));
        dr0_k = - (c0 * r_levinson[n_step] + c(0) * r_levinson[n_step - 1] + c(1) * r_levinson[n_step - 2] + c(2) * r_levinson[n_step - 3]);

        c = A_piv.solve(b2);
        c0 = - (c(0) + c(1) + c(2));
        dr0_k2 = 2 * (c0 * r_levinson[n_step] + c(0) * r_levinson[n_step - 1] + c(1) * r_levinson[n_step - 2] + c(2) * r_levinson[n_step - 3]);

        c = A_piv.solve(b3);
        c0 = - (c(0) + c(1) + c(2));
        dr0_k3 = - 6 * (c0 * r_levinson[n_step] + c(0) * r_levinson[n_step - 1] + c(1) * r_levinson[n_step - 2] + c(2) * r_levinson[n_step - 3]);
        double prev_k = *(k_vals.end() - 2);
        int n_step = r_levinson.size() - 1;
        double extrapol1 = prev_r0 + dk * dr0_k;
        double extrapol2 = extrapol1 + pow(dk, 2) * dr0_k2 / 2;
        double extrapol3 = extrapol2 + pow(dk, 3) * dr0_k3 / 6;
        double extrapol_r0 = extrapol2;

        if (verbose && false) std::cout << "Do jump " << k << ", " << dk << ", " << prev_k << " ?"
            << "\n\t\t r : " << r_levinson[n_step - 2] / sigma[i][j] << ", " << r_levinson[n_step - 1] / sigma[i][j] << ", " << r_levinson[n_step] / sigma[i][j]
            << "\n\t\t r : " << prev_r0 / sigma[i][j] << ", " << r_levinson.back() / sigma[i][j] 
            << "\n\t\t ex : " << extrapol1 / sigma[i][j] << ", " << extrapol2 / sigma[i][j] << ", " << extrapol3 / sigma[i][j] 
            << "\n\t\t e : " << (extrapol1 - r_levinson.back()) / prev_r0 << ", " << (extrapol2 - r_levinson.back()) / prev_r0 << ", " << abs(extrapol3 - r_levinson.back()) / prev_r0
            << "\n\t\t e : " << abs(r_levinson.back() - extrapol_r0) / prev_r0 << " / " << abs(prev_r0 - extrapol_r0) / prev_r0
            << "\n\t\t dr : " << dr0_k / sigma[i][j] << ", " << dr0_k2 / sigma[i][j] << ", " << dr0_k3 / sigma[i][j] 
            << "\n\t\t n : " << prev_n << ", " << n
            << "\n\t\t k : " << *(k_vals.end() - 2) << ", " << k_vals.back() << std::endl;

        if ( abs(r_levinson.back() - extrapol_r0) > abs(prev_r0 - extrapol_r0) 
             && abs(phase_shifts.back() - prev_delta) < PI 
             && r_levinson.back() < 15 * sigma[i][j]
             && r_levinson.back() < prev_r0
             && !worst_is_over
             && (*(k_vals.end() - 2)) - (*(k_vals.end() - 3)) > min_dk
             && l != 0
             && 2 * red_mass[i][j] > 30e-3 / AVOGADRO){
            done_jump = true; n_jump++;
            n += 1;
            if (verbose) std::cout << "Jumping resonance at (" << l << ", " << k << ") : " << prev_r0 / sigma[i][j] << ", " << r_levinson.back() / sigma[i][j] 
                                    << ", " << extrapol_r0 / sigma[i][j] << ", " << dr0_k / sigma[i][j] 
                                    << " / " << prev_n << " => " << n << std::endl;
            for (size_t jmp_idx = 0; jmp_idx < 5; jmp_idx++){
                k += min_dk / 10;
                r_lev = n;
                delta = quantum_phase_shift(i, j, l, E_from_k(k / sigma[i][j]), r_lev);
                k_vals.push_back(k);
                phase_shifts.push_back(delta + n * PI);
                r_levinson.push_back(r_lev);
            }
            if (verbose) std::cout << "Finished resonance jump : " << l << ", " << k << std::endl;
            // extrapol_tol = 1e-5;
            // max_dk = 1e-2;
            // min_dk = 1e-6;
        }

        if (d2pdk2 < 0 && done_jump && !worst_is_over){
            worst_is_over = true;
            extrapol_tol = 1e-2;
            max_dk = 1.0;
            min_dk = 1e-3;
        }

        // std::cout << "Check tail " << k << " : " << d2pdk2 << ", " << dpdk << ", " << d2pdk2 / dpdk << " / " << phase_shifts.back() / PI << ", " << n << ", " << nl0 << std::endl;
        if (( d2pdk2 / dpdk < 0.1 && d2pdk2 < 0 && dpdk < 0) || n < prev_n) {
            if (l == 0 && phase_shifts.back() > 0) continue;
            if (verbose) std::cout << "Entering tail : " << k << ", " << dpdk << ", " << d2pdk2 / dpdk << ", " << phase_shifts.back() / PI << std::endl;
            is_linear = true;
        }
    }
    // std::cout << "Entering tail... " << k << ", " << phase_shifts.back() << std::endl;
    max_dk = 5;
    min_dk = 0.5;
    extrapol_tol = 1e-2;
    std::array<double, 3> prev_dk_vals = dk_vals;
    while (k + max_dk < k_max){
        // std::array<double, 3> extrapol_coeff = quadric_extrapolate_coeff(k_vals, phase_shifts);
        n_step = k_vals.size() - 1;
        dk_vals[0] = k_vals[n_step] - k_vals[n_step - 1];
        dk_vals[1] = k_vals[n_step] - k_vals[n_step - 2];
        dk_vals[2] = k_vals[n_step] - k_vals[n_step - 3];
        if ( (dk_vals[0] != prev_dk_vals[0]) || (dk_vals[1] != prev_dk_vals[1]) || (dk_vals[2] != prev_dk_vals[2]) ) {
            prev_dk_vals = dk_vals;
            for (unsigned int dki = 0; dki < 3; dki++){
                for (unsigned int fi = 0; fi < 3; fi++){
                    A(fi, dki) = pow(dk_vals[dki], fi + 1);
                }
            }
            c = A_piv.compute(A).solve(b3);
            c0 = - (c(0) + c(1) + c(2));
        }
        
        d3pdk3 = - 6 * (c0 * phase_shifts[n_step] + c(0) * phase_shifts[n_step - 1] + c(1) * phase_shifts[n_step - 2] + c(2) * phase_shifts[n_step - 3]);
        double tol_dk = pow(6 * extrapol_tol / abs(d3pdk3), 1. / 3.);

        dk = std::max(min_dk, std::min(tol_dk, max_dk));
        // std::cout << "In Tail : " << k << ", " << n << ", " << dk << " ( " << dk_vals[0] << ", " << dk_vals[1] << ", " << dk_vals[2] << " ) " << d3pdk3 << std::endl;
        k += dk;
        // double extrapol_next_delta = extrapol_coeff[0] * pow(dk, 2) + extrapol_coeff[1] * dk + extrapol_coeff[2];
        fill_absolute_phase_shifts_tail(i, j, l, k, n, k_vals, phase_shifts);
    }

    if (abs(k - k_max) > 1e-8) fill_absolute_phase_shifts(i, j, l, k_max, n, k_vals, phase_shifts, r_levinson);

    std::cout << "\tFinished : " << l << ", " << k_max << std::endl;
}

void Quantum::trace_total_phase_shifts(int i, int j, double T){
    double tol = 1e-8;
    double beta = 1 / (BOLTZMANN * T);
    double pre_exp = - beta * pow(HBAR, 2) / (2 * red_mass[i][j] * pow(sigma[i][j], 2));
    double pre_k = pow(HBAR, 2) / (red_mass[i][j] * sigma[i][j]);
    double k_max = sqrt(log(tol) / pre_exp);

    if (stored_total_phase_shifts.size() > 0){
        if (stored_total_phase_shifts[0].back() >= k_max) return;
    }
    double dk = 0.01;
    size_t n_vals = static_cast<size_t>(k_max / dk);
    vector1d k_vals(n_vals);
    std::vector<double> wts(n_vals);
    
    
    k_vals[0] = 0;
    for (size_t k_idx = 1; k_idx < k_vals.size(); k_idx++){
        k_vals[k_idx] = k_vals[k_idx - 1] + dk;
        wts[k_idx] = exp(pre_exp * pow(k_vals[k_idx], 2)) * pre_k * k_vals[k_idx];
    }

    double tot = 0;
    double part_sum = 0;
    double L = de_broglie_wavelength(i, j, T);
    double B0 = 0;
    const vector2d& Eb = E_bound[i][j];
    for (size_t vib_i = 0; vib_i < Eb.size(); vib_i++){
        for (size_t j = 0; j < Eb[vib_i].size(); j += 2){
            B0 += (2 * j + 1) * pow(L, 3) * AVOGADRO;
        }
    }
    vector1d integrand(n_vals);
    vector1d phase_shifts(n_vals, 0);
    vector1d phase_shift_l(n_vals);
    int l = 0;
    do {
        trace_absolute_phase_shifts(i, j, l, k_max);
        phase_shift_l = interpolate_grid(k_vals, absolute_phase_shift_map[l][0], absolute_phase_shift_map[l][1]);
        for (size_t idx = 0; idx < phase_shifts.size(); idx++){
            phase_shifts[idx] += (2 * l + 1) * phase_shift_l[idx];
            integrand[idx] = (2 * l + 1) * phase_shift_l[idx] * wts[idx];
        }
        part_sum = trapezoid(dk, integrand);
        tot += part_sum;
        // if (l > 100) break;
        double B_th = B0 - tot * pow(L, 3) * AVOGADRO  / (PI * BOLTZMANN * T * sigma[i][j]);
        double B_part = - part_sum * pow(L, 3) * AVOGADRO  / (PI * BOLTZMANN * T * sigma[i][j]);
        double B_th_tot = - tot * pow(L, 3) * AVOGADRO  / (PI * BOLTZMANN * T * sigma[i][j]);
        std::cout << "Converging total : " << B_th_tot << ", " << B_part << ", " 
                    << B_part / B_th_tot << ", " << B0 << ", " << B_th << std::endl;
        l += 2;
        
    } while (abs(part_sum / tot) > 1e-4); // abs(phase_shift_l.back() * (2 * l + 1) / phase_shifts.back()) > 1e-5);
    stored_total_phase_shifts = vector2d({k_vals, phase_shifts});
}

vector2d Quantum::total_phase_shifts(int i, int j, double k_max){
    trace_total_phase_shifts(i, j, k_max);
    return stored_total_phase_shifts;
}

double Quantum::integral_phase_shift(int i, int j, int l, double T){
    double tol = 1e-8;
    double beta = 1 / (BOLTZMANN * T);
    double pre_exp = - beta * pow(HBAR, 2) / (2 * red_mass[i][j] * pow(sigma[i][j], 2));
    double pre_k = pow(HBAR, 2) / (red_mass[i][j] * sigma[i][j]);
    double k_max = sqrt(log(tol) / pre_exp);

    trace_absolute_phase_shifts(i, j, l, k_max);
    vector1d& k_vals = absolute_phase_shift_map[l][0];
    vector1d& phase_shifts = absolute_phase_shift_map[l][1];
    
    const auto integrand = [&](size_t idx){return phase_shifts[idx] * exp(pre_exp * pow(k_vals[idx], 2)) * pre_k * k_vals[idx];};
    const auto integrand_term = [&](size_t idx){return integrand(idx) + integrand(idx - 1);};
    double I = 0;
    double I_term;
    size_t I_idx = 1;
    for (; I_idx < phase_shifts.size(); I_idx++){
        I_term = 0.5 * (k_vals[I_idx] - k_vals[I_idx - 1]) * integrand_term(I_idx);
        I += I_term;
    }
    std::cout << "Initial integral to : " << k_vals.back() << " = " << I << " (" << I_term << ", " << abs(I_term / I) << ")" << std::endl;
    while (abs(I_term / I) > 1e-6) {
        k_max *= 1.1;
        trace_absolute_phase_shifts(i, j, l, k_max);
        for (; I_idx < phase_shifts.size(); I_idx++){
            I_term = 0.5 * (k_vals[I_idx] - k_vals[I_idx - 1]) * integrand_term(I_idx);
            I += I_term;
        }
        std::cout << "Integrate to " << k_vals.back() << ", Result : " << I << " (" << I_term << ", " << abs(I_term / I) << ")" << std::endl;
    } 
    return I / pow(sigma[i][j], 3);
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

double Quantum::cross_section_kernel(int i, int j, int n, int l, double E){
    // Input is reduced energy (E = E(Joule) / eps[i][j])
    switch (n){
    case 0: return (2 * l + 1) * pow(sin(phase_shift(i, j, l, E)), 2);
    case 1: return (l + 1) * pow(sin(phase_shift(i, j, l, E) - phase_shift(i, j, l + 1, E)), 2);
    case 2: return (l + 1) * (l + 2) * pow(sin(phase_shift(i, j, l, E) - phase_shift(i, j, l + 2, E)), 2) / (2 * l + 3);
    default:
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
}

double Quantum::cross_section(int i, int j, const int n, const double E){
    // Input is reduced energy (E = E(Joule) / eps[i][j])
    if (E < 5e-4) return cross_section(i, j, n, 5e-4); // Cross sections approach constant value at small E, but numerical solver fails for very small E
    if (is_singlecomp) i = j;
    
    {
        std::shared_lock lock(cross_section_map_mutex);
        const auto pos = cross_section_map.find({i, j, n, E});
        if (pos != cross_section_map.end()) return pos->second;
    }

    auto [sym_prefactor, anti_prefactor] = get_symmetry_weights(i, j);
    double Q = 0;
    if (sym_prefactor > 0.){
        double Q_even = 0;
        int l_even = 0;
        double q_even;
        do {
            q_even = cross_section_kernel(i, j, n, l_even, E);
            Q_even += q_even;
            l_even += 2;
            // std::cout << "Converging cross section (" << E << ", " << l_even << ") : " << Q_even << ", " << q_even << std::endl;
        } while (abs(q_even) > 1e-6 * Q_even);
        Q += sym_prefactor * Q_even;
    } 
    if (anti_prefactor > 0.){
        double Q_odd = 0;
        int l_odd = 1;
        double q_odd;
        do {
            q_odd = cross_section_kernel(i, j, n, l_odd, E);
            Q_odd += q_odd;
            l_odd += 2;
        } while (abs(q_odd) > 1e-6 * Q_odd);
        Q += anti_prefactor * Q_odd;
    } 

    double k2 = 4. * red_mass[i][j] * E * eps[i][j] / pow(HBAR, 2);
    Q *= 4. * PI / k2;

    {
        std::unique_lock lock(cross_section_map_mutex);
        const CrossSectionPoint point = {i, j, n, E};
        const auto pos = cross_section_map.find(point);
        if (pos != cross_section_map.end()) return pos->second;
        cross_section_map[point] = Q;
    }

    return Q;
}

double Quantum::classical_cross_section(int i, int j, int l, double E){
    return Spherical::cross_section(i, j, l, E);
}

double Quantum::quantum_omega(int i, int j, int n, int s, double T){
    double beta = 1. / (T * BOLTZMANN);
    std::cout << "Integrating : " << n << ", " << s << std::endl;
    const auto kernel = [&](double Eb){
        double Q = cross_section(i, j, n, Eb / (beta * eps[i][j]));
        double pref = exp(- Eb) * pow(Eb, s + 1);
        // std::cout << "Evaluate integrand : " << Eb / (beta * eps[i][j]) << ", " << Q << ", " << pref << std::endl;
        return pref * Q;
    };
    double maxpoint = s + 1.;
    double I = simpson_inf(kernel, 0, maxpoint / 3);
    double val = I * sqrt(2. / (PI * beta * red_mass[i][j]));
    return val;
}

double Quantum::classical_omega(int i, int j, int l, int r, double T){
    double val = Spherical::omega(i, j, l, r, T);
    return val;
}

double Quantum::omega(int i, int j, int l, int r, double T){
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

double Quantum::partial_scattering_volume(int i, int j, double E, int l_max){
    int l, l_step{2};
    if (is_singlecomp) i = j;
    if (i != j) {
        l = 0; l_step = 1;
    }
    else {
        switch (half_spin[i]){
        case 0: l = 0; break;
        case 2: l = 1; break;
        default: throw std::runtime_error("Invalid half-spin!");
        }
    }

    double S = 0;
    double S_term = 0;
    for (; l <= l_max; l += l_step) {
        double delta = absolute_phase_shift(i, j, l, E);
        S_term = (2 * l + 1) * delta;        
        S += S_term;
        // std::cout << "Converging (" << l << ", " << E << ") : " << S << " / " << S_term << std::endl;
        // if (l > 100) break;
    }
    return S;
}

double Quantum::scattering_volume(int i, int j, double E){
    int l, l_step{2};
    if (is_singlecomp) i = j;
    if (i != j) {
        l = 0; l_step = 1;
    }
    else {
        switch (half_spin[i]){
        case 0: l = 0; break;
        case 2: l = 1; break;
        default: throw std::runtime_error("Invalid half-spin!");
        }
    }

    double S = 0;
    double S_term = 0;
    do {
        double delta = absolute_phase_shift(i, j, l, E);
        S_term = (2 * l + 1) * delta;        
        S += S_term;
        // std::cout << "Converging (" << l << ", " << E << ") : " << S << " / " << S_term << std::endl;
        l += l_step;
        // if (l > 100) break;
    } while (abs(S_term / S) > 1e-3);
    return S;
}

vector1d Quantum::second_virial(int i, int j, const vector1d& T){
    const size_t Ncores = 8;
    size_t T_idx = 0;
    vector1d B_vals(T.size(), 0.);
    while (T_idx < T.size()){
        std::vector<std::thread> threads;
        for (size_t thread_idx = 0; thread_idx < Ncores; thread_idx++){
            threads.emplace_back([&](double T_val, double& B){B = second_virial(i, j, T_val);}, T[T_idx], std::ref(B_vals[T_idx]));
            if (++T_idx >= T.size()) break;
        }
        for (auto& thread : threads){
            thread.join();
        }
    }
    return B_vals;
}

double Quantum::second_virial(int i, int j, double T){
    set_internals(0, T, {1.});
    if (!quantum_is_active) return Spherical::second_virial(i, j, T);

    double B = 0;
    std::map<char, double> contribs = second_virial_contribs(i, j, T, "ibt");
    std::cout << "Contribs : " << contribs['i'] << ", " << contribs['b'] << ", " << contribs['t'] << std::endl;
    for (const auto& [key, c] : contribs){
        std::cout << "Contrib (" << key << ") : " << c << std::endl;
        B += c;
    }
    return B;
}

double Quantum::bound_second_virial(int i, int j, double T){
    set_internals(0, T, {1.});
    if (!quantum_is_active) return Spherical::bound_second_virial(i, j, T);
    return second_virial_contribs(i, j, T, "b")['b'];
}

double Quantum::dimer_constant(int i, int j, double T){
    set_internals(0, T, {1.});
    if (!quantum_is_active) return Spherical::dimer_constant(i, j, T);
    
    std::array<double, 2> wts = get_symmetry_weights(i, j);
    const auto [sym_wt, anti_wt] = wts;

    size_t l_start, l_step;
    if (sym_wt == 0){
        l_start = 1; l_step = 2;
    }
    else if (anti_wt == 0){
        l_start = 0; l_step = 2;
    }
    else {
        l_start = 0; l_step = 1;
    }

    double K = 0;
    double beta = 1 / (BOLTZMANN * T);
    const vector2d& Eb = E_bound[i][j];
    for (size_t vib_i = 0; vib_i < Eb.size(); vib_i++){
        for (size_t l = l_start; l < Eb[vib_i].size(); l+=l_step){
            K += wts[l % 2] * (2 * l + 1) * exp(- beta * Eb[vib_i][l]);
        }
    }
    double lamb = de_broglie_wavelength(i, j, T);
    double V_th = pow(lamb, 3);
    return K * V_th * AVOGADRO;
}

std::map<char, double> Quantum::second_virial_contribs(int i, int j, double T, const std::string& contribs){
    // See: Kilpatrick et al. Second virial coefficients of He3 and He4, Physical Review vol. 94, nr. 5 (1954)
    //      DOI: https://doi.org/10.1103/PhysRev.94.1103 
    if (!quantum_is_active) return Spherical::second_virial_contribs(i, j, T, contribs);
    set_internals(0, T, {1.});
    
    std::array<double, 2> wts = get_symmetry_weights(i, j);
    const auto [sym_wt, anti_wt] = wts;

    size_t l_start, l_step;
    if (sym_wt == 0){
        l_start = 1; l_step = 2;
    }
    else if (anti_wt == 0){
        l_start = 0; l_step = 2;
    }
    else {
        l_start = 0; l_step = 1;
    }

    std::map<char, double> B_contribs;
    const auto contribs_contain = [&](const std::string& c){return str_contains(contribs, c);};

    if (contribs_contain("i")){
        B_contribs['i'] = (sym_wt - anti_wt) / 16.;;
    }

    double beta = 1 / (BOLTZMANN * T);
    if (contribs_contain("b")){
        const vector2d& Eb = E_bound[i][j];
        double B_bound = 0;
        for (size_t vib_i = 0; vib_i < Eb.size(); vib_i++){
            for (size_t l = l_start; l < Eb[vib_i].size(); l += l_step){
                B_bound += wts[l % 2] * (2 * l + 1) * exp(- beta * Eb[vib_i][l]);
            }
        }
        B_contribs['b'] = B_bound;
    }

    if (contribs_contain("t")){
        double B_th = 0;
        size_t l = l_start;
        double B_th_term;
        std::cout << "Starting at T = " << T << std::endl;

        double tol = 1e-8;
        double beta = 1 / (BOLTZMANN * T);
        double pre_exp = - beta * pow(HBAR, 2) / (2 * red_mass[i][j] * pow(sigma[i][j], 2));
        double pre_k = pow(HBAR, 2) / (red_mass[i][j] * sigma[i][j]);
        double k_max = sqrt(log(tol) / pre_exp);
        trace_total_phase_shifts(i, j, T);

        const vector1d& k_vals = stored_total_phase_shifts[0];
        const vector1d& phase_shifts = stored_total_phase_shifts[1];

        const auto integrand = [&](size_t idx){return phase_shifts[idx] * exp(pre_exp * pow(k_vals[idx], 2)) * pre_k * k_vals[idx];};
        const auto integrand_term = [&](size_t idx){return integrand(idx) + integrand(idx - 1);};

        double I = 0;
        for (size_t I_idx = 1; I_idx < phase_shifts.size(); I_idx++){
            double I_term = 0.5 * (k_vals[I_idx] - k_vals[I_idx - 1]) * integrand_term(I_idx);
            I += I_term;
        }
        B_th = I / (PI * BOLTZMANN * T * sigma[i][j]);
        std::cout << "Bth : " << B_th << ", ";
        const vector2d& Eb = E_bound[i][j];
        for (size_t vib_i = 0; vib_i < Eb.size(); vib_i++){
            for (size_t l = l_start; l < Eb[vib_i].size(); l += l_step){
                B_th -= wts[l % 2] * (2 * l + 1);
            }
        }
        std::cout << B_th << std::endl;
        B_contribs['t'] = B_th;
    }

    if (contribs_contain("l")){
        B_contribs['l'] = 0;
        const vector2d& Eb = E_bound[i][j];
        for (size_t vib_i = 0; vib_i < Eb.size(); vib_i++){
            for (size_t l = l_start; l < Eb[vib_i].size(); l += l_step){
                B_contribs['l'] -= wts[l % 2] * (2 * l + 1);
                std::cout << "Bl : " << l << ", " << B_contribs['l'] << std::endl;
            }
        }
    }

    double lamb = de_broglie_wavelength(i, j, T);
    double V_th = pow(lamb, 3);
    for (auto& [_, c] : B_contribs){
        c *= - AVOGADRO * V_th;
    }

    return B_contribs;
}

double Quantum::semiclassical_second_virial(int i, int j, double T){
    double tmp_JKWB_E_limit = JKWB_E_limit;
    int tmp_JKWB_l_limit = JKWB_l_limit;
    bool tmp_quantum_active = quantum_is_active;
    set_JKWB_limits(0, 0);
    set_quantum_active(true);
    double B = second_virial(i, j, T);
    set_quantum_active(tmp_quantum_active);
    set_JKWB_limits(tmp_JKWB_E_limit, tmp_JKWB_l_limit);
    return B;
}

void Quantum::set_quantum_active(bool active){
    if (active && !quantum_supported){
        std::cout << "WARNING : Some components do not support Quantum! (I'm letting you activate anyway at your own discretion ...)\n";
    }
    if (active != quantum_is_active){
        clear_all_caches();
    }
    quantum_is_active = active;
}