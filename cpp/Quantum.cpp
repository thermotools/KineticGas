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

Quantum::Quantum(std::string comps_) 
    : Spherical(comps_, true), 
    E_bound(Ncomps, std::vector<vector2d>(Ncomps)),
    half_spin(Ncomps, 0),
    spin(Ncomps, 0),
    rot_ground_state(Ncomps, 0)
{   
    for (size_t i = 0; i < Ncomps; i++){
        const auto qdata = compdata[i]["Quantum"];
        bool has_quantum_params = static_cast<bool>(qdata["Quantum"]);
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

double Quantum::de_broglie_wavelength(int i, double T){
    return PLANCK / sqrt(PI * m[i] * BOLTZMANN * T);
}

double Quantum::de_broglie_wavelength(int i, int j, double T){
    return PLANCK / sqrt(2 * PI * red_mass[i][j] * BOLTZMANN * T);
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
        // std::cout << "shift (" << E / eps[i][j] << ") : " << prev_delta / PI << ", " << new_delta / PI << std::endl;
    } while (abs(new_delta - prev_delta) > 1e-5);
    return new_delta;
}

double Quantum::phase_shift(int i, int j, int l, double E){
    if (E < 5e-4) return 0;
    
    const std::pair<int, double> point(l, E);
    const auto pos = phase_shift_map.find(point);
    if (pos != phase_shift_map.end()) return pos->second;

    double delta = ((E > JKWB_E_limit) || (l > JKWB_l_limit)) ? JKWB_phase_shift(i, j, l, E) : quantum_phase_shift(i, j, l, E);
    phase_shift_map[point] = delta;
    return delta;
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
    double delta = (l <= rot_ground_state[i]) ? PI : 0;
    double Ei = 1e-1;
    if (l <= 2) {
        Ei = 1e-2;
    }
    else if (l > 10) {
        Ei = 1.;
    }
    double dlnE = .25;
    double dE = 1.15; // exp(dlnE);
    while (Ei < E){
        delta = absolute_phase_shift(i, j, l, Ei, delta);
        Ei *= dE;
    }
    delta = absolute_phase_shift(i, j, l, E, delta);
    return delta;
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
    
    auto [sym_prefactor, anti_prefactor] = get_symmetry_weights(i, j);
    // if ((i != j) && (!is_singlecomp)) {
    //     sym_prefactor = anti_prefactor = .5;
    // }
    // else {
    //     double S = half_spin[i] / 2.;
    //     switch (half_spin[i]){
    //     case 0:
    //         sym_prefactor = 1;
    //         anti_prefactor = 0;
    //         break;
    //     case 1:
    //         sym_prefactor = 3. / 4.;
    //         anti_prefactor = 1. / 4.;
    //         break;
    //     case 2:
    //         sym_prefactor = 5. / 9.; // (S + 1) / (2 * S + 1);
    //         anti_prefactor = 4. / 9.; // S / (2 * S + 1);
    //         break;
    //     default:
    //         throw std::runtime_error("Invalid half-spin!");
    //     }
    //     if (half_spin[i] % 2 != 0) { // Odd Half-Integer spin, Fermions: Swap prefactors
    //         double tmp = sym_prefactor;
    //         sym_prefactor = anti_prefactor;
    //         anti_prefactor = tmp;
    //     }
    // }
    double Q = 0;
    if (sym_prefactor > 0.){
        double Q_even = 0;
        int l_even = 0;
        double q_even;
        do {
            q_even = cross_section_kernel(i, j, n, l_even, E);
            Q_even += q_even;
            l_even += 2;
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
    return Q;
}

double Quantum::classical_cross_section(int i, int j, int l, double E){
    return Spherical::cross_section(i, j, l, E);
}

double Quantum::quantum_omega(int i, int j, int n, int s, double T){
    double beta = 1. / (T * BOLTZMANN);
    const auto kernel = [&](double Eb){
        return exp(- Eb) * pow(Eb, s + 1) * cross_section(i, j, n, Eb / (beta * eps[i][j]));
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

double Quantum::second_virial(int i, int j, double T){
    return (quantum_is_active) ? quantum_second_virial(i, j, T) : Spherical::second_virial(i, j, T);
}

double Quantum::quantum_second_virial(int i, int j, double T){
    auto [sym_wt, anti_wt] = get_symmetry_weights(i, j);
    double B = 0;
    double B_s = 0; double B_a = 0;
    if (sym_wt != 0){
        std::vector<double> contribs = second_virial_contribs(i, j, T, StatisticType::BoseEinstein);
        for (const double c : contribs){
            B += sym_wt * c; B_s += sym_wt * c;
            std::cout << sym_wt * c << ", ";
        }
    }
    std::cout << std::endl;
    if (anti_wt != 0){
        std::vector<double> contribs = second_virial_contribs(i, j, T, StatisticType::FermiDirac);
        for (const double c : contribs){
            B += anti_wt * c; B_a += anti_wt * c;
            std::cout << anti_wt * c << ", ";
        }
    }
    std::cout << "\nSym / anti / tot : " << B_s << " / " << B_a << " / " << B << std::endl;
    return B;
}

vector1d Quantum::second_virial_contribs(int i, int j, double T, int istats){
    // See: Kilpatrick et al. Second virial coefficients of He3 and He4, Physical Review vol. 94, nr. 5 (1954)
    //      DOI: https://doi.org/10.1103/PhysRev.94.1103 
    double beta = 1 / (BOLTZMANN * T);
    double lamb = de_broglie_wavelength(i, j, T);
    double V_th = pow(lamb, 3);
    double B_id;
    double B_bound = 0;
    double B_th = 0;
    const vector2d& Eb = E_bound[i][j];
    size_t l_start;
    size_t l_step = 2;
    if (istats == StatisticType::BoseEinstein){
        l_start = 0;
        B_id = 1. / 16.;
    }
    else if (istats == StatisticType::FermiDirac){
        l_start = 1;
        B_id = - 1. / 16.;
    }
    else {
        throw std::runtime_error("Invalid statistic type!");
    }

    for (size_t vib_i = 0; vib_i < Eb.size(); vib_i++){
        for (size_t l = l_start; l < Eb[vib_i].size(); l+=l_step){
            B_bound += (2 * l + 1) * (exp(- beta * Eb[vib_i][l]) - 1.);
        }
    }

    size_t l = l_start;
    double B_th_term;
    do {
        auto integrand = [&](double beta_E){
            double delta = absolute_phase_shift(i, j, l, beta_E / (beta * eps[i][j]));
            return delta * exp(- beta_E);
        };
        double I = simpson(integrand, 0, 0.2, 10) 
                 + simpson(integrand, 0.2, 1, 10) 
                 + simpson(integrand, 1, 3, 20) 
                 + simpson_inf(integrand, 3, 7);

        B_th_term = (2 * l + 1) * I / PI;        
        B_th += B_th_term;
        l += l_step;
    } while (abs(B_th_term / B_th) > 1e-3);

    double B = B_id + B_bound + B_th;
    std::cout << "Contribs : " << B_id << ", " << B_bound << ", " << B_th << std::endl;
    vector1d contribs = {B_id, B_bound, B_th};
    for (double& c : contribs){
        c *= - AVOGADRO * V_th;
    }
    return  contribs;
}

double Quantum::semiclassical_second_virial(int i, int j, double T){
    double tmp_JKWB_E_limit = JKWB_E_limit;
    int tmp_JKWB_l_limit = JKWB_l_limit;
    set_JKWB_limits(0, 0);
    double B = quantum_second_virial(i, j, T);
    set_JKWB_limits(tmp_JKWB_E_limit, tmp_JKWB_l_limit);
    return B;
}

void Quantum::set_quantum_active(bool active){
    if (active != quantum_is_active){
        clear_all_caches();
    }
    quantum_is_active = active;
}