#include "KineticGas.h"
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

class Quantum : public KineticGas {
public:
    Quantum(std::string comps) : KineticGas(comps, false) {
        half_spin = std::vector<size_t>(Ncomps, 0);
        for (size_t i = 0; i < Ncomps; i++){
            half_spin[i] = static_cast<size_t>(static_cast<double>(compdata[i]["spin"]) * 2 + 0.5);
        }
    }

    double potential(int i, int j, double r){
        return 4 * (pow(1. / r, 12) - pow(1. / r, 6));
    }
    double potential_r(int i, int j, double r){
        return 0;
    }
    double potential_rr(int i, int j, double r){
        return 0;
    }

    double phase_shift(int i, int j, int l, double E){
        const double step_size{1e-3};
        double k2 = 2. * (m[i] * m[j] / (m[i] + m[j])) / (pow(HBAR, 2) * E);
        double k = sqrt(k2);
        std::array<double, 3> r = {- step_size, 0., step_size}; // rolls one step at the start of each iteration, so first iteration is {0, step_size, 2 * step_size}
        std::array<double, 3> psi = {0., 0., 0.};
        const auto g_fun = [&](const double r_i){return - (l * (l + 1)) / pow(r_i, 2) + k2 * (E - potential(i, j, r_i));};
        std::array<double, 3> g_vals = {0., 0., g_fun(r[2])}; // only g_vals[1] and g_vals[2] are used. Using double[3] to match indexes with r and psi.
        const auto next_psi = [&](const std::array<double, 3>& r_n, const std::array<double, 3>& psi_n, const std::array<double, 3>& g_n) -> double {
            double s2 = pow(step_size, 2) / 12;
            return 2. * psi_n[1] * ((1. - 5. * s2 * g_n[1]) / (1 + s2 * g_n[2])) - psi_n[0];
        };
        double prev_delta;
        double new_delta = 2. * PI;
        do {
            prev_delta = new_delta;
            r[0] = r[1]; r[1] = r[2]; r[2] = r[1] + step_size;
            g_vals[1] = g_vals[2]; g_vals[2] = g_fun(r[2]); // g_vals[0] = g_vals[1]; Not needed, because g_vals[0] is never used.
            psi[0] = psi[1]; psi[1] = psi[2]; psi[2] = next_psi(r, psi, g_vals);

            double dpsi = (psi[2] - psi[0]) / (2 * step_size);
            double gamma = (dpsi / psi[1]) - (1. / r[1]);
            dual r_dual = r[1];
            auto [jl, djl] = autodiff::derivatives(&spherical_bessel, autodiff::wrt(r_dual), autodiff::at(k * r_dual, l, 1));
            auto [yl, dyl] = autodiff::derivatives(&spherical_bessel, autodiff::wrt(r_dual), autodiff::at(k * r_dual, l, 2));

            new_delta = atan((k * djl - gamma * jl) / (k * dyl - gamma * yl));
        } while (abs(new_delta - prev_delta) > 1e-10);
        return new_delta;
    }

    inline double cross_section_A(int n, int l, size_t k){
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

    double cross_section_kernel(int i, int j, double n, double l, double E){
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

    double cross_section(int i, int j, int n, double E){
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

    double omega(int i, int j, int n, int s, double T){
        double beta = 1. / (T * BOLTZMANN);
        double sfac = 1;
        for (int si = 2; si <= s + 1; si++){
            sfac *= si;
        }
        const auto kernel = [&](double E){return cross_section(i, j, n, E) * exp(- beta * E) * pow(beta * E, s + 2) / sfac;};
        return simpson(kernel, 0, 10, 50);
    }

    std::vector<size_t> half_spin; // Spin of each particle multiplied by two 
protected:
    vector2d model_rdf(double rho, double T, const vector1d& mole_fracs) override {throw std::runtime_error("Method model_rdf not implemented for Quantum!");}
    vector2d model_mtl(double rho, double T, const vector1d& x) override {throw std::runtime_error("Method model_mtl not implemented for Quantum!");}
    vector2d model_etl(double rho, double T, const vector1d& x) override {throw std::runtime_error("Method model_etl not implemented for Quantum!");}

private:
    

};