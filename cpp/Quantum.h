#include "KineticGas.h"
#include "Integration/Integration.h"

class Quantum : public KineticGas {
public:
    Quantum(std::string comps) : KineticGas(comps, false) {

    }

    double potential(int i, int j, double r);
    double potential_r(int i, int j, double r);
    double potential_rr(int i, int j, double r);

    double phase_shift(int i, int j, double E);
    double cross_section_kernel(int i, int j, double n, double l, double E){
        if (n == 0){
            return (2 * l + 1) * pow(sin(phase_shift(i, j, l, E), 2));
        }
        double A{0.}, B{0.};
        double q;
        for (size_t k = 0; (k < 3) && (2 * k + 1 <= n); k++){
            B = (2 * l + 1) * pow(sin(phase_shift(i, j, l, E) - phase_shift(i, j, l + n - 2 * k, E)));
            double B_tmp = 0;
            for (size_t p = 1; p <= n - 2 * k; p++){
                B_tmp += (l + p) / (2 * (l + p) - 1);
            }
            B *= B_tmp;
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
            if (half_spin[i] % != 0) { // Odd Half-Integer spin, Fermions
                double tmp = even_prefactor;
                even_prefactor = odd_prefactor;
                odd_prefactor = tmp;
            }
        }
        double Q = 0;
        if (even_prefactor > 0.){
            double Q_even = 0;
            int l_even = 0;
            do {
                double q_even = cross_section_kernel(i, j, n, l_even, E);
                Q_even += q_even;
                l_even += 2;
            } while (abs(q_even) > 1e-6 * Q_even);
            Q += even_prefactor * Q_even;
        } 
        if (odd_prefactor > 0.){
            double Q_odd = 0;
            int l_odd = 1;
            do {
                double q_odd = cross_section_kernel(i, j, n, l_odd, E);
                Q_odd += q_odd;
                l_odd += 2;
            } while (abs(q_odd) > 1e-6 * Q_odd);
            Q += odd_prefactor * Q_odd;
        } 
        kappa_mul_E_squared = 2. * m[i] * m[j] / ((m[i] + m[j]) * HBAR);
        Q *= 4. * PI / pow(kappa, 2);
        return Q;
    }

    double omega(int i, int j, int n, int s, double T){
        double beta = 1. / (T * BOLTZMANN);
        double sfac = 1;
        for (int si = 2; si <= s + 1; si++){
            sfac *= si;
        }
        const auto kernel = [&](double E){return cross_section(i, j, n, E) * exo(- beta * E) * pow(beta * E, s + 2) / sfac;};
        return simpson_inf(kernel, 0);
    }

protected:

    vector2d model_rdf(double rho, double T, const vector1d& mole_fracs) override {throw std::runtime_error("Method model_rdf not implemented for Quantum!")};
    vector2d model_mtl(double rho, double T, const vector1d& x) override {throw std::runtime_error("Method model_mtl not implemented for Quantum!")};
    vector2d model_etl(double rho, double T, const vector1d& x) override {throw std::runtime_error("Method model_etl not implemented for Quantum!")};

private:
    std::vector<size_t> half_spin; // Spin of each particle multiplied by two 

};