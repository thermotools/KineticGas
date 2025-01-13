#include "multiparam.h"
#include "Factorial.h"
#include <algorithm>

ModTangToennis::ModTangToennis(TangToennisParam param, vector1d mole_weights, bool is_idealgas)
    : Spherical(mole_weights, 
                vector2d(Ncomps, vector1d(Ncomps, param.sigma)), 
                vector2d(Ncomps, vector1d(Ncomps, param.eps_div_k * BOLTZMANN)), 
                is_idealgas, true), param{param}
{
    if (!is_idealgas) throw std::runtime_error("Modified Tang-Toennis only implemented for ideal gas!");
}

ModTangToennis::ModTangToennis(std::string comps, bool is_idealgas, std::string parameter_ref)
    : Spherical(comps, is_idealgas), param()
{
    const auto cdata = compdata[0]["ModTangToennis"][parameter_ref];
    const double A_div_k = cdata["A_div_k"];
    const double b = cdata["b"];
    const double A_tilde = cdata["A_tilde_div_k"];
    const vector1d a = cdata["a"];
    const double a_tilde = cdata["a_tilde"];
    const double eps_div_k = cdata["eps_div_k"];
    const double Re = cdata["Re"];
    const double sigma_ = cdata["sigma"];
    vector1d C = cdata["C"];

    param = TangToennisParam(A_div_k, b, A_tilde, a, a_tilde, eps_div_k, Re, sigma_, C);
    
    eps = vector2d(Ncomps, vector1d(Ncomps, param.eps_div_k * BOLTZMANN));
    sigma = vector2d(Ncomps, vector1d(Ncomps, param.sigma));
}

dual2 ModTangToennis::potential(int i, int j, dual2 r){
    r *= 1e9; // Using nm internally
    if (r < 0.4 * param.Re){
        return (param.A_tilde / r) * exp(- param.a_tilde * r) * BOLTZMANN;
    }
    dual2 u = param.A * exp(param.a1 * r + param.a2 * pow(r, 2) + param.am1 * pow(r, -1) + param.am2 * pow(r, -2));
    dual2 exp_prefactor = exp(- param.b * r);
    for (int n = 3; n <= 8; n++){
        dual2 tmp = 0.;
        int k = 0;
        // Prevent factorial overflow: Compute (power / factorial) in several steps, so that each individual fraction 
        // can be computed without overflow, before multiplying the fractions. Otherwise we get issues with either the
        // numerator, denominator, or both, overflowing, even though the evaluated fraction is a reasonable number.
        for (; k <= std::min(2 * n, 10); k++){
            tmp += pow(param.b * r, k) / partialfactorial(1, k);
        }
        for (; k <= 2 * n; k++){
            tmp += (pow(param.b * r, 10) / partialfactorial(1, 10)) * (pow(param.b * r, k - 10) / partialfactorial(10, k));
        }
        u -= (param.C[n - 3] / pow(r, 2 * n)) * (1 - exp_prefactor * tmp);
    }
    return u * BOLTZMANN;
}

double ModTangToennis::potential(int i, int j, double r){
    r *= 1e9; // Using nm internally
    if (r < 0.4 * param.Re){
        return (param.A_tilde / r) * exp(- param.a_tilde * r) * BOLTZMANN;
    }
    double u = param.A * exp(param.a1 * r + param.a2 * pow(r, 2) + param.am1 * pow(r, -1) + param.am2 * pow(r, -2));
    double exp_prefactor = exp(- param.b * r);
    for (int n = 3; n <= 8; n++){
        double tmp = 0.;
        int k = 0;
        for (; k <= std::min(2 * n, 10); k++){
            tmp += pow(param.b * r, k) / partialfactorial(1, k);
        }
        for (; k <= 2 * n; k++){ // Prevent factorial overflow
            tmp += (pow(param.b * r, 10) / partialfactorial(1, 10)) * (pow(param.b * r, k - 10) / partialfactorial(10, k));
        }
        u -= (param.C[n - 3] / pow(r, 2 * n)) * (1 - exp_prefactor * tmp);
    }
    return u * BOLTZMANN;
}

HFD_B2::HFD_B2(std::string comps) : Quantum(comps) {
    const auto cdata = compdata[0]["HFD-B2"]["default"];
    param.A = cdata["A"];
    param.alpha = cdata["alpha"];
    param.c6 = cdata["c6"];
    param.c8 = cdata["c8"];
    param.c10 = cdata["c10"];

    param.c_vec[0] = param.c6;
    param.c_vec[1] = param.c8;
    param.c_vec[2] = param.c10;

    param.C6 = cdata["C6"];
    param.C8 = cdata["C8"];
    param.C10 = cdata["C10"];
    param.beta_star = cdata["beta_star"];
    param.beta = cdata["beta"];
    param.D = cdata["D"];
    param.eps_div_k = cdata["eps_div_k"];
    param.rm = cdata["rm"];
    param.sigma = cdata["sigma"];

    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = 0; j < Ncomps; j++){
            sigma[i][j] = param.sigma;
            eps[i][j] = param.eps_div_k * BOLTZMANN;
        }
    }

    potential_terms.emplace_back(Polynomial::constant(param.A),
                                 Polynomial(1, 2, {-param.alpha, param.beta_star}));

    potential_terms.emplace_back(Polynomial(-10, -6, {param.c10, param.c8, param.c6}, 2),
                                 Polynomial(-2, 0, {- pow(param.D, 2), 2 * param.D, -1}));
    
    potential_terms.emplace_back(Polynomial(-10, -6, {param.c10, param.c8, param.c6}, 2),
                                 Polynomial::zero());
    
    
}

dual2 HFD_B2::potential(int i, int j, dual2 r){
    r *= 1e9; // working in nm internally
    dual2 x = r / param.rm;
    dual2 V = param.A * exp(- param.alpha * x + param.beta_star * pow(x, 2));

    if (x < 1e-10) return eps[i][j] * V;

    dual2 F = (x > param.D) ? static_cast<dual2>(1.) : exp(-pow((param.D / x) - 1, 2));
    dual2 C = 0;
    for (size_t j = 0; j <= 2; j++){
        C += param.c_vec[j] / pow(x, 2 * j + 6);
    }
    V -= F * C;
    return V * eps[i][j];
}

double HFD_B2::potential(int i, int j, double r){
    r *= 1e9; // working in nm internally
    double x = r / param.rm;
    double V = param.A * exp(- param.alpha * x + param.beta_star * pow(x, 2));

    if (x < 1e-10) return eps[i][j] * V;

    double F = (x > param.D) ? 1. : exp(-pow((param.D / x) - 1, 2));
    double C = 0;
    for (size_t j = 0; j <= 2; j++){
        C += param.c_vec[j] / pow(x, 2 * j + 6);
    }
    V -= F * C;
    return V * param.eps_div_k * BOLTZMANN;
}

double HFD_B2::potential_dn(int i, int j, double r, size_t n){
    r *= 1e9; // Working in nm internally
    double x = r / param.rm;
    double scaling = eps[i][j] * pow(1e9 / param.rm, n);

    double dvdx = potential_terms[0].derivative(x, n);
    if (x < 1e-10) return  dvdx * scaling;

    if (x < param.D) {
        double FC = potential_terms[1].derivative(x, n);
        dvdx -= FC;
    }
    else {
        double FC = potential_terms[2].derivative(x, n);
        dvdx -= FC;
    }

    return dvdx * scaling;
}

Patowski::Patowski(std::string comps)
    : Quantum(comps)
{
    const auto cdata = compdata[0]["Patowski"]["default"];
    param.Rc = cdata["Rc"];
    param.Ac = cdata["Ac"];
    param.Bc = cdata["Bc"];
    param.Cex1 = cdata["Cex1"];
    param.Cex2 = cdata["Cex2"];
    param.Csp1 = cdata["Csp1"];
    param.Csp2 = cdata["Csp2"];
    param.Csp3 = cdata["Csp3"];
    param.Csp4 = cdata["Csp4"];
    param.delta = cdata["delta"];
    param.C6 = cdata["C6"];
    param.C8 = cdata["C8"];
    param.C10 = cdata["C10"];
    param.Cn[0] = param.C6; param.Cn[1] = param.C8; param.Cn[2] = param.C10;
    param.sigma = cdata["sigma"];
    param.r_min = cdata["r_min"];
    param.eps_div_k = cdata["eps_div_k"];

    for (size_t i = 0; i < Ncomps; i++){
        for (size_t j = 0; j < Ncomps; j++){
            sigma[i][j] = param.sigma;
            eps[i][j] = param.eps_div_k * BOLTZMANN;
        }
    }

    potential_terms.emplace_back(Polynomial(0, 3, {param.Csp1, param.Csp2, param.Csp3, param.Csp4}), 
                                 Polynomial(0, 1, {param.Cex1, param.Cex2}));
    potential_terms.emplace_back(Polynomial(-10, -6, {param.C10, param.C8, param.C6}, 2),
                                 Polynomial::zero());

    for (int n = 3; n <= 5; n++){
        std::vector<double> coeff;
        for (int k = 0; k <= 2 * n; k++){
            coeff.push_back(- param.Cn[n - 3] * pow(param.delta, k) / Fac(k).eval_d());
        }
        potential_terms.emplace_back(Polynomial(- 2 * n, 0, coeff),
                                     Polynomial::linear(- param.delta));
    }
}

dual2 Patowski::potential(int i, int j, dual2 r){
    r *= 1e10; // Working in Å internally
    if (r < param.Rc){
        // Potential is only valid for r > Rc, but we need a continuous extrapolation that is well behaved for numerical purposes.
        // In practice, something is probably very wrong if we are ever at distances r < Rc, except when iterating some solver.
        // This extension is made by requiring a contiuous potential and first derivative at r = Rc.
        return param.Ac * exp(param.Bc * (r - param.Rc)) * BOLTZMANN;
    }
    dual2 p = (param.Csp1 + r * param.Csp2 + pow(r, 2) * param.Csp3 + pow(r, 3) * param.Csp4) * exp(param.Cex1 + param.Cex2 * r);
    for (size_t n = 3; n <= 5; n++){
        dual2 tmp = 0;
        for (int k = 0; k <= 2 * n; k++){
            tmp += pow(param.delta * r, k) / Fac(k).eval_d();
        }
        p += (param.Cn[n - 3] / pow(r, 2 * n)) * (1 - tmp * exp(- param.delta * r));
    }
    return p * BOLTZMANN;
}

double Patowski::potential(int i, int j, double r){
    r *= 1e10; // Working in Å internally
    if (r < param.Rc){ 
        // Potential is only valid for r > Rc, but we need a continuous extrapolation that is well behaved for numerical purposes.
        // In practice, something is probably very wrong if we are ever at distances r < Rc, except when iterating some solver.
        // This extension is made by requiring a contiuous potential and first derivative at r = Rc.
        return param.Ac * exp(param.Bc * (r - param.Rc)) * BOLTZMANN;
    }
    double p = (param.Csp1 + r * param.Csp2 + pow(r, 2) * param.Csp3 + pow(r, 3) * param.Csp4) * exp(param.Cex1 + param.Cex2 * r);
    for (size_t n = 3; n <= 5; n++){
        double tmp = 0;
        for (int k = 0; k <= 2 * n; k++){
            tmp += pow(param.delta * r, k) / Fac(k).eval_d();
        }
        p += (param.Cn[n - 3] / pow(r, 2 * n)) * (1 - tmp * exp(- param.delta * r));
    }
    return p * BOLTZMANN;
}

double Patowski::potential_dn(int i, int j, double r, size_t n){
    // std::cout << "potential_dn : " << n << std::endl;
    // if (n == 0) return potential(i, j, r);
    r *= 1e10; // Working in Å internally
    if (r < param.Rc){
        return param.Ac * pow(param.Bc, n) * exp(param.Bc * (r - param.Rc)) * BOLTZMANN;
    }
    double val = 0;
    for (const PolyExp& term : potential_terms){
        double t = term.derivative(r, n);
        val += t;
    }
    return val * BOLTZMANN * pow(1e10, n);
}

PatowskiFH1::PatowskiFH1(std::string comps) 
    : Patowski(comps), D_factors(Ncomps, vector1d(Ncomps, 0.)) 
{
    set_quantum_active(false);
}

double PatowskiFH1::potential(int i, int j, double r){
    dual4th rd = r;
    auto [u, u1, u2, u3, u4] = autodiff::derivatives([&](dual4th r_){return core_potential(i, j, r_);}, autodiff::wrt(rd), autodiff::at(rd));
    return u + D_factors[i][j] * u2;
}

double PatowskiFH1::potential_derivative_r(int i, int j, double r){
    dual4th rd = r;
    auto [u, u1, u2, u3, u4] = autodiff::derivatives([&](dual4th r_){return core_potential(i, j, r_);}, autodiff::wrt(rd), autodiff::at(rd));
    return u1 + D_factors[i][j] * u3;
}

double PatowskiFH1::potential_dblderivative_rr(int i, int j, double r){
    dual4th rd = r;
    auto [u, u1, u2, u3, u4] = autodiff::derivatives([&](dual4th r_){return core_potential(i, j, r_);}, autodiff::wrt(rd), autodiff::at(rd));
    return u2 + D_factors[i][j] * u4;
}

dual4th PatowskiFH1::core_potential(int i, int j, dual4th r){
    r *= 1e10; // Working in Å internally
    if (r < param.Rc){
        // Potential is only valid for r > Rc, but we need a continuous extrapolation that is well behaved for numerical purposes.
        // In practice, something is probably very wrong if we are ever at distances r < Rc, except when iterating some solver.
        // This extension is made by requiring a contiuous potential and first derivative at r = Rc.
        return param.Ac * exp(param.Bc * (r - param.Rc)) * BOLTZMANN;
    }
    dual4th p = (param.Csp1 + r * param.Csp2 + pow(r, 2) * param.Csp3 + pow(r, 3) * param.Csp4) * exp(param.Cex1 + param.Cex2 * r);
    for (size_t n = 3; n <= 5; n++){
        dual4th tmp = 0;
        for (int k = 0; k <= 2 * n; k++){
            tmp += pow(param.delta * r, k) / Fac(k).eval_d();
        }
        p += (param.Cn[n - 3] / pow(r, 2 * n)) * (1 - tmp * exp(- param.delta * r));
    }
    return p * BOLTZMANN;
}

size_t PatowskiFH1::set_internals(double rho, double T, const vector1d& x){
    return set_current_T(T);
}

size_t PatowskiFH1::set_current_T(double T){
    if (T != current_T){
        double beta = 1 / (BOLTZMANN * T);
        for (size_t i = 0; i < Ncomps; i++){
            for (size_t j = 0; j < Ncomps; j++){
                D_factors[i][j] = beta * pow(HBAR, 2) / (24 * red_mass[i][j]);
            }
        }
        current_T = T;
        return 1;
    }
    return 0;
}