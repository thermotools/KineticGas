#include "multiparam.h"
#include "Factorial.h"
#include <algorithm>

ModTangToennis::ModTangToennis(std::string comps, std::string parameter_ref)
    : Quantum(comps), param(compdata[0]["ModTangToennis"][parameter_ref]), 
    short_range_potential(Polynomial(-1, -1, {param.A_tilde}), Polynomial::linear(- param.a_tilde))
{    
    eps = vector2d(Ncomps, vector1d(Ncomps, param.eps_div_k * BOLTZMANN));
    sigma = vector2d(Ncomps, vector1d(Ncomps, param.sigma));

    potential_terms.emplace_back(Polynomial::constant(param.A),
                                 Polynomial(-2, 2, {param.am2, param.am1, 0., param.a1, param.a2}));
    potential_terms.emplace_back(Polynomial(-16, -6, {-param.C[5], -param.C[4], -param.C[3], -param.C[2], -param.C[1], -param.C[0]}, 2),
                                Polynomial::zero());

    // The double sum over n and k in the attractive part of the potential can be compressed to a 
    // single polynomial. That is what is going on here, it saves us some significant runtime.
    vector1d coeff(17, 0.);
    for (int n = 3; n <= 8; n++){
        int k = 0;
        for (; k <= 2 * n; k++){
            coeff[2 * n - k] += param.C[n - 3] * pow(param.b, k) / partialfactorial(1, k);
        }
    }
    vector1d coeff_r(coeff.rbegin(), coeff.rend());
    potential_terms.emplace_back(Polynomial(-16, 0, coeff_r),
                                Polynomial::linear(- param.b));
    
}

dual2 ModTangToennis::potential(int i, int j, dual2 r){
    r /= param.L_unit; // Convert to units used in parameter set (probably Å or nm)
    if (r < param.short_range_lim){
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
    r /= param.L_unit; // Using nm internally
    if (r < param.short_range_lim){
        return (param.A_tilde / r) * exp(- param.a_tilde * r) * BOLTZMANN;
    }
    std::array<double, 17> rpow;
    for (size_t n = 0; n < 17; n++){
        rpow[n] = pow(r, n);
    }

    double u = param.A * exp(param.a1 * rpow[1] + param.a2 * rpow[2] + param.am1 / rpow[1] + param.am2 / rpow[2]);
    double exp_prefactor = exp(- param.b * r);
    for (int n = 3; n <= 8; n++){
        u -= (param.C[n - 3] / rpow[2 * n]);
    }
    for (int n = 0; n <= 16; n++){
        u += param.C_exp[n] * exp_prefactor / rpow[n];
    }
    return u * BOLTZMANN;
}

double ModTangToennis::potential_dn(int i, int j, double r, size_t n){
    r /= param.L_unit;
    double val = 0;
    if (r < param.short_range_lim){
        val = short_range_potential.derivative(r, n);
    }
    else {
        for (const PolyExp& term : potential_terms){
            val += term.derivative(r, n);
        }
    }
    return val * BOLTZMANN / pow(param.L_unit, n);
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