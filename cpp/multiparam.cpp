#include "multiparam.h"
#include <algorithm>

double partialfactorial(int start, int stop){
    double fac = 1.;
    start = (start == 0) ? 1 : start;
    for (int i = start; i <= stop; i++){
        fac *= i;
    }
    return fac;
}

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
    const double A_div_k = compdata[0]["ModTangToennis"][parameter_ref]["A_div_k"];
    const double b = compdata[0]["ModTangToennis"][parameter_ref]["b"];
    const double A_tilde = compdata[0]["ModTangToennis"][parameter_ref]["A_tilde_div_k"];
    const vector1d a = compdata[0]["ModTangToennis"][parameter_ref]["a"];
    const double a_tilde = compdata[0]["ModTangToennis"][parameter_ref]["a_tilde"];
    const double eps_div_k = compdata[0]["ModTangToennis"][parameter_ref]["eps_div_k"];
    const double Re = compdata[0]["ModTangToennis"][parameter_ref]["Re"];
    const double sigma_ = compdata[0]["ModTangToennis"][parameter_ref]["sigma"];
    vector1d C = compdata[0]["ModTangToennis"][parameter_ref]["C"];

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
        for (; k <= std::min(2 * n, 10); k++){
            tmp += pow(param.b * r, k) / partialfactorial(1, k);
        }
        for (; k <= 2 * n; k++){ // Prevent factorial overflow
            tmp += (pow(param.b * r, 10) / partialfactorial(1, 10)) * (pow(param.b * r, k - 10) / partialfactorial(11, k));
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
            tmp += (pow(param.b * r, 10) / partialfactorial(1, 10)) * (pow(param.b * r, k - 10) / partialfactorial(11, k));
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