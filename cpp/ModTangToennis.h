#pragma once
#include "KineticGas.h"
#include "Spherical.h"
#include "global_params.h"
#include <autodiff/forward/dual.hpp>

using namespace autodiff;

struct TangToennisParam{
    double A, b, A_tilde, a_tilde, eps_div_k, Re, sigma;
    // C = {C6, C8, C10, C12, C14, C16}
    vector1d C;
    double a1, a2, am1, am2;

    TangToennisParam(double A, double b, double A_tilde, vector1d a,
                    double a_tilde, double eps_div_k, double Re, double sigma, vector1d C)
                    : A{A}, b{b}, A_tilde{A_tilde}, a_tilde{a_tilde}, eps_div_k{eps_div_k}, Re{Re}, C{C}
        {
        a1 = a[0]; a2 = a[1];
        am1 = a[2]; am2 = a[3];
        }
};

class ModTangToennis : public Spherical {
    public:
    TangToennisParam param;
    ModTangToennis(TangToennisParam param, vector1d mole_weights, vector2d sigma, bool is_idealgas);

    dual2 potential(int i, int j, dual2 r) override;

    vector2d model_rdf(double rho, double T, const vector1d& x){
        throw std::runtime_error("Modified Tang-Toennis only implemented for ideal gas!");
    }

};