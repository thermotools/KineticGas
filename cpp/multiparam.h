#pragma once
#include "KineticGas.h"
#include "Spherical.h"
#include "Quantum.h"
#include "global_params.h"
#include "extentions.h"
#include <autodiff/forward/dual.hpp>

using namespace autodiff;

struct TangToennisParam{
    double A, b, A_tilde, a_tilde, a1, a2, am1, am2, eps_div_k, Re, sigma, L_unit, short_range_lim;
    // C = {C6, C8, C10, C12, C14, C16}
    vector1d C;
    vector1d C_exp;

    TangToennisParam() = default;
    TangToennisParam(double A, double b, double A_tilde, vector1d a,
                    double a_tilde, double eps_div_k, double Re, double sigma, double L_unit, vector1d C)
                    : A{A}, b{b}, A_tilde{A_tilde}, a_tilde{a_tilde}, eps_div_k{eps_div_k}, Re{Re}, sigma{sigma}, L_unit{L_unit}, C{C}
        {
        a1 = a[0]; a2 = a[1];
        am1 = a[2]; am2 = a[3];
        }
    TangToennisParam(const json& param)
        : A{param["A_div_k"]}, b{param["b"]}, A_tilde{param["A_tilde_div_k"]}, a_tilde{param["a_tilde"]},
        a1{param["a1"]}, a2{param["a2"]}, am1{param["am1"]}, am2{param["am2"]},
        eps_div_k{param["eps_div_k"]}, Re{param["Re"]}, sigma{param["sigma"]}, L_unit{param["L_unit"]},
        short_range_lim{param["short_range_lim"]},
        C(param["C"]), C_exp(17, 0.)
    {
        for (int n = 3; n <= 8; n++){
            int k = 0;
            for (; k <= 2 * n; k++){
                C_exp[2 * n - k] += C[n - 3] * pow(b, k) / partialfactorial(1, k);
            }
        }
    }
};

class ModTangToennis : public Quantum {
public:
    TangToennisParam param;

    ModTangToennis(std::string comps, std::string parameter_ref="default");

    dual2 potential(int i, int j, dual2 r) const override;
    double potential(int i, int j, double r) const override;
    double potential_dn(int i, int j, double r, size_t n) const override;
    double potential_derivative_r(int i, int j, double r) const override {return potential_dn(i, j, r, 1);}
    double potential_dblderivative_rr(int i, int j, double r) const override {return potential_dn(i, j, r, 2);}

    vector2d model_rdf(double rho, double T, const vector1d& x) override {
        throw std::runtime_error("Modified Tang-Toennis only implemented for ideal gas!");
    }

    TangToennisParam get_param(){return param;};

private:
    std::vector<PolyExp> potential_terms;
    PolyExp short_range_potential;
};

class FH_ModTangToennies : public FH_Corrected<ModTangToennis> {
public:
    FH_ModTangToennies(std::string comps, size_t FH_order, std::string parameter_ref)
        : FH_Corrected<ModTangToennis>(FH_order, comps, parameter_ref)
    {
        set_quantum_active(false);
    }
};

struct HFD_B2_Param {
    double A, alpha, c6, c8, c10, C6, C8, C10, beta_star, beta, D, eps_div_k, rm, sigma;
    std::array<double, 3> c_vec = {0, 0, 0};
};

class HFD_B2 : public Quantum {
public:
    HFD_B2(std::string comps);
    HFD_B2(HFD_B2_Param param);

    dual2 potential(int i, int j, dual2 r) const override;
    double potential(int i, int j, double r) const override;
    double potential_dn(int i, int j, double r, size_t n) const override;

    HFD_B2_Param get_param(){return param;}
private:
    HFD_B2_Param param;
    std::vector<PolyExp> potential_terms;
};

class FH_HFD_B2 : public FH_Corrected<HFD_B2> {
public:
    FH_HFD_B2(std::string comps, size_t FH_order)
        : FH_Corrected<HFD_B2>(FH_order, comps)
    {
        set_quantum_active(false);
    }
};

struct PatowskiParam{
    double Rc, Ac, Bc, Cex1, Cex2, Csp1, Csp2, Csp3, Csp4, delta, C6, C8, C10, sigma, eps_div_k, r_min;
    std::array<double, 3> Cn;     
};

class Patowski : public Quantum {
public:
    Patowski(std::string comps);
    Patowski(PatowskiParam param);

    using Spherical::potential;
    dual2 potential(int i, int j, dual2 r) const override;
    double potential(int i, int j, double r) const override;
    double potential_dn(int i, int j, double r, size_t n) const override;

    inline PatowskiParam get_param(){return param;}
protected:
    PatowskiParam param;
    std::vector<PolyExp> potential_terms;
};

class PatowskiFH : public FH_Corrected<Patowski> {
public:
    PatowskiFH(std::string comps, size_t FH_order)
        : FH_Corrected<Patowski>(FH_order, comps)
    {
        set_quantum_active(false);
    }
};

using PatowskiTab = Tabulated<Patowski, 1000, 3>;