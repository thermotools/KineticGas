#pragma once
#include "KineticGas.h"
#include "Spherical.h"
#include "Quantum.h"
#include "global_params.h"
#include "extentions.h"
#include <autodiff/forward/dual.hpp>

using namespace autodiff;

struct TangToennisParam{
    double A, b, A_tilde, a_tilde, eps_div_k, Re, sigma;
    // C = {C6, C8, C10, C12, C14, C16}
    vector1d C;
    double a1, a2, am1, am2;

    TangToennisParam() = default;
    TangToennisParam(double A, double b, double A_tilde, vector1d a,
                    double a_tilde, double eps_div_k, double Re, double sigma, vector1d C)
                    : A{A}, b{b}, A_tilde{A_tilde}, a_tilde{a_tilde}, eps_div_k{eps_div_k}, Re{Re}, sigma{sigma}, C{C}
        {
        a1 = a[0]; a2 = a[1];
        am1 = a[2]; am2 = a[3];
        }
};

class ModTangToennis : public Quantum {
    public:
    TangToennisParam param;

    ModTangToennis(std::string comps, bool is_idealgas, std::string parameter_ref="default");

    dual2 potential(int i, int j, dual2 r) override;
    double potential(int i, int j, double r) override;

    vector2d model_rdf(double rho, double T, const vector1d& x) override {
        throw std::runtime_error("Modified Tang-Toennis only implemented for ideal gas!");
    }

    TangToennisParam get_param(){return param;};
};

struct HFD_B2_Param {
    double A, alpha, c6, c8, c10, C6, C8, C10, beta_star, beta, D, eps_div_k, rm, sigma;
    std::array<double, 3> c_vec = {0, 0, 0};
};

class HFD_B2 : public Quantum {
public:
    HFD_B2(std::string comps);
    HFD_B2(HFD_B2_Param param);

    dual2 potential(int i, int j, dual2 r) override;
    double potential(int i, int j, double r) override;
    double potential_dn(int i, int j, double r, size_t n) override;

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
    dual2 potential(int i, int j, dual2 r) override;
    double potential(int i, int j, double r) override;
    double potential_dn(int i, int j, double r, size_t n) override;

    inline PatowskiParam get_param(){return param;}
protected:
    PatowskiParam param;
    std::vector<PolyExp> potential_terms;
};

class PatowskiFH1 : public Patowski {
public:

    PatowskiFH1(std::string comps);

    double potential(int i, int j, double r, double T) {
        set_internals(0., T, {0.});
        return potential(i, j, r);
    }

    double potential(int i, int j, double r) override;
    double potential_derivative_r(int i, int j, double r) override;
    double potential_dblderivative_rr(int i, int j, double r) override;

    dual4th core_potential(int i, int j, dual4th r);

    size_t set_internals(double rho, double T, const vector1d& x) override;

private:
    size_t set_current_T(double T);
    double current_T;
    vector2d D_factors;
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
// using PatowskiFH = FH_Corrected<Patowski>;
// using PatowskiFH1 = FH_Corrected<Patowski, 1>; // Tabulated<PatowskiFH1, 1000, 3>;
// using PatowskiFH2 = FH_Corrected<Patowski, 2>; // FH_Corrected<Splined<PatowskiCore, 1000, 6>, 2>;
// using PatowskiFH3 = FH_Corrected<Patowski, 3>; // FH_Corrected<Splined<PatowskiCore, 1000, 8>, 3>;
// using PatowskiFH5 = FH_Corrected<Splined<PatowskiCore, 500, 10>, 4>; // Spline degree must be at least 2 * FH_order + 2