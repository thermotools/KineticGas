/*
References:
    svrm:
        Accurate statistical associating fluid theory for chain molecules formed from Mie segments
        Thomas Lafitte; Anastasia Apostolakou; Carlos Avendaño; Amparo Galindo; Claire S. Adjiman; Erich A. Müller; George Jackson
        J. Chem. Phys. 139, 154504 (2013)
        https://doi.org/10.1063/1.4819786

-------      TAKE HEED YE WHO COMETH HERE AND DESIRES TO INHERIT THIS CLASS      --------

Because this class implements a potential that is in general dependent on both temperature 
and density, we need to use set_internals actively. Additionally, we want to use hyperduals,
so the class implements `set_effective_params` which is used to propagate derivatives.

Essentially, the class holds a bunch of "effective parameters" that are valid only at the current
density and temperature, in order to not need to re-compute these parameters in every call in a given
call-chain. The setting of these parameters is handled by `set_effective_params`.

If you override this class, and want to implement a potential that is either only temperature 
dependent, or only density dependent, still want to take advantage of caching of collision integrals
and transfer lengths, you will need to override `set_effective_params` as appropriate.

Further, the effective parameters from the previous calculation are used as initial guesses for the 
next calculation. This gives a huge speedup when doing e.g calculations along an isoline, without
a notable cost at other times. The caveat is that you need to initialize the effective parameters
to something reasonable in the constructor. If you are inheriting this class, you should call
the `init_effective_params` method after setting all potential parameters in your constructor.
    - Note: This is not neccesary if you forward the construction call to one of the ExtSutherland
            constructors that takes the potential parameters as input.

----------------------------     YE HATH NOW TAKEN HEED      ----------------------------
*/
#pragma once
#include "KineticGas.h"
#include "Spherical.h"
#include <vector>

using vector1d = std::vector<double>;
using vector2d = std::vector<vector1d>;
using vector3d = std::vector<vector2d>;

using vector1d2 = std::vector<dual2>;
using vector2d2 = std::vector<vector1d2>;
using vector3d2 = std::vector<vector2d2>;

using vector1d1 = std::vector<dual>;
using vector2d1 = std::vector<vector1d1>;
using vector3d1 = std::vector<vector2d1>;

using namespace autodiff;

inline void throw_notimplemented(){
    throw std::runtime_error("This method is not implemented for class ExtSutherland!");
}

vector2d dual_to_double(const vector2d2& vin);

struct RDFConstants{
    double gl_x[20]; // Gauss Legendre points for computing barker henderson diamenter
    double gl_w[20]; // Gauss Legendre points for computing barker henderson diamenter
    double C_coeff_matr[4][4]; // See Eq. (A18) in svrm (https://doi.org/10.1063/1.4819786)
    double phi[7][3]; // Table II in svrm (https://doi.org/10.1063/1.4819786)
};

class ExtSutherland : public Spherical{
public:
    ExtSutherland(std::string comps, size_t nterms, bool is_idealgas=false);
    ExtSutherland(vector1d mole_weights, vector2d sigma, vector2d eps, 
                size_t nterms, bool is_idealgas=false, bool is_singlecomp=false);
    ExtSutherland(vector1d mole_weights, vector2d sigma, vector2d eps, 
                vector3d C, vector3d lambda, vector3d beta_exp, vector3d rho_exp, 
                bool is_idealgas=false, bool is_singlecomp=false);

    dual2 potential(int i, int j, dual2 r, dual2 T, dual2 rho);
    dual2 potential_r(int i, int j, dual2 r, dual2 T, dual2 rho);
    dual2 potential_rr(int i, int j, dual2 r, dual2 T, dual2 rho);
    double potential(int i, int j, double r, double T, double rho);

     // To directly compute the RDF at different pertubation orders.
    vector2d saft_rdf(double rho, double T, const std::vector<double>& x, int order=2, bool g2_correction=true);
    vector3d get_rdf_terms(double rho, double T, const vector1d& x); // Return g0, g1, g2 (no correction), g2 (with correction)
    virtual vector2d2 get_BH_diameters(dual2 rho, dual2 T);
    vector2d get_BH_diameters(double rho, double T);

    std::string potential_params_to_string(size_t i){return potential_params_to_string(i, i);}
    std::string potential_params_to_string(size_t i, size_t j){
        std::stringstream strm;
        strm << "Sigma : " << sigma[i][j] << "\nEpsilon : " << eps[i][j] << "\n";
        strm << "Term   |   C   |  beta_exp   |  rho_exp   |   lambda   \n";
        for (size_t ti = 0; ti < nterms; ti++){
            strm << ti << " | " << C[ti][i][j] << " | " << beta_exp[ti][i][j] << " | " << rho_exp[ti][i][j] << " | " << lambda[ti][i][j] << "\n";
        }
        return strm.str();
    }

    vector2d get_sigma_eff(double rho, double T){
        set_sigma_eff(rho, T);
        return dual_to_double(sigma_eff);
    }

    vector2d get_sigma_min(double rho, double T){
        set_epsilon_eff(rho, T);
        return dual_to_double(r_min);
    }

    vector2d get_epsilon_eff(double rho, double T){
        set_epsilon_eff(rho, T);
        return dual_to_double(eps_eff);
    }

    vector2d get_vdw_alpha(double rho, double T){
        set_effective_params(rho, T);
        return dual_to_double(vdw_alpha);
    }

    // ------------------------------------------------------------------------------------------------------------------- //
    // -------------------------- Sutherland Internals are below here ---------------------------------------------------- //
    // ---------------------- End users should not need to care about anything below ------------------------------------- //
    // ------------------------------------------------------------------------------------------------------------------- //

protected:
    vector3d C;
    vector3d lambda;
    vector3d beta_exp;
    vector3d rho_exp;
    size_t nterms;

    // Because effective parameters from previous calculation are used as initial guesses in next calculation,
    // we need to set them to some value at init. We just set them to the potential parameters.
    // NOTE: Inheriting classes that set their potential parameters after calling the ExtSutherland constructor
    //      must call the init_effective_params method to ensure that effective parameters are initialised to something reasonable.
    void init_effective_params();
    void mix_sigma();
    void mix_epsilon();
    void mix_exponents(vector2d& expo);

    virtual void set_effective_params(dual2 rho, dual2 T);
    void set_internals(double rho, double T, const vector1d& x) override {set_effective_params(rho, T);}
    double current_T{-1.}, current_rho{-1.};
    bool C_set{false}, sigma_set{false}, eps_set{false}, alpha_set{false};
    vector2d2 sigma_eff;
    vector2d2 r_min;
    vector2d2 eps_eff;
    vector2d2 vdw_alpha;
    vector3d2 C_eff;

    OmegaPoint get_omega_point(int i, int j, int l, int r, double T) override {
        if (T != current_T){
            throw std::runtime_error("Something is very wrong ... (in ExtSutherland::get_omega_point)");
        }
        return OmegaPoint(i, j, l, r, T, current_rho);
    }

    StatePoint get_transfer_length_point(double rho, double T, const vector1d& x) override {
        if (T != current_T){
            throw std::runtime_error("Something is even more wrong ... (in ExtSutherland::get_transfer_length_point)");
        }
        return StatePoint(T, rho);
    }

    inline vector2d model_rdf(double rho, double T, const std::vector<double>& x) override {
        return saft_rdf(rho, T, x, 2, true);
    }

    dual2 potential(int i, int j, dual2 r);
    dual2 potential_r(int i, int j, dual2 r);
    dual2 potential_rr(int i, int j, dual2 r);
    double potential(int i, int j, double r) override;
    using Spherical::potential_derivative_r;
    using Spherical::potential_dblderivative_rr;

    void set_C_eff(dual2 rho, dual2 T);
    void set_sigma_eff(dual2 rho, dual2 T);
    void set_epsilon_eff(dual2 rho, dual2 T);
    void set_vdw_alpha(dual2 rho, dual2 T);

    vector2d get_b_max(double T) override;
    vector2d get_BH_diameters(double T);

    vector2d2 rdf_g0_func(dual2 rho, const vector1d& x, const vector2d2& d_BH, const vector2d2& x_eff);
    vector2d1 rdf_g1_func(dual2 rho, dual2 T, const vector1d& x, const vector2d2& d_BH, const vector2d2& x_eff);
    vector2d1 rdf_g2_func(dual2 rho, dual2 T, const vector1d& x, const vector2d2& d_BH, const vector2d2& x_eff, bool g2_correction);
    
    vector2d get_lambda_kl(size_t k, size_t l); // lambda_kl[i, j] = lambda_k[i, j] + lambda_l[i, j]
    vector2d2 get_x0(const vector2d2& d_BH); // x0 = sigma / d_BH
    vector2d2 get_xeff(const vector2d2& d_BH); // x_eff = sigma_eff / d_BH

    vector2d2 gamma_corr(dual2 zeta_x, dual2 T); // Eq. (A37) in svrm (https://doi.org/10.1063/1.4819786)
    dual2 zeta_x_func(dual2 rho, const vector1d& x, const vector2d2& d_BH); // Eq. (A13) in svrm
    dual2 zeta_eff_func(dual2 rho, const vector1d& x, dual2 zeta_x, double lambdaijk); // Eq. (A17) in svrm
    
    dual2 a1_func(int i, int j, dual2 rho, dual2 T, const vector1d& x);
    dual2 a1_func(int i, int j, dual2 rho, const vector1d& x, const vector2d2& d_BH);
    dual2 a2_div_chi_func(int i, int j, dual2 rho, dual2 T, const vector1d& x);
    dual2 a2_div_chi_func(int i, int j, dual2 rho, const vector1d& x, dual2 K_HS, const vector2d2& d_BH, const vector2d2& x_eff); // Eq. (A20) in svrm
    vector2d2 rdf_chi_func(dual2 rho, dual2 T, const vector1d& x); // Eq. (A22) in svrm

    vector1d2 f_corr(dual2 alpha); // Eq. (A26) in svrm

    dual2 a_1s_func(int i, int j, dual2 rho, const vector1d& x, dual2 zeta_x, const vector2d2& d_BH, const vector2d& lambda_k); // Eq. (A16) in svrm (https://doi.org/10.1063/1.4819786)

    dual2 B_func(int i, int j, dual2 rho, const vector1d& x, dual2 zeta_x, const vector2d2& x_eff, const vector2d2& d_BH, const vector2d& lambda_k); // Eq. (A12) in svrm
    dual2 I_func(int i, int j, const vector2d2& xeff, const vector2d& lambda_k); // Eq. (A14) in svrm (https://doi.org/10.1063/1.4819786)
    dual2 J_func(int i, int j, const vector2d2& xeff, const vector2d& lambda_k); // Eq. (A15) in svrm (https://doi.org/10.1063/1.4819786)

    inline dual2 K_HS_func(dual2 zeta_x){ // Eq. (A21) in svrm
        return pow(1 - zeta_x, 4) / (1 + 4 * zeta_x + 4 * pow(zeta_x, 2) - 4 * pow(zeta_x, 3) + pow(zeta_x, 4));
    }

    static constexpr RDFConstants rdf_constants{
                            { // Gauss-Legendre nodes (20 points), see get_BH_diameters
                            -0.9931285991850949, -0.9639719272779139,
                            -0.912234428251326, -0.8391169718222187,
                            -0.7463319064601508, -0.636053680726515,
                            -0.5108670019508271, -0.37370608871541955,
                            -0.22778585114164507, -0.07652652113349737,
                            0.07652652113349737, 0.22778585114164507,
                            0.37370608871541955, 0.5108670019508271,
                            0.636053680726515, 0.7463319064601508,
                            0.8391169718222187, 0.912234428251326,
                            0.9639719272779139, 0.9931285991850949
                            },
                            { // Gauss-Legendre weights (20 points), see get_BH_diameters
                            0.017614007139152742, 0.040601429800386134,
                            0.06267204833410799, 0.08327674157670514,
                            0.1019301198172403, 0.11819453196151845,
                            0.13168863844917675, 0.14209610931838218,
                            0.149172986472604, 0.15275338713072611,
                            0.15275338713072611, 0.149172986472604,
                            0.14209610931838218, 0.13168863844917675,
                            0.11819453196151845, 0.1019301198172403,
                            0.08327674157670514, 0.06267204833410799,
                            0.040601429800386134, 0.017614007139152742
                            },
                            { // C_coeff_matr (See: Eq. (A17-A18) in svrm).
                                // Parameters regressed by T. Maltby, without LRC (2024, unpublished)
                                {0.6382, 0.9540, -31.9965, 83.7914},
                                {2.8596, -18.8315, 150.6853, -468.9992},
                                {-7.0331, 21.9390, -229.9483, 970.6560},
                                {5.5441, -4.0284, 106.4364, -679.1889}

                                // Parameters from Lafitte et al. (Table II in svrm)
                                // {0.81096, 1.7888, -37.578, 92.284},
                                // {1.0205, -19.341, 151.26, -463.50},
                                // {-1.9057, 22.845, -228.14, 973.92},
                                // {1.0885, -6.1962, 106.98, -677.64}

                                // Parameters from T. Maltby, with LRC (2024, unpublished)
                                // {0.9745, -0.1438, -5.2072, -105.4317},
                                // {-1.0980, -10.9314, 252.9727, -274.7130},
                                // {7.3301, -69.6106, -386.7004, 1045.6686},
                                // {-9.7649, 140.1891, -15.4988, -652.1140}

                            },
                            { // phi (See: Table II in svrm)
                                {7.5365557, -359.44 , 1550.9 },
                                {-37.60463, 1825.6  , -5070.1},
                                {71.745953, -3168.0 , 6534.6 },
                                {-46.83552, 1884.2  , -3288.7},
                                {-2.467982, -0.82376, -2.7171},
                                {-0.50272 , -3.1935 , 2.0883 },
                                {8.0956883, 3.7090  , 0.0}
                            }
                            };
};