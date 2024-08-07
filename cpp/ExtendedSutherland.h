/*
References:
    svrm:
        Accurate statistical associating fluid theory for chain molecules formed from Mie segments
        Thomas Lafitte; Anastasia Apostolakou; Carlos Avendaño; Amparo Galindo; Claire S. Adjiman; Erich A. Müller; George Jackson
        J. Chem. Phys. 139, 154504 (2013)
        https://doi.org/10.1063/1.4819786
*/
#pragma once
#include "KineticGas.h"
#include "Spherical.h"
#include "Sutherland.h"
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

class ExtSutherland : public Spherical{
    public:
    ExtSutherland(vector1d mole_weights, vector2d sigma, vector2d eps, 
                vector3d C, vector3d lambda, vector3d beta_exp, vector3d rho_exp, 
                bool is_idealgas=false, bool is_singlecomp=false)
        : Spherical(mole_weights, sigma, eps, is_idealgas, is_singlecomp), C{C}, 
        lambda{lambda}, beta_exp(beta_exp), rho_exp{rho_exp}, nterms{C.size()}
        {}

    // ExtSutherland(vector1d mole_weights, vector2d sigma, vector2d eps, size_t nterms, bool is_idealgas=false, bool is_singlecomp=false)
    //     : Spherical(mole_weights, sigma, eps, is_idealgas, is_singlecomp), nterms{nterms}
    //     {
    //     C = vector3d(nterms, vector2d(Ncomps, vector1d(Ncomps, 0.)));
    //     lambda = vector3d(nterms, vector2d(Ncomps, vector1d(Ncomps, 0.)));
    //     beta_exp = vector3d(nterms, vector2d(Ncomps, vector1d(Ncomps, 0.)));
    //     rho_exp = vector3d(nterms, vector2d(Ncomps, vector1d(Ncomps, 0.)));
    //     }

    dual2 potential(int i, int j, dual2 r, dual2 T, dual2 rho);
    dual2 potential_r(int i, int j, dual2 r, dual2 T, dual2 rho);
    dual2 potential_rr(int i, int j, dual2 r, dual2 T, dual2 rho);
    double potential(int i, int j, double r, double T, double rho);

    inline vector2d model_rdf(double rho, double T, const std::vector<double>& x) override {
        return saft_rdf(rho, T, x, 2);
    }
     // To directly compute the RDF at different pertubation orders. Not used in property computations.
    vector2d saft_rdf(double rho, double T, const std::vector<double>& x, int order=2, bool g2_correction=true);

    vector2d2 rdf_g0_func(dual2 rho, const vector1d& x, const vector2d2& d_BH, const vector2d2& sigma_eff);
    vector2d1 rdf_g1_func(dual2 rho, const vector1d& x, const vector2d2& d_BH, const vector2d2& x_eff);
    vector2d1 rdf_g2_func(dual2 rho, dual2 T, const vector1d& x, const vector2d2& d_BH, const vector2d2& x_eff, const vector2d2& sigma_eff, bool g2_correction);
    virtual vector2d2 get_BH_diameters(dual2 rho, dual2 T);

    vector3d get_rdf_terms(double rho, double T, const vector1d& x); // Return g0, g1, g2 (no correction), g2 (with correction)

    // ------------------------------------------------------------------------------------------------------------------- //
    // -------------------------- Sutherland Internals are below here ---------------------------------------------------- //
    // ---------------------- End users should not need to care about anything below -------------------------------------- //
    // ------------------------------------------------------------------------------------------------------------------- //
    vector2d get_b_max(double rho, double T);
    vector2d get_b_max(double T){throw_notimplemented(); return vector2d();}

    protected:
    vector3d C;
    vector3d lambda;
    vector3d beta_exp;
    vector3d rho_exp;
    size_t nterms;

    dual2 potential(int i, int j, dual2 r) override {
        throw_notimplemented(); return 0.;
    }
    using Spherical::potential_derivative_r;
    using Spherical::potential_dblderivative_rr;

    vector2d2 compute_sigma_eff(dual2 rho, dual2 T);
    vector2d2 compute_epsilon_eff(dual2 rho, dual2 T);
    vector2d2 compute_vdw_alpha(dual2 rho, dual2 T);

    vector2d rdf_g0_func(double rho, const vector1d& x, const vector2d& d_BH);
    vector2d rdf_g1_func(double rho, const vector1d& x, const vector2d& d_BH);
    vector2d rdf_g2_func(double rho, double T, const vector1d& x, const vector2d& d_BH, const vector2d& x_eff, bool g2_correction=true);

    vector2d get_lambda_kl(size_t k, size_t l); // lambda_kl[i, j] = lambda_k[i, j] + lambda_l[i, j]
    vector2d2 get_x0(const vector2d2& d_BH); // x0 = sigma / d_BH
    vector2d2 get_xeff(const vector2d2& d_BH, const vector2d2& sigma_eff); // x_eff = sigma_eff / d_BH

    vector2d2 gamma_corr(dual2 zeta_x, dual2 T, const vector2d2& eps_eff, const vector2d2& vdw_alpha); // Eq. (A37) in svrm (https://doi.org/10.1063/1.4819786)
    dual2 zeta_x_func(dual2 rho, const vector1d& x, const vector2d2& d_BH); // Eq. (A13) in svrm
    // double dzetax_drho_func(const vector1d& x, const vector2d& d_BH); // Derivative of zeta_x wrt. density
    dual2 zeta_eff_func(dual2 rho, const vector1d& x, dual2 zeta_x, double lambdaijk); // Eq. (A17) in svrm
    // double dzeta_eff_drho_func(double rho, const std::vector<double>& x, const vector2d& d_BH, double lambdakij);

    dual2 a1_func(int i, int j, dual2 rho, const vector1d& x, const vector2d2& d_BH, const vector2d2& sigma_eff);
    // virtual vector2d2 a1_func(dual2 rho, dual2 T, const vector1d& x);
    // vector2d2 a1_func(dual2 rho, const vector1d& x, const vector2d2& d_BH);
    // virtual vector2d2 da1_drho_func(dual2 rho, dual2 T, const vector1d& x){
    //     const vector2d2 d_BH = get_BH_diameters(T, rho);
    //     return da1_drho_func(rho, x, d_BH);
    // }
    // vector2d da1_drho_func(double rho, const vector1d& x, const vector2d& d_BH); // Derivative of a1 wrt. density
    // virtual vector2d2 a2ij_div_chi_func(dual2 rho, dual2 T, const vector1d& x){
    //     const vector2d2 d_BH = get_BH_diameters(T, rho);
    //     const vector2d2 sigma_eff = compute_sigma_eff(T, rho);
    //     const vector2d2 x_eff = get_xeff(d_BH, sigma_eff);
    //     const vector2d2 rdf_chi_HS = rdf_chi_func(rho, x);
    //     const dual2 zeta_x = zeta_x_func(rho, x, d_BH);
    //     const dual2 K_HS = K_HS_func(zeta_x);
    //     const vector2d2 a2ij = a2ij_func(rho, x, K_HS, rdf_chi_HS, d_BH, x_eff);
    //     vector2d2 a2ij_div_chi(Ncomps, vector1d2(Ncomps, 0.0));
    //     for (size_t i = 0; i < Ncomps; i++){
    //         for (size_t j = i; j < Ncomps; j++){
    //             a2ij_div_chi[i][j] = a2ij[i][j] / (1 + rdf_chi_HS[i][j]);
    //             a2ij_div_chi[j][i] = a2ij_div_chi[i][j];
    //         }
    //     }
    //     return a2ij_div_chi;
    // }
    dual2 a2_div_chi_func(int i, int j, dual2 rho, const vector1d& x, dual2 K_HS, const vector2d2& d_BH, const vector2d2& x_eff); // Eq. (A20) in svrm
    // vector2d da2ij_div_chi_drho_func(double rho, const vector1d& x, double K_HS, const vector2d& d_BH, const vector2d& x_eff); // Derivative of (a_2 / (1 + chi)) wrt. density
    // virtual vector2d2 da2ij_div_chi_drho_func(dual2 rho, dual2 T, const vector1d& x){
    //     const vector2d2 d_BH = get_BH_diameters(T, rho);
    //     const vector2d2 sigma_eff = compute_sigma_eff(T, rho);
    //     const vector2d2 x_eff = get_xeff(d_BH, sigma_eff);
    //     const dual2 zeta_x = zeta_x_func(rho, x, d_BH);
    //     const dual2 K_HS = K_HS_func(zeta_x);
    //     return da2ij_div_chi_drho_func(rho, x, K_HS, d_BH, x_eff);
    // }
    vector2d2 rdf_chi_func(dual2 rho, dual2 T, const vector1d& x, const vector2d2& sigma_eff); // Eq. (A22) in svrm
    // vector2d drdf_chi_drho_func(double rho, const vector1d& x); // Derivative of chi (from Eq. (22) in svrm)
    vector1d2 f_corr(dual2 alpha); // Eq. (A26) in svrm

    // virtual dual2 a_1s_func(dual2 rho, dual2 T, const vector1d& x, const vector2d& lambda_k); // Forwards call to a_1s_func
    dual2 a_1s_func(int i, int j, dual2 rho, const vector1d& x, dual2 zeta_x, const vector2d2& d_BH, const vector2d& lambda_k); // Eq. (A16) in svrm (https://doi.org/10.1063/1.4819786)
    // vector2d da1s_drho_func(double rho, const vector1d& x, const vector2d& d_BH, const vector2d& lambda_k); // Derivative of a1s wrt. density

    // virtual vector2d2 B_func(dual2 rho, dual2 T, const vector1d& x, const vector2d& lambda); // Forwards call to B_func
    dual2 B_func(int i, int j, dual2 rho, const vector1d& x, dual2 zeta_x, const vector2d2& x_eff, const vector2d2& d_BH, const vector2d& lambda_k); // Eq. (A12) in svrm
    // vector2d2 B_func(dual2 rho, const vector1d& x, const vector2d2& d_BH, const vector2d& lambda_k); // Forwards call to B_func
    // vector2d dBdrho_func(double rho, const vector1d& x, double zeta_x, const vector2d& x_eff, const vector2d& d_BH, const vector2d& lambda_k); // Derivative wrt. density
    // vector2d dBdrho_func(double rho, double T, const vector1d& x, const vector2d& lambda_k);
    dual2 I_func(int i, int j, const vector2d2& xeff, const vector2d& lambda_k); // Eq. (A14) in svrm (https://doi.org/10.1063/1.4819786)
    dual2 J_func(int i, int j, const vector2d2& xeff, const vector2d& lambda_k); // Eq. (A15) in svrm (https://doi.org/10.1063/1.4819786)

    inline dual2 K_HS_func(dual2 zeta_x){ // Eq. (A21) in svrm
        return pow(1 - zeta_x, 4) / (1 + 4 * zeta_x + 4 * pow(zeta_x, 2) - 4 * pow(zeta_x, 3) + pow(zeta_x, 4));
    }
    // inline double dKHS_drho_func(double zeta_x, double dzx_drho){ // derivative of K_HS wrt. density.
    //     return - 4 * dzx_drho * pow(1 - zeta_x, 3) * (2 + 5 * zeta_x - pow(zeta_x, 2) - 2 * pow(zeta_x, 3))
    //             / pow(1 + 4 * zeta_x + 4 * pow(zeta_x, 2) - 4 * pow(zeta_x, 3) - pow(zeta_x, 4), 2);
    // }

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