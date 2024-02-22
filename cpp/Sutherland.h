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
#include <vector>

using vector1d = std::vector<double>;
using vector2d = std::vector<vector1d>;
using vector3d = std::vector<vector2d>;

struct RDFConstants{
    double gl_x[10]; // Gauss Legendre points for computing barker henderson diamenter
    double gl_w[10]; // Gauss Legendre points for computing barker henderson diamenter
    double C_coeff_matr[4][4]; // See Eq. (A18) in svrm (https://doi.org/10.1063/1.4819786)
    double phi[7][3]; // Table II in svrm (https://doi.org/10.1063/1.4819786)
};

class Sutherland : public Spherical{
    public:
    Sutherland(vector1d mole_weights, vector2d sigma, vector2d eps, vector3d C, vector3d lambda, bool is_idealgas=false)
        : Spherical(mole_weights, sigma, is_idealgas), eps{eps}, C{C}, lambda{lambda}, nterms{C.size()}, 
        sigma_eff{sigma}, sigma_min{sigma}, eps_eff(Ncomps, vector1d(Ncomps, 0.)),
        vdw_alpha(Ncomps, vector1d(Ncomps, 0.))
        {compute_sigma_eff(); compute_epsilon_eff(); compute_vdw_alpha();}

    Sutherland(vector1d mole_weights, vector2d sigma, vector2d eps, size_t nterms, bool is_idealgas=false)
        : Spherical(mole_weights, sigma, is_idealgas), eps{eps}, nterms{nterms}, sigma_eff(Ncomps, vector1d(Ncomps, 0.)),
        sigma_min(Ncomps, vector1d(Ncomps, 0.)), eps_eff(Ncomps, vector1d(Ncomps, 0.)), vdw_alpha(Ncomps, vector1d(Ncomps, 0.))
        {
        C = vector3d(nterms, vector2d(Ncomps, vector1d(Ncomps, 0.)));
        lambda = vector3d(nterms, vector2d(Ncomps, vector1d(Ncomps, 0.)));
        }

    virtual double potential(int i, int j, double r) override;
    virtual double potential_derivative_r(int i, int j, double r) override;
    virtual double potential_dblderivative_rr(int i, int j, double r) override;

    inline vector2d get_sigma_eff(){return sigma_eff;} // Tested vs. Mie : OK
    inline vector2d get_sigma_min(){return sigma_min;} // Tested vs. Mie : OK
    inline vector2d get_epsilon_eff(){return eps_eff;} // Tested vs. Mie : OK
    inline vector2d get_vdw_alpha(){return vdw_alpha;} // Tested vs. Mie : OK

    virtual vector2d model_rdf(double rho, double T, const vector1d& mole_fracs) override; // Tested vs. Mie : OK
    vector2d rdf_g0_func(double rho, double T, const vector1d& x); // Tested vs. Mie : OK
    vector2d rdf_g1_func(double rho, double T, const vector1d& x); // Tested vs. Mie : OK
    vector2d rdf_g2_func(double rho, double T, const vector1d& x); // Tested vs. Mie : OK
    virtual vector2d get_BH_diameters(double T); // Tested vs. Mie : OK
    // vector2d get_collision_diameters(double rho, double T, const std::vector<double>& x); // Implemented in Spherical

    // ------------------------------------------------------------------------------------------------------------------- //
    // -------------------------- Sutherland Internals are below here ---------------------------------------------------- //
    // ---------------------- End users should not need to care about anything below -------------------------------------- //
    // ------------------------------------------------------------------------------------------------------------------- //

    vector2d eps;
    vector3d C;
    vector3d lambda;
    size_t nterms;
    vector2d sigma_eff;
    vector2d sigma_min;
    vector2d eps_eff;
    vector2d vdw_alpha;

    virtual vector2d get_b_max(double T) override;

    void compute_sigma_eff(); // Tested vs. Mie : OK
    void compute_epsilon_eff(); // Tested vs. Mie : OK
    void compute_vdw_alpha(); // Tested vs. Mie : OK

    vector2d rdf_g0_func(double rho, const vector1d& x, const vector2d& d_BH); // Tested vs. Mie : OK
    vector2d rdf_g1_func(double rho, const vector1d& x, const vector2d& d_BH); // Tested vs. Mie : OK
    vector2d rdf_g2_func(double rho, double T, const vector1d& x, const vector2d& d_BH, const vector2d& x_eff); // Tested vs. Mie : OK

    vector2d get_lambda_kl(size_t k, size_t l); // lambda_kl[i, j] = lambda_k[i, j] + lambda_l[i, j]
    vector2d get_x0(const vector2d& d_BH); // x0 = sigma / d_BH
    vector2d get_xeff(const vector2d& d_BH); // x_eff = sigma_eff / d_BH

    vector2d gamma_corr(double zeta_x, double T); // Eq. (A37) in svrm (https://doi.org/10.1063/1.4819786)
    double zeta_x_func(double rho, const vector1d& x, const vector2d& d_BH); // Eq. (A13) in svrm
    double dzetax_drho_func(const vector1d& x, const vector2d& d_BH); // Derivative of zeta_x wrt. density
    double zeta_eff_func(double rho, const vector1d& x, double zeta_x, double lambdaijk); // Eq. (A17) in svrm
    double dzeta_eff_drho_func(double rho, const std::vector<double>& x, const vector2d& d_BH, double lambdaijk); // Tested vs. Mie : OK

    vector2d da1_drho_func(double rho, const vector1d& x, const vector2d& d_BH); // Derivative of a1 wrt. density
    vector2d a2ij_func(double rho, const vector1d& x, double K_HS, const vector2d& rdf_chi_HS, const vector2d& d_BH, const vector2d& x_eff); // Eq. (A20) in svrm
    vector2d da2ij_div_chi_drho_func(double rho, const vector1d& x, double K_HS, const vector2d& d_BH, const vector2d& x_eff); // Derivative of (a_2 / (1 + chi)) wrt. density
    vector2d rdf_chi_func(double rho, const vector1d& x, const vector2d& d_BH); // Eq. (A22) in svrm
    vector2d drdf_chi_drho_func(double rho, const vector1d& x, const vector2d& d_BH); // Derivative of chi (from Eq. (22) in svrm)
    vector1d f_corr(double alpha); // Eq. (A26) in svrm

    vector2d a_1s_func(double rho, const vector1d& x, double zeta_x, const vector2d& d_BH, const vector2d& lambda_k); // Eq. (A16) in svrm (https://doi.org/10.1063/1.4819786)
    vector2d a_1s_func(double rho, double T, const vector1d& x, const vector2d& lambda_k); // Forwards call to a_1s_func
    vector2d da1s_drho_func(double rho, const vector1d& x, const vector2d& d_BH, const vector2d& lambda_k); // Derivative of a1s wrt. density
    
    vector2d B_func(double rho, const vector1d& x, double zeta_x, const vector2d& x_eff, const vector2d& d_BH, const vector2d& lambda_k); // Eq. (A12) in svrm
    vector2d B_func(double rho, const vector1d& x, const vector2d& d_BH, const vector2d& lambda_k); // Forwards call to B_func
    vector2d dBdrho_func(double rho, const vector1d& x, double zeta_x, const vector2d& x_eff, const vector2d& d_BH, const vector2d& lambda_k); // Derivative wrt. density
    vector2d I_func(const vector2d& xeff, const vector2d& lambda_k); // Eq. (A14) in svrm (https://doi.org/10.1063/1.4819786)
    vector2d J_func(const vector2d& xeff, const vector2d& lambda_k); // Eq. (A15) in svrm (https://doi.org/10.1063/1.4819786)

    inline double K_HS_func(double zeta_x){ // Eq. (A21) in svrm
        return pow(1 - zeta_x, 4) / (1 + 4 * zeta_x + 4 * pow(zeta_x, 2) - 4 * pow(zeta_x, 3) + pow(zeta_x, 4));
    }
    inline double dKHS_drho_func(double zeta_x, double dzx_drho){ // derivative of K_HS wrt. density.
        return - 4 * dzx_drho * pow(1 - zeta_x, 3) * (2 + 5 * zeta_x - pow(zeta_x, 2) - 2 * pow(zeta_x, 3))
                / pow(1 + 4 * zeta_x + 4 * pow(zeta_x, 2) - 4 * pow(zeta_x, 3) - pow(zeta_x, 4), 2);
    }

    static constexpr RDFConstants rdf_constants{
                            { // gl_x (see: ThermoPack, SAFT-VR Mie memo)
                            -0.973906528517171720078, -0.8650633666889845107321,
                            -0.6794095682990244062343, -0.4333953941292471907993,
                            -0.1488743389816312108848,  0.1488743389816312108848,
                            0.4333953941292471907993,  0.6794095682990244062343,
                            0.8650633666889845107321,  0.973906528517171720078
                            },
                            { // gl_w (see: ThermoPack, SAFT-VR Mie memo)
                            0.0666713443086881375936, 0.149451349150580593146,
                            0.219086362515982043996, 0.2692667193099963550912,
                            0.2955242247147528701739, 0.295524224714752870174,
                            0.269266719309996355091, 0.2190863625159820439955,
                            0.1494513491505805931458, 0.0666713443086881375936
                            },
                            { // C_coeff_matr (See: Eq. (A18) in svrm)
                                {0.81096, 1.7888, -37.578, 92.284},
                                {1.0205, -19.341, 151.26, -463.50},
                                {-1.9057, 22.845, -228.14, 973.92},
                                {1.0885, -6.1962, 106.98, -677.64}
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



namespace noname_go_away{

// Gauss Legendre points for computing barker henderson diamenter (see: ThermoPack, SAFT-VR-Mie docs)
constexpr double gl_x[10] = {-0.973906528517171720078, -0.8650633666889845107321,
                            -0.6794095682990244062343, -0.4333953941292471907993,
                            -0.1488743389816312108848,  0.1488743389816312108848,
                            0.4333953941292471907993,  0.6794095682990244062343,
                            0.8650633666889845107321,  0.973906528517171720078};
constexpr double gl_w[10] = {0.0666713443086881375936, 0.149451349150580593146,
                            0.219086362515982043996, 0.2692667193099963550912,
                            0.2955242247147528701739, 0.295524224714752870174,
                            0.269266719309996355091, 0.2190863625159820439955,
                            0.1494513491505805931458, 0.0666713443086881375936};

constexpr double C_coeff_matr[4][4] // See Eq. A18 of J. Chem. Phys. 139, 154504 (2013); https://doi.org/10.1063/1.4819786
    {
        {0.81096, 1.7888, -37.578, 92.284},
        {1.0205, -19.341, 151.26, -463.50},
        {-1.9057, 22.845, -228.14, 973.92},
        {1.0885, -6.1962, 106.98, -677.64}
    };

constexpr double phi[7][3] // J. Chem. Phys. 139, 154504 (2013); https://doi.org/10.1063/1.4819786, Table II
        {
            {7.5365557, -359.44 , 1550.9 },
            {-37.60463, 1825.6  , -5070.1},
            {71.745953, -3168.0 , 6534.6 },
            {-46.83552, 1884.2  , -3288.7},
            {-2.467982, -0.82376, -2.7171},
            {-0.50272 , -3.1935 , 2.0883 },
            {8.0956883, 3.7090  , 0.0}
        };

};