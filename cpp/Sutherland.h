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
    double gl_x[20]; // Gauss Legendre points for computing barker henderson diamenter
    double gl_w[20]; // Gauss Legendre points for computing barker henderson diamenter
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

    double potential(int i, int j, double r) override;
    double potential_derivative_r(int i, int j, double r) override;
    double potential_dblderivative_rr(int i, int j, double r) override;

    inline vector2d get_sigma_eff(){return sigma_eff;}
    inline vector2d get_sigma_min(){return sigma_min;}
    inline vector2d get_epsilon_eff(){return eps_eff;}
    inline vector2d get_vdw_alpha(){return vdw_alpha;}

    inline vector2d model_rdf(double rho, double T, const std::vector<double>& x) override {
        if (using_LJ_rdf_correlation){
            double rdf = LJ_rdf_correlation(rho, T);
            vector2d rdf_matr(Ncomps, vector1d(Ncomps, rdf));
            return rdf_matr;
        }
        return saft_rdf(rho, T, x, 2);
    }
     // To directly compute the RDF at different pertubation orders. Not used in property computations.
    virtual vector2d saft_rdf(double rho, double T, const std::vector<double>& x, int order=2, bool g2_correction=true);

    vector2d rdf_g0_func(double rho, double T, const vector1d& x);
    vector2d rdf_g1_func(double rho, double T, const vector1d& x);
    vector2d rdf_g2_func(double rho, double T, const vector1d& x, bool g2_correction=true);
    virtual vector2d get_BH_diameters(double T);
    // vector2d get_collision_diameters(double rho, double T, const std::vector<double>& x); // Implemented in Spherical

    double LJ_rdf_correlation(double rho, double T);
    void set_active_LJ_rdf(bool use_LJ_corr){
        using_LJ_rdf_correlation = use_LJ_corr;
    }

    // ------------------------------------------------------------------------------------------------------------------- //
    // -------------------------- Sutherland Internals are below here ---------------------------------------------------- //
    // ---------------------- End users should not need to care about anything below -------------------------------------- //
    // ------------------------------------------------------------------------------------------------------------------- //
    vector2d get_b_max(double T) override;

    protected:
    bool using_LJ_rdf_correlation = false;
    vector2d eps;
    vector3d C;
    vector3d lambda;
    size_t nterms;
    vector2d sigma_eff;
    vector2d sigma_min;
    vector2d eps_eff;
    vector2d vdw_alpha;

    void compute_sigma_eff();
    void compute_epsilon_eff();
    void compute_vdw_alpha();

    vector2d rdf_g0_func(double rho, const vector1d& x, const vector2d& d_BH);
    vector2d rdf_g1_func(double rho, const vector1d& x, const vector2d& d_BH);
    vector2d rdf_g2_func(double rho, double T, const vector1d& x, const vector2d& d_BH, const vector2d& x_eff, bool g2_correction=true);

    vector2d get_lambda_kl(size_t k, size_t l); // lambda_kl[i, j] = lambda_k[i, j] + lambda_l[i, j]
    vector2d get_x0(const vector2d& d_BH); // x0 = sigma / d_BH
    vector2d get_xeff(const vector2d& d_BH); // x_eff = sigma_eff / d_BH

    vector2d gamma_corr(double zeta_x, double T); // Eq. (A37) in svrm (https://doi.org/10.1063/1.4819786)
    double zeta_x_func(double rho, const vector1d& x, const vector2d& d_BH); // Eq. (A13) in svrm
    double dzetax_drho_func(const vector1d& x, const vector2d& d_BH); // Derivative of zeta_x wrt. density
    double zeta_eff_func(double rho, const vector1d& x, double zeta_x, double lambdaijk); // Eq. (A17) in svrm
    double dzeta_eff_drho_func(double rho, const std::vector<double>& x, const vector2d& d_BH, double lambdakij);

    virtual vector2d a1_func(double rho, double T, const vector1d& x);
    vector2d a1_func(double rho, const vector1d& x, const vector2d& d_BH);
    virtual vector2d da1_drho_func(double rho, double T, const vector1d& x){
        const vector2d d_BH = get_BH_diameters(T);
        return da1_drho_func(rho, x, d_BH);
    }
    vector2d da1_drho_func(double rho, const vector1d& x, const vector2d& d_BH); // Derivative of a1 wrt. density
    virtual vector2d a2ij_div_chi_func(double rho, double T, const vector1d& x){
        const vector2d d_BH = get_BH_diameters(T);
        const vector2d x_eff = get_xeff(d_BH);
        const vector2d rdf_chi_HS = rdf_chi_func(rho, x);
        const double zeta_x = zeta_x_func(rho, x, d_BH);
        const double K_HS = K_HS_func(zeta_x);
        const vector2d a2ij = a2ij_func(rho, x, K_HS, rdf_chi_HS, d_BH, x_eff);
        vector2d a2ij_div_chi(Ncomps, vector1d(Ncomps, 0.0));
        std::cout << "Sut : (" << zeta_x << ", " << K_HS << ") : " << a2ij[0][0] << std::endl;
        for (size_t i = 0; i < Ncomps; i++){
            for (size_t j = i; j < Ncomps; j++){
                a2ij_div_chi[i][j] = a2ij[i][j] / (1 + rdf_chi_HS[i][j]);
                a2ij_div_chi[j][i] = a2ij_div_chi[i][j];
            }
        }
        return a2ij_div_chi;
    }
    vector2d a2ij_func(double rho, const vector1d& x, double K_HS, const vector2d& rdf_chi_HS, const vector2d& d_BH, const vector2d& x_eff); // Eq. (A20) in svrm
    vector2d da2ij_div_chi_drho_func(double rho, const vector1d& x, double K_HS, const vector2d& d_BH, const vector2d& x_eff); // Derivative of (a_2 / (1 + chi)) wrt. density
    virtual vector2d da2ij_div_chi_drho_func(double rho, double T, const vector1d& x){
        const vector2d d_BH = get_BH_diameters(T);
        const vector2d x_eff = get_xeff(d_BH);
        const double zeta_x = zeta_x_func(rho, x, d_BH);
        const double K_HS = K_HS_func(zeta_x);
        return da2ij_div_chi_drho_func(rho, x, K_HS, d_BH, x_eff);
    }
    vector2d rdf_chi_func(double rho, const vector1d& x); // Eq. (A22) in svrm
    vector2d drdf_chi_drho_func(double rho, const vector1d& x); // Derivative of chi (from Eq. (22) in svrm)
    vector1d f_corr(double alpha); // Eq. (A26) in svrm

    virtual vector2d a_1s_func(double rho, double T, const vector1d& x, const vector2d& lambda_k); // Forwards call to a_1s_func
    vector2d a_1s_func(double rho, const vector1d& x, double zeta_x, const vector2d& d_BH, const vector2d& lambda_k); // Eq. (A16) in svrm (https://doi.org/10.1063/1.4819786)
    vector2d da1s_drho_func(double rho, const vector1d& x, const vector2d& d_BH, const vector2d& lambda_k); // Derivative of a1s wrt. density

    virtual vector2d B_func(double rho, double T, const vector1d& x, const vector2d& lambda); // Forwards call to B_func
    vector2d B_func(double rho, const vector1d& x, double zeta_x, const vector2d& x_eff, const vector2d& d_BH, const vector2d& lambda_k); // Eq. (A12) in svrm
    vector2d B_func(double rho, const vector1d& x, const vector2d& d_BH, const vector2d& lambda_k); // Forwards call to B_func
    vector2d dBdrho_func(double rho, const vector1d& x, double zeta_x, const vector2d& x_eff, const vector2d& d_BH, const vector2d& lambda_k); // Derivative wrt. density
    vector2d dBdrho_func(double rho, double T, const vector1d& x, const vector2d& lambda_k);
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