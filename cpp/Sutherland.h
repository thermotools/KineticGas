#pragma once
#include "KineticGas.h"
#include "Spherical.h"
#include <vector>

using vector1d = std::vector<double>;
using vector2d = std::vector<vector1d>;
using vector3d = std::vector<vector2d>;


class Sutherland : public Spherical{
    private:
    vector2d eps;
    vector3d C;
    vector3d lambda;
    size_t nterms;
    vector2d sigma_min;
    vector2d sigma_eff;
    vector2d eps_eff;

    public:
    Sutherland(vector1d mole_weights, vector2d sigma, vector2d eps, vector3d C, vector3d lambda, bool is_idealgas=false)
        : Spherical(mole_weights, sigma, is_idealgas), eps{eps}, C{C}, lambda{lambda}, nterms{C[0][0].size()}
        {compute_sigma_eff(), compute_epsilon_eff();}
    
    virtual double potential(int i, int j, double r) override;
    virtual double potential_derivative_r(int i, int j, double r) override;
    virtual double potential_dblderivative_rr(int i, int j, double r) override;
    virtual vector2d get_contact_diameters(double rho, double T, const vector1d& x) override;

    void compute_sigma_eff();
    void compute_epsilon_eff();
    vector2d get_sigma_min();
    vector2d get_sigma_eff();
    vector2d get_epsilon_eff();

    virtual vector2d model_rdf(double rho, double T, const vector1d& mole_fracs) override;

    vector2d rdf_g0_func(double rho, const vector1d& x, const vector2d& d_BH);
    vector2d rdf_g1_func(double rho, const vector1d& x, const vector2d& d_BH);
    vector2d rdf_g1_func(double rho, double T, const vector1d& x);
    vector2d rdf_g2_func(double rho, double T, const vector1d& x, const vector2d& d_BH, const vector2d& x0);
    vector2d rdf_g2_func(double rho, double T, const vector1d& x);
    
    virtual vector2d get_BH_diameters(double T);
    vector2d get_x0(const vector2d& d_BH);
    vector2d get_xeff(const vector2d& d_BH);

    private:
    vector2d a_1s_func(double rho, const vector1d& x, double zeta_x, const vector2d& d_BH, const vector2d& lambda_k);
    vector2d a_1s_func(double rho, double T, const vector1d& x, const vector2d& lambda_k);
    vector2d da1s_drho_func(double rho, const vector1d& x, const vector2d& d_BH, const vector2d& lambda_k);
    vector2d da1_drho_func(double rho, const vector1d& x, const vector2d& d_BH);
    vector2d B_func(double rho, const std::vector<double>& x, double zeta_x, const vector2d& x_eff, const vector2d& d_BH, const vector2d& lambda_k);
    vector2d dBdrho_func(double rho, const vector1d& x, double zeta_x, const vector2d& x_eff, const vector2d& d_BH, const vector2d& lambda_k);
    vector2d I_func(const vector2d& xeff, const vector2d& lambda_k);
    vector2d J_func(const vector2d& xeff, const vector2d& lambda_k);

    double zeta_x_func(double rho, const vector1d& x, const vector2d& d_BH);
    double dzetax_drho_func(const vector1d& x, const vector2d& d_BH);
    double zeta_eff_func(double rho,  const vector1d& x, double zeta_x, double lambdaijk);
    double dzeta_eff_drho_func(double rho, const std::vector<double>& x, const vector2d& d_BH, double lambdaijk);
};

namespace sutherland_rdf{

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