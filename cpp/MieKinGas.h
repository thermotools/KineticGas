/*
Author: Vegard Gjeldvik Jervell
Contains: The MieKinGas class. This class is the model used to evaluate the Enskog solutions for a Mie potential.
            MieKinGas overrides the potential and potential derivative functions in Spherical.

            Also contains functions required to compute the radial distribution function at contact for a Mie-fluid
            and related constant factors that have been previously regressed (namespace mie_rdf_constants).
            See : J. Chem. Phys. 139, 154504 (2013); https://doi.org/10.1063/1.4819786
*/

#pragma once
#include "Spherical.h"

class MieKinGas : public Spherical {
    public:
    vector2d la, lr, C, alpha;

    MieKinGas(vector1d mole_weights, vector2d sigma, vector2d eps, vector2d la, vector2d lr, bool is_idealgas, bool is_singlecomp);
    MieKinGas(std::string comps, bool is_idealgas=false);

    void set_C_alpha();
    void mix_sigma();
    void mix_epsilon();
    void mix_exponents(vector2d& expo);

    dual2 potential(int i, int j, dual2 r) override {
        return C[i][j] * eps[i][j] * (pow(sigma[i][j] / r, lr[i][j]) - pow(sigma[i][j] / r, la[i][j]));
    }

    double potential(int i, int j, double r) override {
        return C[i][j] * eps[i][j] * (pow(sigma[i][j] / r, lr[i][j]) - pow(sigma[i][j] / r, la[i][j]));
    }
    
    double potential_derivative_r(int i, int j, double r) override {
        return C[i][j] * eps[i][j] * ((la[i][j] * pow(sigma[i][j] / r, la[i][j]) / r)
                                        - (lr[i][j] * pow(sigma[i][j] / r, lr[i][j]) / r));
    }

    double potential_dblderivative_rr(int i, int j, double r) override {
        return C[i][j] * eps[i][j] * ((lr[i][j] * (lr[i][j] + 1) * pow(sigma[i][j] / r, lr[i][j]) / pow(r, 2))
                                    - (la[i][j] * (la[i][j] + 1) * pow(sigma[i][j] / r, la[i][j]) / pow(r, 2)));
    }


    double omega(int i, int j, int l, int r, double T) override;
    double omega_correlation(int i, int j, int l, int r, double T_star);
    double omega_recursive_factor(int i, int j, int l, int r, double T);
    // The hard sphere integrals are used as the reducing factor for the correlations.
    // So we need to compute the hard-sphere integrals to convert the reduced collision integrals from the
    // Correlation by Fokin et. al. to the "real" collision integrals.
    inline double omega_hs(int i, int j, int l, int r, double T){
        double w = PI * pow(sigma[i][j], 2) * 0.5 * (r + 1);
        for (int ri = r; ri > 1; ri--) {w *= ri;}
        if (l % 2 == 0){
            w *= (1. - (1.0 / (l + 1.)));
        }
        if (i == j) return sqrt((BOLTZMANN * T) / (PI * m[i])) * w;
        return sqrt(BOLTZMANN * T * (m[i] + m[j]) / (2. * PI * m[i] * m[j])) * w;
    }

    // Contact diameter related methods
    // bmax[i][j] is in units of sigma[i][j]
    // bmax = The maximum value of the impact parameter at which deflection angle (chi) is positive
    std::vector<std::vector<double>> get_b_max(double T) override;
    virtual std::vector<std::vector<double>> get_BH_diameters(double T);
    std::vector<std::vector<double>> get_vdw_alpha(){return alpha;}

    // Methods for computing the radial distribution function at contact
    // Note: A lot of these methods have two overloads: One that takes the temperature and density, and comptes
    //      the BH diameter, packing fractions etc, and another that takes the BH diameters, packing fractions,
    //      and other variables that must be pre-computed directly. I'm not sure which version of these is in primary use
    //      In the "standard" call chain when `model_rdf` is called, but someone should at some point ensure that we
    //      are not computing a bunch of unnessecary BH diameters.
    // Note: For SAFT-type models, we
    inline std::vector<std::vector<double>> model_rdf(double rho, double T, const std::vector<double>& x) override {
        return saft_rdf(rho, T, x, 2);
    }

     // To directly compute the RDF at different pertubation orders. Not used in property computations.
    std::vector<std::vector<double>> saft_rdf(double rho, double T, const std::vector<double>& x, int order=2, bool g2_correction=true);

    std::vector<std::vector<double>> rdf_HS(double rho, const std::vector<double>& x,
                                            const std::vector<std::vector<double>>& d_BH);
    std::vector<std::vector<double>> rdf_HS(double rho, double T, const std::vector<double>& x);
    std::vector<std::vector<double>> rdf_g1_func(double rho, const std::vector<double>& x,
                                                const std::vector<std::vector<double>>& d_BH);
    std::vector<std::vector<double>> rdf_g1_func(double rho, double T, const std::vector<double>& x);
    std::vector<std::vector<double>> rdf_g2_func(double rho, double T, const std::vector<double>& x,
                                                const std::vector<std::vector<double>>& d_BH,
                                                const std::vector<std::vector<double>>& x0,
                                                bool g2_correction=true);
    std::vector<std::vector<double>> rdf_g2_func(double rho, double T, const std::vector<double>& x, bool g2_correction=true);
    std::vector<std::vector<double>> get_x0(const std::vector<std::vector<double>>& d_BH);
    std::vector<std::vector<double>> a_1s_func(double rho, const std::vector<double>& x,
                                                const std::vector<std::vector<double>>& d_BH,
                                                const std::vector<std::vector<double>>& lambda);
    std::vector<std::vector<double>> a_1s_func(double rho, double T, const std::vector<double>& x,
                                                const std::vector<std::vector<double>>& lambda);
    std::vector<std::vector<double>> da1s_drho_func(double rho, const std::vector<double>& x,
                                                    const std::vector<std::vector<double>>& d_BH,
                                                    const std::vector<std::vector<double>>& lambda);
    std::vector<std::vector<double>> da1s_drho_func(double rho, double T, const std::vector<double>& x,
                                                const std::vector<std::vector<double>>& lambda);

    std::vector<std::vector<double>> I_func(const std::vector<std::vector<double>>& x0,
                                            const std::vector<std::vector<double>>& lambda);
    std::vector<std::vector<double>> J_func(const std::vector<std::vector<double>>& x0,
                                            const std::vector<std::vector<double>>& lambda);
    std::vector<std::vector<double>> B_func(double rho, const std::vector<double>& x,
                                            const std::vector<std::vector<double>>& d_BH,
                                            const std::vector<std::vector<double>>& lambda);
    std::vector<std::vector<double>> B_func(double rho, double T, const std::vector<double>& x,
                                                    const std::vector<std::vector<double>>& lambda);
    std::vector<std::vector<double>> dBdrho_func(double rho, const std::vector<double>& x,
                                                const std::vector<std::vector<double>>& d_BH,
                                                const std::vector<std::vector<double>>& lambda);
    std::vector<std::vector<double>> dBdrho_func(double rho, double T, const std::vector<double>& x,
                                                    const std::vector<std::vector<double>>& lambda);

    std::vector<std::vector<double>> a1ij_func(double rho, double T, const std::vector<double>& x);

    std::vector<std::vector<double>> da1ij_drho_func(double rho, const std::vector<double>& x,
                                                    const std::vector<std::vector<double>>& d_BH);
    std::vector<std::vector<double>> da1ij_drho_func(double rho, double T, const std::vector<double>& x);                                                
    
    std::vector<std::vector<double>> a2ij_func(double rho, const std::vector<double>& x, double K_HS,
                                                const std::vector<std::vector<double>>& rdf_chi,
                                                const std::vector<std::vector<double>>& d_BH,
                                                const std::vector<std::vector<double>>& x0);
    std::vector<std::vector<double>> a2ij_func(double rho, double T, const std::vector<double>& x);
    std::vector<std::vector<double>> a2ij_div_chi_func(double rho, double T, const std::vector<double>& x){
        std::vector<std::vector<double>> a2ij = a2ij_func(rho, T, x);
        std::vector<std::vector<double>> d_BH = get_BH_diameters(T);
        std::vector<std::vector<double>> rdf_chi = rdf_chi_func(rho, x);
        std::vector<std::vector<double>> a2ij_div_chi(Ncomps, std::vector<double>(Ncomps, 0.0));
        for (size_t i = 0; i < Ncomps; i++){
            for (size_t j = i; j < Ncomps; j++){
                a2ij_div_chi[i][j] = a2ij[i][j] / (1 + rdf_chi[i][j]);
                a2ij_div_chi[j][i] = a2ij_div_chi[i][j];
            }
        }
        return a2ij_div_chi;
    }

    std::vector<std::vector<double>> da2ij_drho_func(double rho, const std::vector<double>& x, double K_HS, 
                                                    const std::vector<std::vector<double>>& d_BH,
                                                    const std::vector<std::vector<double>>& x0);
    std::vector<std::vector<double>> da2ij_drho_func(double rho, double T, const std::vector<double>& x);

    std::vector<std::vector<double>> da2ij_div_chi_drho_func(double rho, const std::vector<double>& x, double K_HS, 
                                                const std::vector<std::vector<double>>& d_BH,
                                                const std::vector<std::vector<double>>& x0);
    std::vector<std::vector<double>> da2ij_div_chi_drho_func(double rho, double T, const std::vector<double>& x){
        std::vector<std::vector<double>> d_BH = get_BH_diameters(T);
        double zeta_x = zeta_x_func(rho, x, d_BH);
        double K_HS = K_HS_func(zeta_x);
        std::vector<std::vector<double>> x0 = get_x0(d_BH);
        return da2ij_div_chi_drho_func(rho, x, K_HS, d_BH, x0);
    }

    std::vector<std::vector<double>> rdf_chi_func(double rho, const std::vector<double>& x);
    std::vector<std::vector<double>> drdf_chi_drho_func(double rho, const std::vector<double>& x);
    
    std::vector<std::vector<double>> gamma_corr(double zeta_x, double T);
    std::vector<std::vector<double>> gamma_corr(double rho, double T, const std::vector<double>& x){
        std::vector<std::vector<double>> d_BH = get_BH_diameters(T);
        double zeta_x = zeta_x_func(rho, x, d_BH);
        return gamma_corr(zeta_x, T);
    }

    inline double K_HS_func(double zeta_x){
        return pow(1 - zeta_x, 4) / (1 + 4 * zeta_x + 4 * pow(zeta_x, 2) - 4 * pow(zeta_x, 3) + pow(zeta_x, 4));
    }
    inline double dKHS_drho_func(double zeta_x, double dzx_drho){
        return - 4 * dzx_drho * pow(1 - zeta_x, 3) * (2 + 5 * zeta_x - pow(zeta_x, 2) - 2 * pow(zeta_x, 3))
                / pow(1 + 4 * zeta_x + 4 * pow(zeta_x, 2) - 4 * pow(zeta_x, 3) - pow(zeta_x, 4), 2);
    }
    inline double eta_func(double rho, double d_BH){return rho * PI * pow(d_BH, 3.0) / 6.0;}
    double zeta_x_func(double rho,
                    const std::vector<double>& x,
                    const std::vector<std::vector<double>>& d_BH);
    double dzetax_drho_func(const std::vector<double>& x,
                    const std::vector<std::vector<double>>& d_BH);

    double zeta_eff_func(double rho,
                    const std::vector<double>& x,
                    const std::vector<std::vector<double>>& d_BH,
                    double lambdaij);

    double dzeta_eff_drho_func(double rho,
                    const std::vector<double>& x,
                    const std::vector<std::vector<double>>& d_BH,
                    double lambdaij);

    std::vector<double> f_corr(double alpha);

    static constexpr double omega_correlation_factors[2][6][4] =
    {
        {
          {0., -0.145269e1, 0.294682e2, 0.242508e1},
          {0.107782e-1, 0.587725, -0.180714e3, .595694e2},
          {0.546646e-1, -0.651465e1, 0.374457e3, -0.137807e3},
          {0.485352, 0.245523e2, -0.336782e3, 0.814187e2},
          {-0.385355, -0.206868e2, 0.132246e3, 0.},
          {0.847232e-1, 0.521812e1, -0.181140e2, -0.747215e1}
        },
        {
          {0., 0.113086e1, 0.234799e2, 0.310127e1},
          {0., 0.551559e1, -0.137023e3, 0.185848e2},
          {0.325909e-1, -0.292925e2, 0.243761e3, 0.},
          {0.697682, 0.590792e2, -0.143670e3, -0.123518e3},
          {-0.564238, -0.430549e2, 0., 0.137282e3},
          {0.126508, 0.104273e2, 0.150601e2, -0.408911e2}
        }
    };
    
};

namespace mie_rdf_constants{

// Gauss Legendre points for computing barker henderson diamenter (see: ThermoPack, SAFT-VR-Mie docs)
constexpr double gl_x[20] = {-0.9931285991850949, -0.9639719272779139,
                            -0.912234428251326, -0.8391169718222187,
                            -0.7463319064601508, -0.636053680726515,
                            -0.5108670019508271, -0.37370608871541955,
                            -0.22778585114164507, -0.07652652113349737,
                            0.07652652113349737, 0.22778585114164507,
                            0.37370608871541955, 0.5108670019508271,
                            0.636053680726515, 0.7463319064601508,
                            0.8391169718222187, 0.912234428251326,
                            0.9639719272779139, 0.9931285991850949};
constexpr double gl_w[20] = {0.017614007139152742, 0.040601429800386134,
                            0.06267204833410799, 0.08327674157670514,
                            0.1019301198172403, 0.11819453196151845,
                            0.13168863844917675, 0.14209610931838218,
                            0.149172986472604, 0.15275338713072611,
                            0.15275338713072611, 0.149172986472604,
                            0.14209610931838218, 0.13168863844917675,
                            0.11819453196151845, 0.1019301198172403,
                            0.08327674157670514, 0.06267204833410799,
                            0.040601429800386134, 0.017614007139152742};

constexpr double C_coeff_matr[4][4] // See Eq. A17-A18 of J. Chem. Phys. 139, 154504 (2013); https://doi.org/10.1063/1.4819786;
    {
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

}