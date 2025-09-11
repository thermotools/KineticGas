/*Author: Johannes Salomonsen Løken
  Description: Revised Enskog theory for the Lennard-Jones/spline (LJ/s) fluid. The implementation uses the RDF at contact model of Lafitte et al. (https://doi.org/10.1063/1.4819786)
               adapted to suit the LJ/s fluid. For more details, see (Løken, Johannes Salomonsen. Revised Enskog theory and extended corresponding states models for the transport 
               properties of the Lennard-Jones/Spline fluid, MS thesis. NTNU, 2025)
*/

#pragma once
#include "Spherical.h"
#include "cppThermopack/ljs.h"
#include "HardSphere.h"

class LJSpline : public Spherical {
    public:
    double a, b, rs, rc; 
    LJSpline(bool is_idealgas, bool is_singlecomp)
    : Spherical({6.64e-26, 6.64e-26},
                {{3.42e-10, 3.42e-10}, {3.42e-10, 3.42e-10}},
                {{1.71200476e-21, 1.71200476e-21}, {1.71200476e-21, 1.71200476e-21}},
                is_idealgas, is_singlecomp),

                a{-24192. / 3211.},
                b{-387072. / 61009.},
                rs{pow(26./7. , 1./6.)*sigma[0][0]},
                rc{67. * rs / 48.}
    {
        LJs_uv uv_eos{"Default",1.0};
        // bh_eos.set_sigma_eps(sigma[0][0],eps[0][0]);
        GenericEoS ljs_eos{ThermoWrapper(std::move(uv_eos))};
        this -> set_eos(std::move(ljs_eos));
        if ((sigmaij.size() > 2) | (mole_weights.size() > 2) | (eps.size() > 2)) 
        {
            throw std::invalid_argument("The Lennard-Jones/spline is not implemented for multicomponent systems.");
        }
    }   

    //Changing Spherical::theta_integrand, such that a 0 / 0 error is avoided when potential is exactly 0.
    double theta_integrand(int i, int j, double T, double r, double g, double b) override{
        if (g == 0) return 0;
        double t = (pow(r, 4) / pow(b, 2)) * (1.0 - potential(i, j, r) / (BOLTZMANN * T * pow(g, 2))) - pow(r, 2);
        if (t < 0) return 0;
        return pow(t,-0.5); 
    }

    //Potential definitions:
    dual2 potential(int i, int j, dual2 r) override {
        if (r < rs) {
            return 4 * eps[0][0] * (pow(sigma[0][0] / r, 12) - pow(sigma[0][0] / r, 6));
        }
        else if (r < rc) {
            return eps[0][0]*(a * pow(r/rs - rc/rs , 2) + b * pow(r/rs-rc/rs , 3));
        }
        else {
            return static_cast<dual2>(0);
        }
    }

    double potential(int i, int j, double r) override{
        return potential(i, j, static_cast<dual>(r)).val.val;
    };

    double potential_derivative_r(int i, int j, double r) override {
        if (r < rs) {
            return -24*eps[0][0]/r*(2*pow(sigma[0][0] / r, 12) - pow(sigma[0][0] / r, 6));
        }
        else if (r < rc) {
            return eps[0][0]*(2*a*(r/rs - rc/rs)/rs + 3*b*pow(r/rs-rc/rs , 2)/rs);
        }
        else {
            return 0;
        }
    }
    double potential_dblderivative_rr(int i, int j, double r) override {
        if (r < rs) {
            return 24*eps[0][0]/pow(r,2)*(26*pow(sigma[0][0] / r, 12) - 7*pow(sigma[0][0] / r, 6));
        }
        else if (r < rc) {
            return eps[0][0]*(2*a/pow(rs,2) + 6*b*(r/rs-rc/rs)/pow(rs,2));
        }
        else {
            return 0;
        }
    }

    // HS-RDF functions 

    double get_BH_diameter(double T);
    double get_eta(double T, double rho) {return PI*rho*pow(get_BH_diameter(T),3)/6.;}
    double get_x0(double T) {return sigma[0][0]/get_BH_diameter(T);}

    // g at contact functions. Calclulates the g(\sigma) = g_0 + (beta*epsilon)*g_1 + (beta*epsilon)²*g_2

    inline std::vector<std::vector<double>> model_rdf(double rho, double T, const std::vector<double>& x) override {
        return saft_rdf(rho, T, 2);
    }

    double get_g0(double rho, double T);
    double get_g1(double rho, double T);
    double get_g2_MCA(double rho, double T);
    double gamma_corr(double rho, double T);
    double get_g2(double rho, double T) {return (1+gamma_corr(rho,T))*get_g2_MCA(rho,T);}

    double I1(double rho, double T);
    double I1_diff(double rho, double T);
    double I2(double rho, double T);
    double I2_diff(double rho, double T);

    double J1(double rho, double T);
    double J2(double rho, double T);

    double K_HS(double rho, double T);
    double K_HS_deta(double rho, double T);

    std::vector<std::vector<double>> saft_rdf(double rho, double T, int order) {
        double beta_star = eps[0][0]/(BOLTZMANN * T);
        double g;
        double g0 = get_g0(rho, T);
        double g1 = get_g1(rho,T);
        double g2 = get_g2(rho,T);
        if (order == 0) {g = g0;} 
        else if (order == 1) {g = g0 + beta_star*g1;}
        else if (order == 2) {g = g0 + beta_star*g1 + pow(beta_star,2)*g2;}
        else {throw std::invalid_argument("RDF order is invalid");}
        std::vector<std::vector<double>> rdf_vec = {{g,g},{g,g}}; 
        return rdf_vec;
    };

    // Collision integral - parametrizations. 

    double omega(int i, int j, int l, int r, double T) override;
    double omega_correlation(int i, int j, int l, int r, double T_star);
    double omega_recursive_factor(int i, int j, int l, int r, double T);


    // The hard sphere integrals are used as the reducing factor for the correlations.
    // So we need to compute the hard-sphere integrals to convert the reduced collision integrals from 
    // a modified version of the correlation by Fokin et. al. to the "real" collision integrals.
    
    static constexpr double omega_correlation_factors[2][6] =
        {
         {0.12381066, -0.19473855, -0.12605004,  1.42026419, -1.14753069,  0.25563792},
         {0.2525206,  -0.12502955, -1.02633768,  3.31321146, -2.58470784,  0.62745361}
        };

    inline double omega_hs(int i, int j, int l, int r, double T){
        double w = PI * pow(sigma[i][j], 2) * 0.5 * (r + 1);
        for (int ri = r; ri > 1; ri--) {w *= ri;}
        if (l % 2 == 0){
            w *= (1. - (1.0 / (l + 1.)));
        }
        if (i == j) return sqrt((BOLTZMANN * T) / (PI * m[i])) * w;
        return sqrt(BOLTZMANN * T * (m[i] + m[j]) / (2. * PI * m[i] * m[j])) * w;
    }
    double omega_star(int i, int j, double T, int l, int r) {
        //Useful for making/adjusting empirical correlations for the collision integrals. 
        //Returns \Omega^{*(l,r)} = \Omega^{(l,r)} / \Omega^{(l,r)}_{HS}
        return Spherical::omega(i,j,l,r,T) / omega_hs(i,j,l,r,T);}

    double omega_star_approx(int i, int j, double T, int l, int r) {
        return omega(i,j,l,r,T) / omega_hs(i,j,l,r,T);}
};



namespace ljs_rdf_constants{


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

// Constants for integral parametrizations in the RDF at contact model.

constexpr double q_I1[11] = {-0.53572874, -0.37806724,  -0.23016988,   0.34630713,   0.07308556,   5.27289035,
                            -6.43165849, -1.63392187,  16.74650716, -27.36299942,  -0.43778915};

constexpr double q_J1[11] = {1.60718622, -7.09047972e-02, -9.43081376e-01, -5.09423352e+00,  1.36949006e+00,
                            -5.81631691e+01, -1.39795273e+00,  3.66229139e+00, -2.77904034e+02,
                            5.65474349e+01,  1.01557784e+01};

constexpr double q_I2[11] = {0.3592528, 0.27758414, 0.25126487, -0.22023453, -0.11005555, -4.46920946,
                            6.117315, 0.95067716, -13.88450783, 24.30452001, -0.63636868};

constexpr double q_J2[11] = {-0.5388792, 1.16814143e-02,  8.83906862e-02,  2.36081648e+00, -5.31094654e-01,
                            2.62590123e+01, -2.13956077e+00, -1.20815977e+00,  1.18753081e+02,
                            -2.84345337e+01, -2.53564086e+00}; 

constexpr double phi[4] = {13.2795715, 3.84921999, -2.16074213, -1.40190623};
}