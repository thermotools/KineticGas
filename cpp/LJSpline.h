/*Author: Johannes Salomonsen Løken
*/

#pragma once
#include "Spherical.h"
#include "cppThermopack/ljs.h"
#include "HardSphere.h"

class IdealDummy{
    public:
    int VAPPH;
    IdealDummy() {};
    inline std::vector<std::vector<double>> dmudn(double T, double V, const std::vector<double> n) const {return {{GAS_CONSTANT*T/n[0]}};}
    inline double Cp_ideal(double T, int ci) const {return 5. * GAS_CONSTANT / 2.;}
    inline double pressure_tv(double T, double V, const std::vector<double> n) const {throw std::runtime_error("pressure_tv not implemented!");}
    inline double specific_volume(double t, double p, const std::vector<double> n, int phase) const {throw std::runtime_error("specific_volume not implemented!");}
    inline std::vector<double> dvdn(double t, double p, const std::vector<double> n, int phase) const {throw std::runtime_error("dvdn not implemented!");}
};

class LJSpline : public Spherical {
    public:
    double a, b, rs, rc; 
    //
    LJSpline(std::vector<double> mole_weights, std::vector<std::vector<double>> sigmaij, std::vector<std::vector<double>> eps, bool is_idealgas, bool is_singlecomp)
    : Spherical(mole_weights, sigmaij, eps, is_idealgas, is_singlecomp), a{-24192. / 3211.}, b{-387072. / 61009.},
    rs{pow(26./7. , 1./6.)*sigma[0][0]}, rc{67. * rs / 48.} 
    {
        LJs_bh bh_eos{"Default",1.0};
        // bh_eos.set_sigma_eps(sigma[0][0],eps[0][0]);
        GenericEoS ljs_eos{ThermoWrapper(std::move(bh_eos))};
        this -> set_eos(std::move(ljs_eos));
        if ((sigmaij.size() > 2) | (mole_weights.size() > 2) | (eps.size() > 2)) 
        {
            throw std::invalid_argument("The Lennard-Jones/spline is not implemented for multicomponent systems (yet)!");
        }
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

    // Empirical Hard-sphere RDF by H. Liu, https://doi.org/10.1080/00268976.2023.2299250

    autodiff::dual g_hs_dep(const autodiff::dual& r_dual, const autodiff::dual& eta_dual);
    autodiff::dual g_hs_str(const autodiff::dual& r_dual, const autodiff::dual& eta_dual);

    double g_hs(double r, double rho, double T) {
        double d_BH = get_BH_diameter(T);
        r = r / d_BH;
        autodiff::dual r_dual = r;
        autodiff::dual eta_dual = get_eta(T, rho);
        double rm = 1.960492 - 0.240131*eta_dual.val - 3.26102*pow(eta_dual.val,2) + 4.08014*pow(eta_dual.val,3);
        if (r < 1) {
            return 0;
        }
        else if (r < rm) {
            return g_hs_dep(r_dual,eta_dual).val;
        }
        else {
            return g_hs_str(r_dual,eta_dual).val;
        }
    }

    double dg_hs_drho(autodiff::dual r_dual, autodiff::dual rho_dual, autodiff::dual T_dual) {
        const auto ghs_dep_ = [&](autodiff::dual r_dual_, autodiff::dual eta_dual_){return g_hs_dep(r_dual_,eta_dual_);};
        const auto ghs_str_ = [&](autodiff::dual r_dual_, autodiff::dual eta_dual_){return g_hs_str(r_dual_,eta_dual_);};
        double d_BH = get_BH_diameter(T_dual.val);
        r_dual = r_dual /d_BH;
        double r = r_dual.val;
        autodiff::dual eta_dual = get_eta(T_dual.val, rho_dual.val);
        double rm = 1.960492 - 0.240131*eta_dual.val - 3.26102*pow(eta_dual.val,2) + 4.08014*pow(eta_dual.val,3);
        if (r < 1) {
            return 0;
        }
        else if (r < rm) {
            double dg_deta = autodiff::derivative(ghs_dep_,autodiff::wrt(eta_dual), autodiff::at(r_dual,eta_dual));
            return dg_deta*PI*pow(d_BH,3)/6.;
        }
        else {
            double dg_deta = autodiff::derivative(ghs_str_,autodiff::wrt(eta_dual), autodiff::at(r_dual,eta_dual));
            return dg_deta*PI*pow(d_BH,3)/6.;
        }
    }
    
    // g at contact functions. Calclulates the g(\sigma) = g_0 + (beta*epsilon)*g_1 + (beta*epsilon)²*g_2

    inline std::vector<std::vector<double>> model_rdf(double rho, double T, const std::vector<double>& x) override {
        return saft_rdf(rho, T, 2);
    }

    double get_g0(double rho, double T);
    double get_g1(double rho, double T);
    double get_g2_MCA(double rho, double T);
    double gamma_corr(double rho, double T);
    double get_g2(double rho, double T) {return (1+gamma_corr(rho,T))*get_g2_MCA(rho,T);}

    std::vector<std::vector<double>> saft_rdf(double rho, double T, int order) {
        double beta = 1./(BOLTZMANN * T);
        double g = get_g0(rho,T) + beta*eps[0][0]*get_g1(rho,T)+pow((beta*eps[0][0]),2)*get_g2(rho,T); 
        std::vector<std::vector<double>> rdf_vec = {{g,g},{g,g}}; //quick fix, burde generaliseres
        return rdf_vec;
    };

    double integrand_u_g(double r, double rho, double T);
    double integrand_du_g(double r, double rho, double T);
    double integrand_u_du_g(double r,double rho, double T);
    double integrand_u_dg(double r, double rho, double T);
    double integrand_u2_g(double r,double rho, double T);
    double integrand_u2_dg(double r,double rho, double T);
    
    double int_u_g(double rho, double T); // //Retruns integral of U(r)*g_HS(r)*r²dr
    double int_du_g(double rho, double T); //Retruns integral of dU(r)/dr*g_HS(r)*r³dr
    double int_u_du_g(double rho, double T); //Returns integral of dU(r)/dr*U(r)g_HS(r)*r³dr
    double int_u_dg(double rho, double T); //Returns integral of U(r)*dg_HS/d(rho)*r²dr
    double int_u2_g(double rho, double T); //Returns integral of U(r)²*g_HS*r²dr
    double int_u2_dg(double rho, double T); //Returns integral of U(r)²d(g_HS)/d(rho)*r²dr

    autodiff::dual K_HS(autodiff::dual rho_dual, autodiff::dual T);
    double dK_HS_drho(double rho, double T);

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
}