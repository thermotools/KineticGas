#include "LJSpline.h"
using namespace ljs_rdf_constants;

double LJSpline::get_BH_diameter(double T){
    // 20-point Gauss-Legendre from 0.5 sigma to 1 sigma. Using constant value of 1 for integrand within (0, 0.5) sigma.
    double d_BH = 2.;
    double beta = 1. / (BOLTZMANN * T);
    for (int n = 0; n < 20; n++){
        d_BH += gl_w[n] * (1. - exp(- beta * potential(0, 0, sigma[0][0] * (gl_x[n] / 4. + 3. / 4.))));
    }
    d_BH *= sigma[0][0] / 4.;
    return d_BH;
}

autodiff::dual LJSpline::g_hs_dep(const autodiff::dual& r_dual, const autodiff::dual& eta_dual) {
    autodiff::dual alpha = 0.595157 - 2.469402*eta_dual - 31.99603*pow(eta_dual,2) - 0.897749*pow(eta_dual,3);
    autodiff::dual rm = 1.960492 - 0.240131*eta_dual - 3.26102*pow(eta_dual,2) + 4.08014*pow(eta_dual,3);
    autodiff::dual gm = 1.072954 - 1.13263*eta_dual + 3.05087*pow(eta_dual,2) - 5.94853*pow(eta_dual,3);
    autodiff::dual mu = 0.58236 + 8.43296*eta_dual - 34.52268*pow(eta_dual,2) + 58.63440*pow(eta_dual,3) - 35.13525*pow(eta_dual,4);
    autodiff::dual g_sig = (1. - 1./2.*eta_dual + 5./52.*pow(eta_dual,2) - 1./4.*pow(eta_dual,3) + 1./8.*pow(eta_dual,4))/pow((1. - eta_dual),3);
    autodiff::dual B = (rm*gm - g_sig*exp(mu*(rm-1.)))/(exp(alpha*(rm-1.))-exp(mu*(rm-1.)));
    autodiff::dual A = g_sig - B; 
    return (A/r_dual*exp(mu*(r_dual-1.)) + B/r_dual*exp(alpha*(r_dual-1.)));
}


autodiff::dual LJSpline::g_hs_str(const autodiff::dual& r_dual, const autodiff::dual& eta_dual) {
    autodiff::dual omega = 6.56745 + 0.721795*eta_dual; 
    autodiff::dual kappa = 1.75879 - 3.23977*eta_dual + 0.449588*pow(eta_dual,2) - 0.148890*pow(eta_dual,3);
    autodiff::dual rm = 1.960492 - 0.240131*eta_dual - 3.26102*pow(eta_dual,2) + 4.08014*pow(eta_dual,3);
    autodiff::dual gm = 1.072954 - 1.13263*eta_dual + 3.05087*pow(eta_dual,2) - 5.94853*pow(eta_dual,3);
    autodiff::dual delta = -omega*rm - atan((kappa*rm + 1.)/(omega*rm));
    autodiff::dual C = rm*(gm - 1.)*exp(kappa*rm)/cos(omega*rm + delta);
    return 1. + C/r_dual*cos(omega*r_dual+delta)*exp(-kappa*r_dual);
}

double LJSpline::integrand_u_g(double r,double rho, double T) {
    return potential(0,0,r)*g_hs(r,rho,T)*pow(r,2);
}
double LJSpline::integrand_u_dg(double r,double rho, double T) {
    return potential(0,0,r)*dg_hs_drho(r,rho,T)*pow(r,2);
}

double LJSpline::integrand_du_g(double r,double rho, double T) {
    return potential_derivative_r(0,0,r)*g_hs(r,rho,T)*pow(r,3);
}

double LJSpline::integrand_u2_g(double r,double rho, double T) {
    return pow(potential(0,0,r),2)*g_hs(r,rho,T)*pow(r,2);
}

double LJSpline::integrand_u2_dg(double r,double rho, double T) {
    return pow(potential(0,0,r),2)*dg_hs_drho(r,rho,T)*pow(r,2);
}

double LJSpline::integrand_u_du_g(double r,double rho, double T) {
    return potential(0,0,r)*potential_derivative_r(0,0,r)*g_hs(r,rho,T)*pow(r,3);
}

double LJSpline::int_u_g(double rho, double T) {
    std::function<double(double)> func = std::bind(&LJSpline::integrand_u_g, this, std::placeholders::_1, rho, T);
    return simpson(func, sigma[0][0], rc, 1000);
}

double LJSpline::int_u_dg(double rho, double T) {
    std::function<double(double)> func = std::bind(&LJSpline::integrand_u_dg, this, std::placeholders::_1, rho, T);
    return simpson(func, sigma[0][0], rc, 1000);
}

double LJSpline::int_du_g(double rho, double T) {
    std::function<double(double)> func = std::bind(&LJSpline::integrand_du_g, this, std::placeholders::_1, rho, T);
    return simpson(func, sigma[0][0], rc, 1000);
}

double LJSpline::int_u2_g(double rho, double T) {
    std::function<double(double)> func = std::bind(&LJSpline::integrand_u2_g, this, std::placeholders::_1, rho, T);
    return simpson(func, sigma[0][0], rc, 1000);
}

double LJSpline::int_u2_dg(double rho, double T) {
    std::function<double(double)> func = std::bind(&LJSpline::integrand_u2_dg, this, std::placeholders::_1, rho, T);
    return simpson(func, sigma[0][0], rc, 1000);
}

double LJSpline::int_u_du_g(double rho, double T) {
    std::function<double(double)> func = std::bind(&LJSpline::integrand_u_du_g, this, std::placeholders::_1, rho, T);
    return simpson(func, sigma[0][0], rc, 1000);
}


double LJSpline::get_g0(double rho, double T) {
    double eta = get_eta(T, rho);
    double x0 = get_x0(T);
    double k0 = -log(1.-eta)+(42.*eta-39.*pow(eta,2)+9*pow(eta,3)-2.*pow(eta,4))/(6.*pow((1.-eta),3));
    double k1 = (pow(eta,4)+6.*pow(eta,2)-12.*eta)/(2.*pow((1-eta),3));
    double k2 = -3.*pow(eta,2)/(8.*pow((1-eta),2));
    double k3 = (-pow(eta,4)+3.*pow(eta,2)+3.*eta)/(6.*pow((1-eta),3));
    return exp(k0+k1*x0+k2*pow(x0,2)+k3*pow(x0,3));
}

double LJSpline::get_g1(double rho,double T) {
    double I1 = int_du_g(rho, T);
    double J1 = int_u_g(rho, T);
    double dJ1 = int_u_dg(rho,T);
    double d_BH = get_BH_diameter(T);
    return 3/(eps[0][0]*pow(d_BH,3))*(J1+rho*dJ1+I1/3);
}

autodiff::dual LJSpline::K_HS(autodiff::dual rho_dual, autodiff::dual T_dual) {
    autodiff::dual eta_dual = PI*rho_dual*pow(get_BH_diameter(T_dual.val),3)/6.;
    return pow((1-eta_dual),4)/(1+4.*eta_dual+4.*pow(eta_dual,2)-4.*pow(eta_dual,3)+pow(eta_dual,4));
}

double LJSpline::dK_HS_drho(double rho, double T) {
    const auto K_HS_ = [&](autodiff::dual rho_dual_, autodiff::dual T_dual_){return K_HS(rho_dual_,T_dual_);};
    autodiff::dual rho_dual = rho;
    autodiff::dual T_dual = T;
    return autodiff::derivative(K_HS_, autodiff::wrt(rho_dual), autodiff::at(rho_dual,T_dual));
}

double LJSpline::get_g2_MCA(double rho, double T) {
    autodiff::dual rho_dual = rho;
    autodiff::dual T_dual = T;
    double KHS = K_HS(rho_dual,T_dual).val;
    double dKHS = dK_HS_drho(rho,T);
    double I2 = int_u_du_g(rho,T);
    double J2 = int_u2_g(rho,T);
    double dJ2 = int_u2_dg(rho,T);
    double d_BH = get_BH_diameter(T);
    return -3/(2*pow(d_BH,3)*pow(eps[0][0],2))*(KHS*J2+rho*J2*dKHS+rho*KHS*dJ2+2./3.*KHS*I2); 
}

double LJSpline::gamma_corr(double rho, double T) {
    std::vector<double> phi = {16.1284916, -10.15024328, 13.30221007, -0.54347841, 7.21743605};
    double eta = get_eta(T,rho);
    double x0 = get_x0(T);
    double T_reduced = T*BOLTZMANN/eps[0][0];
    double theta = exp(1/T_reduced)-1;
    double gamma = phi[0]*eta*pow(x0,3)*theta*exp(phi[1]*eta*pow(x0,3)+phi[2]*pow(eta,2)*pow(x0,6))*(1+phi[3]*theta+phi[4]*eta*pow(x0,3));
    return gamma;
}

//OMEGA APPROXIMATION FUNCTIONS, SEE CORRESPONDING IMPLEMENTATION FOR MieKinGas CLASS. THESE FUNCTIONS
//ARE ADOPTED FOR THE LJSpline CLASS FROM THE MieKinGas IMPLEMENTATION:

double LJSpline::omega(int i, int j, int l, int r, double T){
    if ((l <= 2) && (r >= l) && (r <= 3)){ 
        // Use modified version of correlation by Fokin, Popov and Kalashnikov, High Temperature, Vol. 37, No. 1 (1999)
        // See [REF master thesis] for details
        // The correlation gives the logarithm of the reduced collision integral.
        // The collision integral is reduced using the hard-sphere value, i.e. lnomega_star = log(omega / omega_hs)
        double T_star = T * BOLTZMANN / eps[i][j];
        if (T_star > 0.4) return omega_correlation(i, j, l, r, T_star) * omega_hs(i, j, l, r, T);
    }
    return Spherical::omega(i, j, l, r, T);
}

double LJSpline::omega_correlation(int i, int j, int l, int r, double T_star){
    // Modified version of correlation by Fokin, Popov and Kalashnikov, High Temperature, Vol. 37, No. 1 (1999)
    // The correlation gives the logarithm of the reduced collision integral.
    // The collision integral is reduced using the hard-sphere value, i.e. lnomega_star = log(omega / omega_hs)
    if (l == r){
        double lnomega_star = - (1. / 6.) * log(T_star);
        double a_m;
        for (int n = 1; n < 7; n++){
            a_m = omega_correlation_factors[l-1][n-1];
            lnomega_star += a_m * pow(T_star, (1. - n) * 0.5);
        }
        return exp(lnomega_star);
    }
    return omega_recursive_factor(i, j, l, r - 1, T_star) * omega_correlation(i, j, l, r - 1, T_star);
}

double LJSpline::omega_recursive_factor(int i, int j, int l, int r, double T_star){
    // Higher order collision integrals can be computed from the derivative of lower order ones using the recursion
    // Given in Fokin, Popov and Kalashnikov, High Temperature, Vol. 37, No. 1 (1999)
    // See also: Hirchfelder, Curtiss & Bird, Molecular Theory of Gases and Liquids.
    // For reduced integrals : omega(l, r + 1) / omega(l, r) = 1 + (d omega(l, r) / d lnT^*) / (r + 2)
    if ((l > 2) || (r > 3)) {throw std::runtime_error("No recursive factor available!");}
    double dlnomega_dlnT = - (1. / 6.);
    double a_m;
    for (int n = 2; n < 7; n++){
        a_m = omega_correlation_factors[l-1][n-1];
        dlnomega_dlnT += a_m * pow(T_star, - (n - 1.) / 2.) * (1. - n) / 2.;
    }
    return 1. + dlnomega_dlnT / (r + 2);
}
