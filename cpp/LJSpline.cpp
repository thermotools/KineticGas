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


double LJSpline::get_g0(double rho, double T) {
    double eta = get_eta(T, rho);
    double x0 = get_x0(T);
    double k0 = -log(1.-eta)+(42.*eta-39.*pow(eta,2)+9*pow(eta,3)-2.*pow(eta,4))/(6.*pow((1.-eta),3));
    double k1 = (pow(eta,4)+6.*pow(eta,2)-12.*eta)/(2.*pow((1-eta),3));
    double k2 = -3.*pow(eta,2)/(8.*pow((1-eta),2));
    double k3 = (-pow(eta,4)+3.*pow(eta,2)+3.*eta)/(6.*pow((1-eta),3));
    return exp(k0+k1*x0+k2*pow(x0,2)+k3*pow(x0,3));
}

double LJSpline::I1(double rho_star, double T) {
    // Uses rho_star = rho*sigma^3*N_avogadro
    double d = get_BH_diameter(T) / sigma[0][0];
    double integral = q_I1[0] + q_I1[1]*rho_star + q_I1[2]*pow(rho_star,2) + q_I1[3]*pow(rho_star,3) + q_I1[4]*pow(rho_star,4) + 
          rho_star*(q_I1[5]*pow(rho_star,2) + q_I1[6]*rho_star + q_I1[7])*(d-1) + 
          rho_star*(q_I1[8]*pow(rho_star,2) + q_I1[9]*rho_star + q_I1[10])*pow((d-1),2);
    return integral;
}

double LJSpline::I1_diff(double rho_star, double T) {
    double d = get_BH_diameter(T) / sigma[0][0];
    double integral = q_I1[1] + 2*q_I1[2]*rho_star + 3*q_I1[3]*pow(rho_star,2) + 4*q_I1[4]*pow(rho_star,3) + 
          (3*q_I1[5]*pow(rho_star,2) + 2*q_I1[6]*rho_star + q_I1[7])*(d-1) + 
          (3*q_I1[8]*pow(rho_star,2) + 2*q_I1[9]*rho_star + q_I1[10])*pow((d-1),2);
    return integral;
}

double LJSpline::I2(double rho_star, double T) {
    // Uses rho_star = rho*sigma^3*N_avogadro
    double d = get_BH_diameter(T) / sigma[0][0];
    double integral = q_I2[0] + q_I2[1]*rho_star + q_I2[2]*pow(rho_star,2) + q_I2[3]*pow(rho_star,3) + q_I2[4]*pow(rho_star,4) + 
          rho_star*(q_I2[5]*pow(rho_star,2) + q_I2[6]*rho_star + q_I2[7])*(d-1) + 
          rho_star*(q_I2[8]*pow(rho_star,2) + q_I2[9]*rho_star + q_I2[10])*pow((d-1),2);
    return integral;
}

double LJSpline::I2_diff(double rho_star, double T) {
    double d = get_BH_diameter(T) / sigma[0][0];
    double integral = q_I2[1] + 2*q_I2[2]*rho_star + 3*q_I2[3]*pow(rho_star,2) + 4*q_I2[4]*pow(rho_star,3) + 
          (3*q_I2[5]*pow(rho_star,2) + 2*q_I2[6]*rho_star + q_I2[7])*(d-1) + 
          (3*q_I2[8]*pow(rho_star,2) + 2*q_I2[9]*rho_star + q_I2[10])*pow((d-1),2);
    return integral;
}

double LJSpline::J1(double rho_star, double T) {
    // Uses rho_star = rho*sigma^3*N_avogadro
    double d = get_BH_diameter(T) / sigma[0][0];
    double integral = q_J1[0] + q_J1[1]*rho_star + q_J1[2]*pow(rho_star,2) + q_J1[3]*pow(rho_star,3) + q_J1[4]*pow(rho_star,4) + 
          rho_star*(q_J1[5]*pow(rho_star,2) + q_J1[6]*rho_star + q_J1[7])*(d-1) + 
          rho_star*(q_J1[8]*pow(rho_star,2) + q_J1[9]*rho_star + q_J1[10])*pow((d-1),2);
    return integral;
}

double LJSpline::J2(double rho_star, double T) {
    // Uses rho_star = rho*sigma^3*N_avogadro
    double d = get_BH_diameter(T) / sigma[0][0];
    double integral = q_J2[0] + q_J2[1]*rho_star + q_J2[2]*pow(rho_star,2) + q_J2[3]*pow(rho_star,3) + q_J2[4]*pow(rho_star,4) + 
          rho_star*(q_J2[5]*pow(rho_star,2) + q_J2[6]*rho_star + q_J2[7])*(d-1) + 
          rho_star*(q_J2[8]*pow(rho_star,2) + q_J2[9]*rho_star + q_J2[10])*pow((d-1),2);
    return integral;
}

double LJSpline::K_HS(double rho, double T) {
    // Isothermal compressibility of HS fluid from Carnahan-Starling EOS
    double eta = get_eta(T, rho);
    return pow((1-eta),4)/(1.+4.*eta+4.*pow(eta,2)-4.*pow(eta,3)+pow(eta,4));
}

double LJSpline::K_HS_deta(double rho, double T) {
    double eta = get_eta(T, rho);
    double temp = pow(eta, 4) - 4. * pow(eta, 3) + 4. * pow(eta, 2) + 4. * eta + 1.;
    return 4. * (pow(eta, 2) - 5. * eta - 2.) * pow(1. - eta, 3) / pow(temp, 2);
}

double LJSpline::get_g1(double rho,double T) {
    double rho_star = rho*pow(sigma[0][0],3);
    double i1 = I1(rho_star, T);
    double di1 = I1_diff(rho_star, T);
    double j1 = J1(rho_star, T);
    double d = get_BH_diameter(T) / sigma[0][0];
    return (3*di1*rho_star + 3*i1 + j1)/pow(d,3);
}

double LJSpline::get_g2_MCA(double rho, double T) {
    double rho_star = rho*pow(sigma[0][0],3);
    double i2 = I2(rho_star, T);
    double di2 = I2_diff(rho_star, T);
    double j2 = J2(rho_star, T);
    double d = get_BH_diameter(T) / sigma[0][0];
    double khs = K_HS(rho, T);
    double dkhs = K_HS_deta(rho,T)*(PI*pow(d,3))/6;
    return (-3/(2*pow(d,3)))*(i2*khs + rho_star*i2*dkhs + rho_star*khs*di2 + 2*khs*j2/3);
}


double LJSpline::gamma_corr(double rho, double T) {
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
        if (T_star > 0.4) {
            return omega_correlation(i, j, l, r, T_star) * omega_hs(i, j, l, r, T);}
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
