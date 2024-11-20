#include "MieKinGas.h"
#include "cppThermopack/saftvrqmie.h"
#include "QuantumMie.h"
#include <iostream>
#include <vector>

int main(){
    std::cout << "Default fluid dir : " << get_fluid_dir() << std::endl;
    set_fluid_dir("../fluids");
    std::cout << "Reading from fluid dir : " << get_fluid_dir() << std::endl;

    double T = 300;
    double Vm = 20e-3;
    std::vector<double> x = {0.5, 0.5};

    double ep = 26. * BOLTZMANN;
    double sig = 3.02e-10;
    double MW = 2 * 1e-3 / AVOGADRO;
    double la = 6.;
    double lr = 9.;
    vector2d sigma(2, vector1d(2, sig));
    vector2d eps(2, vector1d(2, ep));
    vector2d lamb_a(2, vector1d(2, la));
    vector2d lamb_r(2, vector1d(2, lr));
    std::vector<int> fh_orders = {1, 1};

    Saftvrqmie svrqm("H2", 1);
    for (size_t i = 0; i < 1; i++){
        svrqm.set_pure_fluid_param(i + 1, 1., sigma[i][i], eps[i][i] / BOLTZMANN, lamb_a[i][i], lamb_r[i][i]);
    }

    MieKinGas fh0("H2");
    QuantumMie fh1("H2", 1);
    QuantumMie fh1_1({MW, MW}, sigma, eps, lamb_a, lamb_r, fh_orders, false, true);
    fh1_1.set_eos(ThermoWrapper(std::move(svrqm)));

    std::cout << "Params : " << fh1.potential_params_to_string(0) << "\n\n" << fh1_1.potential_params_to_string(0) << std::endl;
    // QuantumMie fh2("H2", 2);
    Units unt = fh0.get_reducing_units(0, 0);
    double rho_r = 0.01;
    while (rho_r < 0.8){
        Vm = 1. / (rho_r * unt.rho);
        double visc = fh1.viscosity(T, Vm, x, 2);
        double tcond = fh1.thermal_conductivity(T, Vm, x, 2);
        Eigen::MatrixXd diff = fh1.interdiffusion(T, Vm, x, 2);
        std::cout << rho_r << " : " << visc << ", " << tcond << ", " << diff(0, 0) << std::endl;
        // visc = fh1_1.viscosity(T, Vm, x, 2);
        // tcond = fh1_1.thermal_conductivity(T, Vm, x, 2);
        // diff = fh1_1.interdiffusion(T, Vm, x, 2);
        // std::cout << "\t" << rho_r << " : " << visc << ", " << tcond << ", " << diff(0, 0) << std::endl;
        rho_r += 0.05;
    }
    return 0;
    
    // MieKinGas mie2({MW, MW}, 
    //                 {{sig, sig}, {sig, sig}}, 
    //                 {{ep, ep}, {ep, ep}}, 
    //                 {{6., 6.}, {6., 6.}}, 
    //                 {{12., 12.}, {12., 12.}},
    //                 false, true);
// 
    // std::cout << "Finished direct init ..." << std::endl;
    // GenericEoS eos(ThermoWrapper(Saftvrmie("AR")));
    // mie2.set_eos(std::move(eos));
    // std::cout << "Finished set_eos ... " << std::endl;
// 
    // std::cout << "Thermal conductivity : " << mie2.thermal_conductivity(T, Vm, x) << std::endl;
    // std::cout << "Viscosity : " << mie2.viscosity(T, Vm, x) << std::endl;
// 
    // MieKinGas mie("AR,KR,NE");
// 
    // std::cout << "Finished init 2 ... " << std::endl;
    // std::vector<double> x2 = {0.5, 0.3, 0.2};
    // std::cout << "Thermal conductivity : " << mie.thermal_conductivity(T, Vm, x2) << std::endl;
    // std::cout << "Viscosity : " << mie.viscosity(T, Vm, x2) << std::endl;
// 
// 
    // std::cout << "CoM => CoV :\n" << mie.CoM_to_CoV_matr(T, Vm, x2) << std::endl;
    return 0;
}