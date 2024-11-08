#include "MieKinGas.h"
#include "cppThermopack/saftvrmie.h"
#include "HardSphere.h"
#include <iostream>
#include <vector>

int main(){
    std::cout << "Default fluid dir : " << get_fluid_dir() << std::endl;
    set_fluid_dir("../fluids");
    std::cout << "Reading from fluid dir : " << get_fluid_dir() << std::endl;

    double T = 300;
    double Vm = 20e-3;
    std::vector<double> x = {0.5, 0.5};

    double ep = 120 * BOLTZMANN;
    double sig = 3.4e-10;
    double MW = 40*1e-3 / AVOGADRO;

    HardSphere hs({MW, MW}, {{sig, sig}, {sig, sig}}, false, false);
    std::cout << hs.interdiffusion(T, Vm, x, 2, FrameOfReference::CoN, 0, 0, true);
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