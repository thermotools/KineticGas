#include "Quantum.h"
#include <iostream>
#include <vector>

int main(){
    Quantum pH2("P-H2,P-H2");
    std::cout << "Finished init, half-spins are : " << pH2.half_spin[0] << std::endl;

    Quantum oH2("O-H2,O-H2");
    std::cout << "Finished init, half-spins are : " << oH2.half_spin[0] << std::endl;

    Quantum nH2("P-H2,O-H2");
    std::cout << "Finished init, half-spins are : " << nH2.half_spin[0] << ", " << nH2.half_spin[1] << std::endl;

    Quantum he("HE,HE3");
    std::cout << "Finished init, half-spins are : " << he.half_spin[0] << ", " << he.half_spin[1] << std::endl;
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