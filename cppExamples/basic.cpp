#include "MieKinGas.h"
#include "HardSphere.h"
#include <vector>

int main(){
    MieKinGas mie("AR,KR,NE"); // MieKinGas for argon/krypton/neon mixture

    double T = 300; // Temperature (K)
    double Vm = 20e-3; // Molar volume (m3 / mol)
    std::vector<double> x = {0.5, 0.3, 0.2}; // Molar composition

    std::cout << "Computing at T = " << T << "K, Vm = " << Vm << " m3 / mol :\n";
    std::cout << "Thermal conductivity : " << mie.thermal_conductivity(T, Vm, x) << std::endl;
    std::cout << "Viscosity : " << mie.viscosity(T, Vm, x) << std::endl;
    
    std::cout << "\nDiffusion coefficients depend on choice of dependent component ...\n";
    std::cout << "With last component (Ne) as dependent : \n" << mie.interdiffusion(T, Vm, x) << "\n\n";
    std::cout << "With second component (Kr) as dependent : \n" << mie.interdiffusion(T, Vm, x, 2, FrameOfReference::CoN, 1) << "\n\n";
    std::cout << "With first component (Ar) as dependent : \n" << mie.interdiffusion(T, Vm, x, 2, FrameOfReference::CoN, 0) << "\n\n";
    return 0;
}