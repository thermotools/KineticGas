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
    return 0;
}