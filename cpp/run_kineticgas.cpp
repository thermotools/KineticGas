#include "MieKinGas.h"
#include <iostream>
#include <vector>

int main(){
    MieKinGas mie("AR,KR,NE");

    std::cout << "Finished init for " << mie.Ncomps << std::endl;
    double T = 300;
    double Vm = 20e-3;
    std::vector<double> x = {0.5, 0.3, 0.2};
    std::cout << "Viscosity : " << mie.viscosity(T, Vm, x) << std::endl;


    std::cout << "CoM => CoV :\n" << mie.CoM_to_CoV_matr(T, Vm, x) << std::endl;
    return 0;
}