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
}