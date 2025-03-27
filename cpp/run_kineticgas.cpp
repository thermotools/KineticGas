#include "multiparam.h"
#include <iostream>
#include <vector>

void print_partitions(const std::vector<std::vector<int>>& partitions){
    for (const auto& p : partitions){
        for (int i : p) {
            std::cout << i << ", ";
        }
        std::cout << std::endl;
    }
}

void do_partitions(){
    std::vector<std::vector<int>> partitions = get_partitions(3, 3);
    std::cout << "Unique : \n";
    print_partitions(partitions);

    std::vector<std::vector<std::vector<int>>> built = build_partitions(6, 3);
    std::cout << "Build : " << std::endl;
    print_partitions(built[3]);
}

int main(){
    // do_partitions();
    // return 0;
    // FH_HFD_B2 kin("HE", 2);
    // FH_ModTangToennies kin("AR", 1, "default");
    ModTangToennis kin("AR");
    // PatowskiFH kin("P-H2", 2);
    // std::vector<double> T_lst = {20};
    // for (const double T : T_lst){
    //     double B = kin.viscosity(T, 1, {0.5, 0.5}, 1);
    //     std::cout << "Computed visc(" << T << ") = " << B << std::endl;
    // }
    // Patowski kin("P-H2");
    kin.set_JKWB_limits(1e9, 100);
    for (int l = 0; l < 101; l += 2){
        kin.trace_absolute_phase_shifts(0, 0, l, 325.);
    }
    // kin.dump_phase_shift_map();
    std::cout << "Second virial : " << kin.second_virial(0, 0, 600.) << std::endl;
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