#pragma once
#include "ModTangToennis.h"

template <typename T>
class AxilrodTellerCorrected: public T {
    public:

    AxilrodTellerCorrected(std::string comps, bool is_idealgas) : T(comps), 
        alpha_pol(Ncomps, vector1d(Ncomps, 0.)), ion_energy(Ncomps, vector1d(Ncomps, 0.)) 
    {
        for (size_t i = 0; i < Ncomps; i++){
            for (size_t j = 0; j < Ncomps; j++){
                alpha_pol[i][j] = sqrt(compdata[i]["polarizability"] * compdata[j]["polarizability"]);
                ion_energy[i][j] = = sqrt(compdata[i]["ionization_energy"] * compdata[j]["ionization_energy"]);
            }
        }
    }

    dual2 potential(int i, int j, dual2 r, dual2 rho){
        dual2 u_dilute = T::potential(i, j, r);
        return u_dilute + ion_energy[i][j] * pow(alpha_pol[i][j] / pow(r, 2), 3) * (PI / 2.) * rho;
    }

    private:
    vector2d alpha_pol; // Polarizability
    vector2d ion_energy; // Ionization energy
};

using AT_TangToennies = AxilrodTellerCorrected<ModTangToennis>;