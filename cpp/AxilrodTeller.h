#pragma once
#include "ModTangToennis.h"

template <typename T>
class AxilrodTellerCorrected: public T {
    public:

    AxilrodTellerCorrected(std::string comps, bool is_idealgas) : T(comps, is_idealgas)
    {
        alpha_pol = vector2d(this->Ncomps, vector1d(this->Ncomps, 0.));
        ion_energy = vector2d(this->Ncomps, vector1d(this->Ncomps, 0.)); 
        for (size_t i = 0; i < this->Ncomps; i++){
            for (size_t j = 0; j < this->Ncomps; j++){
                alpha_pol[i][j] = sqrt(static_cast<double>(this->compdata[i]["polarizability"])
                                        * static_cast<double>(this->compdata[j]["polarizability"]));
                ion_energy[i][j] = sqrt(static_cast<double>(this->compdata[i]["ionization_energy"])
                                        * static_cast<double>(this->compdata[j]["ionization_energy"]));
            }
        }
    }

    double potential(int i, int j, double r, double rho){
        dual2 rd{r}, rhod{rho};
        return potential(i, j, rd, rhod).val.val;
    }

    double potential_derivative_r(int i, int j, double r, double rho){
        dual2 rd = r;
        const auto func = [&](dual2 r_){return potential(i, j, r_, rho);};
        auto [u0, ur, urr] = autodiff::derivatives(func, autodiff::wrt(rd), autodiff::at(rd));
        return ur;
    }
    
    double potential_dblderivative_rr(int i, int j, double r, double rho){
        dual2 rd = r;
        const auto func = [&](dual2 r_){return potential(i, j, r_, rho);};
        auto [u0, ur, urr] = autodiff::derivatives(func, autodiff::wrt(rd, rd), autodiff::at(rd));
        return urr;
    }

    dual2 potential(int i, int j, dual2 r, dual2 rho){
        dual2 u_dilute = T::potential(i, j, r);
        return u_dilute + ion_energy[i][j] * pow(alpha_pol[i][j] / pow(r, 2), 3) * (PI / 2.) * rho;
    }

    private:
    vector2d alpha_pol; // Polarizability
    vector2d ion_energy; // Ionization energy

    using T::potential;
    using T::potential_derivative_r;
    using T::potential_dblderivative_rr;
};

using AT_TangToennies = AxilrodTellerCorrected<ModTangToennis>;