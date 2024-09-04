#pragma once
#include "ModTangToennis.h"

template <typename T>
class AxilrodTellerCorrected: public T {
public:

    AxilrodTellerCorrected(std::string comps, bool is_idealgas, std::string parameter_ref="default") : T(comps, is_idealgas, parameter_ref)
    {   
        alpha_pol = vector2d(this->Ncomps, vector1d(this->Ncomps, 0.));
        ion_energy = vector2d(this->Ncomps, vector1d(this->Ncomps, 0.)); 
        at_energy = vector2d(this->Ncomps, vector1d(this->Ncomps, 0.)); 
        for (size_t i = 0; i < this->Ncomps; i++){
            for (size_t j = 0; j < this->Ncomps; j++){
                alpha_pol[i][j] = sqrt(static_cast<double>(this->compdata[i]["polarizability"])
                                        * static_cast<double>(this->compdata[j]["polarizability"]));
                ion_energy[i][j] = sqrt(static_cast<double>(this->compdata[i]["ionization_energy"])
                                        * static_cast<double>(this->compdata[j]["ionization_energy"]));
                at_energy[i][j] = sqrt(static_cast<double>(this->compdata[i]["AxilrodTellerCorrection"]["default"]["nu"])
                                        * static_cast<double>(this->compdata[j]["AxilrodTellerCorrection"]["default"]["nu"]));
            }
        }
    }

    dual2 potential(int i, int j, dual2 r, dual2 rho){
        dual2 E2 = T::potential(i, j, r);
        dual2 correction_factor = - 0.85 * (at_energy[i][j] * rho) / (T::eps[i][j] * pow(T::sigma[i][j], 6));
        dual2 E3 = correction_factor * E2;
        // const double nu_red = at_energy[i][j] / (T::eps[i][j] * pow(T::sigma[i][j], 9));
        // const double rho_red = rho.val.val * pow(T::sigma[i][j], 3);
        // dual2 E3 = ion_energy[i][j] * pow(alpha_pol[i][j] / pow(r, 2), 3) * (PI / 2.) * rho;
        return E2 + E3;
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
protected:
    void set_internals(double rho, double temp, const vector1d& x) override {
        current_rho = rho;
    }

    StatePoint get_transfer_length_point(double rho, double temp, const vector1d& x) override {
        return StatePoint(rho, temp);
    }

    OmegaPoint get_omega_point(int i, int j, int l, int r, double temp){
        return OmegaPoint(i, j, l, r, temp, current_rho);
    }


private:
    vector2d alpha_pol; // Polarizability
    vector2d ion_energy; // Ionization energy
    vector2d at_energy;

    double current_rho{-1};

    dual2 potential(int i, int j, dual2 r) override {
        return potential(i, j, r, current_rho);
    }

    double potential(int i, int j, double r) override {
        return potential(i, j, r, current_rho);
    }

    double potential_derivative_r(int i, int j, double r) override {
        return potential_derivative_r(i, j, r, current_rho);
    }

    double potential_dblderivative_rr(int i, int j, double r) override {
        return potential_dblderivative_rr(i, j, r, current_rho);
    }
};

using AT_TangToennies = AxilrodTellerCorrected<ModTangToennis>;