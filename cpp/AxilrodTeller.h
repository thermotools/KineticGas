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

    double get_sigma_eff(int i, int j, double rho){
        double sigma_eff = T::sigma[i][j];
        const double tol = 1e-10;
        double f = potential(i, j, sigma_eff, rho);
        while (abs(f / T::eps[i][j]) > tol){
            sigma_eff -= f / potential_derivative_r(i, j, sigma_eff, rho);
            f = potential(i, j, sigma_eff, rho);
        }
        return sigma_eff;
    }

    double get_rmin(int i, int j, double rho){
        double r_min = T::sigma[i][j];
        const double tol = 1e-10;
        double f = potential_derivative_r(i, j, r_min, rho);
        while (abs(f * T::sigma[i][j] / T::eps[i][j]) > tol){
            r_min -= f / potential_dblderivative_rr(i, j, r_min, rho);
            f = potential_derivative_r(i, j, r_min, rho);
        }
        return r_min;
    }

    double get_epsilon_eff(int i, int j, double rho){
        const double r_min = get_rmin(i, j, rho);
        return - potential(i, j, r_min, rho);
    }

    double get_dBH(int i, int j, double temp, double rho){
        double sigma_eff = get_sigma_eff(i, j, rho);
        double beta = 1 / (BOLTZMANN * temp);
        double r_cut = 0.5;
        const auto integrand = [&](double r){return (1 - exp(- beta * potential(i, j, r * sigma_eff, rho)));};
        while (abs(integrand(r_cut) - 1) < 1e-12){
            r_cut += 1e-2;
        }
        r_cut -= 1e-2;
        double integral = r_cut;
        integral += simpson(integrand, r_cut, 1, 50);
        return integral * sigma_eff;
    }

    double get_vdw_alpha_eff(int i, int j, double rho){
        const double sigma_eff = get_sigma_eff(0, 0, rho);
        const double r_min = get_rmin(0, 0, rho) / sigma_eff;
        const double eps_eff = get_epsilon_eff(0, 0, rho);
        const auto integrand = [&](double r){return - potential(i, j, r * sigma_eff, rho) * pow(r, 2) / eps_eff;};
        double r_start = 1;
        double r_stop = r_min;
        double tol = 1e-2;
        double integral = simpson(integrand, r_start, r_stop, 100);
        r_start = r_stop; r_stop = 2 * r_min;
        integral += simpson(integrand, r_start, r_stop, 100);
        r_start = r_stop; r_stop = 3 * r_min;
        for (size_t ti = 0; ti < 8; ti++){
            tol /= 10;
            while (integrand(r_stop) > tol){
                r_stop += 1;
            }
            double part = simpson(integrand, r_start, r_stop, 50);
            integral += part;
            r_start = r_stop;
        }
        return integral;
    }

    double get_vdw_alpha(int i, int j, double rho){
        const double vdw_alpha_eff = get_vdw_alpha_eff(i, j, rho);
        const double eps_eff = get_epsilon_eff(i, j, rho);
        return vdw_alpha_eff * (eps_eff / T::eps[i][j]);
    }

protected:
    void set_internals(double rho, double temp, const vector1d& x) override {
        std::cout << "Setting current rho : " << rho * pow(T::sigma[0][0], 3) << std::endl;
        current_rho = rho;
    }

    StatePoint get_transfer_length_point(double rho, double temp, const vector1d& x) override {
        return StatePoint(temp, rho);
    }

    OmegaPoint get_omega_point(int i, int j, int l, int r, double temp){
        return OmegaPoint(i, j, l, r, temp, current_rho);
    }


private:
    vector2d alpha_pol; // Polarizability
    vector2d ion_energy; // Ionization energy
    vector2d at_energy; // Axilrod-Teller nonadditive energy

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