#pragma once
#include "utils.h"

template<typename T, size_t N, size_t deg>
class Tabulated : public T {
public:

    template<typename... Args>
    Tabulated(Args&&... args)
        : T(std::forward<Args>(args)...), 
          tab_potential(T::Ncomps, std::vector(T::Ncomps, FuncTable<N, deg>())),
          tab_potential_r(T::Ncomps, std::vector(T::Ncomps, FuncTable<N, deg>())),
          tab_potential_rr(T::Ncomps, std::vector(T::Ncomps, FuncTable<N, deg>()))
    {
        update_tables();
    }

    inline void set_using_tabulated(bool use_tab){
        using_tabulated = use_tab; 
        update_tables();
    };

    void update_tables(){
        for (size_t i = 0; i < T::Ncomps; i++){
            for (size_t j = 0; j < T::Ncomps; j++){
                double tab_min = 0.5 * T::sigma[i][j];
                double tab_mid = T::sigma[i][j];
                double tab_mid_n = 100.;
                double tab_max = (tab_mid - tab_min) * (N / tab_mid_n) + tab_min;
                tab_potential[i][j] = FuncTable<N, deg>([&](double r){return T::potential(i, j, r);}, tab_min, tab_max);
                tab_potential_r[i][j] = FuncTable<N, deg>([&](double r){return T::potential_derivative_r(i, j, r);}, tab_min, tab_max);
                tab_potential_rr[i][j] = FuncTable<N, deg>([&](double r){return T::potential_dblderivative_rr(i, j, r);}, tab_min, tab_max);
            }
        }
    }

    size_t set_internals(double rho, double temp, const vector1d& x) override {
        size_t r = T::set_internals(rho, temp, x);
        if (using_tabulated && r) update_tables();
        return r;
    }

    using T::potential;
    double potential(int i, int j, double r) override {
        if (using_tabulated){
            if (j < i) std::swap(i, j);
            if (tab_potential[i][j].in_valid_range(r)){
                return tab_potential[i][j].unsafe_eval(r);
            }
        }
        return T::potential(i, j, r);
    }

    double potential_derivative_r(int i, int j, double r) override {
        if (using_tabulated){
            if (j < i) std::swap(i, j);
            if (tab_potential_r[i][j].in_valid_range(r)){
                return tab_potential_r[i][j].unsafe_eval(r);
            }
        }
        return T::potential_derivative_r(i, j, r);
    }

    double potential_dblderivative_rr(int i, int j, double r) override {
        if (using_tabulated){
            if (j < i) std::swap(i, j);
            if (tab_potential_rr[i][j].in_valid_range(r)){
                return tab_potential_rr[i][j].unsafe_eval(r);
            }
        }
        return T::potential_dblderivative_rr(i, j, r);
    }
    
private:
    bool using_tabulated = true;
    std::vector<std::vector<FuncTable<N, deg>>> tab_potential;
    std::vector<std::vector<FuncTable<N, deg>>> tab_potential_r;
    std::vector<std::vector<FuncTable<N, deg>>> tab_potential_rr;
};


template<typename T, size_t N, size_t deg>
class Splined : public T {
public:

    template<typename... Args>
    Splined(Args&&... args)
        : T(std::forward<Args>(args)...), 
        spl_min{0.5}, spl_max{10.},
        A_head(T::Ncomps, vector1d(T::Ncomps, 0.)),
        B_head(T::Ncomps, vector1d(T::Ncomps, 0.)),
        C_tail(T::Ncomps, vector1d(T::Ncomps, 0.)),
        spl_potential(T::Ncomps, std::vector(T::Ncomps, FuncSpline<N, deg>()))
    {
        update_spline();
    }

    inline void set_using_spline(bool use_spl){
        using_spline = use_spl; 
        update_spline();
    };

    bool get_spline_is_active(){return using_spline;}

    void set_spline_domain(double start, double stop){
        spl_min = start, spl_max = stop;
        update_spline();
    }

    std::pair<double, double> get_spline_domain(){
        double dx = deg * (spl_max - spl_min) / N;
        return std::pair<double, double>(spl_min + dx, spl_max - dx);
    }

    void update_spline(){
        for (size_t i = 0; i < T::Ncomps; i++){
            for (size_t j = 0; j < T::Ncomps; j++){
                spl_potential[i][j] = FuncSpline<N, deg>([&](double r){return T::potential(i, j, r);}, spl_min * T::sigma[i][j], spl_max * T::sigma[i][j]);
                C_tail[i][j] = T::potential(i, j, spl_max * T::sigma[i][j]) * pow(spl_max, 6);
                double u_head = T::potential(i, j, spl_min * T::sigma[i][j]);
                double ur_head = T::potential_derivative_r(i, j, spl_min * T::sigma[i][j]);
                B_head[i][j] = - u_head / ur_head;
                A_head[i][j] = u_head / exp(- B_head[i][j] * spl_min * T::sigma[i][j]);
            }
        }
        
    }

    size_t set_internals(double rho, double temp, const vector1d& x) override {
        size_t r = T::set_internals(rho, temp, x);
        if (using_spline && r) update_spline();
        return r;
    }

    double potential_dn(int i, int j, double r, int n){
        if (n < 0) throw std::out_of_range("Cannot compute derivative < 0.");
        if (using_spline){
            if (j < i) std::swap(i, j);
            if (spl_potential[i][j].in_valid_range(r)){
                return spl_potential[i][j].derivative(r, n);
            }
            if (n > 2) {
                std::pair<double, double> domain = get_spline_domain();
                if (r >= domain.second * T::sigma[i][j]) {
                    return pow(-1, n + 1) * C_tail[i][j] * pow(T::sigma[i][j] / r, 6 + n) / pow(T::sigma[i][j], n);
                }
                return pow( - B_head[i][j], n) * A_head[i][j] * exp( - B_head[i][j] * r);
                throw std::out_of_range("Derivative order (" + std::to_string(n) + " > 2) only available in Splined region (" 
                                        + std::to_string(domain.first) + ", " + std::to_string(domain.second) + "), got " + std::to_string(r / T::sigma[i][j]) );
            }
        }
        switch (n){
            case 0: return T::potential(i, j, r); 
            case 1: return T::potential_derivative_r(i, j, r);
            case 2: return T::potential_dblderivative_rr(i, j, r);
            default: throw std::out_of_range("Derivative order (" + std::to_string(n) + " > 2) only available when Spline is active (use set_using_spline).");
        }
    }

    using T::potential;
    double potential(int i, int j, double r) override {return potential_dn(i, j, r, 0); }
    double potential_derivative_r(int i, int j, double r) override {return potential_dn(i, j, r, 1);}
    double potential_dblderivative_rr(int i, int j, double r) override {return potential_dn(i, j, r, 2);}
    

private:
    bool using_spline = true;
    double spl_min, spl_max;
    vector2d A_head, B_head, C_tail;
    std::vector<std::vector<FuncSpline<N, deg>>> spl_potential;

};

template<typename T>
class FH_Corrected : public T {
public:
    template<typename... Args>
    FH_Corrected(size_t FH_order, Args&&... args)
        : T(std::forward<Args>(args)...), FH_order{FH_order},
        D_factors(T::Ncomps, vector1d(T::Ncomps)),
        current_sigma_eff(T::Ncomps, vector1d(T::Ncomps)),
        current_eps_eff(T::Ncomps, vector1d(T::Ncomps))
    {
        set_D_factors();
    }

    double potential(int i, int j, double r, double temp){
        set_temperature(temp);
        return potential(i, j, r);
    }

    double potential_derivative_r(int i, int j, double r, double temp){
        set_temperature(temp);
        return potential_derivative_r(i, j, r);
    }

    double potential_dblderivative_rr(int i, int j, double r, double temp){
        set_temperature(temp);
        return potential_dblderivative_rr(i, j, r);
    }

    double potential_dn(int i, int j, double r, double temp, int n){
        set_temperature(temp);
        return potential_dn(i, j, r, n);
    }

    using T::potential;
    double potential(int i, int j, double r) override {
        double u = T::potential(i, j, r);
        for (size_t n = 1; n <= FH_order; n++){
            double nfac = 1;
            for (size_t ni = 2; ni <= n; ni++) nfac *= ni;
            u += pow(D_factors[i][j] / (BOLTZMANN * current_T), n) * T::potential_dn(i, j, r, 2 * n) / nfac;
        }
        return u;
    }

    using T::potential_derivative_r;
    double potential_derivative_r(int i, int j, double r) override {
        double ur = 0;
        for (size_t n = 0; n <= FH_order; n++){
            double nfac = 1;
            for (size_t ni = 2; ni <= n; ni++) nfac *= ni;
            ur += pow(D_factors[i][j] / (BOLTZMANN * current_T), n) * T::potential_dn(i, j, r, 2 * n + 1) / nfac;
        }
        return ur;
    }

    using T::potential_dblderivative_rr;
    double potential_dblderivative_rr(int i, int j, double r) override {
        double urr = 0;
        for (size_t n = 0; n <= FH_order; n++){
            double nfac = 1;
            for (size_t ni = 2; ni <= n; ni++) nfac *= ni;
            urr += pow(D_factors[i][j] / (BOLTZMANN * current_T), n) * T::potential_dn(i, j, r, 2 * n + 2) / nfac;
        }
        return urr;
    }

    using T::potential_dn;
    double potential_dn(int i, int j, double r, size_t n) override {
        double un = T::potential_dn(i, j, r, n);
        for (size_t k = 1; k <= FH_order; k++){
            double kfac = 1;
            for (size_t ki = 2; ki <= k; ki++) kfac *= ki;
            un += pow(D_factors[i][j] / (BOLTZMANN * current_T), k) * T::potential_dn(i, j, r, 2 * k + n) / kfac;
        }
        return un;
    }

    using T::get_sigma_eff;
    double get_sigma_eff(int i, int j, double temp) override {
        set_temperature(temp);
        return current_sigma_eff[i][j];
    }

    using T::get_eps_eff;
    double get_eps_eff(int i, int j, double temp) override {
        set_temperature(temp);
        return current_eps_eff[i][j];
    }

    inline size_t set_temperature(double temp) {
        return set_internals(0, temp, {0.});
    }

    size_t set_internals(double rho, double temp, const vector1d& x) override {
        size_t r = T::set_internals(rho, temp, x);
        if (temp != current_T){
            current_T = temp;
            set_effective_params();
            return r + 1;
        }
        return r;
    }

    void set_masses() override {
        T::set_masses();
        set_D_factors();
    }

    void set_D_factors() {
        for (size_t i = 0; i < T::Ncomps; i++){
            for (size_t j = 0; j < T::Ncomps; j++){
                D_factors[i][j] = pow(HBAR, 2) / (24 * T::red_mass[i][j]);
            }
        }
    }

    void set_FH_order(size_t order) {
        if (order != FH_order){
            T::clear_all_caches();
            FH_order = order;
        }
    }

    void set_effective_params(){
        update_sigma_eff();
        update_eps_eff();
    }

    void update_sigma_eff(){
        for (size_t i = 0; i < T::Ncomps; i++){
            for (size_t j = i; j < T::Ncomps; j++){
                double r0 = 3.5 * T::sigma[i][j];
                while (potential(i, j, r0) < 0){
                    r0 -= 0.2 * T::sigma[i][j];
                }
                current_sigma_eff[i][j] 
                    = current_sigma_eff[j][i] 
                    = newton([&](double r){return potential(i, j, r) / T::eps[i][j];},
                             [&](double r){return potential_derivative_r(i, j, r) / T::eps[i][j];},
                             r0);
            }
        }
    }

    void update_eps_eff(){
        for (size_t i = 0; i < T::Ncomps; i++){
            for (size_t j = i; j < T::Ncomps; j++){
                double r0 = newton([&](double r){return potential_derivative_r(i, j, r) * T::sigma[i][j] / T::eps[i][j];},
                                   [&](double r){return potential_dblderivative_rr(i, j, r) * T::sigma[i][j] / T::eps[i][j];},
                                   current_sigma_eff[i][j]);
                current_eps_eff[i][j] = current_eps_eff[j][i] = - potential(i, j, r0);
            }
        }
    }

private:
    double current_T;
    size_t FH_order;
    vector2d D_factors;

    vector2d current_sigma_eff;
    vector2d current_eps_eff;
};