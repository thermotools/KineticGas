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
