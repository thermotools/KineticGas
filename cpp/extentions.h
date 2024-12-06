#pragma once
#include "utils.h"

template<typename T, size_t N, size_t deg>
class Tabulated : public T {
public:

    template<typename... Args>
    Tabulated(Args&&... args)
        : T(std::forward<Args>(args)...)
    {
        for (size_t i = 0; i < T::Ncomps; i++){
            tab_potential.push_back(std::vector<FuncTable<N, deg>>());
            tab_potential_r.push_back(std::vector<FuncTable<N, deg>>());
            tab_potential_rr.push_back(std::vector<FuncTable<N, deg>>());
            for (size_t j = i; j < T::Ncomps; j++){
                double tab_min = 0.5 * T::sigma[i][j];
                double tab_mid = T::sigma[i][j];
                double tab_mid_n = 100.;
                double tab_max = (tab_mid - tab_min) * (N / tab_mid_n) + tab_min;
                tab_potential[i].push_back(FuncTable<N, deg>([&](double r){return T::potential(i, j, r);}, tab_min, tab_max));
                tab_potential_r[i].push_back(FuncTable<N, deg>([&](double r){return T::potential_derivative_r(i, j, r);}, tab_min, tab_max));
                tab_potential_rr[i].push_back(FuncTable<N, deg>([&](double r){return T::potential_dblderivative_rr(i, j, r);}, tab_min, tab_max));
            }
        }
    }

    inline void set_using_tabulated(bool use_tab){using_tabulated = use_tab;};

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
