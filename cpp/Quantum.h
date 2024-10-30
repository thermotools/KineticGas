#pragma once
#include "Spherical.h"

class Quantum : public Spherical {
public:
    Quantum(std::string comps);

    using Spherical::potential;
    dual2 potential(int i, int j, dual2 r) override;
    double potential(int i, int j, double r) override;
    double potential_derivative_r(int i, int j, double r);

    vector2d wave_function(int i, int j, int l, double E, double r_end, double dr);
    double JKWB_phase_shift(int i, int j, int l, double E);
    double phase_shift(int i, int j, int l, double E);
    double cross_section_A(int n, int l, size_t k);
    double cross_section_kernel(int i, int j, double n, double l, double E);
    double cross_section(int i, int j, int n, double E);
    double quantum_omega(int i, int j, int n, int s, double T);

    std::vector<size_t> half_spin; // Spin of each particle multiplied by two 
protected:
    vector2d model_rdf(double rho, double T, const vector1d& mole_fracs) override {throw std::runtime_error("Method model_rdf not implemented for Quantum!");}
private:
    

};