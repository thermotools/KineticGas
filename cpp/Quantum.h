#pragma once
#include "Spherical.h"
#include <mutex>
#include <shared_mutex>

class Quantum : public Spherical {
public:
    Quantum(std::string comps);
    vector2d get_E_bound_from_file(const std::string& comp);

    double omega(int i, int j, int l, int r, double T) override;

    using Spherical::potential;

    vector2d get_de_boer();
    double get_de_boer(int i, int j);
    double get_de_boer(int i){return get_de_boer(i, i);}
    void set_de_boer_mass(int i, double de_boer); // Set particle masses to obtain specified de Boer parameter

    vector2d wave_function(int i, int j, int l, double E, double r_end, double dr);
    double JKWB_phase_shift(int i, int j, int l, double E);
    double JKWB_upper_E_limit(int i, int j);
    double quantum_phase_shift(int i, int j, int l, double E);
    double phase_shift(int i, int j, int l, double E);
    double absolute_phase_shift(int i, int j, int l, double E, double prev_delta);
    double absolute_phase_shift(int i, int j, int l, double E);
    vector2d absolute_phase_shifts(int i, int j, int l, double k_max);

    void trace_absolute_phase_shifts(int i, int j, int l, double k_max);
    void fill_absolute_phase_shifts(int i, int j, int l, double next_k, int& n, vector1d& k_vals, vector1d& phase_shifts);
    void fill_absolute_phase_shifts_tail(int i, int j, int l, double next_k, int& n, vector1d& k_vals, vector1d& phase_shifts);

    int get_levinson_multiple(int i, int j, int l);
    void trace_total_phase_shifts(int i, int j, double k_max);
    vector2d total_phase_shifts(int i, int j, double k_max);
    void dump_phase_shift_map();
    void clear_phase_shift_maps();

    double integral_phase_shift(int i, int j, int l, double T);
    double r_classical_forbidden(int i, int j, int l, double E);

    double cross_section_A(int n, int l, size_t k);
    double cross_section_kernel(int i, int j, int n, int l, double E);
    double JKWB_cross_section(int i, int j, int n, double E);
    double cross_section(int i, int j, int n, double E) override;
    double quantum_omega(int i, int j, int n, int s, double T);
    double classical_omega(int i, int j, int l, int r, double T);

    double classical_cross_section(int i, int j, int l, double E);

    double second_virial(int i, int j, double T) override;
    vector1d second_virial(int i, int j, const vector1d& T);
    double semiclassical_second_virial(int i, int j, double T);
    double bound_second_virial(int i, int j, double T) override;
    double dimer_constant(int i, int j, double T) override;

    double scattering_volume(int i, int j, double E);
    double partial_scattering_volume(int i, int j, double E, int l_max);
    std::map<char, double> second_virial_contribs(int i, int j, double T, const std::string& contribs) override;

    void set_quantum_active(bool active);
    bool get_quantum_active(){return quantum_is_active;}

    void set_JKWB_limits(double E, int l){
        JKWB_E_limit = E; JKWB_l_limit = l;
        phase_shift_map.clear();
    }

    std::pair<double, int> get_JKWB_limits(){
        return std::pair<double, int>(JKWB_E_limit, JKWB_l_limit);
    }

    int get_interaction_statistics(int i, int j);
    std::array<double, 2> get_symmetry_weights(int i, int j);

protected:
    vector2d model_rdf(double rho, double T, const vector1d& mole_fracs) override {throw std::runtime_error("Method model_rdf not implemented for Quantum!");}
    void clear_all_caches() override;

    std::vector<std::vector<vector2d>> E_bound;

private:
    bool quantum_supported = true;
    bool quantum_is_active = true;
    double JKWB_E_limit = 120.;
    int JKWB_l_limit = 1000;
    
    std::vector<unsigned int> half_spin; // Spin of each particle multiplied by two 
    std::vector<double> spin;
    std::vector<unsigned int> rot_ground_state; // Rotational ground state
    std::vector<std::vector<int>> interaction_statistics;
    std::map<std::pair<int, double>, double> phase_shift_map;
    std::map<int, vector2d> absolute_phase_shift_map;
    std::mutex abs_phase_shift_map_mutex;
    vector2d stored_total_phase_shifts;

    std::shared_mutex cross_section_map_mutex;
    std::unordered_map<CrossSectionPoint, double, CrossSectionHash> cross_section_map;

};