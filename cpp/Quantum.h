/*
Author: Vegard G. Jervell
Purpose: Interface to quatum mechanical calulations 

Overview:
The `Quantum` class inherits from Spherical, and implements quantum mechanical evaluation of collision integrals and second
virial coefficients, as well as some utility functionality that is useful when working with quantum systems.

For a full overview of references used for the implementations, see the source file `Quantum.cpp`

- Any sub-class of `Quantum` must implement the interaction potential and derivatives.
- Sub-classes of `Quantum` can activate/deactivate quantum mechanical calculations by using the `set_quantum_active` method.
- If quantum mechanical calculations are deactivated, this class will simply forward calculations to `Spherical`.
*/
#pragma once
#include "Spherical.h"
#include <mutex>
#include <shared_mutex>

class Quantum : public Spherical {
public:
    Quantum(std::string comps);

    // Reroute to one of the below based on whether quantum calculations are active
    double omega(int i, int j, int l, int r, double T) override;

    double quantum_omega(int i, int j, int n, int s, double T);
    double classical_omega(int i, int j, int l, int r, double T);

    // Compute the second virial coefficient, or ceertain contributions to it.
    double second_virial(int i, int j, double T) override;
    vector1d second_virial(int i, int j, const vector1d& T);
    // Compute several contributions: (i, t, b) : (ideal, thermal, bound)
    std::map<char, double> second_virial_contribs(int i, int j, double T, const std::string& contribs) override;
    double semiclassical_second_virial(int i, int j, double T); // Uses JKWB approximation (I think)
    double bound_second_virial(int i, int j, double T) override;
    double bound_second_virial_lim(int i, int j) override; // High temperature limit of the bound contribution
    double dimer_constant(int i, int j, double T) override;

    using Spherical::potential;

    /*
        Utility methods related to the de Boer parameter.
        de Broglie Wavelengths are implemented in `KineticGas`
    */
    vector2d get_de_boer();
    double get_de_boer(int i, int j);
    double get_de_boer(int i){return get_de_boer(i, i);}
    void set_de_boer_mass(int i, double de_boer); // Set particle masses to obtain specified de Boer parameter

    // Set the upper energy limit for when the JKWB approximation is used.
    double JKWB_upper_E_limit(int i, int j);

    /*
        Phase-shift calculations

        This is what makes up the core of the quantum mechanical calculations. The real leg-work of this class is in these methods.

        Transport properties rely only on the `phase_shift`, while second virial coefficients need the absolute phase shift. See
        "The Limits of the Feynman-Hibbs corrections ..." and related references for more details.
    */

    // Solve the Schrödinger equation for the pair (i, j), with angular momentum quantum number l, and energy E (> 0).
    // Terminating at r_end, with step size (dr)
    // The starting point for the integration point is determined by the function itself.
    // Args:
    //     r_end : Dimensionless end position (r / sigma)
    //     dr : Dimensionless Step size (dr / sigma) 
    //     E : Dimensionless energy (E / eps[i][j])
    // Returns: 
    //  [r, psi, phase_shift] : The positions, wave function values, and local phase shifts
    vector2d wave_function(int i, int j, int l, double E, double r_end, double dr);

    // Compute the relative phase shift
    // Either using the JKWB approximation, without the JKWB approximation, or conditionally using the JKWB approximation based on the set limits
    double JKWB_phase_shift(int i, int j, int l, double E);
    double quantum_phase_shift(int i, int j, int l, double E, double& r_lev);
    double phase_shift(int i, int j, int l, double E);

    // Compute the abolute phase shifts...
    // When the absolute phase shift at the previous step is known (relatively cheap)
    double absolute_phase_shift(int i, int j, int l, double E, double prev_delta);

    // When the absolute phase shift at the previous step is not known (much more expensive)
    double absolute_phase_shift(int i, int j, int l, double E);

    // Compute all phase shifts up to the relative wave-vector k_max
    vector2d absolute_phase_shifts(int i, int j, int l, double k_max);

    // These methods do the leg-work
    // They store the computed results in the internal "phase_shift_map"
    void trace_absolute_phase_shifts(int i, int j, int l, double k_max);
    void fill_absolute_phase_shifts(int i, int j, int l, double next_k, int& n, vector1d& k_vals, vector1d& phase_shifts, vector1d& node_count);
    void fill_absolute_phase_shifts_tail(int i, int j, int l, double next_k, int& n, vector1d& k_vals, vector1d& phase_shifts);

    // Get the "Levinson multiple" at E = 0 (count number of vibrational states with angular momentum l)
    int get_levinson_multiple(int i, int j, int l);

    // The "Levinson-radius" is the radial position of the outermost root of the wave function that is inside the classically forbidden region.
    // We track this position to determine whether there are "resonances" in the phase shifts
    vector1d get_levinson_r(int i, int j, int l, const vector1d k_vals);

    // "Total" phase shifts are the sum of absolute phase shifts over the angular momentum quantum number.
    // Result is stored in the internal "total_phase_shift_map" after tracing, and can be retrieved with the total_phase_shifts method.
    void trace_total_phase_shifts(int i, int j, double k_max);
    vector2d total_phase_shifts(int i, int j, double k_max);
    void dump_phase_shift_map(); // Dumps to a custom binary format. Use the below for json.
    void dump_phase_shift_map_to_json(); // Format: {l : [[...k], [...delta]]} or something like that. "k" is probably in "sigma^-1".
    void clear_phase_shift_maps();

    // The "integral phase shift" is the integral of the total phase shift, weighted with a Boltzmann factor.
    double integral_phase_shift(int i, int j, int l, double T);
    double r_classical_forbidden(int i, int j, int l, double E);

    // Compute the collision cross section
    double cross_section_A(int n, int l, size_t k);
    double cross_section_kernel(int i, int j, int n, int l, double E);
    double JKWB_cross_section(int i, int j, int n, double E);
    double cross_section(int i, int j, int n, double E) override;
    double classical_cross_section(int i, int j, int l, double E);

    double scattering_volume(int i, int j, double E);
    double partial_scattering_volume(int i, int j, double E, int l_max);

    void set_quantum_active(bool active);
    bool get_quantum_active(){return quantum_is_active;}

    void set_JKWB_limits(double E, int l){
        JKWB_E_limit = E; JKWB_l_limit = l;
        phase_shift_map.clear();
    }

    std::pair<double, int> get_JKWB_limits(){
        return std::pair<double, int>(JKWB_E_limit, JKWB_l_limit);
    }

    // Whether we have Fermi-Dirac, Bose-Einstein or Boltzmann statistics (uses an enum, see utils.h)
    int get_interaction_statistics(int i, int j);
    std::array<double, 2> get_symmetry_weights(int i, int j);

    vector2d get_E_bound_from_file(const std::string& comp);
    vector2d get_E_bound(int i, int j);

protected:
    vector2d model_rdf(double rho, double T, const vector1d& mole_fracs) override {throw std::runtime_error("Method model_rdf not implemented for Quantum!");}
    void clear_all_caches() override;

    // The bound state energies, organised as E_bound[v][l], where v is the vibrational quantum number, and l is the angular momentum quantum number
    // These are pre-computed values stored in the fluid files, or associated E_bound files.
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