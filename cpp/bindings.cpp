#include "KineticGas.h"
#include "Factorial.h"
#include "MieKinGas.h"
#include "QuantumMie.h"
#include "HardSphere.h"
#include "PseudoHardSphere.h"
#include "multiparam.h"
#include "Quantum.h"
#include "LJSpline.h"
#include "LJTS.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/operators.h"
#include "pybind11/eigen.h"
#include "pybind11/numpy.h"
#include <sstream>

namespace py = pybind11;
using vector1d = std::vector<double>;
using vector2d = std::vector<vector1d>;

PYBIND11_MODULE(libpykingas, handle){
    handle.doc() = "Is this documentation? This is documentation.";
    handle.def("ipow", &ipow);
    handle.def("factorial_tests", &factorial_tests);

    handle.def("get_partitions", &get_partitions);
    handle.def("partition_multiplicity", &partition_multiplicity);

    handle.def("simpson", [](double ulim, int N){return simpson([](double x){return exp(- x) * pow(x, 2);}, 0, ulim, N);});
    handle.def("simpson_inf", [](double ulim){return simpson_inf([](double x){return exp(- x) * pow(x, 2);}, 0, ulim);});

    py::class_<Product>(handle, "Product")
        .def(py::init<int>())
        .def(py::init<double>())
        .def(py::init<Fac>())
        .def("eval", &Product::eval);

    py::class_<Fac>(handle, "Fac")
        .def(py::init<int>())
        .def("eval", &Fac::eval);

    py::class_<Polynomial>(handle, "Polynomial")
        .def(py::init<int, int, vector1d, int>())
        .def("derivative", py::vectorize(&Polynomial::derivative))
        ;

    py::class_<PolyExp>(handle, "PolyExp")
        .def(py::init<Polynomial, Polynomial>())
        .def("derivative", py::vectorize(&PolyExp::derivative))
        ;

    py::class_<OmegaPoint>(handle, "OmegaPoint")
        .def_readwrite("i", &OmegaPoint::i)
        .def_readwrite("j", &OmegaPoint::j)
        .def_readwrite("l", &OmegaPoint::l)
        .def_readwrite("r", &OmegaPoint::r)
        .def_readwrite("T_dK", &OmegaPoint::T_dK)
        .def("__repr__",
        [](const OmegaPoint &p) {
            std::stringstream strm;
            strm << "ij = " << p.i << p.j << ", l = " << p.l << ", r = " << p.r << " T = " << p.T_dK << " dK";
            return strm.str();
        });
    
    py::class_<Units>(handle, "cppUnits")
        .def(py::init<double, double, double>())
        .def_readonly("T", &Units::T)           // Temperature (K)
        .def_readonly("L", &Units::L)           // Length (m)
        .def_readonly("m", &Units::m)           // Mass (kg)
        .def_readonly("E", &Units::E)           // Energy (J)
        .def_readonly("V", &Units::V)           // Volume (m3)
        .def_readonly("t", &Units::t)           // Time (s)
        .def_readonly("F", &Units::F)           // Force (N)
        .def_readonly("speed", &Units::speed)   // Speed (m / s)
        .def_readonly("rho", &Units::rho)       // Density (mol / m3)
        .def_readonly("D", &Units::D)           // Diffusion (m^2 / s)
        .def_readonly("p", &Units::p)           // Pressure (Pa)
        .def_readonly("visc", &Units::visc)     // Shear viscosity (Pa s)
        .def_readonly("kvisc", &Units::kvisc)   // Kinematic viscosity (m^2 / s)
        .def_readonly("tdiff", &Units::tdiff)   // Thermal diffusivity (m^2 / s)
        .def_readonly("tcond", &Units::tcond)   // Thermal conductivity
        ;
    
    py::class_<KineticGas>(handle, "cpp_KineticGas")
        .def("viscosity", &KineticGas::viscosity)
        .def("viscosity_tp", &KineticGas::viscosity_tp)
        .def("kinematic_viscosity", &KineticGas::kinematic_viscosity)
        .def("kinematic_viscosity_tp", &KineticGas::kinematic_viscosity_tp)
        .def("thermal_conductivity", &KineticGas::thermal_conductivity)
        .def("thermal_conductivity_contributions", &KineticGas::thermal_conductivity_contributions)
        .def("thermal_conductivity_tp", &KineticGas::thermal_conductivity_tp)
        .def("thermal_diffusivity", &KineticGas::thermal_diffusivity)
        .def("thermal_diffusivity_tp", &KineticGas::thermal_diffusivity_tp)
        .def("soret_coefficient", &KineticGas::soret_coefficient)
        .def("soret_coefficient_tp", &KineticGas::soret_coefficient_tp)
        .def("interdiffusion", &KineticGas::interdiffusion)
        .def("interdiffusion_tp", &KineticGas::interdiffusion_tp)
        .def("thermal_diffusion_ratio", &KineticGas::thermal_diffusion_ratio)
        .def("thermal_diffusion_coeff", &KineticGas::thermal_diffusion_coeff)
        .def("thermal_diffusion_factor", &KineticGas::thermal_diffusion_factor)
        .def("second_virial", py::overload_cast<int, int, double>(&KineticGas::second_virial))
        .def("bound_second_virial", &KineticGas::bound_second_virial)
        .def("dimer_constant", &KineticGas::dimer_constant)
        .def("get_mtl", &KineticGas::get_mtl)
        .def("get_etl", &KineticGas::get_etl)
        .def("get_rdf", &KineticGas::get_rdf)
        .def("set_eos", py::overload_cast<py::object>(&KineticGas::set_eos))
        .def("frame_of_reference_map", &KineticGas::frame_of_reference_map)
        .def("get_reducing_units", &KineticGas::get_reducing_units)
        .def("de_broglie_wavelength", py::overload_cast<int, double>(&KineticGas::de_broglie_wavelength))

        .def("set_tl_model", &KineticGas::set_transfer_length_model)
        .def("get_tl_model", &KineticGas::get_transfer_length_model)
        .def("get_valid_tl_model", &KineticGas::get_valid_transfer_length_models)
        
        .def("get_conductivity_vector", &KineticGas::get_conductivity_vector)
        .def("get_diffusion_vector", &KineticGas::get_diffusion_vector)
        .def("get_diffusion_matrix", &KineticGas::get_diffusion_matrix)
        .def("get_conductivity_matrix", &KineticGas::get_conductivity_matrix)
        .def("get_viscosity_matrix", &KineticGas::get_viscosity_matrix)
        .def("get_viscosity_vector", &KineticGas::get_viscosity_vector)
        .def("get_bulk_viscosity_matrix", &KineticGas::get_bulk_viscosity_matrix)
        .def("get_bulk_viscosity_vector", &KineticGas::get_bulk_viscosity_vector)
        .def("get_K_factors", &KineticGas::get_K_factors)
        .def("get_K_prime_factors", &KineticGas::get_K_prime_factors)
        .def("get_chemical_potential_factors", &KineticGas::get_chemical_potential_factors)
        .def("get_ksi_factors", &KineticGas::get_ksi_factors)
        
        .def("thermal_conductivity", &KineticGas::thermal_conductivity)
        .def("thermal_conductivity_tp", &KineticGas::thermal_conductivity_tp)
        .def("viscosity", &KineticGas::viscosity)
        .def("interdiffusion", &KineticGas::interdiffusion)
        .def("selfdiffusion", &KineticGas::selfdiffusion)
        ;

    py::class_<Spherical, KineticGas>(handle, "cpp_Spherical")
        .def("potential", py::overload_cast<int, int, double>(&Spherical::potential, py::const_))
        .def("potential_derivative_r", &Spherical::potential_derivative_r)
        .def("potential_dblderivative_rr", &Spherical::potential_dblderivative_rr)
        .def("potential_r", &Spherical::potential_derivative_r)
        .def("potential_rr", &Spherical::potential_dblderivative_rr)

        .def("get_r_min", &Spherical::get_r_min)
        .def("get_sigma_eff", &Spherical::get_sigma_eff)
        .def("get_eps_eff", &Spherical::get_eps_eff)
        .def("get_alpha_eff", &Spherical::get_alpha_eff)

        .def("omega_tester", &Spherical::omega_tester)
        .def("w_integral_tester", &Spherical::w_integral_tester)
        .def("w_integrand", &Spherical::w_integrand)
        ;

    py::class_<ExtSutherland, Spherical>(handle, "cpp_ExtSutherland")
        .def(py::init<vector1d, vector2d, vector2d, 
                        vector3d, vector3d, vector3d, vector3d,
                        bool, bool>())
        .def("potential", py::overload_cast<int, int, double>(&ExtSutherland::potential, py::const_))
        .def("potential_derivative_r", py::overload_cast<int, int, double>(&ExtSutherland::potential_derivative_r, py::const_))
        .def("potential_dblderivative_rr", py::overload_cast<int, int, double>(&ExtSutherland::potential_dblderivative_rr, py::const_))
        .def("saft_rdf", &ExtSutherland::saft_rdf)
        .def("get_rdf_terms", &ExtSutherland::get_rdf_terms)
        .def("get_sigma_eff", &ExtSutherland::get_sigma_eff)
        .def("get_rmin", &ExtSutherland::get_sigma_min)
        .def("get_epsilon_eff", &ExtSutherland::get_epsilon_eff)
        .def("get_dBH", py::overload_cast<double, double>(&ExtSutherland::get_BH_diameters))
        .def("get_vdw_alpha", &ExtSutherland::get_vdw_alpha)
        ;

    py::class_<MieKinGas, Spherical>(handle, "cpp_MieKinGas")
        .def(py::init<vector1d,
                      vector2d,
                      vector2d,
                      vector2d,
                      vector2d,
                      bool,
                      bool
                    >()
            )
        .def(py::init<std::string, bool>())
        .def("set_omega_correlation_active", &MieKinGas::set_omega_correlation_active)
        .def("get_BH_diameters", &MieKinGas::get_BH_diameters)
        .def("get_vdw_alpha", &MieKinGas::get_vdw_alpha)
        .def("saft_rdf", &MieKinGas::saft_rdf)
        .def("rdf_g0", py::overload_cast<double, double, const vector1d&>(&MieKinGas::rdf_HS))
        .def("rdf_g1", py::overload_cast<double, double, const vector1d&>(&MieKinGas::rdf_g1_func))
        .def("rdf_g2", py::overload_cast<double, double, const vector1d&, bool>(&MieKinGas::rdf_g2_func))
        .def("a1", &MieKinGas::a1ij_func)
        .def("da1drho", py::overload_cast<double, double, const vector1d&>(&MieKinGas::da1ij_drho_func))
        .def("a1s", py::overload_cast<double, double, const vector1d&, const vector2d&>(&MieKinGas::a_1s_func))
        .def("da1s_drho", py::overload_cast<double, double, const vector1d&, const vector2d&>(&MieKinGas::da1s_drho_func))
        .def("B_func", py::overload_cast<double, double, const vector1d&, const vector2d&>(&MieKinGas::B_func))
        .def("dBdrho_func", py::overload_cast<double, double, const vector1d&, const vector2d&>(&MieKinGas::dBdrho_func))
        .def("a2", py::overload_cast<double, double, const vector1d&>(&MieKinGas::a2ij_func))
        .def("a2_div_chi", &MieKinGas::a2ij_div_chi_func)
        .def("da2_div_chi_drho", py::overload_cast<double, double, const vector1d&>(&MieKinGas::da2ij_div_chi_drho_func))

        .def("chi", &MieKinGas::chi)
        .def("get_R", &MieKinGas::get_R)
        .def("theta_r", py::overload_cast<int, int, double, double, double, double>(&MieKinGas::theta_r))
        .def("theta_r", py::overload_cast<int, int, double, double, double, double, double>(&MieKinGas::theta_r))
        ;
        // .def("da1_drho", py::overload_cast<double, double, const vector1d&>(&MieKinGas::da1ij_drho_func))
        // .def("da1s_drho", py::overload_cast<double, const vector1d&, const vector2d&, const vector2d&>(&MieKinGas::da1s_drho_func))
        // .def("zeta_eff", &MieKinGas::zeta_eff_func)
        // .def("dzeta_eff_drho", &MieKinGas::dzeta_eff_drho_func)
        // .def("zeta_x", &MieKinGas::zeta_x_func)
        // .def("a1s", py::overload_cast<double, double, const vector1d&, const vector2d&>(&MieKinGas::a_1s_func))
        // .def("B_func", py::overload_cast<double, const vector1d&, const vector2d&, const vector2d&>(&MieKinGas::B_func))
        // .def("get_dBH", &MieKinGas::get_BH_diameters)
        // .def("a1ij", &MieKinGas::a1ij_func)
        // .def("a2ij",  py::overload_cast<double, double, const vector1d&>(&MieKinGas::a2ij_func))
        // .def("da2ij_drho", py::overload_cast<double, double, const vector1d&>(&MieKinGas::da2ij_drho_func))
        // .def("da2_div_chi_drho", py::overload_cast<double, double, const vector1d&>(&MieKinGas::da2ij_div_chi_drho_func))
        // .def("gamma_corr", py::overload_cast<double, double, const vector1d&>(&MieKinGas::gamma_corr))
        // ;

    py::class_<QuantumMie, Spherical>(handle, "cpp_QuantumMie")
        .def(py::init<vector1d, vector2d, vector2d, vector2d, vector2d, std::vector<int>, bool, bool>())
        .def("potential", py::overload_cast<int, int, double, double>(&QuantumMie::potential, py::const_))
        .def("get_sigma_eff", py::overload_cast<double>(&QuantumMie::get_sigma_eff))
        .def("get_sigma_min", py::overload_cast<double>(&QuantumMie::get_sigma_min))
        .def("get_epsilon_eff", py::overload_cast<double>(&QuantumMie::get_epsilon_eff))
        .def("saft_rdf", &QuantumMie::saft_rdf)
        .def("get_rdf_terms", &QuantumMie::get_rdf_terms)
        ;

   py::class_<IntegrationParam>(handle, "IntegrationParam")
        .def(py::init<vector1d, vector1d, double, double, int, int, double>())
        .def("set_end", &IntegrationParam::set_end)
        .def("set_dg", &IntegrationParam::set_dg)
        .def("set_db", &IntegrationParam::set_db)
        .def("set_rlg", &IntegrationParam::set_rlg)
        .def("set_rlb", &IntegrationParam::set_rlb)
        .def("set_dd_lim", &IntegrationParam::set_dd_lim)
        ;

    py::class_<HardSphere, KineticGas>(handle, "cpp_HardSphere")
        .def(py::init<
                        vector1d,
                        vector2d,
                        bool, bool
                    >()
            )
        .def("chi", &HardSphere::chi)
        .def("omega", &HardSphere::omega)
        .def("w_integral", &HardSphere::w_integral)
        .def("cross_section", &HardSphere::cross_section)
        ;
    
    py::class_<Quantum, Spherical>(handle, "cpp_Quantum")
        .def("get_E_bound", &Quantum::get_E_bound)
        .def("cross_section", &Quantum::cross_section)
        .def("classical_cross_section", &Quantum::classical_cross_section)
        .def("reduced_cross_section", &Quantum::reduced_cross_section)
        .def("wave_function", &Quantum::wave_function)
        .def("phase_shift", &Quantum::phase_shift)
        .def("JKWB_phase_shift", &Quantum::JKWB_phase_shift)
        .def("quantum_phase_shift", &Quantum::quantum_phase_shift)
        .def("absolute_phase_shift", py::overload_cast<int, int, int, double, double>(&Quantum::absolute_phase_shift))
        .def("absolute_phase_shift", py::overload_cast<int, int, int, double>(&Quantum::absolute_phase_shift))
        .def("absolute_phase_shifts", &Quantum::absolute_phase_shifts)
        .def("get_levinson_r", &Quantum::get_levinson_r)
        .def("integral_phase_shift", &Quantum::integral_phase_shift)
        .def("total_phase_shifts", &Quantum::total_phase_shifts)
        .def("dump_phase_shift_map", &Quantum::dump_phase_shift_map)
        .def("clear_phase_shift_maps", &Quantum::clear_phase_shift_maps)
        .def("dump_phase_shift_map_to_json", &Quantum::dump_phase_shift_map_to_json)
        .def("r_classical_forbidden", &Quantum::r_classical_forbidden)
        .def("get_de_boer", py::overload_cast<>(&Quantum::get_de_boer))
        .def("get_de_boer", py::overload_cast<int, int>(&Quantum::get_de_boer))
        .def("get_de_boer", py::overload_cast<int>(&Quantum::get_de_boer))
        .def("set_de_boer_mass", &Quantum::set_de_boer_mass)
        .def("JKWB_upper_E_limit", &Quantum::JKWB_upper_E_limit)
        .def("omega", &Quantum::omega)
        .def("quantum_omega", &Quantum::quantum_omega)
        .def("second_virial", py::overload_cast<int, int, double>(&Quantum::second_virial))
        .def("second_virial", py::overload_cast<int, int, const vector1d&>(&Quantum::second_virial))
        .def("second_virial_contribs", &Quantum::second_virial_contribs)
        .def("bound_second_virial_lim", &Quantum::bound_second_virial_lim)
        .def("semiclassical_second_virial", &Quantum::semiclassical_second_virial)
        .def("classical_omega", &Quantum::classical_omega)
        .def("scattering_volume", &Quantum::scattering_volume)
        .def("partial_scattering_volume", &Quantum::partial_scattering_volume)
        .def("set_quantum_active", &Quantum::set_quantum_active)
        .def("get_quantum_active", &Quantum::get_quantum_active)
        .def("set_JKWB_limits", &Quantum::set_JKWB_limits)
        .def("get_JKWB_limits", &Quantum::get_JKWB_limits)
        .def("omega_tester", &Quantum::omega_tester)
        .def("w_integrand", &Quantum::w_integrand)
        .def("theta", &Quantum::theta)
        .def("get_R", &Quantum::get_R)
        ;
    
    py::class_<TangToennisParam>(handle, "cpp_TangToennisParam")
        .def_readwrite("sigma", &TangToennisParam::sigma)
        .def_readwrite("eps_div_k", &TangToennisParam::eps_div_k)
        .def_readwrite("Re", &TangToennisParam::Re)
        .def("__repr__",
        [](const TangToennisParam &t) {
            std::stringstream strm;
            strm << "TangToennisParam\n"
                 << "\tA       : " << t.A << "\n"
                 << "\tb       : " << t.b << "\n"
                 << "\tA_tilde : " << t.A_tilde << "\n"
                 << "\ta_tilde : " << t.a_tilde << "\n"
                 << "\ta       : [ " << t.am2 << ", " << t.am1 << ", " << t.a1 << ", " << t.a2 << " ]\n"
                 << "\tC       : [ ";
            for (size_t i = 0; i < 5; i++){
                strm << t.C[i] << ", ";
            }
            strm << t.C[5] << " ]\n\t-----\n"
                 << "\tsigma   : " << t.sigma << "\n"
                 << "\teps_div_k : " << t.eps_div_k << "\n"
                 << "\tRe      : " << t.Re << std::endl;
            return strm.str();
        })
        ;

    py::class_<ModTangToennis, Quantum>(handle, "cpp_ModTangToennis")
        .def(py::init<std::string, std::string>())
        .def("potential_dn", &ModTangToennis::potential_dn)
        ;
    
    py::class_<FH_ModTangToennies, Quantum>(handle, "cpp_FH_ModTangToennies")
        .def(py::init<std::string, size_t, std::string>())
        .def("potential", py::overload_cast<int, int, double, double>(&FH_ModTangToennies::potential))
        .def("potential_r", py::overload_cast<int, int, double, double>(&FH_ModTangToennies::potential_derivative_r))
        .def("potential_rr", py::overload_cast<int, int, double, double>(&FH_ModTangToennies::potential_dblderivative_rr))
        .def("set_FH_order", &FH_ModTangToennies::set_FH_order)
        ;

    py::class_<HFD_B2, Quantum>(handle, "cpp_HFD_B2")
        .def(py::init<std::string>())
        .def("potential_dn", &HFD_B2::potential_dn)
        ;
    
    py::class_<FH_HFD_B2, Quantum>(handle, "cpp_FH_HFD_B2")
        .def(py::init<std::string, size_t>())
        .def("potential", py::overload_cast<int, int, double, double>(&FH_HFD_B2::potential))
        .def("potential_r", py::overload_cast<int, int, double, double>(&FH_HFD_B2::potential_derivative_r))
        .def("potential_rr", py::overload_cast<int, int, double, double>(&FH_HFD_B2::potential_dblderivative_rr))
        .def("set_FH_order", &FH_HFD_B2::set_FH_order)
        ;

    py::class_<PatowskiParam>(handle, "cpp_PatowskiParam")
        .def_readonly("sigma", &PatowskiParam::sigma)
        .def_readonly("eps_div_k", &PatowskiParam::eps_div_k)
        .def_readonly("r_min", &PatowskiParam::r_min)
        .def("__repr__",
        [](const PatowskiParam &t) {
            std::stringstream strm;
            strm << "PatowskiParam\n"
                 << "\tCex       : " << t.Cex1 << ", " << t.Cex2 << "\n"
                 << "\tCsp       : " << t.Csp1 << ", " << t.Csp2 << ", " << t.Csp3 << ", " << t.Csp4 << "\n"
                 << "\tCn        : " << t.C6 << ", " << t.C8 << ", " << t.C10 << "\n"
                 << "\t----------\n"
                 << "\tsigma     : " << t.sigma << "\n"
                 << "\teps_div_k : " << t.eps_div_k << "\n"
                 << "\tr_min     : " << t.r_min << std::endl;
            return strm.str();
        })
        ;

    py::class_<Patowski, Quantum>(handle, "cpp_Patowski")
        .def(py::init<std::string>())
        .def("get_param", &Patowski::get_param)
        .def("potential_dn", &Patowski::potential_dn)
        ;

    py::class_<PatowskiFH, Quantum>(handle, "cpp_PatowskiFH")
        .def(py::init<std::string, size_t>())
        .def("potential", py::overload_cast<int, int, double, double>(&PatowskiFH::potential))
        .def("potential_r", py::overload_cast<int, int, double, double>(&PatowskiFH::potential_derivative_r))
        .def("potential_rr", py::overload_cast<int, int, double, double>(&PatowskiFH::potential_dblderivative_rr))
        .def("set_FH_order", &PatowskiFH::set_FH_order)
        .def("get_param", &PatowskiFH::get_param)
        ;

    py::class_<PseudoHardSphere, Spherical>(handle, "cpp_PseudoHardSphere")
        .def(py::init<
                        vector1d,
                        vector2d,
                        bool, bool
                    >()
            )
        ;

    py::class_<LJSpline, Spherical>(handle, "cpp_LJSpline")
        .def(py::init<bool, bool>()
            )
        .def("omega", &LJSpline::omega)
        .def("omega_star", &LJSpline::omega_star)
        .def("omega_star_approx", &LJSpline::omega_star_approx)
        ;

    py::class_<LJTS, Spherical>(handle, "cpp_LJTS")
        .def(py::init<
                        vector1d,
                        vector2d,
                        vector2d,
                        bool, bool
                    >()
            )
        ;
}