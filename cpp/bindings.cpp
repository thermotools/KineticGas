#include "KineticGas.h"
#include "Factorial.h"
#include "MieKinGas.h"
#include "HardSphere.h"
#include "PseudoHardSphere.h"
#include "multiparam.h"
#include "Quantum.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/operators.h"
#include "pybind11/eigen.h"
#include <sstream>

namespace py = pybind11;
using vector1d = std::vector<double>;
using vector2d = std::vector<vector1d>;

#define KineticGas_bindings(Model) \
        .def("viscosity", &Model::viscosity) \
        .def("viscosity_tp", &Model::viscosity_tp) \
        .def("kinematic_viscosity", &Model::kinematic_viscosity) \
        .def("kinematic_viscosity_tp", &Model::kinematic_viscosity_tp) \
        .def("thermal_conductivity", &Model::thermal_conductivity) \
        .def("thermal_conductivity_contributions", &Model::thermal_conductivity_contributions) \
        .def("thermal_conductivity_tp", &Model::thermal_conductivity_tp) \
        .def("thermal_diffusivity", &Model::thermal_diffusivity) \
        .def("thermal_diffusivity_tp", &Model::thermal_diffusivity_tp) \
        .def("soret_coefficient", &Model::soret_coefficient) \
        .def("soret_coefficient_tp", &Model::soret_coefficient_tp) \
        .def("interdiffusion", &Model::interdiffusion) \
        .def("interdiffusion_tp", &Model::interdiffusion_tp) \
        .def("thermal_diffusion_ratio", &Model::thermal_diffusion_ratio) \
        .def("thermal_diffusion_coeff", &Model::thermal_diffusion_coeff) \
        .def("thermal_diffusion_factor", &Model::thermal_diffusion_factor) \
        .def("get_mtl", &Model::get_mtl) \
        .def("get_etl", &Model::get_etl) \
        .def("get_rdf", &Model::get_rdf) \
        .def("set_eos", py::overload_cast<py::object>(&Model::set_eos)) \
        .def("frame_of_reference_map", &Model::frame_of_reference_map) \
        .def("get_reducing_units", &Model::get_reducing_units) \
        \
        .def("set_tl_model", &Model::set_transfer_length_model) \
        .def("get_tl_model", &Model::get_transfer_length_model) \
        .def("get_valid_tl_models", &Model::get_valid_transfer_length_models) \
        \
        .def("get_conductivity_vector", &Model::get_conductivity_vector) \
        .def("get_diffusion_vector", &Model::get_diffusion_vector) \
        .def("get_diffusion_matrix", &Model::get_diffusion_matrix) \
        .def("get_conductivity_matrix", &Model::get_conductivity_matrix) \
        .def("get_viscosity_matrix", &Model::get_viscosity_matrix)\
        .def("get_viscosity_vector", &Model::get_viscosity_vector)\
        .def("get_K_factors", &Model::get_K_factors) \
        .def("get_K_prime_factors", &Model::get_K_prime_factors) \
        .def("get_chemical_potential_factors", &Model::get_chemical_potential_factors) \
        .def("get_ksi_factors", &Model::get_ksi_factors)
        \
        
#define Spherical_potential_bindings(Model) \
        .def("potential", py::overload_cast<int, int, double>(&Model::potential)) \
        .def("potential_derivative_r", &Model::potential_derivative_r) \
        .def("potential_dblderivative_rr", &Model::potential_dblderivative_rr) \
        .def("potential_r", &Model::potential_derivative_r) \
        .def("potential_rr", &Model::potential_dblderivative_rr) \
        .def("get_reducing_units", &Model::get_reducing_units) \

PYBIND11_MODULE(libpykingas, handle){
    handle.doc() = "Is this documentation? This is documentation.";
    handle.def("ipow", &ipow);
    handle.def("factorial_tests", &factorial_tests);

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

    py::class_<MieKinGas>(handle, "cpp_MieKinGas")
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
        KineticGas_bindings(MieKinGas)
        Spherical_potential_bindings(MieKinGas)
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

   py::class_<IntegrationParam>(handle, "IntegrationParam")
        .def(py::init<vector1d, vector1d, double, double, int, int, double>())
        .def("set_end", &IntegrationParam::set_end)
        .def("set_dg", &IntegrationParam::set_dg)
        .def("set_db", &IntegrationParam::set_db)
        .def("set_rlg", &IntegrationParam::set_rlg)
        .def("set_rlb", &IntegrationParam::set_rlb)
        .def("set_dd_lim", &IntegrationParam::set_dd_lim)
        ;

    py::class_<TangToennisParam>(handle, "cpp_TangToennisParam")
        .def(py::init<double, double, double, vector1d,
                        double, double, double, double,
                        vector1d>()
             )
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

    py::class_<ModTangToennis>(handle, "cpp_ModTangToennis")
        .def(py::init<TangToennisParam, vector1d, bool>())
        KineticGas_bindings(ModTangToennis)
        .def("potential", py::overload_cast<int, int, double>(&ModTangToennis::potential))
        .def("potential_r", &ModTangToennis::potential_derivative_r)
        .def("potential_rr", &ModTangToennis::potential_dblderivative_rr)
        ;

    py::class_<HardSphere>(handle, "cpp_HardSphere")
        .def(py::init<
                        vector1d,
                        vector2d,
                        bool, bool
                    >()
            )
        KineticGas_bindings(HardSphere)
        .def("chi", &HardSphere::chi)
        .def("omega", &HardSphere::omega)
        .def("w_integral", &HardSphere::w_integral)
        .def("cross_section", &HardSphere::cross_section)
        ;
    
    py::class_<Quantum>(handle, "cpp_Quantum")
        // .def(py::init<std::string>())
        KineticGas_bindings(Quantum)
        Spherical_potential_bindings(Quantum)
        .def("cross_section", &Quantum::cross_section)
        .def("classical_cross_section", &Quantum::classical_cross_section)
        .def("reduced_cross_section", &Quantum::reduced_cross_section)
        .def("wave_function", &Quantum::wave_function)
        .def("phase_shift", &Quantum::phase_shift)
        .def("JKWB_phase_shift", &Quantum::JKWB_phase_shift)
        .def("get_de_boer", py::overload_cast<>(&Quantum::get_de_boer))
        .def("get_de_boer", py::overload_cast<int, int>(&Quantum::get_de_boer))
        .def("get_de_boer", py::overload_cast<int>(&Quantum::get_de_boer))
        .def("set_de_boer_mass", &Quantum::set_de_boer_mass)
        .def("de_broglie_wavelength", &Quantum::de_broglie_wavelength)
        .def("JKWB_upper_E_limit", &Quantum::JKWB_upper_E_limit)
        .def("omega", &Quantum::omega)
        .def("quantum_omega", &Quantum::quantum_omega)
        .def("classical_omega", &Quantum::classical_omega)
        .def("set_quantum_active", &Quantum::set_quantum_active)
        .def("get_quantum_active", &Quantum::get_quantum_active)
        .def("set_JKWB_limits", &Quantum::set_JKWB_limits)
        .def("get_JKWB_limits", &Quantum::get_JKWB_limits)
        ;
    
    py::class_<HFD_B2, Quantum>(handle, "cpp_HFD_B2")
        .def(py::init<std::string>())
        KineticGas_bindings(HFD_B2)
        Spherical_potential_bindings(HFD_B2)
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
        KineticGas_bindings(Patowski)
        Spherical_potential_bindings(Patowski)
        .def("get_param", &Patowski::get_param)
        ;

    py::class_<PseudoHardSphere>(handle, "cpp_PseudoHardSphere")
        .def(py::init<
                        vector1d,
                        vector2d,
                        bool, bool
                    >()
            )
        KineticGas_bindings(PseudoHardSphere)
        ;
}