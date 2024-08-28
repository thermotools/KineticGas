#include "KineticGas.h"
#include "Factorial.h"
#include "MieKinGas.h"
#include "QuantumMie.h"
#include "HardSphere.h"
#include "PseudoHardSphere.h"
#include "Sutherland.h"
#include "ModTangToennis.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/operators.h"
#include <sstream>

namespace py = pybind11;
using vector1d = std::vector<double>;
using vector2d = std::vector<vector1d>;

#define KineticGas_bindings(Model) \
        .def("get_conductivity_vector", &Model::get_conductivity_vector) \
        .def("get_diffusion_vector", &Model::get_diffusion_vector) \
        .def("get_diffusion_matrix", &Model::get_diffusion_matrix) \
        .def("get_conductivity_matrix", &Model::get_conductivity_matrix) \
        .def("get_viscosity_matrix", &Model::get_viscosity_matrix)\
        .def("get_viscosity_vector", &Model::get_viscosity_vector)\
        .def("get_mtl", &Model::get_mtl) \
        .def("get_etl", &Model::get_etl) \
        .def("get_rdf", &Model::get_rdf) \
        .def("get_K_factors", &Model::get_K_factors) \
        .def("get_K_prime_factors", &Model::get_K_prime_factors) \
        .def("set_eos", &Model::set_eos) \
        \
        .def("thermal_conductivity", &Model::thermal_conductivity) \
        .def("thermal_conductivity_tp", &Model::thermal_conductivity_tp) \


#define Spherical_potential_bindings(Model) \
        .def("potential", py::overload_cast<int, int, double>(&Model::potential)) \
        .def("potential_derivative_r", &Model::potential_derivative_r) \
        .def("potential_dblderivative_rr", &Model::potential_dblderivative_rr) \

#define Spherical_bindings(Model) \
        .def("omega", &Model::omega) \
        .def("chi", &Model::chi) \
        .def("get_R", &Model::get_R) \
        .def("theta_r", py::overload_cast<int, int, double, double, double, double>(&Model::theta_r)) \
        .def("set_tl_model", &Sutherland::set_transfer_length_model)


PYBIND11_MODULE(libpykingas, handle){
    handle.doc() = "Is this documentation? This is documentation.";
    handle.def("ipow", &ipow);
    handle.def("factorial_tests", &factorial_tests);

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

    py::class_<Sutherland>(handle, "cpp_Sutherland")
        .def(py::init<vector1d, vector2d, vector2d, vector3d, vector3d, bool, bool>())
        KineticGas_bindings(Sutherland)
        Spherical_potential_bindings(Sutherland)
        Spherical_bindings(Sutherland)
        .def("set_active_LJ_rdf", &Sutherland::set_active_LJ_rdf)
        .def("get_sigma_eff", &Sutherland::get_sigma_eff)
        .def("get_sigma_min", &Sutherland::get_sigma_min)
        .def("get_epsilon_eff", &Sutherland::get_epsilon_eff)
        .def("get_vdw_alpha", &Sutherland::get_vdw_alpha)
        .def("get_BH_diameters", &Sutherland::get_BH_diameters)
        .def("rdf_g0", py::overload_cast<double, double, const vector1d&>(&Sutherland::rdf_g0_func))
        .def("rdf_g1", py::overload_cast<double, double, const vector1d&>(&Sutherland::rdf_g1_func))
        .def("rdf_g2", py::overload_cast<double, double, const vector1d&, bool>(&Sutherland::rdf_g2_func))
        ;
        // Functions below this comment are only exposed for testing purposes
        // .def("da1_drho", &Sutherland::da1_drho_func)
        // .def("da1s_drho", &Sutherland::da1s_drho_func)
        // .def("dzeta_eff_drho", py::overload_cast<double, const vector1d&, const vector2d&, double>(&Sutherland::dzeta_eff_drho_func))
        // .def("zeta_x", &Sutherland::zeta_x_func)
        // .def("a1s", py::overload_cast<double, double, const vector1d&, const vector2d&>(&Sutherland::a_1s_func))
        // .def("B_func", py::overload_cast<double, const vector1d&, const vector2d&, const vector2d&>(&Sutherland::B_func))

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
        Spherical_bindings(MieKinGas)
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

   py::class_<QuantumMie>(handle, "cpp_QuantumMie")
        .def(py::init<vector1d, vector2d, vector2d, vector2d, vector2d, std::vector<int>, bool, bool>())
        KineticGas_bindings(QuantumMie)
        Spherical_bindings(QuantumMie)
        .def("omega_tester", &QuantumMie::omega_tester)
        .def("potential", py::overload_cast<int, int, double, double>(&QuantumMie::potential))
        .def("potential_derivative_r", py::overload_cast<int, int, double, double>(&QuantumMie::potential_derivative_r))
        .def("potential_dblderivative_rr", py::overload_cast<int, int, double, double>(&QuantumMie::potential_dblderivative_rr))
        .def("get_sigma_eff", py::overload_cast<double>(&QuantumMie::get_sigma_eff))
        .def("get_sigma_min", py::overload_cast<double>(&QuantumMie::get_sigma_min))
        .def("get_epsilon_eff", py::overload_cast<double>(&QuantumMie::get_epsilon_eff))
        .def("get_BH_diameters", &QuantumMie::get_BH_diameters)
        .def("saft_rdf", &QuantumMie::saft_rdf)
        ;

    py::class_<TangToennisParam>(handle, "cpp_TangToennisParam")
        .def(py::init<double, double, double, vector1d,
                        double, double, double, double,
                        vector1d>()
             )
        ;

    py::class_<ModTangToennis>(handle, "cpp_ModTangToennis")
        .def(py::init<TangToennisParam,
                        vector1d,
                        vector2d,
                        bool>()
             )
        KineticGas_bindings(ModTangToennis)
        .def("potential", py::overload_cast<int, int, double>(&ModTangToennis::potential))
        .def("potential_r", &ModTangToennis::potential_derivative_r)
        .def("potential_rr", &ModTangToennis::potential_dblderivative_rr)
        .def("set_tl_model", &ModTangToennis::set_transfer_length_model)
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
        .def("w_integral", &HardSphere::w_integral);
    
    py::class_<PseudoHardSphere>(handle, "cpp_PseudoHardSphere")
        .def(py::init<
                        vector1d,
                        vector2d,
                        bool, bool
                    >()
            )
        KineticGas_bindings(PseudoHardSphere)
        Spherical_bindings(PseudoHardSphere);
}