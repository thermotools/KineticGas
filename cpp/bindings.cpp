#include "KineticGas.h"
#include "Factorial.h"
#include "MieKinGas.h"
#include "HardSphere.h"
#include "PseudoHardSphere.h"
#include "Sutherland.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/operators.h"
#include <sstream>

namespace py = pybind11;

#define KineticGas_bindings(Model) \
        .def("A", &Model::A) \
        .def("A_prime", &Model::A_prime) \
        \
        .def("H_ij", &Model::H_ij) \
        .def("H_i", &Model::H_i) \
        .def("get_conductivity_vector", &Model::get_conductivity_vector) \
        .def("get_diffusion_vector", &Model::get_diffusion_vector) \
        .def("get_diffusion_matrix", &Model::get_diffusion_matrix) \
        .def("get_conductivity_matrix", &Model::get_conductivity_matrix) \
        .def("get_viscosity_matrix", &Model::get_viscosity_matrix)\
        .def("get_viscosity_vector", &Model::get_viscosity_vector)\
        \
        .def("B_prime", &Model::B_prime) \
        .def("B_dblprime", &Model::B_dblprime) \
        \
        .def("L_ij", &Model::L_ij) \
        .def("L_i", &Model::L_i) \
        \
        .def("get_rdf", &Model::get_rdf) \
        .def("get_K_factors", &Model::get_K_factors) \
        .def("get_K_prime_factors", &Model::get_K_prime_factors)\
        .def("get_contact_diameters", &Model::get_contact_diameters) \
        \
        .def_readwrite("omega_map", &Model::omega_map)

#define Spherical_potential_bindings(Model) \
        .def("potential", &Model::potential) \
        .def("potential_derivative_r", &Model::potential_derivative_r) \
        .def("potential_dblderivative_rr", &Model::potential_dblderivative_rr) \

#define Spherical_bindings(Model) \
        .def("chi", &Model::chi) \
        .def("get_R", &Model::get_R) \
        .def("omega", &Model::omega) \
        \
        .def("get_R_rootfunc", &Model::get_R_rootfunc) \
        .def("get_R_rootfunc_derivative", &Model::get_R_rootfunc_derivative) \
        \
        .def("theta", &Model::theta) \
        .def("theta_lim", &Model::theta_lim) \
        .def("theta_integral", &Model::theta_integral) \
        .def("theta_integrand", &Model::theta_integrand) \
        .def("theta_integrand_dblderivative", &Model::theta_integrand_dblderivative) \
         \
        .def("w_integrand", &Model::w_integrand) \
        .def("w_integral", &Model::w_integral)


#ifndef DEBUG
PYBIND11_MODULE(KineticGas_r, handle){
#else
PYBIND11_MODULE(KineticGas_d, handle){
#endif
    handle.doc() = "Is this documentation? This is documentation.";
    handle.def("ipow", &ipow);
    handle.def("factorial_tests", &factorial_tests);

    // handle.def("zeta_x", &mie_rdf::zeta_x_func);
    // handle.def("dzeta_x_drho", &mie_rdf::dzetax_drho_func);
    // handle.def("zeta_eff", &mie_rdf::zeta_eff_func);
    // handle.def("dzeta_eff_drho", &mie_rdf::dzeta_eff_drho_func);

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
        .def(py::init<vector1d, vector2d, vector2d, vector3d, vector3d, bool>())
        KineticGas_bindings(Sutherland)
        Spherical_potential_bindings(Sutherland)
        Spherical_bindings(Sutherland)
        .def("get_sigma_eff", &Sutherland::get_sigma_eff)
        .def("get_sigma_min", &Sutherland::get_sigma_min)
        .def("get_epsilon_eff", &Sutherland::get_epsilon_eff)
        .def("get_BH_diameters", &Sutherland::get_BH_diameters)
        .def("rdf_g0", py::overload_cast<double, double, const vector1d&>(&Sutherland::rdf_g0_func))
        .def("rdf_g1", py::overload_cast<double, double, const vector1d&>(&Sutherland::rdf_g1_func))
        .def("da1_drho", &Sutherland::da1_drho_func)
        .def("da1s_drho", &Sutherland::da1s_drho_func)
        .def("zeta_eff", py::overload_cast<double, const vector1d&, const vector2d&, double>(&Sutherland::zeta_eff_func))
        .def("dzeta_eff_drho", py::overload_cast<double, const vector1d&, const vector2d&, double>(&Sutherland::dzeta_eff_drho_func))
        .def("zeta_x", &Sutherland::zeta_x_func)
        .def("a1s", py::overload_cast<double, double, const vector1d&, const vector2d&>(&Sutherland::a_1s_func))
        .def("B_func", py::overload_cast<double, const vector1d&, const vector2d&, const vector2d&>(&Sutherland::B_func))
        // Functions above this comment have been tested to reproduce MieKinGas
        
        ;
        // .def("rdf_g1", &Sutherland::rdf_g1_func);

    py::class_<MieKinGas>(handle, "cpp_MieKinGas")
        .def(py::init<vector1d,
                      vector2d,
                      vector2d,
                      vector2d,
                      vector2d,
                      bool
                    >()
            )
        KineticGas_bindings(MieKinGas)
        Spherical_potential_bindings(MieKinGas)
        Spherical_bindings(MieKinGas)
        .def("get_BH_diameters", &MieKinGas::get_BH_diameters)
        .def("rdf_g0", &MieKinGas::rdf_HS)
        .def("rdf_g1", py::overload_cast<double, double, const vector1d&>(&MieKinGas::rdf_g1_func))
        .def("da1_drho", py::overload_cast<double, double, const vector1d&>(&MieKinGas::da1ij_drho_func))
        .def("da1s_drho", py::overload_cast<double, const vector1d&, const vector2d&, const vector2d&>(&MieKinGas::da1s_drho_func))
        .def("zeta_eff", &MieKinGas::zeta_eff_func)
        .def("dzeta_eff_drho", &MieKinGas::dzeta_eff_drho_func)
        .def("zeta_x", &MieKinGas::zeta_x_func)
        .def("a1s", py::overload_cast<double, double, const vector1d&, const vector2d&>(&MieKinGas::a_1s_func))
        .def("B_func", py::overload_cast<double, const vector1d&, const vector2d&, const vector2d&>(&MieKinGas::B_func))

        .def("g2", py::overload_cast<double, double, const vector1d&>(&MieKinGas::rdf_g2_func))
        .def("get_dBH", &MieKinGas::get_BH_diameters)
        .def("a1ij", &MieKinGas::a1ij_func)
        .def("a2ij",  py::overload_cast<double, double, const vector1d&>(&MieKinGas::a2ij_func))
        .def("da2ij_drho", py::overload_cast<double, double, const vector1d&>(&MieKinGas::da2ij_drho_func));
        
    
    py::class_<HardSphere>(handle, "cpp_HardSphere")
        .def(py::init<
                        vector1d,
                        vector2d,
                        bool
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
                        bool
                    >()
            )
        KineticGas_bindings(PseudoHardSphere)
        Spherical_bindings(PseudoHardSphere);
}