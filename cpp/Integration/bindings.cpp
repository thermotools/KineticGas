#include "Integration.h"
#include "integrator_tests.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/operators.h"
#include "pybind11/functional.h"

namespace py = pybind11;

#ifndef DEBUG
PYBIND11_MODULE(Integration_r, handle){
#else
PYBIND11_MODULE(Integration_d, handle){
#endif
    handle.doc() = "Integration module";
    handle.def("get_plane", &get_plane);
    handle.def("get_line", &get_line);
    handle.def("integrate_plane", &integrate_plane_py);
    handle.def("integrate2d", &integrate2d);
    handle.def("mesh2d", &mesh2d);
    handle.def("trisurf", &trisurf);

    handle.def("integrator_test", &integrator_test);
    handle.def("integrator_test_linear", &integrator_test_linear);
    handle.def("mesh_test", &mesh_test);
    handle.def("exp_testfunc", &testfun);
    handle.def("lin_testfunc", &testfun_linear);

    py::class_<Point>(handle, "Point")
        .def(py::init<double, double, double>())
        .def(py::init<double, double>())
        .def_readonly("x", &Point::x)
        .def_readonly("y", &Point::y)
        .def_readonly("z", &Point::z);

    py::class_<Plane>(handle, "Plane")
        .def(py::init<double, double, double>())
        .def_readonly("A", &Plane::A)
        .def_readonly("B", &Plane::B)
        .def_readonly("C", &Plane::C);

    py::class_<Line>(handle, "Line")
        .def(py::init<double, double>())
        .def_readonly("a", &Line::a)
        .def_readonly("b", &Line::b);
}