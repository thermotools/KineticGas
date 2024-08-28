/*
Contains: 
    GenericEoS : Interface describing the requirements of an EoS that is held by a KineticGas object.
    ThermoWrapper : Thin wrapper for the thermopack Thermo class that conforms to the GenericEoS interface.
*/
#pragma once
#include <iostream>
#include <functional>
#include <memory>

/*
KineticGas objects hold an internal std::unique_ptr<GenericEoS>, which is used to compute various properties that we need an EoS for.

The GenericEoS class is written such that any class implementing methods with signatures equivalent to the public methods of GenericEoS
can be held by an std::unique_ptr<GenericEoS>. Note that an "equivalent" signature also includes signatures for which implicit conversions
exist which handle neccessary conversions.

Also note: Upon constructing a GenericEoS from some other type, that type is likely to be moved, such that the recommended construction pattern is
std::unique_ptr<GenericEoS> eos = std::make_unique<GenericEoS>(MyEosClass(arg1, arg2, ...)),
to ensure that the moved-from object is not handled at any point.

Finally: The signatures of the required methods are defined as macros (sorry) because they are reused several times in the wrapping class (sorry)
        If anyone has a better solution for this I'm all ears.
*/

#define dmudn_signature std::vector<std::vector<double>> dmudn(double t, double V, const std::vector<double> n) const
#define Cpid_signature double Cp_ideal(double t, int c) const
#define pressure_signature double pressure_tv(double t, double V, const std::vector<double> n) const
#define volume_signature double specific_volume(double t, double p, const std::vector<double> n, int phase) const
#define dvdn_signature std::vector<double> dvdn(double t, double p, const std::vector<double> n, int phase) const

class GenericEoS {
public:
    int VAPPH;
    template<typename T>
    GenericEoS(T eos_) : eos(std::make_unique<WrappedEoS<T>>(std::move(eos_))) {VAPPH = eos->VAPPH;}

    /*
    All inputs are:
        t : Temperature (K)
        V : Volume (m^33)
        p : Pressure (Pa)
        n : Mole numbers (mol)
    */
    // Derivative of chemical potentials wrt. mole numbers at constant (T, V), symmetric (J / mol^2)
    dmudn_signature {return eos->dmudn(t, V, n);}

    // Ideal gas heat capacity of component c (J / mol)
    Cpid_signature {return eos->Cp_ideal(t, c);}

    // Pressure (Pa)
    pressure_signature {return eos->pressure_tv(t, V, n);}

    // Specific volume (m^3 / mol). When (phase == this->VAPPH), this method should return the vapour phase specific volume.
    volume_signature {return eos->specific_volume(t, p, n, phase);}
    
    // Derivative of specific volume wrt. mole numbers at constant (T, p) (m^3 / mol^2). When (phase == this->VAPPH), this method should return the derivative of the vapour phase specific volume.
    dvdn_signature {return eos->dvdn(t, p, n, phase);}

private:
    struct EoS_Interface {
        int VAPPH;
        virtual ~EoS_Interface() = default;
        virtual dmudn_signature = 0;
        virtual Cpid_signature = 0;
        virtual pressure_signature = 0;
        virtual volume_signature = 0;
        virtual dvdn_signature = 0;
    };

    template<typename T>
    struct WrappedEoS : EoS_Interface {
        T wrapped;
        WrappedEoS(T eos) : wrapped(std::move(eos)) {VAPPH = wrapped.VAPPH;}

        inline dmudn_signature override {return wrapped.dmudn(t, V, n);}
        inline Cpid_signature override {return wrapped.Cp_ideal(t, c);}
        inline pressure_signature override {return wrapped.pressure_tv(t, V, n);}
        inline volume_signature override {return wrapped.specific_volume(t, p, n, phase);}
        inline dvdn_signature override {return wrapped.dvdn(t, p, n, phase);}
    };

    std::unique_ptr<const EoS_Interface> eos;
};

// #ifdef NOPYTHON
#include <cppThermoPack/thermo.h>
class ThermoWrapper{
    public:
    int VAPPH;
    ThermoWrapper(Thermo&& eos_) : eos(std::make_unique<Thermo>(std::move(eos_))) {VAPPH = eos->VAPPH;}

    inline std::vector<std::vector<double>> dmudn(double T, double V, const std::vector<double> n) const {return eos->chemical_potential_tv(T, V, n, false, false, true).dn();}
    inline double Cp_ideal(double T, int ci) const {return eos->idealenthalpysingle(T, ci, true).dt();}
    inline double pressure_tv(double T, double V, const std::vector<double> n) const {return eos->pressure_tv(T, V, n);}
    inline double specific_volume(double t, double p, const std::vector<double> n, int phase) const {return eos->specific_volume(t, p, n, phase);}
    inline std::vector<double> dvdn(double t, double p, const std::vector<double> n, int phase) const {return eos->specific_volume(t, p, n, phase, false, false, true).dn();}

    private:
    std::unique_ptr<Thermo> eos;
};
// #endif

#ifndef NOPYTHON
#include <pybind11/pybind11.h>
namespace py = pybind11;

class PyWrapper{
    public:
    int VAPPH;
    PyWrapper(py::object eos) : VAPPH{eos.attr("VAPPH").cast<int>()}, eos{eos} {
        std::vector<std::string> required_attr = {"chemical_potential_tv", "idealenthalpysingle", "pressure_tv", "specific_volume"};
        for (const auto& attr : required_attr){
            py::object tmp = eos.attr(attr.data()); // Trigger a python-side AttributeError if any of the required attributes are missing
        }
    }

    inline std::vector<std::vector<double>> dmudn(double T, double V, const std::vector<double> n) const {
        py::tuple rt = eos.attr("chemical_potential_tv")(T, V, n, false, false, true);
        std::vector<std::vector<double>> r = rt[1].cast<std::vector<std::vector<double>>>();
        return r;
    }
    inline double Cp_ideal(double T, int ci) const {
        py::tuple rt = eos.attr("idealenthalpysingle")(T, ci, true);
        return rt[1].cast<double>();
    }
    inline double pressure_tv(double T, double V, const std::vector<double> n) const {
        py::tuple rt = eos.attr("pressure_tv")(T, V, n);
        return rt[0].cast<double>();
    }
    inline double specific_volume(double T, double p, const std::vector<double> n, int phase) const {
        py::tuple rt = eos.attr("specific_volume")(T, p, n, phase);
        return rt[0].cast<double>();
    }
    inline std::vector<double> dvdn(double T, double p, const std::vector<double> n, int phase) const {
        py::tuple rt = eos.attr("specific_volume")(T, p, n, phase, false, false, true);
        return rt[1].cast<std::vector<double>>();
    }

    private:
    py::object eos;
};

#endif