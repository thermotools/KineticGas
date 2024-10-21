/*
Contains: 
    GenericEoS : Interface describing the requirements of an EoS that is held by a KineticGas object.
    ThermoWrapper : Thin wrapper for the thermopack Thermo class that conforms to the GenericEoS interface.
    PyWrapper : Wrapper for python-side EoS objects with signatures equivalent to those used in ThermoPack v2
*/
#pragma once
#include <iostream>
#include <functional>
#include <memory>
#include <cppThermopack/thermo.h>

/*
KineticGas objects hold an internal std::unique_ptr<GenericEoS>, which is used to compute various properties that we need an EoS for.

If you are implementing or wrapping an EoS to use with KineticGas, note that not all methods must be implemented in order to use a given method
in KineticGas. Refer to the list below for an overview of what transport properties that require a given method.

    * Cp_ideal - Ideal gas isobaric heat capacity : Used for thermal conductivity of multi-atomic molecules
    * dmudn - Derivative of chemical potential wrt. mole numbers at constant (T, V) : Used for diffusion coefficients
    * specific_volume - Molar volume given T, p : Used for TP-interface methods (compute volume, feed call to the corresponding TV-interface method)
    * pressure_tv - Pressure : Used for Center of Volume diffusion coefficients and thermal diffusion coefficients
    * dvdn - Partial molar volume : Used for Center of Volume diffusion coefficients and thermal diffusion coefficients
    * Cp_real - Isobaric heat capacity : Used for thermal diffusivity

The below list indicates which methods must be implemented in an EoS to compute a given property. If a method is not listed, it can be
used without the EoS implementing anything at all.

    * thermal_conductivity : Cp_ideal
    * interdiffusion : dmudn
    * thermal_diffusion_coeff : dmudn
    * thermal_diffusion_ratio : dmudn
    * thermal_diffusion_factor : dmudn
    * interdiffusion / thermal_diffusion_coeff (with frame_of_reference=CoV) : dmudn, pressure_tv, dvdn
    * thermal_diffusivity : Cp_real
    * <method>_tp : specific_volume

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
#define Cp_signature double Cp_real(double t, double V, const std::vector<double> n) const
#define pressure_signature double pressure_tv(double t, double V, const std::vector<double> n) const
#define volume_signature double specific_volume(double t, double p, const std::vector<double> n, int phase) const
#define dvdn_signature std::vector<double> dvdn(double t, double p, const std::vector<double> n, int phase) const

class GenericEoS {
public:
    int VAPPH;
    template<typename T>
    GenericEoS(T eos_) : eos(std::make_unique<WrappedEoS<T>>(std::move(eos_))) {VAPPH = eos->VAPPH;}
    
    GenericEoS(GenericEoS&& other){
        eos.swap(other.eos);
        VAPPH = other.VAPPH;
    }
    /*
    All inputs are:
        t : Temperature (K)
        V : Volume (m^33)
        p : Pressure (Pa)
        n : Mole numbers (mol)
        phase : Phase flag (this->VAPPH should indicate vapour)
    */
    // Derivative of chemical potentials wrt. mole numbers at constant (T, V), symmetric (J / mol^2)
    dmudn_signature {return eos->dmudn(t, V, n);}

    // Ideal gas heat capacity of component c (J / mol K)
    Cpid_signature {return eos->Cp_ideal(t, c);}

    // Isobaric heat capacity (J / mol K)
    Cp_signature {return eos->Cp_real(t, V, n);}

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
        virtual Cp_signature = 0;
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
        inline Cp_signature override {return wrapped.Cp_real(t, V, n);}
        inline pressure_signature override {return wrapped.pressure_tv(t, V, n);}
        inline volume_signature override {return wrapped.specific_volume(t, p, n, phase);}
        inline dvdn_signature override {return wrapped.dvdn(t, p, n, phase);}
    };

    std::unique_ptr<const EoS_Interface> eos;
};

class ThermoWrapper{
    public:
    int VAPPH;
    ThermoWrapper(Thermo&& eos_) : eos(std::make_unique<Thermo>(std::move(eos_))) {VAPPH = eos->VAPPH;}

    inline std::vector<std::vector<double>> dmudn(double T, double V, const std::vector<double> n) const {return eos->chemical_potential_tv(T, V, n, false, false, true).dn();}
    inline double Cp_ideal(double T, int ci) const {return eos->idealenthalpysingle(T, ci, true).dt();}
    inline double Cp_real(double T, double V, std::vector<double> n) const {return eos->enthalpy_tvp(T, V, n, true).dt();}
    inline double pressure_tv(double T, double V, const std::vector<double> n) const {return eos->pressure_tv(T, V, n);}
    inline double specific_volume(double t, double p, const std::vector<double> n, int phase) const {return eos->specific_volume(t, p, n, phase);}
    inline std::vector<double> dvdn(double t, double p, const std::vector<double> n, int phase) const {return eos->specific_volume(t, p, n, phase, false, false, true).dn();}

    private:
    std::unique_ptr<Thermo> eos;
};

#ifdef PYLIB
#include <pybind11/pybind11.h>
namespace py = pybind11;

class PyWrapper{
    public:
    int VAPPH;
    PyWrapper(py::object eos) : VAPPH{eos.attr("VAPPH").cast<int>()}, eos{eos} {
        std::vector<std::string> required_attr = {"chemical_potential_tv", "idealenthalpysingle", "enthalpy_tvp", "pressure_tv", "specific_volume"};
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
    inline double Cp_real(double T, double V, std::vector<double> n) const {
        py::tuple rt = eos.attr("enthalpy_tvp")(T, V, n, true);
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