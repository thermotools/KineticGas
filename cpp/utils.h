/*
Various utility structures and enums

    * StatePoint and OmegaPoint: Use for lookup of collision integrals and transfer lengths

    * FrameOfReference: Enum for frames of reference
    
    * TransferLengthModel: Enum for transfer length models

    * Units: Struct to hold reducing units (i.e. Lennard-Jones units) of a species.
            See: KineticGas::get_reducing_units
    
    * StatisticType: Enum for Bose-Einstein/Fermi-Dirac/Boltzmann statistics

    * FuncTable: Callable class for creating a tabulated function from some single-variate function
                See: extensions.h::Tabulated for example usage.
    * FuncSpline: Same as FuncTable, but uses splines

    * KineticGasCache: There's currently a bunch of really messy caching going on. I'm in the process of isolating all this this class, so that management of thread-safe 
        access to the caches can be centralized. That will also make it a lot easier to manage a `const`/`mutable` structure, since we can just have a `mutable` cache

References:
    (I) The Kinetic Gas theory of Mie fluids, V. G. Jervell, Norwegian university of science and technology (2022)
    (II) Revised Enskog theory for Mie fluids: Prediction of diffusion coefficients, thermal diffusion coefficients,
    viscosities and thermal conductivities, V. G. Jervell and Ø. Wilhelmsen, J. Chem. Phys. (2023)
*/

#pragma once
#include "global_params.h"
#include "math.h"
#include <string>
#include <Eigen/Dense>
#include <array>
#include <stdexcept>
#include <shared_mutex>
#include <optional>
#include <map>


struct StatePoint{
    int T_dK;
    double rho;
    StatePoint(double T) : T_dK{static_cast<int>((T * 100.) + 0.5)}, rho{-1} {}
    StatePoint(double T, double rho) : T_dK{static_cast<int>((T * 100.) + 0.5)}, rho{rho}{}

    bool operator<(const StatePoint& other) const {
        if (T_dK < other.T_dK) return true;
        else if (T_dK == other.T_dK){
            if (rho < other.rho) return true;
        }
        return false;
    }
};

/*
   To avoid unneccesary evaluations of the collision integrals, this struct is used to represent a point in 
   The five-dimensional (i, j, l, r, T)-space where the collision integral has been evaluated.

   NB: Resolution along the T-axis is 0.1 K, as (by experience) the collision integrals are a weak enough
   function of T to justify using T-values rounded to the nearest .1 K to improve speed, with little cost to precision.
*/
struct OmegaPoint{
    int i, j, l, r, T_dK;
    double rho;
    OmegaPoint(int i, int j, int l, int r, double T, double rho) : i{i}, j{j}, l{l}, r{r}, rho{rho} {
         T_dK = (int) ((T * 10.0) + 0.5); // Temperature in dK (10^-1 K)
    };

    OmegaPoint(int i, int j, int l, int r, double T) : OmegaPoint(i, j, l, r, T, 0){}

    bool operator<(const OmegaPoint& other) const {
        if (i < other.i) return true;
        else if (i == other.i){
            if (j < other.j) return true;
            else if (j == other.j){
                if (l < other.l) return true;
                else if (l == other.l){
                    if (r < other.r) return true;
                    else if (r == other.r){
                        if (T_dK < other.T_dK) return true;
                        else if (T_dK == other.T_dK){
                            if (rho < other.rho) return true;
                        }
                    }
                }
            }
        }
        return false;
    }
};

struct CrossSectionPoint{
    int i, j, l, E;
    CrossSectionPoint(int i, int j, int l, double E)
        : i{i}, j{j}, l{l}, E{static_cast<int>(std::round(100 * E))} 
    {}

    bool operator==(const CrossSectionPoint& other) const{
        return (i == other.i) && (j == other.j) && (l == other.l) && (E == other.E);
    }
};

struct CrossSectionHash {
    size_t operator()(const CrossSectionPoint& s) const {
        return std::hash<int>()(s.i) ^ (std::hash<int>()(s.j) << 1) ^
               (std::hash<int>()(s.l) << 2) ^ (std::hash<double>()(s.E) << 3);
    }
};

class KineticGasCache {

    bool has_omega(const OmegaPoint& point) const;
    bool has_mtl(const StatePoint& stp) const;
    bool has_etl(const StatePoint& stp) const;

    std::optional<double> get_omega(const OmegaPoint& point) const;
    std::optional<vector2d> get_mtl(const StatePoint& stp) const;
    std::optional<vector2d> get_etl(const StatePoint& stp) const;

    void store_omega(const OmegaPoint& point, const double omega);
    void store_mtl(const StatePoint& stp, const vector2d& mtl);
    void store_etl(const StatePoint& stp, const vector2d& etl);

    void clear();

private:
    std::map<OmegaPoint, double> omega_map;
    mutable std::shared_mutex omega_mutex;

    std::map<StatePoint, vector2d> mtl_map;
    mutable std::shared_mutex mtl_mutex;

    std::map<StatePoint, vector2d> etl_map;
    mutable std::shared_mutex etl_mutex;
};

enum FrameOfReference{
    CoM,
    CoN,
    CoV,
    solvent,
    zarate,
    zarate_x,
    zarate_w
};

enum TransferLengthModel{
    DEFAULT = -1, // Default model
    collision_diameter, // Model presented in Refs. (I) and (II) (see top of file)
    EWCA, // Unpublished model for momentum- and energy transfer lengths (MTL and ETL)
    correlation, // Unpublished correlation for Argon transfer lengths
    INVALID // Used to terminate loops over all model identifiers
};

struct Units {
    const double 
     m      // Mass (kg)
    ,L      // Length (m)
    ,E      // Energy (J)
    ,T      // Temperature (K)
    ,V      // Volume (m3)
    ,t      // Time (s)
    ,F      // Force (N)
    ,speed  // Speed (m / s)
    ,rho    // Density (mol / m3)
    ,D      // Diffusion (m^2 / s)
    ,p      // Pressure (Pa)
    ,visc   // Shear viscosity (Pa s)
    ,kvisc  // Kinematic viscosity (m^2 / s)
    ,tdiff  // Thermal diffusivity (m^2 / s)
    ,tcond  // Thermal conductivity
    ;

    Units(double m_unit, double L_unit, double E_unit) 
        : m{m_unit}
        , L{L_unit}
        , E{E_unit}
        , T{E_unit / BOLTZMANN}
        , V{pow(L, 3)}
        , t{sqrt(m / E) * L}
        , F{E / L}
        , speed{L / t}
        , rho{1 / (AVOGADRO * V)}
        , D{pow(L, 2) / t}
        , p{E / V}
        , visc{p * t}
        , kvisc{D}
        , tdiff{kvisc}
        , tcond{E / (t * T * L)}
    {}
};

void set_fluid_dir(const std::string path);
std::string get_fluid_dir();

enum StatisticType {
    Boltzmann,
    BoseEinstein,
    FermiDirac
};

inline bool str_contains(const std::string& str, const std::string& sub){
    return (str.find(sub) != std::string::npos);
}

/*
    Create a tabulated version of a function

    Tabulation uses N points, and when the tabulated function is called, uses an interpolation of degree "deg" between the points.

    The interpolation degree "deg" must be odd.
*/
template<size_t N, size_t deg>
class FuncTable{
public:
    
    FuncTable(std::function<double(double)> fun, double x_min, double x_max)
        : x_min{x_min}, x_max{x_max}, dx{(x_max - x_min) / (N - 1)}
        {
            std::array<double, N> x_vals;
            std::array<double, N> f_vals;
            double x = x_min;
            for (size_t i = 0; i < N; i++){
                x_vals[i] = x;
                f_vals[i] = fun(x);
                x += dx;
            }
            for (size_t i = 0; i < N; i++){
                C[i] = get_interpolant(i, x_vals, f_vals);
            }
        }

    FuncTable() = default;
    FuncTable(const FuncTable<N, deg>&) = default;
    FuncTable(FuncTable<N, deg>&&) = default;
    FuncTable<N, deg>& operator=(const FuncTable<N, deg>&) = default;
    FuncTable<N, deg>& operator=(FuncTable<N, deg>&&) = default;


    double unsafe_eval(double x){
        size_t x0_idx = static_cast<size_t>((x - x_min) / dx);
        double v = 0.;
        for (size_t n = 0; n <= deg; n++){
            v += C[x0_idx][n] * pow(x, deg - n);
        }
        return v;
    }

    double operator()(double x){
        if (!(in_valid_range(x))) throw std::out_of_range("Invalid region for tabulated function!");
        return unsafe_eval(x);
    }

    inline bool in_valid_range(double x){return ((x > x_min + ((deg + 1) / 2) * dx) && (x + ((deg + 1) / 2) * dx < x_max));}

private:
    double x_min, x_max, dx;    
    static constexpr size_t N_intervals = N - 1;
    std::array<std::array<double, deg + 1>, N> C;

    std::array<double, deg + 1> get_interpolant(size_t x0_idx, const std::array<double, N>& x_vals, const std::array<double, N>& f_vals) {
        static_assert(deg % 2 == 1, "Interpolation degree must be odd.");
        constexpr size_t n_points = deg + 1;

        if (x0_idx < ((deg - 1) / 2)) {
            return get_interpolant((deg - 1) / 2, x_vals, f_vals);
        }
        else if (x0_idx > N - (n_points / 2)){
            return get_interpolant(N - (n_points / 2), x_vals, f_vals);
        }

        size_t start_idx = x0_idx - (deg - 1) / 2;

        Eigen::MatrixXd V(n_points, n_points);
        Eigen::VectorXd f(n_points);

        for (size_t i = 0; i < n_points; ++i) {
            size_t idx = start_idx + i; 
            double x = x_vals[idx];
            for (size_t j = 0; j < n_points; ++j) {
                V(i, j) = pow(x, n_points - j - 1); 
            }
            f(i) = f_vals[idx];
        }

        Eigen::VectorXd coeff;
        if (deg < 5){
            coeff = V.partialPivLu().solve(f);
        } 
        else{
            coeff = V.colPivHouseholderQr().solve(f);
        }
        std::array<double, deg + 1> coeff_arr;
        for (size_t i = 0; i < deg + 1; i++){
            coeff_arr[i] = coeff(i);
        }
        return coeff_arr;
    }

};

/*
    Create a segmented spline from a function

    Basically the FuncTable, except it uses N points, and creates splines of degree "deg".
*/
template<size_t N, size_t deg>
class FuncSpline {
public:
    FuncSpline(std::function<double(double)> fun, double x_min, double x_max)
        : x_min{x_min}, x_max{x_max}, dx{(x_max - x_min) / (N - 1)} 
    {
        std::array<double, N> x_vals;
        std::array<double, N> f_vals;
        // Generate tabulated points and function values
        for (size_t i = 0; i < N; ++i) {
            double x = x_min + i * dx;
            x_vals[i] = x;
            f_vals[i] = fun(x);
        }
        // Fit a spline to the tabulated data
        fit_spline(x_vals, f_vals);
    }

    FuncSpline() = default;
    FuncSpline(const FuncSpline&) = default;
    FuncSpline(FuncSpline&&) = default;
    FuncSpline<N, deg>& operator=(const FuncSpline&) = default;
    FuncSpline<N, deg>& operator=(FuncSpline&&) = default;

    double operator()(double x) const {
        // Check if x is in range
        if ( !(in_valid_range(x)) ) throw std::out_of_range("Input x is out of the interpolation range!");
        double val = 0.0;
        for (size_t i = 0; i < coefficients.size(); ++i) {
            val += coefficients[i] * b_spline_basis(x, i, deg);
        }
        return val;
    }

    // Method to evaluate the n'th derivative of the spline at point x
    double derivative(double x, int n) const {
        // Check if x is in range
        if (n < 0) throw std::out_of_range("Cannot compute derivative < 0");
        if (n > deg) throw std::out_of_range("Too high derivative order (" + std::to_string(n) + ") for FuncSpline with degree " + std::to_string(deg));
        if ( !(in_valid_range(x)) ) throw std::out_of_range("Input x is out of the spline range!");
        if (n == 0) return this->operator()(x);

        double val = 0.0;
        for (size_t i = 0; i < coefficients.size(); ++i) {
            val += coefficients[i] * b_spline_basis_derivative(x, i, deg, n);
        }
        return val;
    }

    bool in_valid_range(double x) const {
        return ( (x_min + deg * dx < x ) && (x < x_max - deg * dx) );
    }

private:
    double x_min, x_max, dx;

    static constexpr size_t num_knots = N + deg + 1;
    static constexpr size_t num_coefficients = N;

    std::array<double, num_knots> knots{};
    std::array<double, num_coefficients> coefficients{};

    // Fit the spline to the tabulated data
    void fit_spline(const std::array<double, N>& x_vals, const std::array<double, N>& f_vals) {
        // Construct the knot vector (uniform or clamped)
        for (size_t i = 0; i < num_knots; ++i) {
            if (i < deg + 1){
                knots[i] = x_min;  // Clamped knots at the start
            } else if (i >= N) {
                knots[i] = x_max;  // Clamped knots at the end
            } else {
                knots[i] = x_min + (i - deg) * dx;
            }
        }

        Eigen::MatrixXd basis(N, N);
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                basis(i, j) = b_spline_basis(x_vals[i], j, deg);
            }
        }

        Eigen::VectorXd values(N);
        for (size_t i = 0; i < N; ++i) {
            values(i) = f_vals[i];
        }

        Eigen::VectorXd coeffs = basis.colPivHouseholderQr().solve(values);
        for (size_t i = 0; i < N; ++i) {
            coefficients[i] = coeffs(i);
        }
    }

    // B-spline basis function
    double b_spline_basis(double x, size_t i, size_t k) const {
        if (k == 0) {
            return (x >= knots[i] && x < knots[i + 1]) ? 1.0 : 0.0;
        }

        double term1 = 0.0, term2 = 0.0;
        if (knots[i + k] != knots[i]) {
            term1 = (x - knots[i]) / (knots[i + k] - knots[i]) * b_spline_basis(x, i, k - 1);
        }
        if (knots[i + k + 1] != knots[i + 1]) {
            term2 = (knots[i + k + 1] - x) / (knots[i + k + 1] - knots[i + 1]) * b_spline_basis(x, i + 1, k - 1);
        }

        return term1 + term2;
    }

    // Recursive method to compute the n'th derivative of a B-spline basis function
    double b_spline_basis_derivative(double x, size_t i, size_t k, int n) const {
        if (n == 0) {
            return b_spline_basis(x, i, k);
        }

        double term1 = 0.0, term2 = 0.0;
        if ( (k > 0) && (knots[i + k] != knots[i]) ) {
            term1 = k * b_spline_basis_derivative(x, i, k - 1, n - 1) / (knots[i + k] - knots[i]);
        }

        if ( (k > 0) && (knots[i + k + 1] != knots[i + 1]) ) {
            term2 = k * b_spline_basis_derivative(x, i + 1, k - 1, n - 1) / (knots[i + k + 1] - knots[i + 1]);
        }

        return term1 - term2;
    }
};