/*
Author: Vegard Gjeldvik Jervell
Contains: Implementation of the adaptive 2D integration scheme described in "The Kinetic Gas theory of Mie fluids",
            Vegard Gjeldvik Jervell, Masters Thesis, Norwegian University of Science and Technology (2022).
            (https://ntnuopen.ntnu.no/ntnu-xmlui/handle/11250/3029213)

            If you are using this, the function integrate2d is most likely the function you want to use.

            The integration is based on triangulating a function on a square grid, and refining the size of
            the triangles if the function has large second derivatives. The algorithm is recursive, where
            
            integrate_adaptive() integrates a certain region by calling integration_step(). If integration_step()
            detects large second derivatives it calls integrate_adaptive() on a subregion.

            integrate2d() is simply an interface to integrate_adaptive() for easier use.

            The integrated function is an std::function<double(double, double)> 
            where the large number of arguments is a result of me not bothering to generalise this.
*/
#pragma once
#include <vector>
#include <algorithm>
#include <memory>
#include <utility>
#include <map>
#include <iostream>
#include <functional>

/*
    Some utility structures for holding geometric information
*/
class Point{
    public:
    double x, y, z;
    Point(double x, double y) : x{x}, y{y}, z{0} {};
    Point(double x, double y, double z) : x{x}, y{y}, z{z} {};

    Point operator+(const Point& rhs){
        return Point{x + rhs.x, y + rhs.y, z + rhs.z};
    }

    void operator+=(const Point& rhs){
        x += rhs.x; y += rhs.y; z += rhs.z;
    }
};

struct Plane{ // The plane given by z = Ax + By + C
    const double A, B, C;
    Plane(double A, double B, double C) : A{A}, B{B}, C{C} {};
};

struct Line{ // The line given by y = ax + b
    const double a, b;
    Line(double a, double b) : a{a}, b{b} {};
};

Plane get_plane(const Point& p1, const Point& p2, const Point& p3); // Return vector of [A, B, C] where these are the coefficients of z = Ax + By + C
Line get_line(const Point& p1, const Point& p2); // Return vector of [A, B] where y = Ax + B

// Integrate the plane interpolating (p1, p2, p3) analytically
double integrate_plane_py(const Point& p1, const Point& p2, const Point& p3);
double integrate_plane(std::shared_ptr<const Point> p1, std::shared_ptr<const Point> p2, std::shared_ptr<const Point> p3);

/*
    func : The function to be evaluated
    p : The point at which to evaluate the function
    Nx, Ny : The indexes of the point on some pre-defined grid
    evaluated_points : Map of (Nx, Ny) -> Function value

    Check whether the point indexed (Nx, Ny) has already been evaluated. If not, evaluate the function 
    and store the value before returning it.
*/
double eval_function(std::shared_ptr<Point> p, int Nx, int Ny,
                        const std::function<double(double, double)>& func,
                        std::map<std::pair<int, int>, const double>& evaluated_points);

double eval_function(Point p, int Nx, int Ny,
                        const std::function<double(double, double)>& func,
                        std::map<std::pair<int, int>, const double>& evaluated_points);

/*
    p1, p2, p3 : Three points forming the previous triangular integration domain
    Nx, Ny : Indexes of the point p3
    integral : The value of the integral before adding this domain
    dx, dy : Minimum integration step lengths
    Nxsteps, Nysteps : Multiple of integration step lengths
    subdomain_dblder_limit : If (absolute) second derivative is larger than this, the integration domain will be 
                            subdivided and the integration will be carried out by integrate_adaptive.
    evaluated_points : Map of (Nx, Ny) -> function value, for lazy evaluation
    func : The function to integrate 

    Carry out an integration step, moving (dx * Nxsteps) in the x-direction and (dy * Nysteps) in the y-direction.
*/
void integration_step(std::shared_ptr<Point>& p1, std::shared_ptr<Point>& p2, std::shared_ptr<Point>& p3,
                        int& Nx, int& Ny, double& integral,
                        double dx, double dy,
                        int& Nxsteps, int Nysteps,
                        double subdomain_dblder_limit,
                        std::map<std::pair<int, int>, const double>& evaluated_points,
                        const std::function<double(double, double)>& func);

/*
    See integrate2d for practical usage.

    Nx_origin, Ny_origin : The indexes of the point "origin" in evaluated_points
    Nx_end, Ny_end : The indexes of the last point in evaluated_points
    dx, dy : Minimum integration step
    Nxsteps, Nysteps : Multiple of dx, dy to use for integration step
    subdomain_dblder_limit : Tolerance for refining the domain
    evaluated_points : Map of (Nx, Ny) -> function value, for lazy evaluation

    Integrates the function "func" over the square domain starting at origin, and extending to 
    (dx * Nx_end) in the x-direction, and (dy * Ny_end) in the y-direction. The initial 
    integration step length is dx * Nx_steps, dy * Nysteps.

    Integration steps are carried out by integration_step, which may in turn call this function, 
    such that the two together form a recursive call structure.
*/
double integrate_adaptive(const Point& origin,
                            int Nx_origin, int Ny_origin,
                            int Nx_end, int Ny_end,
                            double dx, double dy,
                            int& Nxsteps, int Nysteps,
                            double subdomain_dblder_limit,
                            std::map<std::pair<int, int>, const double>& evaluated_points,
                            const std::function<double(double, double)>& func);

/*
    dx, dy : The minimum step sizes in the x- and y-directions
    refinement_levels_[x/y] : Number of times integration domain can be refined
    subdomain_dblder_limit : Sensitivity to refinement (lower value gives more refinement)

    Integrate the square with origin in the "lower left" corner and end in the "upper right" corner.
    The initial integration steps are (dx * refinement_levels_x) and (dy * refinement_levels_y),
    if a second derivative exceeds subdomain_dblder_limit, the integration steps will be halved.

    Thus, refinement_levels_[x/y] should be a power of 2, where refinement_level_[x/y] = 1, will 
    give no refinement, refinement_level_[x/y] = 2 gives the option of refining once, and
    refinement_level_[x/y] = pow(2, N) gives the possibility of refining N times.
*/
double integrate2d(const Point& origin, const Point& end,
                    double dx, double dy,
                    int refinement_levels_x, int refinement_levels_y,
                    double subdomain_dblder_limit,
                    const std::function<double(double, double)>& func);

/*
    Standard Simpson-rule integration, with N_intervals between x0 and xN, where N_intervals >= 2
*/
double simpson(std::function<double(double)> func, double x0, double xN, int N_intervals);
/*
    Simpson rule, starting with 10 points between x0 and init_end, then continuing outwards from init_end until infinite integral has converged
    For users: This method is intended to integrate the "long tail" of some function. You likely want to handle integration close to the origin
                yourself, before using this method to capture the tail.

    Example:
    auto func = [](double r){return r * sin(2 * PI * r) * exp(-r);}; // Want to integrate this from 0 to inf
    double I1 = simpson(func, 0, 2, 50); // High resolution at start to capture details
    double I2 = simpson_inf(func, 2, 3); // This will continue outward to capture the tail of the function
    double I = I1 + I2; // Total integral
*/
double simpson_inf(std::function<double(double)> func, double x0, double init_end, double tol=1e-8, double I=0);


/*
    Simpson integration of the function g(x) = w(x)f(x), and w(x) simultaneously.
    Useful if w(x) is expensive to evaluate, and you want to evaluate the weighted average
        a = I[w(x) * f(x)] / I[w(x)],
    which is then evaluated as
        std::pair<double, double> FW = weighted_simpson(func, wt, x0, xN, N);
        double a = FW.first / FW.second;
    Returns : {F, W}, where F is the integral of f(x) * w(x), and W is the integral of w(x)
*/
std::pair<double, double> weighted_simpson(std::function<double(double)> func, std::function<double(double)> wt, double x0, double xN, int N_intervals);

/*
    For use when nesting std::pair<double, double> weighted_simpson(double(double), double(double), double, double, int)
    Takes the function "FW" as an argument, and integrates both outputs using Simpsons rule.

    Returns : {F, W}, where F is the integral of the first output, and W is the integral of the second output.
*/
std::pair<double, double> weighted_simpson(std::function<std::pair<double, double>(double)> FW, double x0, double xN, int N_intervals);


/*
    One-sided tanh-sinh integration scheme.
    Integrates func(u) from 0 to 1, with step size dh in tanh-space,
    giving arbitrarily high resolution near u = 1
*/
double tanh_sinh(std::function<double(double)> func, double dh, double tol=1e-5);

/*
    Standard Newton solver

    * newton_usafe: Will fail silently, signalling an error with the ierr flag.
    * newton: Will throw an error upon failure.

    Use the former if you plan to handle errors yourself. Use the latter if you're doing something that should never fail.
*/
double newton_usafe(const std::function<double(double)>& fun, const std::function<double(double)>& df, double x0, double ftol, double dtol, int& ierr) noexcept;
inline double newton(const std::function<double(double)>& fun, const std::function<double(double)>& df, double x0, double ftol=1e-10, double dtol=-1){
    int ierr = 0;
    double r = newton_usafe(fun, df, x0, ftol, dtol, ierr);
    if (ierr == 0) return r;
    if (ierr == 1) throw std::runtime_error("Newton reached max iter!");
    if (ierr == 2) throw std::runtime_error("Newton encountered NAN!");
    if (ierr == 3) throw std::runtime_error("Newton encountered INF!");
    throw std::runtime_error("Unknown Newton error : " + std::to_string(ierr));
}


/*
    Basic bracket solvers

    * bracet_positive: Will always return the positive side of the bracket
    * bracket_root_usafe: Will never throw, and will set x0 to the value closest to the root (determined from abs(f0) < abs(f1))
    * bracket_root: Will do some sanity checking, and may throw.
*/
double bracket_positive(const std::function<double(double)>& fun, double x0, double x1, double tol=1e-10);
void bracket_root_usafe(const std::function<double(double)>& fun, double& x0, double& x1, double& f0, double& f1, double dtol, double ftol) noexcept;
inline double bracket_root(const std::function<double(double)>& fun, double x0, double x1, double dtol=1e-5, double ftol=1e-10){
    double f0 = fun(x0);
    double f1 = fun(x1);
    if (f0 * f1 > 0){
        throw std::runtime_error("Bracket solver: Initial values have same sign!");
    }
    bracket_root_usafe(fun, x0, x1, f0, f1, dtol, ftol);
    return x0;
}

/*
    Simple trapezoidal integration, either using arrays, or using a constant grid spacing
*/
double trapezoid(const std::vector<double>& x, const std::vector<double>& y);
double trapezoid(const double dx, const std::vector<double>& y);

// Fit the quadric f(x) = ax^2 + bx + c to the (x, y) data 
std::array<double, 3> fit_quadric(const std::array<double, 3>& x, const std::array<double, 3>& y);

// Same as above, but fit to the last three elements in the vectors.
std::array<double, 3> quadric_extrapolate_coeff(const std::vector<double>& x, const std::vector<double>& y);

/*
    Do a linear interpolation on (x, y) to get the value at x_val
    Will throw if x_val is outside the range of the x-values
*/
double interpolate_grid(const double x_val, const std::vector<double>& x, const std::vector<double>& y);

/*
    Get the y-values at every new_x position by linear interpolation of (old_x, y)

    Useful when we solve something on one grid (old_x) which is e.g. non-uniform, and want to extract the values on some new grid (new_x) from the solution.
*/
std::vector<double> interpolate_grid(const std::vector<double>& new_x, const std::vector<double>& old_x, const std::vector<double>& y);