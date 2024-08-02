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
    One-sided tanh-sinh integration scheme.
    Integrates func(u) from 0 to 1, with step size dh in tanh-space,
    giving arbitrarily high resolution near u = 1
*/
double tanh_sinh(std::function<double(double)> func, double dh, double tol=1e-5);

