/*
Author: Vegard Gjeldvik Jervell
Contains: Implementation of Integration.h
          
*/

#include "Integration.h"
#include "../global_params.h"
#include <vector>
#include <algorithm>
#include <memory>
#include <utility>
#include <map>
#include <math.h>

#define FLTEPS 1e-12

Plane get_plane(const Point& p1, const Point& p2, const Point& p3){
    // a, b, c are components of vector normal to plane
    const double a = - p1.y * p3.z - p2.y * p1.z + p2.y * p3.z + p1.y * p2.z + p3.y * p1.z - p3.y * p2.z;
    const double b = p1.x * p3.z + p2.x * p1.z - p2.x * p3.z - p1.x * p2.z - p3.x * p1.z + p3.x * p2.z;
    const double c = - p1.x * p3.y - p2.x * p1.y + p2.x * p3.y + p1.x * p2.y + p3.x * p1.y - p3.x * p2.y;
    const double d = - (a * p1.x + b * p1.y + c * p1.z);
    return Plane{- a / c, - b / c, - d / c};
}

Line get_line(const Point& p1, const Point& p2){
    const double A = (p2.y - p1.y) / (p2.x - p1.x);
    const double B = p2.y - A * p2.x;
    return Line{A, B};
}

double integrate_plane_py(const Point& p1_in, const Point& p2_in, const Point& p3_in){ // For unittests on python-side
    std::shared_ptr<const Point> p1, p2, p3;
    p1 = std::make_shared<const Point>(p1_in);
    p2 = std::make_shared<const Point>(p2_in);
    p3 = std::make_shared<const Point>(p3_in);
    return integrate_plane(p1, p2, p3);
}

//Analytic integral of the triangle with corners (p1, p2, p3)
double integrate_plane(std::shared_ptr<const Point> p1, std::shared_ptr<const Point> p2, std::shared_ptr<const Point> p3){
    
    if ((p1->x == p2->x) && (p1->x == p3->x)){
        return 0;
    }
    const Plane plane = get_plane(*p1, *p2, *p3);
    if ((p1->x <= p2->x) && (p1->x <= p3->x)){ // Hvis p1_in.x er minst ...
        if (p2->x > p3->x){
            std::swap(p2, p3);
        }
    }
    else if ((p1->x >= p2->x) && (p1->x >= p3->x)){ // Hvis p1.x er størst ...
        std::swap(p1, p3); // Nå er p3 størst
        if (p1->x > p2->x){
            std::swap(p1, p2);
        }
    }
    else{ // p1_in.x er midterst
        if (p2->x < p3->x){
            std::swap(p1, p2);
        }
        else{
            std::swap(p1, p3);
            std::swap(p2, p3);
        }
    }


    std::shared_ptr<Line> l_upper_1{new Line{get_line(*p1, *p3)} }; // y    p1 * * * * * * p3   // Assuming this configuration
    std::shared_ptr<Line> l_lower_1{new Line{get_line(*p1, *p2)} }; // ^        *       *       // Swapping if y2 > y3
    std::shared_ptr<Line> l_upper_2{l_upper_1};                     // |          *   *         //
    std::shared_ptr<Line> l_lower_2{new Line{get_line(*p2, *p3)} }; //  => x     p2 *           //
    
    if (p3->y < p2->y){
        std::swap(l_lower_1, l_upper_1);
        std::swap(l_lower_2, l_upper_2);
    }
    
    double integral = 0;
    if (abs(p1->x - p2->x) > FLTEPS){ // Integral from x1 to x2

        const double A31_star = l_upper_1->a - l_lower_1->a;
        const double B31_star = l_upper_1->b - l_lower_1->b;
        const double A31_dblstar = pow(l_upper_1->a, 2) - pow(l_lower_1->a, 2);
        const double B31_dblstar = 2 * (l_upper_1->a * l_upper_1->b - l_lower_1->a * l_lower_1->b);
        const double C31_dblstar = pow(l_upper_1->b, 2) - pow(l_lower_1->b, 2);

        const double A12_tilde = plane.A * A31_star + plane.B * A31_dblstar / 2;
        const double B12_tilde = plane.A * B31_star + plane.C * A31_star + plane.B * B31_dblstar / 2;
        const double C12_tilde = plane.C * B31_star + plane.B * C31_dblstar / 2;

        integral += (A12_tilde / 3) * (pow(p2->x, 3) - pow(p1->x, 3)) + (B12_tilde / 2) * (pow(p2->x, 2) - pow(p1->x, 2)) + C12_tilde * (p2->x - p1->x);
    }
    if (abs(p2->x - p3->x) > FLTEPS){ // Integral from x2 to x3

        const double A23_star = l_upper_2->a - l_lower_2->a;
        const double B23_star = l_upper_2->b - l_lower_2->b;
        const double A23_dblstar = pow(l_upper_2->a, 2) - pow(l_lower_2->a, 2);
        const double B23_dblstar = 2 * (l_upper_2->a * l_upper_2->b - l_lower_2->a * l_lower_2->b);
        const double C23_dblstar = pow(l_upper_2->b, 2) - pow(l_lower_2->b, 2);

        const double A23_tilde = plane.A * A23_star + plane.B * A23_dblstar / 2;
        const double B23_tilde = plane.A * B23_star + plane.C * A23_star + plane.B * B23_dblstar / 2;
        const double C23_tilde = plane.C * B23_star + plane.B * C23_dblstar / 2;

        integral += (A23_tilde / 3) * (pow(p3->x, 3) - pow(p2->x, 3)) + (B23_tilde / 2) * (pow(p3->x, 2) - pow(p2->x, 2)) + C23_tilde * (p3->x - p2->x);
    }
    return integral;

}

// Checks the map to see if the function has already been evaluated, in which case value is retrieved from the map.
// Returns the function value, and also stores the value in the z-component of the point.
double eval_function(std::shared_ptr<Point> p, int Nx, int Ny,
                        const std::function<double(double, double)>& func,
                        std::map<std::pair<int, int>, const double>& evaluated_points){
        std::pair<int, int> pos{Nx, Ny};
        if (evaluated_points.find(pos) == evaluated_points.end()){
            double val = func(p->x, p->y);
            evaluated_points.insert(std::pair<std::pair<int, int>, const double>(pos, val));
            p->z = val;
            return val;
        }
        p->z = evaluated_points[pos];
        return evaluated_points[pos];
}

// Evaluate function without storing in point
double eval_function(Point p, int Nx, int Ny,
                        const std::function<double(double, double)>& func,
                        std::map<std::pair<int, int>, const double>& evaluated_points){
        std::pair<int, int> pos{Nx, Ny};
        if (evaluated_points.find(pos) == evaluated_points.end()){
            double val = func(p.x, p.y);
            evaluated_points.insert(std::pair<std::pair<int, int>, const double>(pos, val));
            return val;
        }
        return evaluated_points[pos];
}

// Conducts a single integration step (x += dx) if second derivative of integrand exceeds a given limit,
// integrate_adaptive() is called on the subdomain (x, x + dx)X(y, y + dy) with increased refinement
void integration_step(std::shared_ptr<Point>& p1, std::shared_ptr<Point>& p2, std::shared_ptr<Point>& p3,
                        int& Nx, int& Ny, double& integral,
                        double dx, double dy,
                        int& Nxsteps, int Nysteps,
                        double subdomain_dblder_limit,
                        std::map<std::pair<int, int>, const double>& evaluated_points,
                        const std::function<double(double, double)>& func){

    Point ystep{0, dy * Nysteps};
    Point xstep{dx * Nxsteps, - dy * Nysteps};

    double hx = abs(dx * Nxsteps);
    double hy = abs(dy * Nysteps);

    double f = eval_function(*p3, Nx, Ny, func, evaluated_points);
    double f1x = eval_function(*p3 + Point{hx, 0}, Nx + abs(Nxsteps), Ny,  func, evaluated_points);
    double f2x = eval_function(*p3 + Point{2 * hx, 0}, Nx + 2 * abs(Nxsteps), Ny, func, evaluated_points);
    double f1y = eval_function(*p3 + Point{0, hy}, Nx, Ny + abs(Nysteps), func, evaluated_points);
    double f2y = eval_function(*p3 + Point{0, 2 * hy}, Nx, Ny + 2 * abs(Nysteps), func, evaluated_points);
    double f1x1y = eval_function(*p3 + Point{hx, hy}, Nx + abs(Nxsteps), Ny + abs(Nysteps), func, evaluated_points);
    double f2x2y = eval_function(*p3 + Point{2 * hx, 2 * hy}, Nx + 2 * abs(Nxsteps), Ny + 2 * abs(Nysteps), func, evaluated_points);
    double fxx = (f - 2 * f1x + f2x) / pow(dx * Nxsteps, 2);
    double fyy = (f - 2 * f1y + f2y) / pow(dy * Nysteps, 2);
    double fxy = (f2x2y - 2 * f1x1y + f - pow(hx, 2) * fxx - pow(hy, 2) * fyy) / (2 * hx * hy);

    double dblder = abs(fxx) + abs(fyy) + abs(fxy);

    if ((dblder * pow(hx + hy, 3) > subdomain_dblder_limit) && ((abs(Nxsteps) > 1) || (abs(Nysteps) > 1))){ // Increase refinement
        Point sub_origin = *p3;
        int sub_Nx_origin = Nx;
        int sub_Ny_origin = Ny;
        int sub_Nxsteps = Nxsteps;
        int sub_Nysteps = Nysteps;
        int sub_Nx_end = Nx + Nxsteps;
        int sub_Ny_end = Ny + Nysteps;
        double sub_subdomain_dblder_limit = subdomain_dblder_limit;
        if (abs(fxx) * pow(hx, 3) > subdomain_dblder_limit && (abs(fyy) * pow(hy, 3) < subdomain_dblder_limit) && (abs(fxy) < subdomain_dblder_limit) && (abs(Nxsteps) > 1)){ // Only need to refine x
            sub_Nxsteps = sub_Nxsteps / 2;
        }
        else if (abs(fyy) * pow(hy, 3) > sub_subdomain_dblder_limit && (abs(fxx)  * pow(hx, 3) < subdomain_dblder_limit) && (abs(fxy) < subdomain_dblder_limit) && (abs(Nysteps) > 1)){
            sub_Nysteps = sub_Nysteps / 2;
        }
        else { // Refine x if possible and y if possible
            if (abs(Nxsteps) > 1) {
                sub_Nxsteps = sub_Nxsteps / 2;
            }
            if (abs(Nysteps) > 1){
                sub_Nysteps = sub_Nysteps / 2;
            }
        }
        integral += integrate_adaptive(sub_origin,
                                       sub_Nx_origin, sub_Ny_origin,
                                       sub_Nx_end, sub_Ny_end,
                                       dx, dy,
                                       sub_Nxsteps, sub_Nysteps,
                                       sub_subdomain_dblder_limit,
                                       evaluated_points,
                                       func);
        // Set all points to the gridpoint at the lower right corner of the subdomain that was just integrated (if Nxsteps is positive, otherwise to the lower left corner)
        p1 = std::move(p2);
        p2 = std::move(p3);
        p3 = std::shared_ptr<Point>{new Point(*p2 + ystep)};
        Ny += Nysteps;
        eval_function(p3, Nx, Ny, func, evaluated_points);
        integral += integrate_plane(p1, p2, p3);

        p1 = p3;
        p2 = p3;
        *p3 += xstep; // Set all points to the point following the refined region
        Nx += Nxsteps;
        Ny -= Nysteps;
        eval_function(p3, Nx, Ny, func, evaluated_points);
    }
    else{
        p1 = std::move(p2);
        p2 = std::move(p3);
        p3 = std::make_shared<Point>(*p2 + ystep);
        Ny += Nysteps;
        eval_function(p3, Nx, Ny, func, evaluated_points);

        integral += integrate_plane(p1, p2, p3);

        p1 = std::move(p2);
        p2 = std::move(p3);
        p3 = std::make_shared<Point>(*p2 + xstep);
        Nx += Nxsteps;
        Ny -= Nysteps;
        eval_function(p3, Nx, Ny, func, evaluated_points);
        integral += integrate_plane(p1, p2, p3);
    }

}

// Integrates the square domain bounded by (origin, end) in the xy-plane by calling integration_step()
// until the end of the domain is reached.
double integrate_adaptive(const Point& origin,
                            int Nx_origin, int Ny_origin,
                            int Nx_end, int Ny_end,
                            double dx, double dy,
                            int& Nxsteps, int Nysteps,
                            double subdomain_dblder_limit,
                            std::map<std::pair<int, int>, const double>& evaluated_points,
                            const std::function<double(double, double)>& func){

    double integral = 0;
    Point ystep{0, dy * Nysteps};
    Point xstep{dx * Nxsteps, - dy * Nysteps};

    std::shared_ptr<Point> p1, p2, p3;
    p1 = std::make_shared<Point>(origin);
    p2 = p1;
    p3 = p2;
    int Nx, Ny;
    Nx = Nx_origin;
    Ny = Ny_origin;
    eval_function(p3, Nx, Ny, func, evaluated_points);

    integration_step(p1, p2, p3, Nx, Ny, integral, dx, dy, Nxsteps, Nysteps, subdomain_dblder_limit, evaluated_points, func);
    while (Ny < Ny_end){
        while (std::min(Nx_origin, Nx_end) < Nx && Nx < std::max(Nx_origin, Nx_end)){
            integration_step(p1, p2, p3, Nx, Ny, integral, dx, dy, Nxsteps, Nysteps, subdomain_dblder_limit, evaluated_points, func);
        }
        p1 = std::move(p2);
        p2 = std::move(p3);
        p3 = std::make_shared<Point>(*p2 + ystep);
        Ny += Nysteps;
        eval_function(p3, Nx, Ny, func, evaluated_points);
        integral += integrate_plane(p1, p2, p3);
        if (Ny < Ny_end){
            Nxsteps *= -1;
            integration_step(p1, p2, p3, Nx, Ny, integral, dx, dy, Nxsteps, Nysteps, subdomain_dblder_limit, evaluated_points, func);
        }
    }

    return integral;
}

// Interface to integrate_adaptive() for more simple use
double integrate2d(const Point& origin, const Point& end,
                    double dx, double dy,
                    int refinement_levels_x, int refinement_levels_y,
                    double subdomain_dblder_limit,
                    const std::function<double(double, double)>& func){

    int Nx_origin{0}, Ny_origin{0};
    double delta_x = end.x - origin.x;
    double delta_y = end.y - origin.y;
    int Nx_end = (int) (delta_x / dx + 0.5);
    int Ny_end = (int) (delta_y / dy + 0.5);
    int Nxsteps = refinement_levels_x;
    int Nysteps = refinement_levels_y;
    std::map<std::pair<int, int>, const double> evaluated_points;

    double val = integrate_adaptive(origin,
                                     Nx_origin, Ny_origin,
                                     Nx_end, Ny_end,
                                     dx, dy,
                                     Nxsteps, Nysteps,
                                     subdomain_dblder_limit,
                                     evaluated_points,
                                     func);
    return val;
}


double simpson(std::function<double(double)> func, double x0, double xN, int N_intervals){
    double dx = (xN - x0) / N_intervals;
    double val = func(x0) + func(xN);
    for (size_t i = 1; i <= (N_intervals / 2) - 1; i++){
        val += 4. * func(x0 + (2 * i - 1) * dx) + 2. * func(x0 + (2 * i) * dx);
    }
    val += 4. * func(xN - dx);
    val *= dx / 3.;
    return val;
}

std::pair<double, double> weighted_simpson(std::function<double(double)> func, std::function<double(double)> wt, double x0, double xN, int N_intervals){
    double dx = (xN - x0) / N_intervals;
    double F{0}, W{0}, w1{0}, w2{0};

    w1 = wt(x0); w2 = wt(xN);
    W += w1 + w2;
    F += func(x0) * w1 + func(xN) * w2;
    for (size_t i = 1; i <= (N_intervals / 2) - 1; i++){
        w1 = wt(x0 + (2 * i - 1) * dx);
        w2 = wt(x0 + (2 * i) * dx);
        W += 4. * w1 + 2. * w2;
        F += 4. * func(x0 + (2 * i - 1) * dx) * w1 + 2. * func(x0 + (2 * i) * dx) * w2;
    }
    w1 = wt(xN - dx);
    W += 4. * w1;
    F += 4. * func(xN - dx) * w1;
    W *= dx / 3.;
    F *= dx / 3.;
    return std::pair<double, double>(F, W);
}

std::pair<double, double> weighted_simpson(std::function<std::pair<double, double>(double)> FW, double x0, double xN, int N_intervals){
    double dx = (xN - x0) / N_intervals;
    double F{0}, W{0};
    std::pair<double, double> FW_1, FW_2;
    FW_1 = FW(x0); FW_2 = FW(xN);
    F += FW_1.first + FW_2.first;
    W += FW_1.second + FW_2.second;
    for (size_t i = 1; i <= (N_intervals / 2) - 1; i++){
        FW_1 = FW(x0 + (2 * i - 1) * dx); 
        FW_2 = FW(x0 + (2 * i) * dx);
        F += 4. * FW_1.first + 2. * FW_2.first;
        W += 4. * FW_1.second + 2. * FW_2.second;
    }
    FW_1 = FW(xN - dx);
    F += 4. * FW_1.first;
    W += 4. * FW_1.second;
    F *= dx / 3;
    W *= dx / 3;
    return std::pair<double, double>(F, W);
}

double tanh_sinh(std::function<double(double)> func, double h, double tol){
    constexpr int n_max_intervals = 10000;
    double I{0.0};
    int k{1};
    double u = tanh(PI * sinh(k * h) / 2.);
    double w = (PI / 2.) * h * cosh(k * h) / pow(cosh(PI * sinh(k * h * 1e-3) / 2.), 2);
    double f = func(u);
    do {
        I += w * f;
        k += 1;
        u = tanh(PI * sinh(k * h) / 2.);
        w = (PI / 2) * h * cosh(k * h) / pow(cosh(PI * sinh(k * h) / 2.), 2);
        f = func(u);
        if (isnan(f) || isnan(w) || isinf(f) || isinf(w)) break;
    } while ((abs((f * w) / I) > tol) && (k < n_max_intervals));
    if (k >= n_max_intervals){
        std::cout << "INTEGRATION WARNING : tanh-sinh integration reached maximum number of intervals (" << n_max_intervals
                    << "). Relative value of last integration step was " << abs((f * w) / I) << ", tolerance is : " << tol << std::endl;
    }
    return I;
}