#include "Integration.h"
#include "../global_params.h"
#include <vector>
#include <algorithm>
#include <memory>
#include <utility>
#include <map>
#include <math.h>

constexpr double FLT_EPS = 1e-12;

double eval_function(std::shared_ptr<Point> p, int Nx, int Ny,
                        std::function<double(double, double)> func,
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

double eval_function(Point p, int Nx, int Ny,
                        std::function<double(double, double)> func,
                        std::map<std::pair<int, int>, const double>& evaluated_points){
        std::pair<int, int> pos{Nx, Ny};
        if (evaluated_points.find(pos) == evaluated_points.end()){
            double val = func(p.x, p.y);
            evaluated_points.insert(std::pair<std::pair<int, int>, const double>(pos, val));
            return val;
        }
        return evaluated_points[pos];
}

void integration_step(std::shared_ptr<Point>& p1, std::shared_ptr<Point>& p2, std::shared_ptr<Point>& p3,
                        int& Nx, int& Ny, double& integral,
                        const double& dx, const double& dy,
                        int& Nxsteps, const int& Nysteps,
                        const double subdomain_dblder_limit,
                        std::map<std::pair<int, int>, const double>& evaluated_points,
                        std::function<double(double, double)> func){

    Point ystep{0, dy * Nysteps};
    Point xstep{dx * Nxsteps, - dy * Nysteps};

    double hx = abs(dx * Nxsteps);
    double hy = abs(dy * Nysteps);

    double f = eval_function(*p3, Nx, Ny, func, evaluated_points);
    double f1x = eval_function(*p3 + Point{hx, 0}, Nx + abs(Nxsteps), Ny, func, evaluated_points);
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
        p3 = std::make_shared<Point>(*p2 + ystep);
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
        p3 = std::shared_ptr<Point>{new Point(*p2 + ystep)};
        Ny += Nysteps;
        eval_function(p3, Nx, Ny, func, evaluated_points);

        integral += integrate_plane(p1, p2, p3);

        p1 = std::move(p2);
        p2 = std::move(p3);
        p3 = std::shared_ptr<Point>{new Point(*p2 + xstep)};
        Nx += Nxsteps;
        Ny -= Nysteps;
        eval_function(p3, Nx, Ny, func, evaluated_points);
        integral += integrate_plane(p1, p2, p3);
    }

}

double integrate2d(const Point& origin, const Point& end,
                    const double& dx, const double& dy,
                    const int& refinement_levels_x, const int& refinement_levels_y,
                    const double& subdomain_dblder_limit,
                    std::function<double(double, double)> func){
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

double integrate_adaptive(const Point& origin,
                            const int& Nx_origin, const int& Ny_origin,
                            const int& Nx_end, const int& Ny_end,
                            const double& dx, const double& dy,
                            int& Nxsteps, const int& Nysteps,
                            const double& subdomain_dblder_limit,
                            std::map<std::pair<int, int>, const double>& evaluated_points,
                            std::function<double(double, double)> func){

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