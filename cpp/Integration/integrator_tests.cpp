#include "integrator_tests.h"
#include "Integration.h"

void mesh_step(std::shared_ptr<Point>& p1, std::shared_ptr<Point>& p2, std::shared_ptr<Point>& p3,
                        int& Nx, int& Ny,
                        const double& dx, const double& dy,
                        int& Nxsteps, const int& Nysteps,
                        const double subdomain_dblder_limit,
                        std::map<std::pair<int, int>, const double>& evaluated_points,
                        const int& arg_i, const int& arg_j, const double& arg_T, const int& arg_l, const int& arg_r,
                        std::function<double(int, int, double, double, double, int, int)> func, std::vector<Point>& points){

    Point ystep{0, dy * Nysteps};
    Point xstep{dx * Nxsteps, - dy * Nysteps};

    double hx = abs(dx * Nxsteps);
    double hy = abs(dy * Nysteps);

    double f = eval_function(*p3, Nx, Ny, arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
    double f1x = eval_function(*p3 + Point{hx, 0}, Nx + abs(Nxsteps), Ny, arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
    double f2x = eval_function(*p3 + Point{2 * hx, 0}, Nx + 2 * abs(Nxsteps), Ny, arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
    double f1y = eval_function(*p3 + Point{0, hy}, Nx, Ny + abs(Nysteps), arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
    double f2y = eval_function(*p3 + Point{0, 2 * hy}, Nx, Ny + 2 * abs(Nysteps), arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
    double f1x1y = eval_function(*p3 + Point{hx, hy}, Nx + abs(Nxsteps), Ny + abs(Nysteps), arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
    double f2x2y = eval_function(*p3 + Point{2 * hx, 2 * hy}, Nx + 2 * abs(Nxsteps), Ny + 2 * abs(Nysteps), arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
    double fxx = (f - 2 * f1x + f2x) / pow(dx * Nxsteps, 2);
    double fyy = (f - 2 * f1y + f2y) / pow(dy * Nysteps, 2);
    double fxy = (f2x2y - 2 * f1x1y + f - pow(hx, 2) * fxx - pow(hy, 2) * fyy) / (2 * hx * hy);

    double dblder = abs(fxx) + abs(fyy) + abs(fxy);

    if ((dblder * pow(hx + hy, 3)) > subdomain_dblder_limit && (abs(Nxsteps) > 1) && (abs(Nysteps) > 1)){ // Increase refinement
        Point sub_origin = *p3;
        int sub_Nx_origin = Nx;
        int sub_Ny_origin = Ny;
        int sub_Nxsteps = Nxsteps;
        int sub_Nysteps = Nysteps;
        int sub_Nx_end = Nx + Nxsteps;
        int sub_Ny_end = Ny + Nysteps;
        double sub_subdomain_dblder_limit = subdomain_dblder_limit;
        if (abs(fxx) * pow(hx, 3) > subdomain_dblder_limit  && (abs(fyy) * pow(hy, 3)) < subdomain_dblder_limit && (abs(Nxsteps) > 1)){ // Only need to refine x
            sub_Nxsteps = sub_Nxsteps / 2;
            // sub_subdomain_dblder_limit = sub_subdomain_dblder_limit / pow(hx, 3);
        }
        else if (abs(fyy)  * pow(hy, 3) > sub_subdomain_dblder_limit && (abs(fxx) * pow(hx, 3)) < subdomain_dblder_limit && (abs(Nysteps) > 1)){
            sub_Nysteps = sub_Nysteps / 2;
            // sub_subdomain_dblder_limit = sub_subdomain_dblder_limit / pow(hy, 3);
        }
        else { // Refine x if possible and y if possible
            if (abs(Nxsteps) > 1) {
                sub_Nxsteps = sub_Nxsteps / 2;
                // sub_subdomain_dblder_limit = sub_subdomain_dblder_limit / pow(hx, 3);
            }
            if (abs(Nysteps) > 1){
                sub_Nysteps = sub_Nysteps / 2;
                // sub_subdomain_dblder_limit = sub_subdomain_dblder_limit / pow(hy, 3);
            }
        }

        mesh_adaptive(sub_origin,
                       sub_Nx_origin, sub_Ny_origin,
                       sub_Nx_end, sub_Ny_end,
                       dx, dy,
                       sub_Nxsteps, sub_Nysteps,
                       sub_subdomain_dblder_limit,
                       evaluated_points,
                       arg_i, arg_j, arg_T, arg_l, arg_r,
                       func, points);

        p1 = std::move(p2);
        p2 = std::move(p3);
        p3 = std::shared_ptr<Point>{new Point(*p2 + ystep)};
        Ny += Nysteps;
        eval_function(p3, Nx, Ny, arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
        points.push_back(Point(*p3));

        p1 = p3;
        p2 = p3;
        *p3 += xstep; // Set all points to the point following the refined region
        Nx += Nxsteps;
        Ny -= Nysteps;
        eval_function(p3, Nx, Ny, arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
        points.push_back(Point(*p3));

    }
    else{
        p1 = std::move(p2);
        p2 = std::move(p3);
        p3 = std::shared_ptr<Point>{new Point(*p2 + ystep)};
        Ny += Nysteps;
        eval_function(p3, Nx, Ny, arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);

        points.push_back(Point(*p3));

        p1 = std::move(p2);
        p2 = std::move(p3);
        p3 = std::shared_ptr<Point>{new Point(*p2 + xstep)};
        Nx += Nxsteps;
        Ny -= Nysteps;
        eval_function(p3, Nx, Ny, arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
        points.push_back(Point(*p3));
    }

}

void mesh_adaptive(const Point& origin,
                    const int& Nx_origin, const int& Ny_origin,
                    const int& Nx_end, const int& Ny_end,
                    const double& dx, const double& dy,
                    int& Nxsteps, const int& Nysteps,
                    const double& subdomain_dblder_limit,
                    std::map<std::pair<int, int>, const double>& evaluated_points,
                    const int& arg_i, const int& arg_j, const double& arg_T, const int& arg_l, const int& arg_r,
                    std::function<double(int, int, double, double, double, int, int)> func, std::vector<Point>& points){

    Point ystep{0, dy * Nysteps};
    Point xstep{dx * Nxsteps, - dy * Nysteps};

    std::shared_ptr<Point> p1, p2, p3;
    p1 = std::make_shared<Point>(origin);
    p2 = p1;
    p3 = p2;
    int Nx, Ny;
    Nx = Nx_origin;
    Ny = Ny_origin;
    eval_function(p3, Nx, Ny, arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
    points.push_back(Point(*p3));

    mesh_step(p1, p2, p3, Nx, Ny, dx, dy, Nxsteps, Nysteps, subdomain_dblder_limit, evaluated_points, arg_i, arg_j, arg_T, arg_l, arg_r, func, points);
    //points.push_back(Point(*p3));
    while (Ny < Ny_end){
        while (std::min(Nx_origin, Nx_end) < Nx && Nx < std::max(Nx_origin, Nx_end)){
            mesh_step(p1, p2, p3, Nx, Ny, dx, dy, Nxsteps, Nysteps, subdomain_dblder_limit, evaluated_points, arg_i, arg_j, arg_T, arg_l, arg_r, func, points);
        }
        p1 = std::move(p2);
        p2 = std::move(p3);
        p3 = std::make_shared<Point>(*p2 + ystep);
        Ny += Nysteps;
        eval_function(p3, Nx, Ny, arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
        points.push_back(Point(*p3));
        if (Ny < Ny_end){
            Nxsteps *= -1;
            mesh_step(p1, p2, p3, Nx, Ny, dx, dy, Nxsteps, Nysteps, subdomain_dblder_limit, evaluated_points, arg_i, arg_j, arg_T, arg_l, arg_r, func, points);
            //points.push_back(Point(*p3));
        }
    }
}

std::vector<std::vector<double>> mesh2d(const Point& origin, const Point& end,
                                        const double& dx, const double& dy,
                                        const int& refinement_levels,
                                        const double& subdomain_dblder_limit,
                                        const int& arg_i, const int& arg_j, const double& arg_T, const int& arg_l, const int& arg_r,
                                        std::function<double(int, int, double, double, double, int, int)> func){

    int Nx_origin{0}, Ny_origin{0};
    double delta_x = end.x - origin.x;
    double delta_y = end.y - origin.y;
    int Nx_end = (int) (delta_x / dx + 0.5);
    int Ny_end = (int) (delta_y / dy + 0.5);
    int Nxsteps = refinement_levels;
    int Nysteps = refinement_levels;
    std::map<std::pair<int, int>, const double> evaluated_points;
    std::vector<Point> points;
    mesh_adaptive(origin,
                 Nx_origin, Ny_origin,
                 Nx_end, Ny_end,
                 dx, dy,
                 Nxsteps, Nysteps,
                 subdomain_dblder_limit,
                 evaluated_points,
                 arg_i, arg_j, arg_T, arg_l, arg_r,
                 func, points);

    std::vector<double> x, y, z;

    for (std::vector<Point>::iterator it = points.begin(); it != points.end(); it++){
        x.push_back(it->x);
        y.push_back(it->y);
        z.push_back(it->z);
    }
    return std::vector<std::vector<double>> {x, y, z};
}

void trisurf_step(std::shared_ptr<Point>& p1, std::shared_ptr<Point>& p2, std::shared_ptr<Point>& p3,
                    int& Nx, int& Ny, std::vector<std::vector<Point>>& surf,
                    const double& dx, const double& dy,
                    int& Nxsteps, const int& Nysteps,
                    const double subdomain_dblder_limit,
                    std::map<std::pair<int, int>, const double>& evaluated_points,
                    const int& arg_i, const int& arg_j, const double& arg_T, const int& arg_l, const int& arg_r,
                    std::function<double(int, int, double, double, double, int, int)> func){

    Point ystep{0, dy * Nysteps};
    Point xstep{dx * Nxsteps, - dy * Nysteps};

    double hx = abs(dx * Nxsteps);
    double hy = abs(dy * Nysteps);

    double f = eval_function(*p3, Nx, Ny, arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
    double f1x = eval_function(*p3 + Point{hx, 0}, Nx + abs(Nxsteps), Ny, arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
    double f2x = eval_function(*p3 + Point{2 * hx, 0}, Nx + 2 * abs(Nxsteps), Ny, arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
    double f1y = eval_function(*p3 + Point{0, hy}, Nx, Ny + abs(Nysteps), arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
    double f2y = eval_function(*p3 + Point{0, 2 * hy}, Nx, Ny + 2 * abs(Nysteps), arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
    double f1x1y = eval_function(*p3 + Point{hx, hy}, Nx + abs(Nxsteps), Ny + abs(Nysteps), arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
    double f2x2y = eval_function(*p3 + Point{2 * hx, 2 * hy}, Nx + 2 * abs(Nxsteps), Ny + 2 * abs(Nysteps), arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
    double fxx = (f - 2 * f1x + f2x) / pow(dx * Nxsteps, 2);
    double fyy = (f - 2 * f1y + f2y) / pow(dy * Nysteps, 2);
    double fxy = (f2x2y - 2 * f1x1y + f - pow(hx, 2) * fxx - pow(hy, 2) * fyy) / (2 * hx * hy);

    double dblder = abs(fxx) + abs(fyy) + abs(fxy);

    if ((dblder > subdomain_dblder_limit) && (abs(Nxsteps) > 1) && (abs(Nysteps) > 1)){ // Increase refinement
        Point sub_origin = *p3;
        int sub_Nx_origin = Nx;
        int sub_Ny_origin = Ny;
        int sub_Nxsteps = Nxsteps;
        int sub_Nysteps = Nysteps;
        int sub_Nx_end = Nx + Nxsteps;
        int sub_Ny_end = Ny + Nysteps;
        double sub_subdomain_dblder_limit = subdomain_dblder_limit;
        if (abs(fxx) > subdomain_dblder_limit && (abs(fyy) < subdomain_dblder_limit) && (abs(Nxsteps) > 1)){ // Only need to refine x
            sub_Nxsteps = sub_Nxsteps / 2;
            sub_subdomain_dblder_limit = sub_subdomain_dblder_limit / pow(hx, 3);
        }
        else if (abs(fyy) > sub_subdomain_dblder_limit && (abs(fxx) < subdomain_dblder_limit) && (abs(Nysteps) > 1)){
            sub_Nysteps = sub_Nysteps / 2;
            sub_subdomain_dblder_limit = sub_subdomain_dblder_limit / pow(hy, 3);
        }
        else { // Refine x if possible and y if possible
            if (abs(Nxsteps) > 1) {
                sub_Nxsteps = sub_Nxsteps / 2;
                sub_subdomain_dblder_limit = sub_subdomain_dblder_limit / pow(hx, 3);
            }
            if (abs(Nysteps) > 1){
                sub_Nysteps = sub_Nysteps / 2;
                sub_subdomain_dblder_limit = sub_subdomain_dblder_limit / pow(hy, 3);
            }
        }

        trisurf_adaptive(sub_origin,
                           sub_Nx_origin, sub_Ny_origin,
                           sub_Nx_end, sub_Ny_end,
                           dx, dy,
                           sub_Nxsteps, sub_Nysteps,
                           sub_subdomain_dblder_limit,
                           evaluated_points,
                           arg_i, arg_j, arg_T, arg_l, arg_r,
                           func, surf);
        // Set all points to the gridpoint at the lower right corner of the subdomain that was just integrated (if Nxsteps is positive, otherwise to the lower left corner)
        p1 = std::move(p2);
        p2 = std::move(p3);
        p3 = std::shared_ptr<Point>{new Point(*p2 + ystep)};
        Ny += Nysteps;
        eval_function(p3, Nx, Ny, arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
        surf.push_back(std::vector<Point>{Point(*p1), Point(*p2), Point(*p3)});

        p1 = p3;
        p2 = p3;
        *p3 += xstep; // Set all points to the point following the refined region
        Nx += Nxsteps;
        Ny -= Nysteps;
        eval_function(p3, Nx, Ny, arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
    }
    else{
        p1 = std::move(p2);
        p2 = std::move(p3);
        p3 = std::shared_ptr<Point>{new Point(*p2 + ystep)};
        Ny += Nysteps;
        eval_function(p3, Nx, Ny, arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);

        surf.push_back(std::vector<Point>{Point(*p1), Point(*p2), Point(*p3)});

        p1 = std::move(p2);
        p2 = std::move(p3);
        p3 = std::shared_ptr<Point>{new Point(*p2 + xstep)};
        Nx += Nxsteps;
        Ny -= Nysteps;
        eval_function(p3, Nx, Ny, arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
        surf.push_back(std::vector<Point>{Point(*p1), Point(*p2), Point(*p3)});
    }

}

void trisurf_adaptive(const Point& origin,
                        const int& Nx_origin, const int& Ny_origin,
                        const int& Nx_end, const int& Ny_end,
                        const double& dx, const double& dy,
                        int& Nxsteps, const int& Nysteps,
                        const double& subdomain_dblder_limit,
                        std::map<std::pair<int, int>, const double>& evaluated_points,
                        const int& arg_i, const int& arg_j, const double& arg_T, const int& arg_l, const int& arg_r,
                        std::function<double(int, int, double, double, double, int, int)> func,
                        std::vector<std::vector<Point>>& surf){
    Point ystep{0, dy * Nysteps};
    Point xstep{dx * Nxsteps, - dy * Nysteps};

    std::shared_ptr<Point> p1, p2, p3;
    p1 = std::make_shared<Point>(origin);
    p2 = p1;
    p3 = p2;
    int Nx, Ny;
    Nx = Nx_origin;
    Ny = Ny_origin;
    eval_function(p3, Nx, Ny, arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);

    trisurf_step(p1, p2, p3, Nx, Ny, surf, dx, dy, Nxsteps, Nysteps, subdomain_dblder_limit, evaluated_points, arg_i, arg_j, arg_T, arg_l, arg_r, func);
    while (Ny < Ny_end){
        while (std::min(Nx_origin, Nx_end) < Nx && Nx < std::max(Nx_origin, Nx_end)){
            trisurf_step(p1, p2, p3, Nx, Ny, surf, dx, dy, Nxsteps, Nysteps, subdomain_dblder_limit, evaluated_points, arg_i, arg_j, arg_T, arg_l, arg_r, func);
        }
        p1 = std::move(p2);
        p2 = std::move(p3);
        p3 = std::make_shared<Point>(*p2 + ystep);
        Ny += Nysteps;
        eval_function(p3, Nx, Ny, arg_i, arg_j, arg_T, arg_l, arg_r, func, evaluated_points);
        surf.push_back(std::vector<Point>{Point(*p1), Point(*p2), Point(*p3)});
        if (Ny < Ny_end){
            Nxsteps *= -1;
            trisurf_step(p1, p2, p3, Nx, Ny, surf, dx, dy, Nxsteps, Nysteps, subdomain_dblder_limit, evaluated_points, arg_i, arg_j, arg_T, arg_l, arg_r, func);
        }
    }
}

std::vector<std::vector<std::vector<double>>> trisurf(const Point& origin, const Point& end,
                                                        const double& dx, const double& dy,
                                                        const int& refinement_levels_x, const int& refinement_levels_y,
                                                        const double& subdomain_dblder_limit,
                                                        const int& arg_i, const int& arg_j, const double& arg_T, const int& arg_l, const int& arg_r,
                                                        std::function<double(int, int, double, double, double, int, int)> func){

    int Nx_origin{0}, Ny_origin{0};
    double delta_x = end.x - origin.x;
    double delta_y = end.y - origin.y;
    int Nx_end = (int) (delta_x / dx + 0.5);
    int Ny_end = (int) (delta_y / dy + 0.5);
    int Nxsteps = refinement_levels_x;
    int Nysteps = refinement_levels_y;
    std::map<std::pair<int, int>, const double> evaluated_points;
    std::vector<std::vector<Point>> surf;
    trisurf_adaptive(origin,
                     Nx_origin, Ny_origin,
                     Nx_end, Ny_end,
                     dx, dy,
                     Nxsteps, Nysteps,
                     subdomain_dblder_limit,
                     evaluated_points,
                     arg_i, arg_j, arg_T, arg_l, arg_r,
                     func, surf);

    std::vector<std::vector<std::vector<double>>> verts;
    for (std::vector<std::vector<Point>>::iterator it = surf.begin(); it != surf.end(); it++){
        verts.push_back(std::vector<std::vector<double>>{{it->at(0).x, it->at(0).y, it->at(0).z},
                                                         {it->at(1).x, it->at(1).y, it->at(1).z},
                                                         {it->at(2).x, it->at(2).y, it->at(2).z}});
    }

    return verts;
}

double testfun(const int i, const int j, const double T, const double x, const double y, const int r, const int l){
    return exp(- (pow(x - 5, 2) + pow(y - 5, 2)));
}

double testfun_linear(const int i, const int j, const double T, const double x, const double y, const int r, const int l){
    return x + y;
}

double integrator_test(double origin_x, double origin_y, double end_x, double end_y,
                       double dx, double dy, int refinement_levels, double subdomain_dblder_limit){
    int i{1}, j{1}, r{1}, l{1}; // Dummy values
    double T{1}; // Dummy values
    Point origin{origin_x, origin_y}, end{end_x, end_y};
    double val = integrate2d(origin, end, dx, dy, refinement_levels, refinement_levels, subdomain_dblder_limit, i, j, T, r, l, &testfun);
    double a = -log(testfun(i, j, T, 1, 0, r, l));
    return val; // Integral of testfun on (-inf, inf) is pi / a
}

double integrator_test_linear(double origin_x, double origin_y, double end_x, double end_y,
                       double dx, double dy, int refinement_levels, double subdomain_dblder_limit){
    int i{1}, j{1}, r{1}, l{1}; // Dummy values
    double T{1}; // Dummy values
    Point origin{origin_x, origin_y}, end{end_x, end_y};
    double val = integrate2d(origin, end, dx, dy, refinement_levels, refinement_levels, subdomain_dblder_limit, i, j, T, r, l, &testfun_linear);
    return val;
}

std::vector<std::vector<double>> mesh_test(double origin_x, double origin_y, double end_x, double end_y,
                                            double dx, double dy, int refinement_levels, double subdomain_dblder_limit){
    int i{1}, j{1}, r{1}, l{1}; // Dummy values
    double T{1}; // Dummy values
    Point origin{origin_x, origin_y}, end{end_x, end_y};
    std::vector<std::vector<double>> mesh = mesh2d(origin, end, dx, dy, refinement_levels, subdomain_dblder_limit, i, j, T, r, l, &testfun);
    return mesh;
}

std::vector<std::vector<double>> mesh_test_linear(double origin_x, double origin_y, double end_x, double end_y,
                       double dx, double dy, int refinement_levels, double subdomain_dblder_limit){
    int i{1}, j{1}, r{1}, l{1}; // Dummy values
    double T{1}; // Dummy values
    Point origin{origin_x, origin_y}, end{end_x, end_y};
    std::vector<std::vector<double>> mesh = mesh2d(origin, end, dx, dy, refinement_levels, subdomain_dblder_limit, i, j, T, r, l, &testfun_linear);
    return mesh;
}

