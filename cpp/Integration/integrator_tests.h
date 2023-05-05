#pragma once
#include "Integration.h"
#include <vector>
#include <algorithm>
#include <memory>
#include <utility>
#include <map>
#include <math.h>
#include <functional>

// mesh_step and mesh_adaptive are together the exact same algorithm as integration_step and integrate_adaptive
// Defined in Integrator.h. Rather than compute the integral, mesh_adaptive returns a vector of all the integration
// points, in the order they are evaluated.
// mesh2d is the corresponding "mirror function" for integrate2d.
void mesh_step(std::shared_ptr<Point>& p1, std::shared_ptr<Point>& p2, std::shared_ptr<Point>& p3,
                        int& Nx, int& Ny,
                        const double& dx, const double& dy,
                        int& Nxsteps, const int& Nysteps,
                        const double subdomain_dblder_limit,
                        std::map<std::pair<int, int>, const double>& evaluated_points,
                        const int& arg_i, const int& arg_j, const double& arg_T, const int& arg_l, const int& arg_r,
                        std::function<double(int, int, double, double, double, int, int)> func, std::vector<Point>& points);

void mesh_adaptive(const Point& origin,
                    const int& Nx_origin, const int& Ny_origin,
                    const int& Nx_end, const int& Ny_end,
                    const double& dx, const double& dy,
                    int& Nxsteps, const int& Nysteps,
                    const double& subdomain_dblder_limit,
                    std::map<std::pair<int, int>, const double>& evaluated_points,
                    const int& arg_i, const int& arg_j, const double& arg_T, const int& arg_l, const int& arg_r,
                    std::function<double(int, int, double, double, double, int, int)> func, std::vector<Point>& points);

std::vector<std::vector<double>> mesh2d(const Point& origin, const Point& end,
                                        const double& dx, const double& dy,
                                        const int& refinement_levels,
                                        const double& subdomain_dblder_limit,
                                        const int& arg_i, const int& arg_j, const double& arg_T, const int& arg_l, const int& arg_r,
                                        std::function<double(int, int, double, double, double, int, int)> func);

/* trisurf_step :
    Mirror function for integration_step. See : trisurf_adaptive.
*/
void trisurf_step(std::shared_ptr<Point>& p1, std::shared_ptr<Point>& p2, std::shared_ptr<Point>& p3,
                        int& Nx, int& Ny, std::vector<std::vector<Point>>& surf,
                        const double& dx, const double& dy,
                        int& Nxsteps, const int& Nysteps,
                        const double subdomain_dblder_limit,
                        std::map<std::pair<int, int>, const double>& evaluated_points,
                        const int& arg_i, const int& arg_j, const double& arg_T, const int& arg_l, const int& arg_r,
                        std::function<double(int, int, double, double, double, int, int)> func);

/*  trisurf_adaptive :
    Follows the exact same algorithm as integrate_adaptive, but instead of integrating three points at a time, appends them
    to the vector "surf".
    This becomes an (N, 3) vector corresponding to each triangle that was integrated.
    The dimentions are then (Triangle, point) such that surf[0, 1] is the 2nd point of the first triangle that was integrated.
*/
void trisurf_adaptive(const Point& origin,
                        const int& Nx_origin, const int& Ny_origin,
                        const int& Nx_end, const int& Ny_end,
                        const double& dx, const double& dy,
                        int& Nxsteps, const int& Nysteps,
                        const double& subdomain_dblder_limit,
                        std::map<std::pair<int, int>, const double>& evaluated_points,
                        const int& arg_i, const int& arg_j, const double& arg_T, const int& arg_l, const int& arg_r,
                        std::function<double(int, int, double, double, double, int, int)> func,
                        std::vector<std::vector<Point>>& surf);

/* trisurf :
    Interface to trisurf_adaptive. Takes the (N, 3) vector of points supplied by trisurf_adaptive and converts it
    to an (N, 3, 3) vector of doubles. The dimentions are then (Triangle, Point, coordinate) such that
    verts[0, 1, 2] gives the 3rd cooridinate (z) of the 2nd point in the 1st triangle that was integrated.
*/
std::vector<std::vector<std::vector<double>>> trisurf(const Point& origin, const Point& end,
                                                        const double& dx, const double& dy,
                                                        const int& refinement_levels_x, const int& refinement_levels_y,
                                                        const double& subdomain_dblder_limit,
                                                        const int& arg_i, const int& arg_j, const double& arg_T, const int& arg_l, const int& arg_r,
                                                        std::function<double(int, int, double, double, double, int, int)> func);


double testfun(int i, int j, double T, double x, double y, int r, int l);
double testfun_linear(int i, int j, double T, double x, double y, int r, int l);
double integrator_test(double origin_x, double origin_y, double end_x, double end_y,
                       double dx, double dy, int refinement_levels, double subdomain_dblder_limit);
double integrator_test_linear(double origin_x, double origin_y, double end_x, double end_y,
                       double dx, double dy, int refinement_levels, double subdomain_dblder_limit);
std::vector<std::vector<double>> mesh_test(double origin_x, double origin_y, double end_x, double end_y,
                                            double dx, double dy, int refinement_levels, double subdomain_dblder_limit);