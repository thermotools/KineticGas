#pragma once
#include <string>

#include <vector>

constexpr double BOLTZMANN = 1.38064852e-23;
constexpr double GAS_CONSTANT = 8.31446261815324;
constexpr double AVOGADRO = GAS_CONSTANT / BOLTZMANN;
constexpr double PI = 3.14159265359;
constexpr double PLANCK = 6.62607015e-34;
constexpr double HBAR = PLANCK / (2.0 * PI);
constexpr double FLTEPS = 1e-10;

using vector1d = std::vector<double>;
using vector2d = std::vector<vector1d>;
using vector3d = std::vector<vector2d>;