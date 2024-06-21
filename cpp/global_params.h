#pragma once
#include <string>

constexpr double BOLTZMANN = 1.38064852e-23;
constexpr double GAS_CONSTANT = 8.31446261815324;
constexpr double AVOGADRO = GAS_CONSTANT / BOLTZMANN;
constexpr double PI = 3.14159265359;
constexpr double PLANCK = 6.62607015e-34;
constexpr double HBAR = PLANCK / (2.0 * PI);
constexpr double FLTEPS = 1e-10;

extern std::string fluid_dir;
void set_fluid_dir(const std::string& path);