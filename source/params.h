#pragma once
#include <cmath>

#define DRAW true
#define VERBOSE true
#define GET_ERROR true
#define eps 0.000001

#define Nx 40
#define Ny 40
#define Nz 40

// #define Nx 120
// #define Ny 60
// #define Nz 20

#define delta_x 1./(Nx - 1)
#define delta_y 1./(Ny - 1)
#define delta_z 1./(Nz - 1)

#define dx 0.25
#define dy 0.15
#define dz 0.1

#define delta_t 0.9/(2 * (dx / std::pow(delta_x, 2) + dy / std::pow(delta_y, 2) + dz / std::pow(delta_z, 2)))