#pragma once
#include <cmath>

bool draw = true;
bool verbose = false;
double eps = 0.000001;

const int Nx = 40;
const int Ny = 40;
const int Nz = 40;

// const int Nx = 120;
// const int Ny = 60;
// const int Nz = 20;

double delta_x = 1. / (Nx - 1);
double delta_y = 1. / (Ny - 1);
double delta_z = 1. / (Nz - 1);
double dx = 0.25, dy = 0.15, dz = 0.1;
double delta_t = 0.9 / (2 * (dx / std::pow(delta_x, 2) + dy / std::pow(delta_y, 2) + dz / std::pow(delta_z, 2)));
int i, j, k;