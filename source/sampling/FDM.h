#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "../mesh/creating.h"

/* Finite Difference Method */

double delta(int Nx, int Ny, int Nz, std::string variable);

double delta_t(int Nx, int Ny, int Nz);

double razn_sh(MeshArray & mesh);