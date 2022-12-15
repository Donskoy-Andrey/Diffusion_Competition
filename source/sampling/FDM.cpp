#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "mesh/creating.h"

/* Finite Difference Method */

double delta(int Nx, int Ny, int Nz, std::string variable){
    double value = 0;
    if (variable == "x") {
        value = 1*1.0 / (Nx-1);
    } else if (variable == "y") {
        value = 1*1.0 / (Ny-1);
    } else if (variable == "z") {
        value = 1*1.0 / (Nz-1);
    }
    return value;
}

double delta_t(int Nx, int Ny, int Nz)
{
    double dx = 0.25, dy = 0.15, dz = 0.1;
    double delta_x = delta(Nx, Ny, Nz, "x");
    double delta_y = delta(Nx, Ny, Nz, "y");
    double delta_z = delta(Nx, Ny, Nz, "z");
    double delta_t = 0.9 / ( 2 * (dx / pow(delta_x, 2) + dy / pow(delta_y, 2) + dz / pow(delta_z, 2)));
    return delta_t;
}

double diff(MeshArray & mesh, int i, int j, int k)
{
    int Nx = mesh.get_Nx();
    int Ny = mesh.get_Ny();
    int Nz = mesh.get_Nz();

    double delta_x = delta(Nx, Ny, Nz, "x");
    double delta_y = delta(Nx, Ny, Nz, "y");
    double delta_z = delta(Nx, Ny, Nz, "z");

    double LxU = (mesh(i-1, j, k) - 2 * mesh(i, j, k) + mesh(i + 1, j, k)) / pow(delta_x, 2);
    double LyU = (mesh(i, j - 1, k) - 2 * mesh(i, j, k) + mesh(i, j + 1, k)) / pow(delta_y, 2);
    double LzU = (mesh(i, j, k + 1) - 2 * mesh(i, j, k) + mesh(i, j, k + 1)) / pow(delta_z, 2);

    return 0;
}