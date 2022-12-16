// #include <iostream>
// #include <vector>
// #include <cmath>
// #include <fstream>
// #include "creating.h"
// #include "FDM.h"


// double diff(MeshArray & mesh, int i, int j, int k, 
//     double delta_x, double delta_y, double delta_z, double delta_t,
//     double dx, double dy, double dz)
// {
//     double LxU = (mesh(i - 1, j, k) - 2 * mesh(i, j, k) + mesh(i + 1, j, k)) / pow(delta_x, 2);
//     double LyU = (mesh(i, j - 1, k) - 2 * mesh(i, j, k) + mesh(i, j + 1, k)) / pow(delta_y, 2);
//     double LzU = (mesh(i, j, k + 1) - 2 * mesh(i, j, k) + mesh(i, j, k + 1)) / pow(delta_z, 2);

//     double x = delta_x * i;
//     double y = delta_y * j;
//     double z = delta_z * k;
//     double f = (dx + dy + dz) * M_PI * M_PI * std::sin(M_PI*x) * std::sin(M_PI*y) * std::sin(M_PI*z);
//     double U_next = mesh(i, j, k) + delta_t * (f + dx*LxU + dy*LyU + dz*LzU);
//     return U_next;
// }

// MeshArray get_final_solution(MeshArray & mesh){
//     double dx = 0.25, dy = 0.15, dz = 0.1;
//     int Nx = mesh.get_Nx();
//     int Ny = mesh.get_Ny();
//     int Nz = mesh.get_Nz();

//     double delta_x = 1. / (Nx - 1);
//     double delta_y = 1. / (Ny - 1);
//     double delta_z = 1. / (Nz - 1);
//     double delta_t = 0.9 / ( 2 * (dx / pow(delta_x, 2) + dy / pow(delta_y, 2) + dz / pow(delta_z, 2)));

//     double t = 0;
//     while (t <= 1){
//         // mesh = mesh.next_solver(delta_x, delta_y, delta_z, delta_t, dx, dy, dz);
//         t += delta_t;
//     }
//     return mesh;
// }
