#include <iostream>
#include "mesh/creating.h"

int main(int argc, char** argv) {
    int Nx = 0;
    int Ny = 0;
    int Nz = 0;
    if (argc == 4) {
        Nx = std::stoi(argv[1]);
        Ny = std::stoi(argv[2]);
        Nz = std::stoi(argv[3]);
    } else {
        std::cout << "WARNING:\tMesh doesn't get size params. Initialise mesh with size = 0x0x0" << std::endl;
    }
    MeshArray mesh(Nx, Ny, Nz);
    mesh.print_projection();
    return 0;
}