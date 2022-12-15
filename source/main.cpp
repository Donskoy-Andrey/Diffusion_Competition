#include <iostream>
#include "mesh/creating.h"

/* 
    START CODE:
    
    pip install -r requirements.txt
    mkdir build
    cd build
    cmake ..
    make
    ./main 5 5 5
*/
int main(int argc, char** argv) {
    int Nx = 0;
    int Ny = 0;
    int Nz = 0;

    if (argc == 4) {
        Nx = std::stoi(argv[1]);
        Ny = std::stoi(argv[2]);
        Nz = std::stoi(argv[3]);
    } else {
        std::cout << "WARNING:\tMesh doesn't get enough size params. Initialise mesh with size = 0x0x0." << std::endl;
        return -1;
    }

    MeshArray mesh(Nx, Ny, Nz);
    MeshArray real_mesh = mesh.real_solution(true);

    /*
    // Examples to get params of mesh.
    std::cout << real_mesh.get_Nx() << std::endl; // x-dimension (amount of nodes on x-axis)
    std::cout << real_mesh.get_Ny() << std::endl; // y-dimension (amount of nodes on y-axis) 
    std::cout << real_mesh.get_Nz() << std::endl; // z-dimension (amount of nodes on z-axis)
    std::cout << real_mesh.get_array() << std::endl; // values of all nodes
    std::cout << real_mesh(3,3,3) << std::endl;  // element with index x=3, y=3, z=3
    */

    return 0;
}