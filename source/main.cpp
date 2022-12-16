#include <iostream>
#include "creating.h"
#include "params.h"

/* 
    START CODE:
    
    pip install -r requirements.txt
    mkdir build
    cd build
    cmake ..
    make
    ./main
*/

int main(int argc, char** argv) {

    /*
    // Parsing of params
    
    if (argc == 4) {
        Nx = std::stoi(argv[1]);
        Ny = std::stoi(argv[2]);
        Nz = std::stoi(argv[3]);
    } else {
        std::cout << "WARNING:\tMesh doesn't get enough size params." << std::endl;
        return -1;
    }
    */

    MeshArray mesh(Nx, Ny, Nz);
    MeshArray real_mesh = mesh.real_solution();
    
    const clock_t begin_time = std::clock();
    MeshArray our_mesh = mesh.get_final_solution();
    std::cout << "RESULT:\t" << float(std::clock () - begin_time ) / CLOCKS_PER_SEC << "s" << std::endl;
    return 0;
}