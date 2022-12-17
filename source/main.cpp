#include <iostream>
#include "mesh/creating.h"
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
    MeshArray mesh2;

    mesh2.real_solution();
    
    MeshArray mesh;
    const clock_t begin_time = std::clock();
        mesh.get_final_solution();
    std::cout << "RESULT:\t" << float(std::clock () - begin_time ) / CLOCKS_PER_SEC << "s" << std::endl;
    return 0;
}