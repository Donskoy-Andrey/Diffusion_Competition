#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>


class MeshArray {
    private:
        std::vector <double> array;
        int capacity = 0;
        int size = 0;
        int Nx = 0;
        int Ny = 0;
        int Nz = 0;

    public:
        MeshArray() = default;
        MeshArray(int Nx, int Ny, int Nz);
        inline double operator() (int i, int j, int k);
        inline void print_projection();
        inline MeshArray real_solution();
};

inline double MeshArray::operator()(int i, int j, int k){
    return MeshArray::array[i + Nx*j + Nx*Ny*k];
};

inline void MeshArray::print_projection(){
    int length;
    if (this->Nx >= 10){
        length = 10;
    }
    else {
        length = this->Nx;
    }

    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < length; ++j) {
            std::cout << std::setw(5) << MeshArray::operator()(i, j, 0) << " ";
        }
        std::cout << std::endl;
    }
};


MeshArray::MeshArray(int Nx, int Ny, int Nz){
    MeshArray::Nx = Nx;
    MeshArray::Ny = Ny;
    MeshArray::Nz = Nz;
    std::cout << "INFO:\tSize of mesh: " << Nx << "x" << Ny << "x" << Nz << std::endl;

    for (double z = 0; z <= 1; z += 1./Nz) {
        for (double y = 0; y <= 1; y += 1./Ny) {
            for (double x = 0; x <= 1; x += 1./Nx) {
                MeshArray::array.push_back(0.);
            }
        }
    }

    MeshArray::capacity = MeshArray::array.capacity();
    MeshArray::size = MeshArray::array.size();
};

double real_function(double x, double y, double z, double t = 1){
    double dx = 0.25;
    double dy = 0.15;
    double dz = 0.1;

    return std::sin(M_PI*x) * std::sin(M_PI*y) * std::sin(M_PI*z) *
        (1 - std::exp(-(dx+dy+dz)*M_PI_2*t));

};


inline MeshArray MeshArray::real_solution(){
    Nx, Ny, Nz = MeshArray::Nx, MeshArray::Ny, MeshArray::Nz;
    std::cout << Nx << " " << Ny << " " << Nz << std::endl;
    MeshArray realmesh;
    realmesh.Nx = Nx;
    realmesh.Ny = Ny;
    realmesh.Nz = Nz;

    std::ofstream file;
    file.open("../data/files/analytic_mesh.txt");

    for (double z = 0; z <= 1; z += 1./Nz) {
        for (double y = 0; y <= 1; y += 1./Ny) {
            for (double x = 0; x <= 1; x += 1./Nx) {
                double value = real_function(x, y, z);
                realmesh.array.push_back(value);
                file << std::to_string(value) + " ";
            }
        }
    }
    
    file.close();
    realmesh.capacity = realmesh.array.capacity();
    realmesh.size = realmesh.array.size();
    std::string command = "python3 ../source/mesh/drawing.py ";
    command += std::to_string(Nx) + " " + std::to_string(Ny) + " " + std::to_string(Nz);
    std::system(command.c_str());
    return realmesh;
};