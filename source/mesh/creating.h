#pragma once

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
        ~MeshArray() = default;
        inline MeshArray(int Nx, int Ny, int Nz);
        inline double operator() (int i, int j, int k);
        inline void print_projection();
        inline MeshArray real_solution(bool draw);
        const inline int get_Nx();
        const inline int get_Ny();
        const inline int get_Nz();
        const inline std::vector <double> get_array();
};

inline double MeshArray::operator()(int i, int j, int k){
    return MeshArray::array[i + Nx*j + Nx*Ny*k];
};

inline void MeshArray::print_projection(){
    if ((MeshArray::Nx != 0) and (MeshArray::Ny != 0) and (MeshArray::Nz != 0))
    {
        int lengthX, lengthY;
        if (this->Nx >= 10){
            lengthX = 10;
        }
        else {
            lengthX = this->Nx+1;
        }
        
        if (this->Ny >= 10){
            lengthY = 10;
        }
        else {
            lengthY = this->Ny+1;
        }

        for (int i = 0; i < lengthY; ++i) {
            for (int j = 0; j < lengthX; ++j) {
                std::cout << std::setw(5) << MeshArray::operator()(i, j, 0) << " ";
            }
            std::cout << std::endl;
        }
    }
};

const inline int MeshArray::get_Nx(){
    return this->Nx;
};

const inline int MeshArray::get_Ny(){
    return this->Ny;
};

const inline int MeshArray::get_Nz(){
    return this->Nz;
};

const inline std::vector <double> MeshArray::get_array(){
    return this->array;
};

inline MeshArray::MeshArray(int Nx, int Ny, int Nz){
    MeshArray::Nx = Nx;
    MeshArray::Ny = Ny;
    MeshArray::Nz = Nz;
    std::cout << "INFO:\tSize of mesh: " << Nx << "x" << Ny << "x" << Nz << "." << std::endl;

    for (double z = 0; z <= 1; z += 1./(Nz-1)) {
        for (double y = 0; y <= 1; y += 1./(Ny-1)) {
            for (double x = 0; x <= 1; x += 1./(Nx-1)) {
                MeshArray::array.push_back(0.);
            }
        }
    }

    MeshArray::capacity = MeshArray::array.capacity();
    MeshArray::size = MeshArray::array.size();
};

inline double real_function(double x, double y, double z, double t = 1){
    double dx = 0.25;
    double dy = 0.15;
    double dz = 0.1;
    return std::sin(M_PI*x) * std::sin(M_PI*y) * std::sin(M_PI*z) *
        (1 - std::exp(-(dx+dy+dz)*M_PI_2*t));

};

inline void data_to_vtu(std::string filename, int Nx, int Ny, int Nz){
    std::string command = "python3 ../source/mesh/drawing.py ";
    std::cout << "INFO:\tRun python script." << std::endl;
    command += std::to_string(Nx) + " " + std::to_string(Ny) + " " + std::to_string(Nz) + " " + filename;
    std::cout << "\t" << command << std::endl;
    std::system(command.c_str());
}

inline MeshArray MeshArray::real_solution(bool draw){
    Nx = this->Nx, Ny = this->Ny, Nz = this->Nz;
    std::cout << "INFO:\tCreate real mesh with size: " << Nx << "x" << Ny << "x" << Nz << "." <<  std::endl;
    this->array.clear();
    std::ofstream file;
    std::string filename = "../data/files/analytical_mesh.txt";
    file.open(filename);

    for (double z = 0; z <= 1; z += 1./(Nz-1)) {
        for (double y = 0; y <= 1; y += 1./(Ny-1)) {
            for (double x = 0; x <= 1; x += 1./(Nx-1)) {
                double value = real_function(x, y, z);
                this->array.push_back(value);
                if (draw) {
                    file << std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + " " + std::to_string(value) + "\n";
                }
            }
        }
    }
    
    file.close();

    this->capacity = this->array.capacity();
    this->size = this->array.size();

    if (draw) {
        /* Run py-script to create mesh */
        data_to_vtu(filename, Nx, Ny, Nz);
    }

    std::string commandDelete = "rm " + filename;
    std::system(commandDelete.c_str());
    return *this;
};