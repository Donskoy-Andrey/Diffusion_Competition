#pragma once

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "FDM.h"

class MeshArray {
    private:
        std::vector <double> array;
        int Nx = 0;
        int Ny = 0;
        int Nz = 0;

    public:
        MeshArray() = default;
        ~MeshArray() = default;
        inline MeshArray(int Nx, int Ny, int Nz);
        inline double operator() (int i, int j, int k);
        inline MeshArray real_solution(bool draw);
        inline MeshArray next_solver(
            double delta_x, double delta_y, double delta_z, double delta_t,
            double dx, double dy, double dz
        );
        
        inline double diff(int i, int j, int k, 
            double delta_x, double delta_y, double delta_z, double delta_t,
            double dx, double dy, double dz);
            
        inline MeshArray get_final_solution(bool draw);
        const inline int get_Nx();
        const inline int get_Ny();
        const inline int get_Nz();
        const inline std::vector <double> get_array();
        inline void get_image(std::string filename);
};

inline double MeshArray::operator()(int i, int j, int k){
    int Nx = MeshArray::Nx; 
    int Ny = MeshArray::Ny; 
    int Nz = MeshArray::Nz;
    return MeshArray::array[i + (Nx-1)*j + (Nx-1)*(Ny-1)*k];
}

const inline int MeshArray::get_Nx(){
    return this->Nx;
}

const inline int MeshArray::get_Ny(){
    return this->Ny;
}

const inline int MeshArray::get_Nz(){
    return this->Nz;
}

const inline std::vector <double> MeshArray::get_array(){
    return this->array;
}

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
}

inline double real_function(double x, double y, double z, double t = 1){
    double dx = 0.25, dy = 0.15, dz = 0.1;
    return std::sin(M_PI*x) * std::sin(M_PI*y) * std::sin(M_PI*z) *
        (1 - std::exp(-(dx+dy+dz)*M_PI*M_PI*t));

}

inline void data_to_vtu(std::string filename, int Nx, int Ny, int Nz){
    std::string command = "python3 ../source/mesh/drawing.py ";
    std::cout << "INFO:\tRun python script." << std::endl;
    command += std::to_string(Nx) + " " + std::to_string(Ny) + " " + std::to_string(Nz) + " " + filename;
    std::cout << "\t" << command << std::endl;
    std::system(command.c_str());
}

inline void MeshArray::get_image(std::string filename){
    int Nx = MeshArray::Nx; 
    int Ny = MeshArray::Ny; 
    int Nz = MeshArray::Nz;

    std::ofstream file;
    file.open(filename);
    for (double z = 0; z <= 1; z += 1. / (Nz - 1)) {
        for (double y = 0; y <= 1; y += 1. / (Ny - 1)) {
            for (double x = 0; x <= 1; x += 1. / (Nx - 1)) {

                    int i = x / (1. / (Nx - 1));
                    int j = y / (1. / (Ny - 1));
                    int k = z / (1. / (Nz - 1));

                double value = MeshArray::operator()(i, j, k);
                file << std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + " " + std::to_string(value) + "\n";
            }
        }
    }

    file.close();

    data_to_vtu(filename, Nx, Ny, Nz);

    std::string commandDelete = "rm " + filename;
    std::system(commandDelete.c_str());
}

inline MeshArray MeshArray::real_solution(bool draw){
    int Nx = MeshArray::Nx; 
    int Ny = MeshArray::Ny; 
    int Nz = MeshArray::Nz;
    MeshArray::array.clear();
    std::cout << "INFO:\tCreate real mesh with size: " << Nx << "x" << Ny << "x" << Nz << "." <<  std::endl;

    for (double z = 0; z <= 1; z += 1. / (Nz - 1)) {
        for (double y = 0; y <= 1; y += 1. / (Ny - 1)) {
            for (double x = 0; x <= 1; x += 1. / (Nx - 1)) {
                double value = real_function(x, y, z);
                MeshArray::array.push_back(value);

                    int i = x / (1. / (Nx - 1));
                    int j = y / (1. / (Ny - 1));
                    int k = z / (1. / (Nz - 1));
            }
        }
    }

    if (draw) {
        std::string filename = "../data/files/analytical_mesh.txt";
        MeshArray::get_image(filename);
    }

    return *this;
}


inline MeshArray MeshArray::next_solver(
    double delta_x, double delta_y, double delta_z, double delta_t,
    double dx, double dy, double dz)
{
    int Nx = this->Nx; 
    int Ny = this->Ny; 
    int Nz = this->Nz;
    double value;

    std::vector <double> tmp;

    int i, j, k;

    for (double z = 0; z <= 1; z += 1. / (Nz-1)) {
        for (double y = 0; y <= 1; y += 1. / (Ny-1)) {
            for (double x = 0; x <= 1; x += 1. / (Nx-1)) {
                i = x / (1. / (Nx - 1));
                j = y / (1. / (Ny - 1));
                k = z / (1. / (Nz - 1));
                if (((i % Nx) * (j%Ny) * (k%Nz)) == 0){
                    value = 0;   
                } else {
                    value = this->diff(i, j, k, delta_x, delta_y, delta_z, delta_t, dx, dy, dz);
                }
                tmp.push_back(value);
            }
        }
    }
    this->array = tmp;
    return *this;
}


inline double MeshArray::diff(int i, int j, int k, 
    double delta_x, double delta_y, double delta_z, double delta_t,
    double dx, double dy, double dz)
{
    double LxU = (this->operator()(i - 1, j, k) - 2 * this->operator()(i, j, k) + this->operator()(i + 1, j, k)) / std::pow(delta_x, 2);
    double LyU = (this->operator()(i, j - 1, k) - 2 * this->operator()(i, j, k) + this->operator()(i, j + 1, k)) / std::pow(delta_y, 2);
    double LzU = (this->operator()(i, j, k + 1) - 2 * this->operator()(i, j, k) + this->operator()(i, j, k + 1)) / std::pow(delta_z, 2);

    double x = delta_x * i;
    double y = delta_y * j;
    double z = delta_z * k;
    double f = (dx + dy + dz) * M_PI * M_PI * std::sin(M_PI*x) * std::sin(M_PI*y) * std::sin(M_PI*z);
    double U_next = this->operator()(i, j, k) + delta_t * (f + dx*LxU + dy*LyU + dz*LzU);
    return U_next;
}

inline MeshArray MeshArray::get_final_solution(bool draw){
    double dx = 0.25, dy = 0.15, dz = 0.1;
    int Nx = this->Nx; 
    int Ny = this->Ny; 
    int Nz = this->Nz;

    double delta_x = 1. / (Nx - 1);
    double delta_y = 1. / (Ny - 1);
    double delta_z = 1. / (Nz - 1);
    double delta_t = 0.9 / (2 * (dx / std::pow(delta_x, 2) + dy / std::pow(delta_y, 2) + dz / std::pow(delta_z, 2)));

    double t = 0;
    MeshArray mesh(Nx, Ny, Nz);
    std::cout << "INFO:\tIteration by time: "  << std::endl << "\t";
    while (t <= 1){
        mesh = mesh.next_solver(delta_x, delta_y, delta_z, delta_t, dx, dy, dz);
        t += delta_t;
        std::cout << t  << " | ";
    }
    std::cout << std::endl;

    if (draw) {
        std::string filename = "../data/files/our_mesh.txt";
        MeshArray::get_image(filename);
    }
    return mesh;
}
