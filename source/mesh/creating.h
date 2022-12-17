#pragma once

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "../params.h"

class MeshArray {
    private:
        /* Size and data of mesh*/
        double array[Nx*Ny*Nz] = {0.0};

    public:
        /* Create mesh s*/
        MeshArray() = default;
        ~MeshArray() = default;
        inline double operator() (int i, int j, int k);

        /* Compute useful params */
        inline double diff(int i, int j, int k);
        inline void next_solver();

        /* Get and draw solutions */
        inline void real_solution();
        inline void get_final_solution();
        inline void get_image(std::string & filename);
};

inline double MeshArray::operator()(int i, int j, int k){ 
    return MeshArray::array[i + Nx * j + Nx * Ny * k];
}

inline double real_function(double x, double y, double z, double t = 1){
    return std::sin(M_PI*x) * std::sin(M_PI*y) * std::sin(M_PI*z) *
        (1 - std::exp(-(dx+dy+dz)*M_PI*M_PI*t));
}

inline void data_to_vtu(std::string & filename){
    std::string command = "python3 ../source/mesh/drawing.py ";
    command += std::to_string(Nx) + " " + std::to_string(Ny) + " " + std::to_string(Nz) + " " + filename;
    #if VERBOSE == true
        std::cout << "INFO:\tRun python script." << std::endl;
        std::cout << "\t" << command << std::endl;
    #endif
    std::system(command.c_str());
}

inline void MeshArray::get_image(std::string & filename){
    int i, j, k;
    std::ofstream file;
    file.open(filename);
    for (double z = 0; z <= 1. + eps; z += 1. / (Nz - 1)) {
        for (double y = 0; y <= 1. + eps; y += 1. / (Ny - 1)) {
            for (double x = 0; x <= 1. + eps; x += 1. / (Nx - 1)) {

                    i = x / (1. / (Nx - 1));
                    j = y / (1. / (Ny - 1));
                    k = z / (1. / (Nz - 1));

                double value = MeshArray::operator()(i, j, k);
                file << std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + " " + std::to_string(value) + "\n";
            }
        }
    }
    file.close();
    data_to_vtu(filename);
    std::string commandDelete = "rm " + filename;
    std::system(commandDelete.c_str());
}

inline void MeshArray::real_solution(){
    #if VERBOSE == true
        std::cout << "INFO:\tCreate real mesh with size: " << Nx << "x" << Ny << "x" << Nz << "." <<  std::endl;
    #endif

    int counter = 0;
    for (double z = 0; z <= 1. + eps; z += 1. / (Nz - 1)) {
        for (double y = 0; y <= 1. + eps; y += 1. / (Ny - 1)) {
            for (double x = 0; x <= 1. + eps; x += 1. / (Nx - 1)) {
                double value = real_function(x, y, z);
                MeshArray::array[counter] = value;
                ++counter;
            }
        }
    }
    
    #if DRAW == true
        std::string filename = "../data/files/analytical_mesh.txt";
        MeshArray::get_image(filename);
    #endif
}

inline void MeshArray::next_solver()
{
    int i, j, k;
    double value;
    int counter = 0;

    for (double z = 0; z <= 1. + eps; z += 1. / (Nz - 1)) {
        for (double y = 0; y <= 1. + eps; y += 1. / (Ny - 1)) {
            for (double x = 0; x <= 1. + eps; x += 1. / (Nx - 1)) {
                    i = x / (1. / (Nx - 1));
                    j = y / (1. / (Ny - 1));
                    k = z / (1. / (Nz - 1));
                if ((x == 0) or (y == 0) or (z == 0) or (x == 1) or (y == 1) or (z == 1)){
                    value = 0;   
                } else {
                    value = this->diff(i, j, k);
                }
                this->array[counter] = value;
                ++counter;
            }
        }
    }
}

inline double MeshArray::diff(int i, int j, int k)
{
    
    double LxU = (MeshArray::operator()(i - 1, j, k) - 
        2 * MeshArray::operator()(i, j, k) + MeshArray::operator()(i + 1, j, k)) / std::pow(delta_x, 2);
    double LyU = (MeshArray::operator()(i, j - 1, k) - 
        2 * MeshArray::operator()(i, j, k) + MeshArray::operator()(i, j + 1, k)) / std::pow(delta_y, 2);
    double LzU = (MeshArray::operator()(i, j, k + 1) - 
        2 * MeshArray::operator()(i, j, k) + MeshArray::operator()(i, j, k + 1)) / std::pow(delta_z, 2);

    double x = delta_x * i;
    double y = delta_y * j;
    double z = delta_z * k;
    double f = (dx + dy + dz) * M_PI * M_PI * std::sin(M_PI*x) * std::sin(M_PI*y) * std::sin(M_PI*z);
    double U_next = this->operator()(i, j, k) + delta_t * (f + dx*LxU + dy*LyU + dz*LzU);
    return U_next;
}

inline void MeshArray::get_final_solution() {
    double t = 0;
    #if VERBOSE == true
        std::cout << "INFO:\tIteration by time: "  << std::endl << "\t";        
    #endif

    while (t <= 1){
        this->next_solver();

        t += delta_t;
        #if VERBOSE == true
            std::cout << t << " | ";
        #endif
    }
    
    #if VERBOSE == true
        std::cout << std::endl;
    #endif

    #if DRAW == true
        std::string filename = "../data/files/our_mesh.txt";
        MeshArray::get_image(filename);
    #endif

    #if GET_ERROR == true
        double sum = 0;
        double counter = 0.;
        double max_value = 0;
        MeshArray real_mesh;
        real_mesh.real_solution();
        for (int i = 0; i < Nx*Ny*Nz; ++i){
            std::cout << this->array[i] << " " << real_mesh.array[i] << std::endl;
            double value = std::fabs(this->array[i] - real_mesh.array[i]);
            sum += value;
            if (value > max_value){
                max_value = value;
            }
            ++counter;
        }
        std::cout << sum << " " << counter << std::endl;
        std::cout << "ERROR:\t" << "MEAN: " << sum / counter << "\tMAX: " << max_value << std::endl;
    #endif
}