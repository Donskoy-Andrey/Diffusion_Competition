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
        std::vector <double> array;

    public:
        /* Create mesh s*/
        MeshArray();
        ~MeshArray() = default;
        inline double const operator()(int i, int j, int k);

        /* Compute useful params */
        inline double const diff(int i, int j, int k);
        inline void next_solver();

        /* Get and draw solutions */
        inline void real_solution();
        inline void get_final_solution();
        inline void get_image(std::string  & filename);
};

MeshArray::MeshArray(){
    int counter = 0;
    for (int i = 0; i < Nx*Ny*Nz; ++i){
        this->array.push_back(0.0);
    }
}

inline double const MeshArray::operator()(int i, int j, int k){
    return this->array[i + Nx*j + Nx*Ny*k];
}

inline double const real_function(double x, double y, double z, double t = 1){
    return std::sin(M_PI*x) * std::sin(M_PI*y) * std::sin(M_PI*z) *
        (1 - std::exp(-(dx+dy+dz)*M_PI*M_PI*t));
}

inline void data_to_vtu(std::string & filename){
    std::string command = "python3 ../source/mesh/drawing.py ";
    command += std::to_string(Nx) + " " + std::to_string(Ny) + " " + std::to_string(Nz) + " " + filename;
    #if VERBOSE
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

                double value = this->operator()(i, j, k);
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
    #if VERBOSE
        std::cout << "INFO:\tCreate real mesh with size: " << Nx << "x" << Ny << "x" << Nz << "." <<  std::endl;
    #endif
    
    int counter = 0;
    for (double z = 0; z <= 1. + eps; z += 1. / (Nz - 1)) {
        for (double y = 0; y <= 1. + eps; y += 1. / (Ny - 1)) {
            for (double x = 0; x <= 1. + eps; x += 1. / (Nx - 1)) {
                double value = real_function(x, y, z);
                this->array[counter] = value;
                ++counter;
            }
        }
    }

    #if DRAW
        std::string filename = "../data/files/analytical_mesh.txt";
        MeshArray::get_image(filename);
    #endif
}

inline void MeshArray::next_solver()
{
    double value = 0;
    std::vector <double> tmp;
    tmp.clear();
    for (double z = 0; z <= 1. + eps; z += 1. / (Nz-1)) {
        for (double y = 0; y <= 1. + eps; y += 1. / (Ny-1)) {
            for (double x = 0; x <= 1. + eps; x += 1. / (Nx-1)) {
                    int i = x / (1. / (Nx - 1));
                    int j = y / (1. / (Ny - 1));
                    int k = z / (1. / (Nz - 1));
                if ((std::fabs(x - 0) < eps) or 
                    (std::fabs(y - 0) < eps) or
                    (std::fabs(z - 0) < eps) or
                    (std::fabs(x - 1) < eps) or
                    (std::fabs(y - 1) < eps) or
                    (std::fabs(z - 1) < eps)){
                    value = 0;   
                } else {
                    value = this->diff(i, j, k);
                }
                tmp.push_back(value);
            }
        }
    }
    this->array.clear();
    for (int i = 0; i < Nx*Ny*Nz; ++i) {
        this->array.push_back(tmp[i]);
    }
}

inline double const MeshArray::diff(int i, int j, int k)
{
    double LxU = (this->operator()(i - 1, j, k) - 
        2 * this->operator()(i, j, k) + this->operator()(i + 1, j, k)) / (delta_x * delta_x);
    double LyU = (this->operator()(i, j - 1, k) - 
        2 * this->operator()(i, j, k) + this->operator()(i, j + 1, k)) / (delta_y * delta_y);
    double LzU = (this->operator()(i, j, k - 1) - 
        2 * this->operator()(i, j, k) + this->operator()(i, j, k + 1)) / (delta_z * delta_z);

    double x = delta_x * i;
    double y = delta_y * j;
    double z = delta_z * k;
    double f = (dx + dy + dz) * M_PI * M_PI * std::sin(M_PI*x) * std::sin(M_PI*y) * std::sin(M_PI*z);
    double U_next = this->operator()(i, j, k) + delta_t * (f + dx*LxU + dy*LyU + dz*LzU);
    return U_next;
}

inline void MeshArray::get_final_solution(){
    
    const clock_t begin_time = std::clock();
    double t = 0;

    #if VERBOSE
        std::cout << "INFO:\tIteration by time: "  << std::endl << "\t";        
    #endif
    int counter = 0;
    while (t < 1){
        this->next_solver();
        t += delta_t;
        #if VERBOSE
            std::cout << t  << " | ";
        #endif
        ++counter;
    }

    #if VERBOSE
        std::cout << std::endl;
    #endif

    std::cout << "---> RESULT:\t" << float(std::clock () - begin_time ) / CLOCKS_PER_SEC << "s" << std::endl;

    #if DRAW
        std::string filename = "../data/files/our_mesh.txt";
        MeshArray::get_image(filename);
    #endif

    #if GET_ERROR
        double max_error = 0.0;
        MeshArray meshnew;
        meshnew.real_solution();
        for (double z = 0; z <= 1. + eps; z += 1. / (Nz-1)) {
            for (double y = 0; y <= 1. + eps; y += 1. / (Ny-1)) {
                for (double x = 0; x <= 1. + eps; x += 1. / (Nx-1)) {
                        int i = x / (1. / (Nx - 1));
                        int j = y / (1. / (Ny - 1));
                        int k = z / (1. / (Nz - 1));
                    // std::cout << meshnew(i,j,k) << " " << this->operator()(i,j,k) << std::endl;
                    double value = std::fabs(meshnew(i,j,k) - this->operator()(i,j,k));
                    if (value > max_error) {
                        max_error = value;
                    }
                }
            }
        }
        std::cout << "ERROR:\t" << "\tMAX: " << max_error << std::endl;
    #endif
}