#pragma once

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>

bool const DRAW = false;

int const Nx = 122;
int const Ny = 62;
int const Nz = 22;

// int const Nx = 20;
// int const Ny = 20;
// int const Nz = 20;

double const delta_x = (1./(Nx - 1));
double const delta_y = (1./(Ny - 1));
double const delta_z = (1./(Nz - 1));

double const delta_x_MPI = delta_x * M_PI;
double const delta_y_MPI = delta_y * M_PI;
double const delta_z_MPI = delta_z * M_PI;

double const delta_x_2 = delta_x * delta_x;
double const delta_y_2 = delta_y * delta_y;
double const delta_z_2 = delta_z * delta_z;

double const dx = 0.25;
double const dy = 0.15;
double const dz = 0.1;

double const delta_t = (0.9/(2 * (dx / std::pow(delta_x, 2) + dy / std::pow(delta_y, 2) + dz / std::pow(delta_z, 2))));
double const my_const = (dx + dy + dz) * M_PI * M_PI;


class MeshArray {
    private:
        /* Size and data of mesh*/
        double array[Nx*Ny*Nz] = {0.0};

    public:
        /* Create mesh s*/
        MeshArray() = default;
        ~MeshArray() = default;
        inline double const operator()(int i, int j, int k);

        /* Compute useful params */
        inline double const diff(int i, int j, int k);
        inline void next_solver();

        /* Get and draw solutions */
        inline void real_solution();
        inline void get_final_solution();
        inline void get_image(std::string & filename);
};

inline double const MeshArray::operator()(int i, int j, int k){
    return this->array[i + Nx*j + Nx*Ny*k];
}

inline double const real_function(double x, double y, double z, double t = 1){
    return std::sin(M_PI*x) * std::sin(M_PI*y) * std::sin(M_PI*z) *
        (1 - std::exp(-(dx+dy+dz)*M_PI*M_PI*t));
}

inline void MeshArray::real_solution(){
    int counter = 0;
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double x = delta_x * i;    
                double y = delta_y * j;    
                double z = delta_z * k;  
                this->array[counter] = real_function(x, y, z);;
                ++counter;
            }
        }
    }
}

inline void MeshArray::next_solver()
{
    double tmp[Nx*Ny*Nz] = {0.0};

    int counter = 0;
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                if ((i != 0) and (j != 0) and (k != 0) and
                    (i != (Nx - 1)) and (j != (Ny - 1)) and (k != (Nz - 1))){
                        tmp[counter] = this->diff(i, j, k);
                }
                ++counter;
            }
        }
    }
    std::swap(this->array, tmp);
}

inline double const MeshArray::diff(int i, int j, int k)
{
    double current = this->operator()(i, j, k);
    double LxU = (this->operator()(i - 1, j, k) - 
        2 * current + this->operator()(i + 1, j, k)) / delta_x_2;
    double LyU = (this->operator()(i, j - 1, k) - 
        2 * current + this->operator()(i, j + 1, k)) / delta_y_2;
    double LzU = (this->operator()(i, j, k - 1) - 
        2 * current + this->operator()(i, j, k + 1)) / delta_z_2;

    double x = delta_x_MPI * i;
    double y = delta_y_MPI * j;
    double z = delta_z_MPI * k;
    double f = my_const * std::sin(x) * std::sin(y) * std::sin(z);
    return current + delta_t * (f + dx*LxU + dy*LyU + dz*LzU);
}

inline void MeshArray::get_final_solution(){
    
    const clock_t begin_time = std::clock();
    double t = 0;
    while (t < 1){
        this->next_solver();
        t += delta_t;
    }

    std::cout << "---> RESULT:\t" << float(std::clock () - begin_time) / CLOCKS_PER_SEC << "s" << std::endl;

    if DRAW {
        std::string filename = "../data/files/our_mesh.txt";
        MeshArray::get_image(filename);
    }

    double max_error = 0.0;
    MeshArray meshnew;
    meshnew.real_solution();
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double x = delta_x * i;    
                double y = delta_y * j;    
                double z = delta_z * k;  
                double value = std::fabs(meshnew(i,j,k) - this->operator()(i,j,k));
                if (value > max_error) {
                    max_error = value;
                }
            }
        }
    }
    std::cout << "ERROR:\t" << "\tMAX: " << max_error << std::endl;
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
    std::ofstream file;
    file.open(filename);
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double x = delta_x * i;    
                double y = delta_y * j;    
                double z = delta_z * k;    

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