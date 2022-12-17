#pragma once

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <algorithm>

#define my_const (dx + dy + dz) * M_PI * M_PI

// double const eps = 1e-6;

int const Nx = 120;
int const Ny = 60;
int const Nz = 20;

// int const Nx = 20;
// int const Ny = 20;
// int const Nz = 20;

double const delta_x = (1./(Nx - 1));
double const delta_y = (1./(Ny - 1));
double const delta_z = (1./(Nz - 1));

double const dx = 0.25;
double const dy = 0.15;
double const dz = 0.1;

double const delta_t = (0.9/(2 * (dx / std::pow(delta_x, 2) + dy / std::pow(delta_y, 2) + dz / std::pow(delta_z, 2))));

class MeshArray {
    private:
        /* Size and data of mesh*/
        // std::vector <double> array;
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
                double value = real_function(x, y, z);
                this->array[counter] = value;
                ++counter;
            }
        }
    }
}

inline void MeshArray::next_solver()
{
    double value = 0;
    double tmp[Nx*Ny*Nz] = {0.0};

    int counter = 0;
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                if ((i == 0) or 
                    (j == 0) or
                    (k == 0) or
                    (i == (Nx - 1)) or
                    (j == (Ny - 1)) or
                    (k == (Nz - 1))){
                    tmp[counter] = 0;
                } else {
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
        2 * current + this->operator()(i + 1, j, k)) / (delta_x * delta_x);
    double LyU = (this->operator()(i, j - 1, k) - 
        2 * current + this->operator()(i, j + 1, k)) / (delta_y * delta_y);
    double LzU = (this->operator()(i, j, k - 1) - 
        2 * current + this->operator()(i, j, k + 1)) / (delta_z * delta_z);

    double x = delta_x * i;
    double y = delta_y * j;
    double z = delta_z * k;
    double f = my_const * std::sin(M_PI*x) * std::sin(M_PI*y) * std::sin(M_PI*z);
    return current + delta_t * (f + dx*LxU + dy*LyU + dz*LzU);
}

inline void MeshArray::get_final_solution(){
    
    const clock_t begin_time = std::clock();
    double t = 0;
    while (t < 1){
        this->next_solver();
        t += delta_t;
    }

    std::cout << "---> RESULT:\t" << float(std::clock () - begin_time ) / CLOCKS_PER_SEC << "s" << std::endl;

    double max_error = 0.0;
    MeshArray meshnew;
    meshnew.real_solution();
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double x = delta_x * i;    
                double y = delta_y * j;    
                double z = delta_z * k;  
                // std::cout << meshnew(i,j,k) << " " << this->operator()(i,j,k) << std::endl;
                double value = std::fabs(meshnew(i,j,k) - this->operator()(i,j,k));
                if (value > max_error) {
                    max_error = value;
                }
            }
        }
    }
    std::cout << "ERROR:\t" << "\tMAX: " << max_error << std::endl;
}