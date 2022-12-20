#pragma once

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <mpi.h>

int const processor_count = 1;
bool const GET_ERROR = true;
bool const DRAW = false;

std::string const analytical_path = "./data/files/analytical_";
std::string const output_path = "./data/files/output_";

int const Nx = 122;
int const Ny = 62;
int const Nz = 22;

// int const Nx = 14;
// int const Ny = 10;
// int const Nz = 10;

double const delta_x = (1. / (Nx - 1));
double const delta_y = (1. / (Ny - 1));
double const delta_z = (1. / (Nz - 1));

double const delta_x_MPI = delta_x * M_PI;
double const delta_y_MPI = delta_y * M_PI;
double const delta_z_MPI = delta_z * M_PI;

double const delta_x_2 = delta_x * delta_x;
double const delta_y_2 = delta_y * delta_y;
double const delta_z_2 = delta_z * delta_z;

double const dx = 0.25;
double const dy = 0.15;
double const dz = 0.1;

double const delta_t = (0.9 / (2 * (dx / std::pow(delta_x, 2) + dy / std::pow(delta_y, 2) + dz / std::pow(delta_z, 2))));
double const my_const = (dx + dy + dz) * M_PI * M_PI;

class MeshArray {
    private:
        /* Size and data of mesh*/
        double array[((Nx - 2) / processor_count + 2) * Ny * Nz] = {};

    public:
        /* Create mesh s*/
        MeshArray() = default;
        ~MeshArray() = default;
        inline double const operator()(int i, int j, int k, int P);

        /* Compute useful params */
        inline double const diff(int i, int j, int k, int P, int myID);
        inline void next_solver(int P, int myID);

        /* Get and draw solutions */
        inline void real_solution(int P, int myID);
        inline void get_final_solution(int P, int myID);
        inline void get_image(std::string &filename, int P, int myID);

        /* Parallel process */
        void calculate_plate(int P, int myID, double * loc_array, std::string i_type);
};

inline double const MeshArray::operator()(int i, int j, int k, int P) {
    return this->array[i + ((Nx - 2) / P + 2) * j + ((Nx - 2) / P + 2) * Ny * k];
}

inline double const real_function(double x, double y, double z, double t = 1) {
    return std::sin(M_PI * x) * std::sin(M_PI * y) * std::sin(M_PI * z) *
           (1 - std::exp(-(dx + dy + dz) * M_PI * M_PI * t));
}

inline void MeshArray::real_solution(int P, int myID) {
    int counter = 0;
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < ((Nx - 2) / P + 2); ++i) {

                double x = delta_x * (i + myID * (((Nx - 2) / P + 2) - 2));
                double y = delta_y * j;
                double z = delta_z * k;

                this->array[counter] = real_function(x, y, z);
                ++counter;
            }
        }
    }

    if (DRAW) {
        std::string filename = analytical_path + std::to_string(myID) + ".txt";
        this->get_image(filename, P, myID);
    }
}

void MeshArray::calculate_plate(int P, int myID, double * loc_array, std::string i_type) {
    int i = 1;
    if (i_type == "right") { 
        i = ((Nx - 2) / P + 2) - 2; // n - 2
    }
    int counter = 0;
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            if ((k != 0) and (k != Nz-1) and (j != 0) and (j != Ny - 1)) {
                loc_array[counter] = this->diff(i, j, k, P, myID); // local i
            }
            ++counter;
        }
    }
}

inline void MeshArray::next_solver(int P, int myID) {
    double tmp[((Nx - 2) / processor_count + 2) * Ny * Nz] = {};

    MPI_Request request[4] = {};
    MPI_Status statuses[4] = {};

    int num_request = 0;

    double loc_array_right_send[Nz*Ny] = {};
    double loc_array_left_load[Nz*Ny] = {};
    double loc_array_left_send[Nz*Ny] = {};
    double loc_array_right_load[Nz*Ny] = {};

    if (myID != (P - 1)) {
        /* Calculate our right edge */
        this->calculate_plate(P, myID, loc_array_right_send, "right");

        /* Send our right edge as neighbour's left edge */
        MPI_Isend(loc_array_right_send, Nz*Ny, MPI_DOUBLE, myID + 1, 0, MPI_COMM_WORLD, &request[num_request++]);
        
        /* Get neighbour's left edge as our right edge */
        MPI_Irecv(loc_array_left_load, Nz*Ny, MPI_DOUBLE, myID + 1, 0, MPI_COMM_WORLD, &request[num_request++]);
    }

    if (myID != 0) {
        /* Calculate our left edge */
        this->calculate_plate(P, myID, loc_array_left_send, "left");

        /* Send our left edge as neighbour's right edge */
        MPI_Isend(loc_array_left_send, Nz * Ny, MPI_DOUBLE, myID - 1, 0, MPI_COMM_WORLD, &request[num_request++]);

        /* Get neighbour's right edge as our left edge */
        MPI_Irecv(loc_array_right_load, Nz * Ny, MPI_DOUBLE, myID - 1, 0, MPI_COMM_WORLD, &request[num_request++]);
    }

    MPI_Waitall(num_request, request, statuses);
    
    /* Write neighbour's left edge as our right edge */
    if (myID != (P - 1)) {
        int i = ((Nx - 2) / P) + 2 - 1; // n - 1 
        int loc_counter = 0;
        for (int k = 0; k < Nz; ++k) {
            for (int j = 0; j < Ny; ++j) {
                tmp[
                    i + j * ((Nx - 2) / P + 2) + k * Ny * ((Nx - 2) / P + 2)
                ] = loc_array_left_load[loc_counter];
                ++loc_counter;
            }
        }
    }

    /* Write neighbour's right edge as our left edge */
    if (myID != 0) {
        int i = 0;
        int loc_counter = 0;
        for (int k = 0; k < Nz; ++k) {
            for (int j = 0; j < Ny; ++j) {
                tmp[
                    i + j * ((Nx - 2) / P + 2) + k * Ny * ((Nx - 2) / P + 2)
                ] = loc_array_right_load[loc_counter];
                ++loc_counter;
            }
        }    
    }

    /* Write kernel's elements */
    for (int k = 1; k < Nz - 1; ++k) {
        for (int j = 1; j < Ny - 1; ++j) {
            for (int i = 1; i < ((Nx - 2) / P + 2) - 1; ++i) {
                tmp[
                    i + j * ((Nx - 2) / P + 2) + k * Ny * ((Nx - 2) / P + 2)
                ] = this->diff(i, j, k, P, myID);
            }
        }
    }

    std::swap(this->array, tmp);
}

inline double const MeshArray::diff(int i, int j, int k, int P, int myID) {

    double current = this->operator()(i, j, k, P);

    double LxU = (this->operator()(i - 1, j, k, P) -
                  2 * current + this->operator()(i + 1, j, k, P)) /
                 delta_x_2;
    double LyU = (this->operator()(i, j - 1, k, P) -
                  2 * current + this->operator()(i, j + 1, k, P)) /
                 delta_y_2;
    double LzU = (this->operator()(i, j, k - 1, P) -
                  2 * current + this->operator()(i, j, k + 1, P)) /
                 delta_z_2;

    double x = delta_x_MPI * (i + myID * (((Nx - 2) / P + 2) - 2));
    double y = delta_y_MPI * j;
    double z = delta_z_MPI * k;
    double f = my_const * std::sin(x) * std::sin(y) * std::sin(z);

    return current + delta_t * (f + dx * LxU + dy * LyU + dz * LzU);
}

inline void MeshArray::get_final_solution(int P, int myID) {
    const double begin_time = MPI_Wtime();

    double t = 0;
    while (t < 1) {
        this->next_solver(P, myID);
        t += delta_t;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (myID == 0) std::cout << "---> RESULT:\t" << MPI_Wtime() - begin_time << "s" << std::endl;
   
    if (DRAW) {
        std::string filename = output_path + std::to_string(myID) + ".txt";
        MeshArray::get_image(filename, P, myID);
    }

    if (GET_ERROR) {
        double max_error = 0.0;
        MeshArray meshnew;
        meshnew.real_solution(P, myID);
        for (int k = 0; k < Nz; ++k) {
            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < ((Nx - 2) / P + 2); ++i) {
                    double value = std::fabs(meshnew(i, j, k, P) - this->operator()(i, j, k, P));
                    if (value > max_error) {
                        max_error = value;
                    }
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout << "ERROR:\t"
                << "\tMAX: " << max_error << std::endl;
    }
}

inline void data_to_vtu(std::string & filename) {
    std::string command = "python3 ./source/mesh/drawing.py ";
    command += std::to_string(Nx) + " " + std::to_string(Ny) + " " + std::to_string(Nz) + " " + filename;
    #if VERBOSE
        std::cout << "INFO:\tRun python script." << std::endl;
        std::cout << "\t" << command << std::endl;
    #endif
    std::system(command.c_str());
}

inline void MeshArray::get_image(std::string & filename, int P, int myID) {
    std::ofstream file;
    file.open(filename);
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i <  ((Nx - 2) / P + 2); ++i) {
                double x = delta_x * (i + myID * (((Nx - 2) / P + 2) - 2)); 
                double y = delta_y * j;
                double z = delta_z * k;

                double value = this->operator()(i, j, k, P);
                file << std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + " " + std::to_string(value) + "\n";
            }
        }
    }
    file.close();
    data_to_vtu(filename);
    std::string commandDelete = "rm " + filename;
    std::system(commandDelete.c_str());
}