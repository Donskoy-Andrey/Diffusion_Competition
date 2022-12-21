#pragma once

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <mpi.h>

// int const processor_count = 12;
// int const processor_count_x = 4;
// int const processor_count_y = 3;

int const processor_count = 4;
int const processor_count_x = 2;
int const processor_count_y = 2;

bool const GET_ERROR = false;
bool const DRAW = true;

std::string const analytical_path = "./data/files/analytical_";
std::string const output_path = "./data/files/output_";

// int const Nx = 122;
// int const Ny = 62;
// int const Nz = 22;

int const Nx = 14;
int const Ny = 14;
int const Nz = 10;

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
        double array[((Nx - 2) / processor_count_x + 2) * ((Ny - 2) / processor_count_y + 2) * Nz] = {};

    public:
        /* Create mesh s*/
        MeshArray() = default;
        ~MeshArray() = default;
        inline double const operator()(int i, int j, int k, int P);

        /* Compute useful params */
        inline double const diff(int i, int j, int k, int P, int myID, int col, int row);
        inline void next_solver(int P, int myID, int col, int row);

        /* Get and draw solutions */
        inline void real_solution(int P, int myID);
        inline void get_final_solution(int P, int myID);
        inline void get_image(std::string &filename, int P, int myID, int col, int row);

        /* Parallel process */
        void calculate_plate(int P, int myID, double * loc_array, std::string send_type, int col, int row);
};

inline double const MeshArray::operator()(int i, int j, int k, int P) {
    return this->array[i + ((Nx - 2) / processor_count_x + 2) * j + ((Nx - 2) / processor_count_x + 2) * ((Ny - 2) / processor_count_y + 2) * k];
}

inline double const real_function(double x, double y, double z, double t = 1) {
    return std::sin(M_PI * x) * std::sin(M_PI * y) * std::sin(M_PI * z) *
           (1 - std::exp(-(dx + dy + dz) * M_PI * M_PI * t));
}

inline void MeshArray::real_solution(int P, int myID) {
    int counter = 0;

    int const col = myID % (processor_count / processor_count_y); 
    int const row = myID / (processor_count / processor_count_y); 

    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < ((Ny - 2) / processor_count_y + 2); ++j) {
            for (int i = 0; i <  ((Nx - 2) / processor_count_x + 2); ++i) {
                double x = delta_x * (i +  col * (((Nx - 2) / processor_count_x + 2) - 2)); 
                double y = delta_y * (j +  row * (((Ny - 2) / processor_count_y + 2) - 2));
                double z = delta_z * k;

                this->array[counter] = real_function(x, y, z);
                ++counter;
            }
        }
    }

    if (DRAW) {
        std::string filename = analytical_path + std::to_string(myID) + ".txt";
        this->get_image(filename, P, myID, col, row);
    }
}

void MeshArray::calculate_plate(int P, int myID, double * loc_array, std::string send_type, int col, int row) {
    if (send_type == "left") { 
        std::cout << myID << " " << "LEFT" << std::endl;
        int i = 1;

        int counter = 0;
        for (int k = 0; k < Nz; ++k) {
            for (int j = 0; j < ((Ny - 2) / P + 2); ++j) {
                if ((k != 0) and (k != Nz - 1) and (j != 0) and (j != Ny - 1)) {
                    loc_array[counter] = this->diff(i, j, k, P, myID, col, row); // local i
                }
                ++counter;
            }
        }

    } else if (send_type == "right") { 
        std::cout << myID << " " << "RIGHT" << std::endl;
        int i = ((Nx - 2) / P + 2) - 2; // n - 2

        int counter = 0;
        for (int k = 0; k < Nz; ++k) {
            for (int j = 0; j < ((Ny - 2) / P + 2); ++j) {
                if ((k != 0) and (k != Nz - 1) and (j != 0) and (j != Ny - 1)) {
                    loc_array[counter] = this->diff(i, j, k, P, myID, col, row); // local i
                }
                ++counter;
            }
        }

    } else if (send_type == "up") {
        std::cout << myID << " " << "UP" << std::endl;
        int j = ((Ny - 2) / P + 2) - 2; // n - 2

        int counter = 0;
        for (int k = 0; k < Nz; ++k) {
            for (int i = 0; i < ((Nx - 2) / P + 2); ++i) {
                if ((k != 0) and (k != Nz - 1) and (j != 0) and (j != Ny - 1)) {
                    std::cout << "UP : " << counter << std::endl;
                    loc_array[counter] = this->diff(i, j, k, P, myID, col, row); // local j
                }
                ++counter;
            }
        }

    } else if (send_type == "down") {
        std::cout << myID << " " << "DOWN" << std::endl;
        int j = 1;

        int counter = 0;
        for (int k = 0; k < Nz; ++k) {
            for (int i = 0; i < ((Nx - 2) / P + 2); ++i) {
                if ((k != 0) and (k != Nz - 1) and (j != 0) and (j != Ny - 1)) {
                    loc_array[counter] = this->diff(i, j, k, P, myID, col, row); // local j
                }
                ++counter;
            }
        }
    }
}

inline void MeshArray::next_solver(int P, int myID, int col, int row) {
    double tmp[((Nx - 2) / processor_count_x + 2) * ((Ny - 2) / processor_count_y + 2) * Nz] = {};

    MPI_Request request[processor_count] = {};
    MPI_Status statuses[processor_count] = {};

    int num_request = 0;

    double loc_array_up_send[Nz*Ny] = {};
    double loc_array_down_load[Nz*Ny] = {};
    double loc_array_down_send[Nz*Ny] = {};
    double loc_array_up_load[Nz*Ny] = {};

    double loc_array_right_send[Nz*Ny] = {};
    double loc_array_left_load[Nz*Ny] = {};
    double loc_array_left_send[Nz*Ny] = {};
    double loc_array_right_load[Nz*Ny] = {};

    std::cout << myID << " ["<< col << ", " << row << "]"<< std::endl;
    if (col != processor_count_x - 1) {
        /* Calculate our right edge */
        this->calculate_plate(P, myID, loc_array_right_send, "right", col, row);

        /* Send our right edge as neighbour's left edge */
        MPI_Isend(loc_array_right_send, Nz * ((Ny - 2) / P + 2), MPI_DOUBLE, myID + 1, 0, MPI_COMM_WORLD, &request[num_request++]);
        std::cout << myID << " Loc_array_right : " << sizeof(loc_array_right_send) / sizeof(double) << std::endl;
        
        /* Get neighbour's left edge as our right edge */
        MPI_Irecv(loc_array_left_load, Nz * ((Ny - 2) / P + 2), MPI_DOUBLE, myID + 1, 0, MPI_COMM_WORLD, &request[num_request++]);
    }

    if (col != 0) {
        /* Calculate our left edge */
        this->calculate_plate(P, myID, loc_array_left_send, "left", col, row);

        /* Send our left edge as neighbour's right edge */
        MPI_Isend(loc_array_left_send, Nz * ((Ny - 2) / P + 2), MPI_DOUBLE, myID - 1, 0, MPI_COMM_WORLD, &request[num_request++]);
        std::cout << myID << " Loc_array_left : " << sizeof(loc_array_left_send) / sizeof(double) << std::endl;

        /* Get neighbour's right edge as our left edge */
        MPI_Irecv(loc_array_right_load, Nz * ((Ny - 2) / P + 2), MPI_DOUBLE, myID - 1, 0, MPI_COMM_WORLD, &request[num_request++]);
    }

    // if (row != processor_count_y - 1) {
    //     /* Calculate our up edge */
    //     this->calculate_plate(P, myID, loc_array_up_send, "up", col, row);

    //     /* Send our up edge as neighbour's down edge */
    //     MPI_Isend(loc_array_up_send, Nz * ((Nx - 2) / P + 2), MPI_DOUBLE, myID - processor_count_y, 0, MPI_COMM_WORLD, &request[num_request++]);
    //     std::cout << myID << " Loc_array_up : " << sizeof(loc_array_up_send) / sizeof(double) << std::endl;

    //     /* Get neighbour's down edge as our up edge */
    //     MPI_Irecv(loc_array_down_load, Nz * ((Nx - 2) / P + 2), MPI_DOUBLE, myID - processor_count_y, 0, MPI_COMM_WORLD, &request[num_request++]);
    // }

    if (row != 0) {
        /* Calculate our down edge */
        this->calculate_plate(P, myID, loc_array_down_send, "down", col, row);

        /* Send our down edge as neighbour's up edge */
        MPI_Isend(loc_array_down_send, Nz * ((Nx - 2) / P + 2), MPI_DOUBLE, myID + processor_count_y, 0, MPI_COMM_WORLD, &request[num_request++]);
        std::cout << myID << " Loc_array_down : " << sizeof(loc_array_down_send) / sizeof(double) << std::endl;
        
        /* Get neighbour's up edge as our down edge */
        MPI_Irecv(loc_array_up_load, Nz * ((Nx - 2) / P + 2), MPI_DOUBLE, myID + processor_count_y, 0, MPI_COMM_WORLD, &request[num_request++]);
    }

    MPI_Waitall(num_request, request, statuses);
    std::cout << "YES!" << std::endl;
    /* Write neighbour's left edge as our right edge */
    if (col != processor_count_x - 1) {
        int i = ((Nx - 2) / P) + 2 - 1; // n - 1 
        int loc_counter = 0;
        for (int k = 0; k < Nz; ++k) {
            for (int j = 0; j < ((Ny - 2) / P + 2); ++j) {
                tmp[
                    i + j * ((Nx - 2) / P + 2) + k * ((Ny - 2) / P + 2) * ((Nx - 2) / P + 2)
                ] = loc_array_left_load[loc_counter];
                ++loc_counter;
            }
        }
    }

    /* Write neighbour's right edge as our left edge */
    if (col != 0) {
        int i = 0;
        int loc_counter = 0;
        for (int k = 0; k < Nz; ++k) {
            for (int j = 0; j < ((Ny - 2) / P + 2); ++j) {
                tmp[
                    i + j * ((Nx - 2) / P + 2) + k * ((Ny - 2) / P + 2) * ((Nx - 2) / P + 2)
                ] = loc_array_right_load[loc_counter];
                ++loc_counter;
            }
        }    
    }

    /* Write neighbour's up edge as our down edge */
    if (row != processor_count_y - 1) {
        int j = 0;
        int loc_counter = 0;
        for (int k = 0; k < Nz; ++k) {
            for (int i = 0; i < ((Nx - 2) / P + 2); ++i) {
                tmp[
                    i + j * ((Nx - 2) / P + 2) + k * ((Ny - 2) / P + 2) * ((Nx - 2) / P + 2)
                ] = loc_array_up_load[loc_counter];
                ++loc_counter;
            }
        }    
    }
    
    /* Write neighbour's down edge as our up edge */
    if (row != 0) {
        int j = ((Ny - 2) / P) + 2 - 1; // n - 1 
        int loc_counter = 0;
        for (int k = 0; k < Nz; ++k) {
            for (int i = 0; i < ((Nx - 2) / P + 2); ++j) {
                tmp[
                    i + j * ((Nx - 2) / P + 2) + k * ((Ny - 2) / P + 2) * ((Nx - 2) / P + 2)
                ] = loc_array_down_load[loc_counter];
                ++loc_counter;
            }
        }    
    }

    /* Write kernel's elements */
    for (int k = 1; k < Nz - 1; ++k) {
        for (int j = 1; j < ((Ny - 2) / P + 2) - 1; ++j) {
            for (int i = 1; i < ((Nx - 2) / P + 2) - 1; ++i) {
                tmp[
                    i + j * ((Nx - 2) / P + 2) + k * ((Ny - 2) / P + 2) * ((Nx - 2) / P + 2)
                ] = this->diff(i, j, k, P, myID, col, row);
            }
        }
    }

    std::swap(this->array, tmp);
}

inline double const MeshArray::diff(int i, int j, int k, int P, int myID, int col, int row) {

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
    double y = delta_y_MPI * (j + myID * (((Ny - 2) / P + 2) - 2));
    double z = delta_z_MPI * k;
    double f = my_const * std::sin(x) * std::sin(y) * std::sin(z);

    return current + delta_t * (f + dx * LxU + dy * LyU + dz * LzU);
}

inline void MeshArray::get_final_solution(int P, int myID) {
    const double begin_time = MPI_Wtime();

    int const col = myID % (processor_count / processor_count_y); 
    int const row = myID / (processor_count / processor_count_y); 

    double t = 0;
    while (t < 1) {
        this->next_solver(P, myID, col, row);
        t += delta_t;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (myID == 0) std::cout << "---> RESULT:\t" << MPI_Wtime() - begin_time << "s" << std::endl;
   
    // if (DRAW) {
    //     std::string filename = output_path + std::to_string(myID) + ".txt";
    //     MeshArray::get_image(filename, P, myID);
    // }

    // if (GET_ERROR) {
    //     double max_error = 0.0;
    //     MeshArray meshnew;
    //     meshnew.real_solution(P, myID);
    //     for (int k = 0; k < Nz; ++k) {
    //         for (int j = 0; j < Ny; ++j) {
    //             for (int i = 0; i < ((Nx - 2) / P + 2); ++i) {
    //                 double value = std::fabs(meshnew(i, j, k, P) - this->operator()(i, j, k, P));
    //                 if (value > max_error) {
    //                     max_error = value;
    //                 }
    //             }
    //         }
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     std::cout << "ERROR:\t"
    //             << "\tMAX: " << max_error << std::endl;
    // }
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

inline void MeshArray::get_image(std::string & filename, int P, int myID, int col, int row) {
    std::ofstream file;
    file.open(filename);
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < ((Ny - 2) / processor_count_y + 2); ++j) {
            for (int i = 0; i <  ((Nx - 2) / processor_count_x + 2); ++i) {
                double x = delta_x * (i +  col * (((Nx - 2) / processor_count_x + 2) - 2)); 
                double y = delta_y * (j +  row * (((Ny - 2) / processor_count_y + 2) - 2));
                double z = delta_z * k;

                double value = this->operator()(i, j, k, P);
                file << std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + " " + std::to_string(value) + "\n";
            }
        }
    }
    file.close();
    data_to_vtu(filename);
    // std::string commandDelete = "rm " + filename;
    // std::system(commandDelete.c_str());
}