#include <iostream>
#include <vector>
#include <iomanip>


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
        inline void print_XtoY();

};

inline double MeshArray::operator()(int i, int j, int k){
    return MeshArray::array[i + Nx*j + Nx*Ny*k];
};

inline void MeshArray::print_XtoY(){
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

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz; ++k) {
                MeshArray::array.push_back(0.);
                // std::cout << "|" << std::endl;
            }
        }
    }

    MeshArray::capacity = MeshArray::array.capacity();
    MeshArray::size = MeshArray::array.size();
};