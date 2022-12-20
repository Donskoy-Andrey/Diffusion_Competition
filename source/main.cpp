#include "mesh/creating.h"

/*
    mpic++ source/main.cpp -O2 && mpirun -np 4 ./a.out
*/

int main(int argc, char** argv) {
    int P, myID; 
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &P);
    MPI_Comm_rank (MPI_COMM_WORLD, &myID);
    MPI_Barrier(MPI_COMM_WORLD);
        MeshArray mesh;
        mesh.get_final_solution(P, myID);
    MPI_Finalize();
    return 0;
}