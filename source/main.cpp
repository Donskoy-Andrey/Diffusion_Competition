#include <iostream>
#include "mesh/creating.h"

int main() {
    MeshArray mesh(7, 7, 7);
    mesh.print_XtoY();
    return 0;
}