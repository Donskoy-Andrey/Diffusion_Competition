import sys
import numpy as np
from pyevtk.hl import gridToVTK


def numpyToVTK(data, Nx, Ny, Nz):
    x = np.linspace(0, 1, Nx)
    y = np.linspace(0, 1, Ny)
    z = np.linspace(0, 1, Nz)
    noSlices = 5

    gridToVTK("../data/files/analytic_mesh", x, y, z, cellData = {'U': data})

def main():
    Nx = int(sys.argv[1]) + 1
    Ny = int(sys.argv[2]) + 1
    Nz = int(sys.argv[3]) + 1
    with open(r"../data/files/analytic_mesh.txt", 'r') as file:
        values = np.array([float(i) for i in file.readlines()[0].split()])
    values = values.reshape((Nx, Ny, Nz))
    print(values[:,:,0])
    print(values[:,:,1])
    print(values[:,:,2])
    print(values[:,:,3])
    print(values[:,:,4])
    print(values[:,:,5])

    numpyToVTK(values, Nx, Ny, Nz)
        

if __name__ == "__main__":
    main()