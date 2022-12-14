import sys
import numpy as np
from pyevtk.hl import pointsToVTK


def numpyToVTK(xs, ys, zs, values):
    xs = np.array(xs)
    ys = np.array(ys)
    zs = np.array(zs)
    values = np.array(values)

    pointsToVTK(
        "../data/files/analytic_mesh", 
        xs, ys, zs, 
        data = {'U': values}
    )

def main():
    Nx = int(sys.argv[1]) + 1
    Ny = int(sys.argv[2]) + 1
    Nz = int(sys.argv[3]) + 1
    xs, ys, zs, values = [], [], [], []
    with open(r"../data/files/analytic_mesh.txt", 'r') as file:
        for line in file.readlines():
            x, y, z, value = [float(i) for i in line.split()]
            xs.append(x)
            ys.append(y)
            zs.append(z)
            values.append(value)

    numpyToVTK(xs, ys, zs, values)
        

if __name__ == "__main__":
    main()