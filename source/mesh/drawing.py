import sys
import numpy as np
from pyevtk.hl import pointsToVTK

verbose = False

def numpyToVTK(xs: list, ys: list, zs: list, values: list, filename: str) -> None:
    xs = np.array(xs)
    ys = np.array(ys)
    zs = np.array(zs)
    values = np.array(values)
    output = filename.split('/')[-1].split('.')[-2]
    pointsToVTK(
        f"../data/files/{output}", 
        xs, ys, zs, 
        data = {'U': values}
    )
    if verbose:
        print(f"INFO:\tImage saved ({values.shape[0]} points).")

def main():
    Nx = int(sys.argv[1])
    Ny = int(sys.argv[2])
    Nz = int(sys.argv[3])
    filename = sys.argv[4]
    xs, ys, zs, values = [], [], [], []
    with open(filename, 'r') as file:
        for line in file.readlines():
            x, y, z, value = [float(i) for i in line.split()]
            xs.append(x)
            ys.append(y)
            zs.append(z)
            values.append(value)

    numpyToVTK(xs, ys, zs, values, filename)
        

if __name__ == "__main__":
    main()