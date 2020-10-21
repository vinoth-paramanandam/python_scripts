import numpy as np
import numba
import matplotlib.pyplot as plt
import time
import os 


def readtecplotblckfile(filename):
    #opening the output file as infile
    with open(filename, 'r') as infile:
        #reading the header data and using it
        infile.readline()
        value = infile.readline().split()
        nx = int(value[3])
        ny = int(value[6])
#reading the x and y data values
        xval = infile.readline().split()
        yval = infile.readline().split()
#converting the list to numpy arrays
        x = np.asarray(xval, dtype=np.float64).reshape(nx, ny, order='F')
        y = np.asarray(yval, dtype=np.float64).reshape(nx, ny, order='F')

#reading the cell center variables
        rhoval = infile.readline().split()
        uval = infile.readline().split()
        vval = infile.readline().split()
        pval = infile.readline().split()
#converting the list to numpy array
        rho = np.asarray(rhoval, dtype=np.float64).reshape(
            nx-1, ny-1, order='F')
        u = np.asarray(uval, dtype=np.float64).reshape(nx-1, ny-1, order='F')
        v = np.asarray(vval, dtype=np.float64).reshape(nx-1, ny-1, order='F')
        p = np.asarray(pval, dtype=np.float64).reshape(nx-1, ny-1, order='F')

        return x, y, rho, u, v, p


def cellcenter2node(var):
    nx = var.shape[0] + 1
    ny = var.shape[1] + 1

    varnode = np.zeros((nx, ny), order='F')

    for j in range(1, ny-1):
        for i in range(1, nx-1):
            varnode[i, j] = (var[i, j] + var[i-1, j] +
                             var[i, j-1] + var[i-1, j-1])*0.25

    varnode[0, :] = 2*varnode[1, :] - varnode[2, :]
    varnode[nx-1, :] = 2*varnode[nx-2, :] - varnode[nx-3, :]
    varnode[:, 0] = 2*varnode[:, 1] - varnode[:, 2]
    varnode[:, ny-1] = 2*varnode[:, ny-2] - varnode[:, ny-3]

    return varnode


@numba.jit(nopython=True, cache=True)
def find_lsq(x, y, var, gradx, grady):
    nx = var.shape[0]
    ny = var.shape[1]

    for j in range(1, ny-1):
        for i in range(1, nx-1):
            a11 = x[i+1, j] - x[i, j]
            a12 = x[i-1, j] - x[i, j]
            a13 = x[i, j+1] - x[i, j]
            a14 = x[i, j-1] - x[i, j]
            b11 = y[i+1, j] - y[i, j]
            b12 = y[i-1, j] - y[i, j]
            b13 = y[i, j+1] - y[i, j]
            b14 = y[i, j-1] - y[i, j]
            c1 = var[i+1, j] - var[i, j]
            c2 = var[i-1, j] - var[i, j]
            c3 = var[i, j+1] - var[i, j]
            c4 = var[i, j-1] - var[i, j]
            A = np.array([[a11, b11], [a12, b12], [a13, b13], [a14, b14]])
            B = np.array([[c1], [c2], [c3], [c4]])

            grad, a, b, c = np.linalg.lstsq(A, B, rcond=-1)
            gradx[i, j] = grad[0][0]
            grady[i, j] = grad[1][0]

    return gradx, grady


def lst_gradient(x, y, var):
        nx = var.shape[0]
        ny = var.shape[1]

        gradnodex = np.zeros((nx, ny), order='F')
        gradnodey = np.zeros((nx, ny), order='F')

        gradnodex, gradnodey = find_lsq(x, y, var, gradnodex, gradnodey)

        gradnodex[0, :] = gradnodex[1, :]
        gradnodey[0, :] = gradnodey[1, :]

        gradnodex[nx-1, :] = gradnodex[nx-2, :]
        gradnodey[nx-1, :] = gradnodey[nx-2, :]

        gradnodex[:, 0] = gradnodex[:, 1]
        gradnodey[:, 0] = gradnodey[:, 1]

        gradnodex[:, ny-1] = gradnodex[:, ny-2]
        gradnodey[:, ny-1] = gradnodey[:, ny-2]

        return gradnodex, gradnodey


if __name__ == "__main__":
    files = [f for f in os.listdir(os.curdir) if os.path.isfile(f) and f.split('.')[1] == 'dat']
    print(files)
    # #starting the time counter
    start = time.time()

    for inputfile in files:
        print(inputfile)
        #filename to be analysed
        filename = inputfile

        #Output filename
        outputfile = filename.split('.')[0]

        # # extracting the data from the file
        x, y, rho, u, v, p = readtecplotblckfile(filename)

        # # converting the cell data to node data
        rhonode = cellcenter2node(rho)

        # #finding the least sqaure gradient for density
        rgradx, rgrady = lst_gradient(x, y, rhonode)

        # # applying the Quirk method of calculating the gradient
        rgrad = np.sqrt(np.square(rgradx) + np.square(rgrady))
        rgradmin = np.amin(rgrad)
        rgradmax = np.amax(rgrad)
        K = 15./(rgradmax - rgradmin)
        rgradfinal = np.exp(-K*(rgrad - rgradmin))

        # #saving the output as numpy file
        np.save(outputfile, rgradfinal)
#         # np.savetxt('test.out', x)

    # #end of the time counter
    end = time.time()

    # #printing the time for the computation
    print(end-start)

    pass