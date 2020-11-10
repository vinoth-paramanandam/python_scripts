import numpy as np
 

if __name__ == "__main__":
    
    x = np.zeros(3)
    y = np.zeros(3)
    y_interpolating = np.zeros(3)

    x[0] = 0.0
    x[1] = 2.0
    x[2] = 3.0

    y[0] = 1.0
    y[1] = 2.0
    y[2] = 4.0

    x_interpolating = 1.0

    for i in range(x.shape[0]):
        print(i)
    pass