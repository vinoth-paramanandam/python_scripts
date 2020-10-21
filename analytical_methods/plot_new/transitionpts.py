import numpy as np


def vonNuemann(I0, *data):

    M, g = data

    f = np.empty((10))

    M0 = I0[0]
    M1 = I0[1]
    beta = I0[2]
    pr1 = I0[3]
    theta = I0[4]

    M10 = I0[5]
    M2 = I0[6]
    beta1 = I0[7]
    pr2 = I0[8]
    theta1 = I0[9]

    f[0] = M0 - M*np.sin(beta)
    f[1] = M1 - (1/np.sin(beta - theta)) * \
        np.sqrt((1 + (g-1)*0.5*(M0**2))/(g*M0**2 - (g-1)*0.5))
    f[2] = pr1 - (1 + (2*g*(M0**2-1.))/(g+1))
    f[3] = 0.5*np.tan(theta) - (M0**2 - 1) / \
        (np.tan(beta)*(M**2*(g + np.cos(2*beta))+2))

    f[4] = M10 - M1*np.sin(beta1)
    f[5] = M2 - (1/np.sin(beta1 - theta1)) * \
        np.sqrt((1 + (g-1)*0.5*(M10**2))/(g*M10**2 - (g-1)*0.5))
    f[6] = pr2 - (1 + (2*g*(M10**2-1.))/(g+1))
    f[7] = 0.5*np.tan(theta1) - (M10**2 - 1) / \
        (np.tan(beta1)*(M1**2*(g + np.cos(2*beta1))+2))

    f[8] = pr1*pr2 - (1 + (2*g*(M**2-1.))/(g+1))
    f[9] = theta1 - theta

    return f


def detachment(I0, *data):

    M, g = data

    f = np.empty((10))

    M0 = I0[0]
    M1 = I0[1]
    beta = I0[2]
    pr1 = I0[3]
    theta = I0[4]

    M10 = I0[5]
    M2 = I0[6]
    beta1 = I0[7]
    pr2 = I0[8]
    theta1 = I0[9]

    f[0] = M0 - M*np.sin(beta)
    f[1] = M1 - (1/np.sin(beta - theta)) * \
        np.sqrt((1 + (g-1)*0.5*(M0**2))/(g*M0**2 - (g-1)*0.5))
    f[2] = pr1 - (1 + (2*g*(M0**2-1.))/(g+1))
    f[3] = 0.5*np.tan(theta) - (M0**2 - 1) / \
        (np.tan(beta)*(M**2*(g + np.cos(2*beta))+2))

    f[4] = M10 - M1*np.sin(beta1)
    f[5] = M2 - (1/np.sin(beta1 - theta1)) * \
        np.sqrt((1 + (g-1)*0.5*(M10**2))/(g*M10**2 - (g-1)*0.5))
    f[6] = pr2 - (1 + (2*g*(M10**2-1.))/(g+1))
    f[7] = 0.5*np.tan(theta1) - (M10**2 - 1) / \
        (np.tan(beta1)*(M1**2*(g + np.cos(2*beta1))+2))

    f[8] = np.sin(beta1) - np.sqrt((1/(g*M1**2))*((g+1)*0.25*M1**2
                                                  - 1 + np.sqrt((g+1)*((g+1)*0.0625*M1**4 + (g-1)*0.5*M1**2 + 1))))
    f[9] = theta1 - theta

    return f


def sonic(I0, *data):

    M, g = data

    f = np.empty((10))

    M0 = I0[0]
    M1 = I0[1]
    beta = I0[2]
    pr1 = I0[3]
    theta = I0[4]

    M10 = I0[5]
    M2 = I0[6]
    beta1 = I0[7]
    pr2 = I0[8]
    theta1 = I0[9]

    f[0] = M0 - M*np.sin(beta)
    f[1] = M1 - (1/np.sin(beta - theta)) * \
        np.sqrt((1 + (g-1)*0.5*(M0**2))/(g*M0**2 - (g-1)*0.5))
    f[2] = pr1 - (1 + (2*g*(M0**2-1.))/(g+1))
    f[3] = 0.5*np.tan(theta) - (M0**2 - 1) / \
        (np.tan(beta)*(M**2*(g + np.cos(2*beta))+2))

    f[4] = M10 - M1*np.sin(beta1)
    f[5] = M2 - (1/np.sin(beta1 - theta1)) * \
        np.sqrt((1 + (g-1)*0.5*(M10**2))/(g*M10**2 - (g-1)*0.5))
    f[6] = pr2 - (1 + (2*g*(M10**2-1.))/(g+1))
    f[7] = 0.5*np.tan(theta1) - (M10**2 - 1) / \
        (np.tan(beta1)*(M1**2*(g + np.cos(2*beta1))+2))

    f[8] = np.sin(beta1) - np.sqrt((1/(g*M1**2))*((g+1)*0.25*M1**2 - (3 - g)*0.25
                                                  + np.sqrt((g+1)*((g+1)*0.0625*M1**4 + (g-3)*0.125*M1**2 + (g+9)/16))))
    f[9] = theta1 - theta

    return f
