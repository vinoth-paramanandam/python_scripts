import numpy as np
from scipy.optimize import fsolve

def theta_beta_m(*data):
    M, g, beta = data

    f = np.arctan(2.0*((M**2*np.sin(beta)**2 - 1) / \
        (np.tan(beta)*(M**2*(g + np.cos(2*beta))+2))))

    return f

if __name__ == "__main__":
    #specific heat capacity
    g = 1.4
    
    #Gas constant
    R = 287.14

    #Input the Mach desired Mach number
    M = 1.5

    beta = np.linspace(np.arcsin(1./M), np.pi/2, endpoint=True, num=10000)
    theta = theta_beta_m(*(M, g, beta))
    print(np.max(theta)*180.0/np.pi)

    # for i in range(beta.size):
    #     print(beta[i]*180.0/np.pi, theta[i]*180.0/np.pi)
    # pass

    # beta1 = np.arcsin(1/M)
    # theta1 = theta_beta_m(*(M, g, beta1))
    # print(beta1*180/np.pi, theta1*180/np.pi)