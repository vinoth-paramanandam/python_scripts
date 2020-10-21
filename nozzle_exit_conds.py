import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

def third_critical_pr(M, g):
    return (1 + 0.5*(g-1)*M**2)**(g/(g-1))

def get_oblique_pr(M, g, beta):
    return (2*g*M**2*np.sin(beta)**2 - (g-1))/(g+1)

def get_beta(M, g, NPR):
    pr = get_shock_jump(M, g, NPR)
    return np.arcsin(np.sqrt((pr*(g+1) + (g-1))/(2*g*M**2)))

def get_shock_jump(M, g, NPR):
    return third_critical_pr(M, g)/NPR

def get_theta(M, g, beta):
    nume = (((g+1)*M**2)/(2*(M**2*np.sin(beta)**2 -1)) -1)*np.tan(beta)
    return np.arctan(1/nume)

def get_M(M, g, beta): 
    nume = (g+1)**2*M**4*np.sin(beta)**2 - 4*(M**2*np.sin(beta)**2 -1)*(g*M**2*np.sin(beta)**2 + 1)
    deno = (2*g*M**2*np.sin(beta)**2 - (g-1))*((g-1)*M**2*np.sin(beta)**2 + 2)
    return np.sqrt(nume/deno)

def get_npr(M, g, beta):
    return third_critical_pr(M, g)/get_oblique_pr(M, g, beta)

if __name__ == "__main__":
    global g, M

    M = 5.0
    g = 1.4
    beta = 39.3*np.pi/180.0

    theta = get_theta(M, g, beta)
    NPR = get_npr(M, g, beta)
    print(NPR, beta*180.0/np.pi, theta*180.0/np.pi)

    NPR = 45
    beta1 = get_beta(M, g, NPR)
    theta1 = get_theta(M, g, beta1)

    print(NPR, beta1*180.0/np.pi, theta1*180.0/np.pi)

    NPR = 44
    beta1 = get_beta(M, g, NPR)
    theta1 = get_theta(M, g, beta1)

    print(NPR, beta1*180.0/np.pi, theta1*180.0/np.pi)

    NPR = 43
    beta2 = get_beta(M, g, NPR)
    theta2 = get_theta(M, g, beta2)

    print(NPR, beta2*180.0/np.pi, theta2*180.0/np.pi)

    NPR = 43.5
    beta2 = get_beta(M, g, NPR)
    theta2 = get_theta(M, g, beta2)

    print(NPR, beta2*180.0/np.pi, theta2*180.0/np.pi)

    NPR = 50.79401850
    betaval = get_beta(5.0, 1.4, NPR)
    thetaval = get_theta(5.0, 1.4, betaval)

    print(NPR, np.rad2deg(betaval), np.rad2deg(thetaval))

    # NPR = np.arange(40, 71, 0.5)

    # nprfinal = np.zeros((NPR.shape[0], 4), dtype=float)
    # beta = get_beta(M, g, NPR)*180.0/np.pi
    # theta = get_theta(M, g, beta*np.pi/180.0)*180.0/np.pi
    # pr = get_shock_jump(M, g, NPR)
    # pj = 101325.0/pr
    # tj = 300.0/(1 + 0.5*(g-1)*M**2)
    
    # nprfinal[:, 0] = NPR
    # nprfinal[:, 1] = pj
    # nprfinal[:, 2] = beta
    # nprfinal[:, 3] = theta
    # np.savetxt('data.csv', nprfinal, delimiter=',')
    pass