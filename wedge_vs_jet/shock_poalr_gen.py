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

if __name__ == "__main__":
    global g, M

    M = 5.0
    g = 1.4
    NPR = 41
    
    npr = np.arange(41, 60, 0.25, dtype=np.float)

    # beta_init = np.arcsin(1/M)
    # beta_final = np.pi/2

    # beta = np.linspace(beta_init, beta_final, num=1000, endpoint=True)
    # theta = get_theta(M, g, beta)
    # pr = get_oblique_pr(M, g, beta)

    # beta1 = get_beta(M, g, NPR)
    # M1 = get_M(M, g, beta1)
    # thetaadd = get_theta(M, g, beta1)
    # pradd = get_oblique_pr(M, g, beta1)
    
    # beta_init = np.arcsin(1./M)
    # beta_final = np.pi/2
    # beta = np.linspace(beta_init, beta_final, num=1000, endpoint=True)

    # theta1 = get_theta(M1, g, beta)
    # pr1 = get_oblique_pr(M1, g, beta1)  
    pass