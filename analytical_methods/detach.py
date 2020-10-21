import numpy as np
from scipy.optimize import fsolve
from transitionpts import detachment, sonic, vonNuemann

"""
Function to calculate the vonNuemann, detachment and sonic condtions. The function uses calls 
from the transitionpts.py various events. It takes the input of Mach number M and specific heat
capcity g and generates the pressure ratios in term of p0/pb for the above metion condtions
"""
def flowpts(*data):
    M, g = data

    p0pj = (1 + 0.5*(g-1)*M**2)**(g/(g-1))
   
    Mguess = 0.6
    pguess = 2.
    I0 = ([2*Mguess, 1.8*Mguess, np.pi/4, pguess, np.pi/8, 1.6*Mguess, 1.5*Mguess, np.pi/4, pguess/2, np.pi/8])
    
    zv = fsolve(vonNuemann, I0, args=data, xtol=1.e-9)

    zd = fsolve(detachment, I0, args=data)

    zs = fsolve(sonic, I0, args=data)

    return p0pj/zv[3], p0pj/zd[3], p0pj/zs[3]

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
    #specific heat capacity
    g = 1.4
    #Nozzle exit height. 
    #H = 1 defines the non dimensional height 
    H = 1
    #Gas constant
    R = 287.14

    #Input the Mach desired Mach number
    M = 1.5
    #Input the desired Nozzle pressure ratio
    npr = 2.1

    #Calculation of vonNuemann and Detachment condtion
    """
    if the npr is above the vonNuemann condtion then the 
    mach reflection is impossible. Hence the vonNuemann
    condition will be computed to check the upper bound of
    the Mach reflection phenomenon.
    """
    npr_vn, npr_dt, npr_sonic = flowpts(*(M, g))
    print(npr_sonic, npr_dt, npr_vn)
    
    beta_vn = get_beta(M, g, npr_vn)
    beta_dt = get_beta(M, g, npr_dt)
    beta_sonic = get_beta(M, g, npr_sonic)

    theta_vn = get_theta(M, g, beta_vn)
    theta_dt = get_theta(M, g, beta_dt)
    theta_sonic = get_theta(M, g, beta_sonic)

    print (theta_vn*180.0/np.pi, theta_dt*180.0/np.pi, theta_sonic*180.0/np.pi)

    print(npr_vn, theta_vn*180.0/np.pi)
    print(npr_dt, theta_dt*180.0/np.pi)

    # M = 1.5   NPR_detach = 2.711