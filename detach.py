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
   
    Mguess = 0.8
    pguess = 2.
    I0 = ([2*Mguess, 1.8*Mguess, np.pi/4, pguess, np.pi/8, 1.6*Mguess, 1.5*Mguess, np.pi/4, pguess/2, np.pi/8])
    
    zv = fsolve(vonNuemann, I0, args=data, xtol=1.e-9)

    zd = fsolve(detachment, I0, args=data)

    zs = fsolve(sonic, I0, args=data)

    return p0pj/zv[3], p0pj/zd[3], p0pj/zs[3]


if __name__ == "__main__":
    #specific heat capacity
    g = 1.4
    #Nozzle exit height. 
    #H = 1 defines the non dimensional height 
    H = 1
    #Gas constant
    R = 287.14

    #Input the Mach desired Mach number
    M = 2.0

    #Calculation of vonNuemann and Detachment condtion
    """
    if the npr is above the vonNuemann condtion then the 
    mach reflection is impossible. Hence the vonNuemann
    condition will be computed to check the upper bound of
    the Mach reflection phenomenon.
    """
    
    Mguess = 1.0
    pguess = 2.
    I0 = ([2*Mguess, 1.8*Mguess, np.pi/4, pguess, np.pi/8, 1.6*Mguess, 1.5*Mguess, np.pi/4, pguess/2, np.pi/8])
    
    for i in range(8):
        data = (M, g)
        p0pj = (1 + 0.5*(g-1)*M**2)**(g/(g-1))
        zd = fsolve(detachment, I0, args=data)
        print(i+1, M, p0pj/zd[3])
        M = M - 0.1
        I0 = zd