import numpy as np
from scipy.optimize import fsolve
from sys import exit
from transitionpts import vonNuemann, detachment, sonic
from shocktheories import threeshock
from geometry import LinBendor, mLinBendor, Mouton

"""
Prandtl-Meyer expansion relation to calculate the flow PM function for given mach number
"""
def pm(*data):
    M, g = data

    a = np.sqrt((g+1)/(g-1))
    b = np.arctan(np.sqrt((M**2 -1)/a**2))
    c = np.arctan(np.sqrt(M**2-1))

    return a*b - c


"""
Prandtl-Meyer expansion relation to calculate the flow mach number for given M
"""
def PMfunction(I0, *data):
    nu_M, g = data

    M = I0

    a = np.sqrt((g+1)/(g-1))
    b = np.arctan(np.sqrt((M**2 - 1)/a**2))
    c = np.arctan(np.sqrt(M**2-1))

    return nu_M - (a*b - c)

"""
Function to calculate the vonNuemann, detachment and sonic condtions. The function uses calls 
from the transitionpts.py various events. It takes the input of Mach number M and specific heat
capcity g and generates the pressure ratios in term of p0/pb for the above metion condtions
"""
def flowpts(*data):
    M, g = data

    p0pj = (1 + 0.5*(g-1)*M**2)**(g/(g-1))
   
    Mguess = 1.0
    pguess = 3.
    I0 = ([2*Mguess, 1.8*Mguess, np.pi/4, pguess, np.pi/8, 1.6*Mguess, 1.5*Mguess, np.pi/4, pguess/2, np.pi/8])
    
    zv = fsolve(vonNuemann, I0, args=data, xtol=1.e-9)

    zd = fsolve(detachment, I0, args=data)

    zs = fsolve(sonic, I0, args=data)

    return p0pj/zv[3], p0pj/zd[3], p0pj/zs[3]


"""
Function to calculate the three shock theory solution with appropriate boundary conditions. 
The function uses calls from the shocktheories.py. It takes the input of Mach number M and 
specific heat capacity g with npr and generates the solution of flow field close to triple point
"""

def triplepoint(*data):
    M, g, npr = data

    p0pj = (1 + 0.5*(g-1)*M**2)**(g/(g-1))
    p1 = p0pj/npr
   
    #calculation of normal shock relations for initial conditions
    p3 = (2*g*M**2 - g +1)/(g+1)
    r3 = ((g+1)*M**2)/(2 + (g-1)*M**2)
    t3 = p3/r3
    beta3 = np.pi/2
    theta3 = np.pi/10
    m3 = M/10

    #calcualtion of oblique shock conditions
    beta1 = np.pi/4
    m1 = M*0.8
    r1 = ((g+1)*m1**2*np.sin(beta1)**2)/((2+(g-1)*m1**2*np.sin(beta1)**2))
    t1 = p1/r1

    beta2 = np.pi/4
    theta2 = np.pi/8
    m2 = M*0.6
    r2 = ((g+1)*m2**2*np.sin(beta2)**2)/((2+(g-1)*m2**2*np.sin(beta2)**2))
    p2 = (1 + (2*g*(m2**2*np.sin(beta2)**2-1.))/(g+1))
    t2 = p2/r2

    data = (M, g, p1)
    # print(data)
    # I0 = [m1, m2, 0.5, beta1, beta2, beta3, theta2, theta2, np.pi/12, 15, 12, 2, 2, 3, 4, 5, 6]
    I0 = [m1, m2, m3, beta1, beta2, beta3, theta2, theta2, theta3, p2, p3, r1, r2, r3, t1, t2, t3]
    zthreeshock = fsolve(threeshock, I0, args=data)
    return zthreeshock


"""
Function to calculate the flow properties across the expansion fan in the nozzle flow.
This function takes the values of g and M2 and resolves the flow field behind the reflected
shock. It gives the flow features in the expansion region so that the flow will be parallel 
to the horizontal. The characteristic for which this happens is calculated in this routine
"""

def expansionregion(*data):
    M2, g, theta3 = data

    # Calculation of first expansion wave angle
    mu_M2 = np.arcsin(1/M2)

    # Calculation of Prandtl-Meyer function for the desired Machnumber
    nu_M2 = pm(*(M2, g))

    # Calculation of PM function at the point where flow is parallel to horizontal(D)
    nu_MD = theta3 + nu_M2

    # Calculation of Mach number from the PM fuction
    dataexp = (nu_MD, g) 
    MD = fsolve(PMfunction, M2, args=dataexp)[0]

    mu_MD = np.arcsin(1/MD)

    return nu_M2, mu_M2, MD, nu_MD, mu_MD


"""
Function to calculate the mass averaged flow in the subsonic pocket dwonstream of the Mach stem
This function takes the values of density, temperature and theta3 close to the triple point and
behind the Mach stem and calculated the mass averaged Mach number based on the normal shock
Mach number Mn from free stream Mach number
"""

def subsonicpocket(*data):
    M, g, R, M3, r3, t3 = data

    #calculation of the normal shock relations 
    Mn = np.sqrt((1 + (g-1)*M**2 + ((0.25*(g+1)**2)-g)*M**4) /
              ((g*M**2 - 0.5*(g-1))*(0.5*(g-1)*M**2 + 1)))

    pn = (2*g*M**2 -g +1)/(g+1)
    rn = ((g+1)*M**2)/(2+(g-1)*M**2)
    tn = pn/rn
    an = np.sqrt(g*R*tn)

    a3 = np.sqrt(g*R*t3)

    un = Mn*an
    u3 = M3*a3

    Mavg = 2*(r3*u3*np.cos(theta3) + rn*un)/((rn+r3)*(an+a3))

    return Mavg


"""
Function to calculate the machstem height using the geometric relationship
It uses the values calculted by three shock theory and geometric relations
given by Li and Bendor and Modified Li and Bendor method
"""
def machstem(*data):
    
    #setting up the initial conditions
    I0 = [1,1,1,1,1,1,1,1,1,1]
    zmachstem = fsolve(LinBendor, I0, args=data)
    zmachstemmodified = fsolve(mLinBendor, I0[0:8], args=data)
    zmachstemmouton = fsolve(Mouton, I0[0:6], args=data)
    return zmachstem, zmachstemmodified, zmachstemmouton
""" 
Definition for the main function. All the variable will be provided here for the computation.
The main function takes the values of the Mach number , gamma and nozzle pressure ratio to 
return the machstem height calculation using Li and Bendor's method and modified Li and Bendor method
"""
if __name__ == "__main__":
    #specific heat capacity
    g = 1.4
    #Nozzle exit height. 
    #H = 1 defines the non dimensional height 
    H = 1
    #Gas constant
    R = 287.14

    #Input the Mach desired Mach number
    M = 5
    #Input the desired Nozzle pressure ratio
    npr = 40.0

    #Calculation of vonNuemann and Detachment condtion
    """
    if the npr is above the vonNuemann condtion then the 
    mach reflection is impossible. Hence the vonNuemann
    condition will be computed to check the upper bound of
    the Mach reflection phenomenon.
    """
    npr_vn, npr_dt, npr_sonic = flowpts(*(M, g))
    
    # if npr > npr_vn :
    #     print("The nozzle pressure ratio is beyond the value of the von-Nuemann Condition")
    #     print("Please enter the value of npr below ", npr_vn)
    #     exit
    counter = 0
    # nprval = np.linspace(npr, npr_vn, num=50, endpoint=True)
    nprval = np.arange(49, 51, 1.0)
    # f1 = open('newM5_LnBMach.csv', 'a')
    # f2 = open('newM5_mLnB.csv', 'a')
    # f3  = open('newM5_Mou.csv', 'a')
    # #Calculation of solution to three shock theory for the given npr and Mach number
    for npr in nprval:
        zthreeshock = triplepoint(*(M, g, npr))

        #Set the values of the variables needed
        theta1 = zthreeshock[6]
        theta2 = zthreeshock[7]
        theta3 = zthreeshock[8]

        beta1 = zthreeshock[3]
        beta2 = zthreeshock[4]
        beta3 = zthreeshock[5]

        M1 = zthreeshock[0]
        M2 = zthreeshock[1]
        M3 = zthreeshock[2]

        p30 = zthreeshock[10]
        r30 = zthreeshock[13]
        t30 = zthreeshock[16]

        # Calculation of expansion fan interaction with the slip line
        nu_M2, mu_M2, MD, nu_MD, mu_MD = expansionregion(*(M2, g, theta3))

        # Calculation of subsonic pocket interaction
        Mbar = subsonicpocket(*(M,g,R,M3,r30,t30))

        # Calculation of Mach stem height
        zLinBendor, zmLinBendor, zMouton = machstem(*(H, mu_MD, beta1, theta1, beta2, mu_M2, theta3, g, Mbar))
        
        # valLnB = np.append(np.append([M], [npr]),
        #                    np.append(zLinBendor, [theta3]))
        # valmLnB = np.append(np.append([M], [npr]),
        #                    np.append(zmLinBendor, [theta3]))
        # valMou = np.append(np.append([M], [npr]),
        #                    np.append(zMouton, [theta3]))
        # newvalMou = np.append(valMou, zthreeshock)
        # valmLnB = np.append([M], [npr], zmLinBendor, [theta3])

        # np.savetxt(f1, [valLnB], delimiter=',', fmt='%f')
        # np.savetxt(f2, [valmLnB], delimiter=',', fmt='%f')
        # np.savetxt(f3, [newvalMou], delimiter=',', fmt='%f')
        print(npr, zLinBendor[1])
    # f1.close
    # f2.close
    # f3.close
    pass


# if __name__ == "__main__":
#     #specific heat capacity
#     g = 1.4
#     #Nozzle exit height.
#     #H = 1 defines the non dimensional height
#     H = 1
#     #Gas constant
#     R = 287.14

#     #Input the Mach desired Mach number
#     M = 3
#     #Input the desired Nozzle pressure ratio
#     npr = 7

#     #Calculation of vonNuemann and Detachment condtion
#     """
#     if the npr is above the vonNuemann condtion then the 
#     mach reflection is impossible. Hence the vonNuemann
#     condition will be computed to check the upper bound of
#     the Mach reflection phenomenon.
#     """
#     npr_vn, npr_dt, npr_sonic = flowpts(*(M, g))

#     if npr > npr_vn:
#         print("The nozzle pressure ratio is beyond the value of the von-Nuemann Condition")
#         print("Please enter the value of npr below ", npr_vn)
#         exit

#     # #Calculation of solution to three shock theory for the given npr and Mach number
#     zthreeshock = triplepoint(*(M, g, npr))

#     #Set the values of the variables needed
#     theta1 = zthreeshock[6]
#     theta2 = zthreeshock[7]
#     theta3 = zthreeshock[8]

#     beta1 = zthreeshock[3]
#     beta2 = zthreeshock[4]
#     beta3 = zthreeshock[5]

#     M1 = zthreeshock[0]
#     M2 = zthreeshock[1]
#     M3 = zthreeshock[2]

#     p30 = zthreeshock[10]
#     r30 = zthreeshock[13]
#     t30 = zthreeshock[16]

#     # Calculation of expansion fan interaction with the slip line
#     nu_M2, mu_M2, MD, nu_MD, mu_MD = expansionregion(*(M2, g, theta3))

#     # Calculation of subsonic pocket interaction
#     Mbar = subsonicpocket(*(M, g, R, M3, r30, t30))

#     # Calculation of Mach stem height
#     zLinBendor, zmLinBendor = machstem(
#         *(H, mu_MD, beta1, theta1, beta2, mu_M2, theta3, g, Mbar))

#     f1 = open('M5_LnB.csv', 'a')
#     np.savetxt(f1, [zLinBendor], delimiter=',', fmt='%f',
#                header='A Sample 2D Numpy Array :: Header', footer='This is footer')

#     f2 = open('M5_mLnB.csv', 'a')
#     np.savetxt(f2, [zmLinBendor], delimiter=',', fmt='%f',
#                header='A Sample 2D Numpy Array :: Header', footer='This is footer')
#     f1.close
#     f2.close
#     pass
