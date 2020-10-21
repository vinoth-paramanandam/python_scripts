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

def threeshock_sonic(I0, *data):
    M0, g = data
   
    M1 = I0[0]
    
    Phi1 = I0[1]
    Phi2 = I0[2]
    Phim = I0[3]

    M2 = I0[4]
    
    f = np.zeros((5),dtype=float)
    
    f[0] =M1-np.sqrt(((1+((g-1)*M0**2*np.sin(Phi1)*np.sin(Phi1)))+((((g+1)**2/4)-g*np.sin(Phi1)*np.sin(Phi1))*M0**4*(np.sin(Phi1))**2))/(((g*M0**2*(np.sin(Phi1))**2)-((g-1)/2))*(((M0**2*(g-1)/2)*(np.sin(Phi1))**2)+1)))
    f[1] =P10-((2/(g+1))*((g*M0**2*(np.sin(Phi1))**2)-(g-1)/2))
    f[2] =M2-np.sqrt(((1+((g-1)*M1**2*np.sin(Phi2)*np.sin(Phi2)))+((((g+1)**2/4)-g*np.sin(Phi2)*np.sin(Phi2))*M1**4*(np.sin(Phi2))**2))/(((g*M1**2*(np.sin(Phi2))**2)-((g-1)/2))*(((M1**2*(g-1)/2)*(np.sin(Phi2))**2)+1)))
    f[3] =P21-((2/(g+1))*((g*M1**2*(np.sin(Phi2))**2)-(g-1)/2))
    f[4] =M3-np.sqrt(((1+((g-1)*M0**2*np.sin(Phi3)**2))+((((g+1)**2/4)-g*np.sin(Phi3)**2)*M0**4*(np.sin(Phi3))**2))/(((g*M0**2*(np.sin(Phi3))**2)-((g-1)/2))*(((M0**2*(g-1)/2)*(np.sin(Phi3))**2)+1)))
    f[5] =P30-((2/(g+1))*((g*M0**2*(np.sin(Phi3))**2)-(g-1)/2))
    f[6] =P30-P21*P10
    f[7] = np.sin(Phi2) - np.sqrt((1/(g*M1**2))*((g+1)*0.25*M1**2 - (3 - g)*0.25
                                                  + np.sqrt((g+1)*((g+1)*0.0625*M1**4 + (g-3)*0.125*M1**2 + (g+9)/16))))

    return f


def refsonicthree(I0, data):
    M0, g = data
   
    M1 = I0[0]
    M2 = I0[1]
    M3 = I0[2]
    
    Phi1 = I0[3]    
    Phi2 = I0[4]    
    Phi3 = I0[5]
    
    Theta1 = I0[6]
    Theta2 = I0[7]    
    Theta3 = I0[8]

    P21 = I0[9]
    P30 = I0[10]
    
    Rho10 = I0[11]
    Rho21 = I0[12]
    Rho30 = I0[13]
    
    T10 = I0[14]
    T21 = I0[15]
    T30 = I0[16]
    
    f = np.zeros((17),dtype=float)
    
    f[0] =M1-np.sqrt(((1+((g-1)*M0**2*np.sin(Phi1)*np.sin(Phi1)))+((((g+1)**2/4)-g*np.sin(Phi1)*np.sin(Phi1))*M0**4*(np.sin(Phi1))**2))/(((g*M0**2*(np.sin(Phi1))**2)-((g-1)/2))*(((M0**2*(g-1)/2)*(np.sin(Phi1))**2)+1)))
    f[1] =Theta1-np.arctan((2*(1/np.tan(Phi1))*(M0**2*(np.sin(Phi1))**2-1))/(M0**2*(g+np.cos(2*Phi1))+2))
    f[2] =P10-((2/(g+1))*((g*M0**2*(np.sin(Phi1))**2)-(g-1)/2))
    f[3] =Rho10-(((g+1)*M0**2*(np.sin(Phi1))**2)/(((g-1)*M0**2*(np.sin(Phi1))**2)+2))
    f[4] =T10-(P10/(Rho10))
    f[5] =M2-np.sqrt(((1+((g-1)*M1**2*np.sin(Phi2)*np.sin(Phi2)))+((((g+1)**2/4)-g*np.sin(Phi2)*np.sin(Phi2))*M1**4*(np.sin(Phi2))**2))/(((g*M1**2*(np.sin(Phi2))**2)-((g-1)/2))*(((M1**2*(g-1)/2)*(np.sin(Phi2))**2)+1)))
    f[6] =Theta2-np.arctan((2*(1/np.tan(Phi2))*(M1**2*(np.sin(Phi2))**2-1))/(M1**2*(g+np.cos(2*Phi2))+2))
    f[7] =P21-((2/(g+1))*((g*M1**2*(np.sin(Phi2))**2)-(g-1)/2))
    f[8] =Rho21-(((g+1)*M1**2*(np.sin(Phi2))**2)/(((g-1)*M1**2*(np.sin(Phi2))**2)+2))
    f[9] =T21-(P21/(Rho21))
    f[10]=M3-np.sqrt(((1+((g-1)*M0**2*np.sin(Phi3)**2))+((((g+1)**2/4)-g*np.sin(Phi3)**2)*M0**4*(np.sin(Phi3))**2))/(((g*M0**2*(np.sin(Phi3))**2)-((g-1)/2))*(((M0**2*(g-1)/2)*(np.sin(Phi3))**2)+1)))
    f[11]=Theta3-np.arctan((2*(1/np.tan(Phi3))*(M0**2*(np.sin(Phi3))**2-1))/(M0**2*(g+np.cos(2*Phi3))+2))
    f[12]=P30-((2/(g+1))*((g*M0**2*(np.sin(Phi3))**2)-(g-1)/2))
    f[13]=Rho30-(((g+1)*M0**2*(np.sin(Phi3))**2)/(((g-1)*M0**2*(np.sin(Phi3))**2)+2))
    f[14]=T30-(P30/(Rho30))
    f[15]=Theta3-Theta1+Theta2
    f[16]=P30-P21*P10

    return f