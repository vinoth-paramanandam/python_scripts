import numpy as np
from numpy import sin, cos, tan, pi, arctan

def threeshock(I0, *data):
    M0, g, P10 = data
   
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
    