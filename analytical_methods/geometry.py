import numpy as np

def LinBendor(I0, *data):
    H, mu_MD, beta1, theta1, beta2, mu_M2, theta3, g, Mavg = data

    xt = I0[0]
    Hm = I0[1]
    xd = I0[2]
    yd = I0[3]
    xe = I0[4]
    Hs = I0[5]
    xf = I0[6]
    yf = I0[7]
    ye = I0[8]
    yt = I0[9]
    
    f = np.empty(10)
    
    
    f[0] = xt*np.tan(beta1) - (H - yt)
    f[1] = xd*np.tan(theta1) - (H - yd)
    f[2] = (xd - xt)*np.tan(beta2 - theta1) - (yd - yt)
    f[3] = (xf - xd)*np.tan(mu_M2 + theta3) + (yf - yd) # check this eqn, i think the present version that only length matters works very well here
    f[4] = (xe - xd)*np.tan(mu_MD) - (yd - ye)
    f[5] = (xf - xt)*np.tan(-theta3) - (yf - yt)
    f[6] = (xf - xe)*(np.tan(-theta3)/(2+ np.tan(theta3)**2)) - (yf - Hs)
    f[7] = yt - (ye/Mavg)*((2/(g+1))*(1 + (g-1)*0.5*Mavg**2))**((g+1)/(2*(g-1)))
    f[8] = yt - Hm
    f[9] = ye - Hs

    return f


def mLinBendor(I0, *data):
    H, mu_MD, beta1, theta1, beta2, mu_M2, theta3, g, Mavg = data

    xt = I0[0]
    Hm = I0[1]
    xd = I0[2]
    yd = I0[3]
    Hs = I0[4]
    xf = I0[5]
    yf = I0[6]
    yt = I0[7]

    f = np.empty(8)

    f[0] = xt*np.tan(beta1) - (H - yt)
    f[1] = xd*np.tan(theta1) - (H - yd)
    f[2] = (xd - xt)*np.tan(beta2 - theta1) - (yd - yt)
    f[3] = (xf - xd)*np.tan(mu_M2 + theta3) + (yf - yd)
    f[4] = (xf - xt)*np.tan(-theta3) - (yf - yt)
    f[5] = yt - (yf/Mavg)*((2/(g+1)) *
                           (1 + (g-1)*0.5*Mavg**2))**((g+1)/(2*(g-1)))
    f[6] = yt - Hm
    f[7] = yf - Hs

    return f


def Mouton(I0, *data):
    H, mu_MD, beta1, theta1, beta2, mu_M2, theta3, g, Mavg = data

    ot = I0[0]
    et = I0[1]
    ef = I0[2]
    oe = I0[3]
    sh = I0[4]
    sstar = I0[5]
    
    f = np.empty(6)
    
    
    f[0] = ot*np.sin(beta1) + sh - 1
    f[1] = et*np.sin(beta2- theta1) + sh - ef*np.sin(mu_MD) - sstar
    f[2] = oe*np.cos(theta1) + ef*np.cos(mu_MD) - ot*np.cos(beta1) - (sh - sstar)/np.tan(theta3)
    f[3] = ot*np.sin(beta1) - oe*np.sin(theta1) - et*np.sin(beta2 - theta1)
    f[4] = ot*np.cos(beta1) + et*np.cos(beta2 - theta1) - oe*np.cos(theta1)
    f[5] = sh - (sstar/Mavg)*((2/(g+1))*(1 + (g-1)*0.5*Mavg**2))**((g+1)/(2*(g-1)))

    return f