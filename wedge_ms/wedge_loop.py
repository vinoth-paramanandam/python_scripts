import numpy as np
from numpy import tan, sin, cos, arcsin, arctan, arccos, pi, sqrt
from scipy.optimize import fsolve
import pandas as pd

def cot(x):
    return 1.0/tan(x)

# Mach number relation across the oblique shock
def F(M, beta, g):
    term1_nume = M*M + (2.0/(g-1))
    term1_deno = (2*g/(g-1))*M*M*sin(beta)*sin(beta) - 1.0

    term2_nume = M*M*cos(beta)*cos(beta)
    term2_deno = (g-1)*0.5*M*M*sin(beta)*sin(beta) + 1

    term1 = term1_nume/term1_deno
    term2 = term2_nume/term2_deno
    return sqrt(term1 + term2)

# pressure relation across the oblique shock
def G(M, beta, g):
    return (2*g/(g+1))*M**2*sin(beta)**2 - (g-1)/(g+1)

# Density relation across the oblique shock
def H(M, beta, g):
    return ((g+1)*M*M*sin(beta)*sin(beta))/((g-1)*M*M*sin(beta)*sin(beta) + 2)

# Sound relation across the oblique shock
def A(M, beta, g):
    return sqrt(G(M, beta, g)/H(M, beta, g))

# ø-ß-M relation
def S(M, beta, g):
    nume = (M*M*sin(beta)*sin(beta) - 1)
    deno = M*M*(g + cos(2*beta)) + 2

    return 2*cot(beta)*(nume/deno)

def calc_beta(M, theta, g):
    return fsolve(lambda beta: S(M, beta, g) - tan(theta), pi/4)

def calc_theta(M, beta, g):
    return arctan(S(M, beta, g))

def nu(M, g):
    
    a = np.sqrt((g+1)/(g-1))
    b = np.arctan(np.sqrt((M**2 -1)/a**2))
    c = np.arctan(np.sqrt(M**2-1))

    return a*b - c

def threeshock(I0, *data):
    M0, p0, r0, a0, theta, g = data

    # Solve the flow field across the incident shockwave
    beta1 = calc_beta(M0, theta, g)
    theta1 = theta
    M1 = F(M0, beta1, g)
    pr1 = G(M0, beta1, g)
    
    M2 = I0[0]
    beta2 = I0[1]
    theta2 = I0[2]
    ar2 = I0[3]
    rr2 = I0[4]
    pr2 = I0[5]

    M3 = I0[6]
    beta3 = I0[7]
    theta3 = I0[8]
    ar3 = I0[9]
    rr3 = I0[10]
    pr3 = I0[11]

    f = np.empty([12])

    f[0] = M2 - F(M1, beta2, g)
    f[1] = pr2 - G(M1, beta2, g)
    f[2] = rr2 - H(M1, beta2, g)
    f[3] = ar2 - A(M1, beta2, g)
    f[4] = tan(theta2) - S(M1, beta2, g)

    f[5] = M3 - F(M0, beta3, g)
    f[6] = pr3 - G(M0, beta3, g)
    f[7] = rr3 - H(M0, beta3, g)
    f[8] = ar3 - A(M0, beta3, g)
    f[9] = tan(theta3) - S(M0, beta3, g)

    f[10] = theta3 - theta1 + theta2
    f[11] = pr2*pr1 - pr3

    return f

def shock_expansion(I0, *data):
    M1, theta1, p1, M2, p2, theta3, g = data

    Mc = I0[0]
    Mcd = I0[1]
    Md = I0[2]

    alpha = I0[3]
    betac = I0[4]
    thetacd = I0[5]
    
    pc = I0[6]
    pcd = I0[7]
    pd = I0[8]
   
    f = np.empty([9])

    f[0] = nu(Mc, g) - nu(M1, g) - theta1 + alpha
    f[1] = (pc/p1) - ((2 + (g-1)*M1*M1)/(2 + (g-1)*Mc*Mc))**(g/(g-1))
    f[2] = Mcd - F(Mc, betac, g)
    f[3] = pcd - pc*G(Mc, betac, g)
    f[4] = tan(thetacd) - S(Mc, betac, g)
    f[5] = alpha - thetacd
    f[6] = pcd - pd
    f[7] = nu(Md, g) - nu(M2, g) - theta3
    f[8] = (pd/p2) - ((2 + (g-1)*M2*M2)/(2 + (g-1)*Md*Md))**(g/(g-1))

    return f

def subsonic_portion(M, r0, a0, g, M3, r3, a3, theta3):
    rg = r0*H(M, pi/2, g)
    ag = a0*A(M, pi/2, g)

    u3 = M3*a3
    ug = F(M, pi/2, g)*ag

    Mbar = (2.0*(r3*u3*cos(theta3) + rg*ug))/((r3 + rg)*(a3 + ag))

    return Mbar

def delta(d1, d2):
    nume = 2.0*tan(d1) + tan(d2 - d1)
    deno = 2.0 - tan(d1)*tan(d2 - d1)
    return arctan(nume/deno)

def geomentry(I0, *data):

    M1, M2, Mc, Mcd, Md, w, H_, theta1, theta3, beta1, beta2, betac, alpha, Mbar = data
    
    mu_b = arcsin(1.0/M1)
    mu_c = arcsin(1.0/Mc)
    mu_d = arcsin(1.0/Md)
    mu_cd = arcsin(1.0/Mcd)
    mu_f = arcsin(1.0/M2)
    xr = w*cos(theta1)
    yr = H_ - w*sin(theta1)

    xb = I0[0]
    yb = I0[1]
    xc = I0[2]
    yc = I0[3]
    xd = I0[4]
    yd = I0[5]
    xe = I0[6]
    ye = I0[7]
    xf = I0[8]
    yf = I0[9]
    xt = I0[10]
    yt = I0[11]
    Hs = I0[12]
    Hm = I0[13]

    f = np.empty([14])
    
    f[0] = yb - yr + (xb - xr)*tan(mu_b + theta1)
    f[1] = yc - yr + (xc - xr)*tan(mu_c + alpha)
    f[2] = yf - yb + (xf - xb)*tan(mu_f + theta3)
    f[3] = ye - yd + (xe - xd)*tan(mu_d)
    f[4] = yb - yt - (xb - xt)*tan(beta2 - theta1)
    f[5] = ye - Hs
    f[6] = yt - Hm
    f[7] = xt - (H_ - Hm)*cot(beta1)
    f[8] = yf - yt + (xf - xt)*tan(theta3)
    f[9] = yb - yc - (xb - xc)*tan(delta(beta2-theta1, betac-alpha))
    f[10] = yc - yd - (xc - xd)*tan(delta(-mu_cd, -mu_d))
    f[11] = yb - yd - (xb - xd)*tan(delta(theta3, 0))
    f[12] = yf - ye - (xf - xe)*tan(delta(-theta3, 0))
    f[13] = Hm - (Hs/Mbar)*((2/(g+1))*(1 + (g-1)*0.5*Mbar**2))**((g+1)/(2*(g-1)))

    return f

def calc_mach_stem(g, M, theta, R, H_, w):
    
    # Initial state
    M0 = M
    p0 = 101325.0
    t0 = 300.0
    r0 = p0/(R*t0)
    a0 = sqrt(g*R*t0)

    # Behind the incident shock wave
    beta1 = calc_beta(M0, theta, g)[0]
    theta1 = theta
    M1 = F(M0, beta1, g)
    pr1 = G(M0, beta1, g)
    ar1 = A(M0, beta1, g)
    rr1 = H(M0, beta1, g)

    print ('Solved the incident shock parameters')
    data = (M0, p0, r0, a0, theta, g)

    I0 = [1.5, pi/4.0, pi/8.0, 1.2, 2, 3, 0.5, pi/2, pi/4, 2, 3, 5]
    zthreeshock = fsolve(threeshock, I0, args=data)

    print ('Solved the three shock theory close to triple point')
    M2 = zthreeshock[0]
    M3 = zthreeshock[6]

    beta2 = zthreeshock[1]
    beta3 = zthreeshock[7]
    
    theta2 = zthreeshock[2]
    theta3 = zthreeshock[8]

    ar2 = zthreeshock[3]
    ar3 = zthreeshock[9]
    
    rr2 = zthreeshock[4]
    rr3 = zthreeshock[10]

    pr2 = zthreeshock[5]
    pr3 = zthreeshock[11]

    p1 = pr1*p0
    p2 = pr2*p1
    p3 = pr3*p0

    r1 = rr1*r0
    r2 = rr2*r1
    r3 = rr3*r0

    a1 = ar1*a0
    a2 = ar2*a1
    a3 = ar3*a0
    
    data = (M1, theta1, p1, M2, p2, theta3, g)
    # I0 = [1.2*M2, 0.8*M2, M2, pi/8, pi/4, pi/8, 0.7*p1, 1.2*p1, p1]
    I0 = [M2, M2, M2, pi/8, pi/4, pi/8, p1, p1, p1]
    zshockexpansion = fsolve(shock_expansion, I0, args=data)
    
    print ('Solved the Shockwave Expansion fan interaction')

    Mc = zshockexpansion[0]
    Mcd = zshockexpansion[1]
    Md = zshockexpansion[2]

    alpha = zshockexpansion[3]
    betac = zshockexpansion[4]
    thetacd = zshockexpansion[5]
    
    pc = zshockexpansion[6]
    pcd = zshockexpansion[7]
    pd = zshockexpansion[8]

    Mbar = subsonic_portion(M0, r0, a0, g, M3, r3, a3, theta3)
    
    print ('Solved the Subsonic Portion behind the Mach stem')
    
    data = (M1, M2, Mc, Mcd, Md, w, H_, theta1, theta3, beta1, beta2, betac, alpha, Mbar)
    I0 = [w, H_, w, H_, w, H_, w, H_, w, H_, w, H_, w, H_]
    zMs = fsolve(geomentry, I0, args=data)
    
    print ('Solved the Geometry to obtain the Mach stem estimate')

    muB = arcsin(1.0/M1)
    numeHtmax = w*sin(muB + theta1)*sin(beta1 - theta1)
    denoHtmax = sin(muB + theta1 - beta1)

    numeHtmin = w*sin(beta2 - theta1)*sin(beta1 - theta1)
    denoHtmin = sin(beta1 + beta2 - theta1)

    Htmax = zMs[13] + numeHtmax/denoHtmax
    Htmin = zMs[13] + numeHtmin/denoHtmin
    return zMs, Htmax, Htmin

if __name__ == "__main__":
    g = 1.4
    M = 4.96
    # theta = 30.0*pi/180.0
    R = 8314.4621/28.9647
    # H_ = 1.0
    # w = 0.77
    weg_length = pd.read_csv('wedlen_4m_jet.csv')
    file  = open('wedge_LnB.csv', 'a')
    zMs = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

    for i in range(weg_length['M'].shape[0]):
        theta = weg_length['theta'][i]
        H_ = 1.0
        w = weg_length['w_h'][i]
        M = weg_length['M'][i]
        
        # if i == 0:
        #     Ival = [w, H_, w, H_, w, H_, w, H_, w, H_, w, H_, w, H_]
        # else:
        #     Ival = zMs
        Ival = [w, H_, w, H_, w, H_, w, H_, w, H_, w, H_, w, H_]
        zMs, Htmax, Htmin = calc_mach_stem(g, M, theta, R, H_, w)
        Hm = zMs[13]

        print(M, theta*180.0/pi, w, Hm)
        thetaval = theta*180/pi
        valLnBweged = np.append(np.append(thetaval, Hm), np.append(Htmax, Htmin))
    
        np.savetxt(file, [valLnBweged], delimiter=',', fmt='%f')
    pass
