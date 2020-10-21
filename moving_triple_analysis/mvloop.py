import numpy as np
from numpy import tan, sin, cos, arcsin, arctan, arccos, pi
from scipy.optimize import fsolve

def LinBendor(I0, *data):
    H, mu_MD, beta1, theta1, beta2, mu_M2, theta3, g, Mavg = data

    # print (mu_MD)
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
    
    fval = np.empty(10)
    
    
    fval[0] = xt*tan(beta1) - (H - yt)#xt*tan(beta1) - (H - yt)
    fval[1] = xd*tan(theta1) - (H - yd)
    fval[2] = (xd - xt)*tan(beta2 - theta1) - (yd - yt)
    fval[3] = (xf - xd)*tan(mu_M2 + theta3) + (yf - yd) # check this eqn, i think the present version that only length matters works very well here
    fval[4] = (xe - xd)*tan(mu_MD) - (yd - ye)
    fval[5] = (xf - xt)*tan(-theta3) - (yf - yt)
    fval[6] = (xf - xe)*(tan(-theta3)/(2+ tan(theta3)**2)) - (yf - Hs)
    fval[7] = yt - (ye/Mavg)*((2/(g+1))*(1 + (g-1)*0.5*Mavg**2))**((g+1)/(2*(g-1)))
    fval[8] = yt - Hm
    fval[9] = ye - Hs

    return fval


def cot(x):
    return 1.0/tan(x)

def F(M, beta, g):
    term1_nume = M*M + (2.0/(g-1))
    term1_deno = (2*g/(g-1))*M*M*sin(beta)*sin(beta) - 1.0

    term2_nume = M*M*cos(beta)*cos(beta)
    term2_deno = (g-1)*0.5*M*M*sin(beta)*sin(beta) + 1

    term1 = term1_nume/term1_deno
    term2 = term2_nume/term2_deno
    return np.sqrt(term1 + term2)

def G(M, beta, g):
    return (2*g/(g+1))*M**2*sin(beta)**2 - (g-1)/(g+1)

def H(M, beta, g):
    return ((g+1)*M*M*sin(beta)*sin(beta))/((g-1)*M*M*sin(beta)*sin(beta) + 2)

def A(M, beta, g):
    return np.sqrt(G(M, beta, g)/H(M, beta, g))

def S(M, beta, g):
    nume = (M*M*sin(beta)*sin(beta) - 1)
    deno = M*M*(g + cos(2*beta)) + 2

    return 2*cot(beta)*(nume/deno)

def PSI(M, a, Utp, phi, psi):
    vx = M*a*cos(phi) - Utp*cos(psi)
    vy = M*a*sin(phi) - Utp*sin(psi)
    return np.sqrt(vx*vx + vy*vy)/a

def PHIu(M, a, Utp, phi, psi, theta):
    vx = M*a*cos(phi+theta) - Utp*cos(psi)
    vy = M*a*sin(phi+theta) - Utp*sin(psi)
    return np.sqrt(vx*vx + vy*vy)/a

def PHId(M, a, Utp, phi, psi, theta):
    vx = M*a*cos(phi-theta) - Utp*cos(psi)
    vy = M*a*sin(phi-theta) - Utp*sin(psi)
    return np.sqrt(vx*vx + vy*vy)/a

def Bu(M, a, Utp, phi, psi, beta):
    vx = M*a*cos(phi) - Utp*cos(psi)
    vy = M*a*sin(phi) - Utp*sin(psi)

    return beta + (phi - arctan(vy/vx))

def Bd(M, a, Utp, phi, psi, beta):
    vx = M*a*cos(phi) - Utp*cos(psi)
    vy = M*a*sin(phi) - Utp*sin(psi)

    return beta - (phi - arctan(vy/vx))

def Qu(Mu, au, Md, ad, Utp, phi, psi, theta):
    vx = Mu*au*cos(phi) - Utp*cos(psi)
    vy = Mu*au*sin(phi) - Utp*sin(psi)

    vx1 = Md*ad*cos(phi+theta) - Utp*cos(psi)
    vy1 = Md*ad*sin(phi+theta) - Utp*sin(psi)
    
    return arctan(vy1/vx1) - arctan(vy/vx)

def Qd(Mu, au, Md, ad, Utp, phi, psi, theta):
    vx = Mu*au*cos(phi) - Utp*cos(psi)
    vy = Mu*au*sin(phi) - Utp*sin(psi)

    vx1 = Md*ad*cos(phi-theta) - Utp*cos(psi)
    vy1 = Md*ad*sin(phi-theta) - Utp*sin(psi)
    
    return -(arctan(vy1/vx1) - arctan(vy/vx))

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

def tp_threeshock(I0, *data):
    M0, p0, r0, a0, beta, U_tp = data

    # Values after the incident shock wave
    M1 = F(M0, beta, g)
    beta1 = beta
    theta1 = arctan(S(M0, beta, g))
    pr1 = G(M0, beta, g)
    ar1 = A(M0, beta, g)
    rr1 = H(M0, beta, g)


    M2star = I0[0]
    pr2star = I0[1]
    beta2star = I0[2]
    rr2star = I0[3]
    ar2star = I0[4]
    theta2star = I0[5]

    M3star = I0[6]
    pr3star = I0[7]
    beta3star = I0[8]
    rr3star = I0[9]
    ar3star = I0[10]
    theta3star = I0[11]

    M1star = I0[12]
    M0star = I0[13]
    
    theta2 = I0[14]
    theta3 = I0[15]
    
    beta2 = I0[16]
    beta3 = I0[17]

    M2 = I0[18]
    M3 = I0[19]

    f = np.empty([20])

    f[0] = M2star - F(M1star, beta2star, g)
    f[1] = pr2star - G(M1star, beta2star, g)
    f[2] = rr2star - H(M1star, beta2star, g)
    f[3] = ar2star - A(M1star, beta2star, g)
    f[4] = tan(theta2star) - S(M1star, beta2star, g)

    f[5] = M3star - F(M0star, beta3star, g)
    f[6] = pr3star - G(M0star, beta3star, g)
    f[7] = rr3star - H(M0star, beta3star, g)
    f[8] = ar3star - A(M0star, beta3star, g)
    f[9] = tan(theta3star) - S(M0star, beta3star, g)

    f[10] = M1star - PSI(M1, ar1*a0, U_tp, -theta1, pi-beta1)
    f[11] = M2star - PHIu(M2, A(M1, beta2, g)*ar1*a0, U_tp, -theta1, pi-beta1, theta2)
    f[12] = beta2star - Bu(M1, ar1*a0, U_tp, -theta1, pi - beta1, beta2)
    f[13] = theta2star - Qu(M1, ar1*a0, M2, A(M1, beta2, g)*ar1*a0, U_tp, -theta1, pi - beta1, theta2)

    f[14] = M0star - PSI(M0, a0, U_tp, 0, pi - beta1)
    f[15] = M3star - PHId(M3, A(M0, beta3, g)*a0, U_tp, 0, pi - beta1, theta3)
    f[16] = beta3star - Bd(M0, a0, U_tp, 0, pi-beta1, beta3)
    f[17] = theta3star - Qd(M0, a0, M3, A(M0, beta3, g)*a0, U_tp, 0, pi-beta1, theta3)

    f[18] = theta3 - theta1 + theta2
    f[19] = G(M1, beta2, g)*pr1 - G(M0, beta3, g)

    return f

def nu(*data):
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

def expansionregion(*data):
    M2, g, theta3 = data

    # Calculation of first expansion wave angle
    mu_M2 = arcsin(1/M2)

    # Calculation of Prandtl-Meyer function for the desired Machnumber
    nu_M2 = nu(*(M2, g))

    # Calculation of PM function at the point where flow is parallel to horizontal(D)
    nu_MD = theta3 + nu_M2

    # Calculation of Mach number from the PM fuction
    dataexp = (nu_MD, g) 
    MD = fsolve(PMfunction, M2, args=dataexp)[0]

    mu_MD = arcsin(1/MD)

    return nu_M2, mu_M2, MD, nu_MD, mu_MD

def subsonicpocket(*data):
    M, g, M3, r3, a3, a0, r0, p0 = data

    #calculation of the normal shock relations 
    Mn = np.sqrt((1 + (g-1)*M**2 + ((0.25*(g+1)**2)-g)*M**4) /
              ((g*M**2 - 0.5*(g-1))*(0.5*(g-1)*M**2 + 1)))

    pn = (2*g*M**2 -g +1)/(g+1)
    rn = ((g+1)*M**2)/(2+(g-1)*M**2)
    tn = pn/rn

    rn = rn*r0
    an = np.sqrt(tn)*a0

    r3 = r3*r0
    a3 = a3*a0

    un = Mn*an
    u3 = M3*a3

    Mavg = 2*(r3*u3*np.cos(theta3) + rn*un)/((rn+r3)*(an+a3))

    return Mavg

def machstem(*data):
    
    #setting up the initial conditions
    I0 = np.ones(10)#[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]
    zmachstem = fsolve(LinBendor, I0, args=data)
    # zmachstemmodified = fsolve(mLinBendor, I0[0:8], args=data)
    # zmachstemmouton = fsolve(Mouton, I0[0:6], args=data)
    return zmachstem#, zmachstemmodified, zmachstemmouton

if __name__ == "__main__":
    g = 1.4
    M = 5.0
    npr = 50.7940185
    R = 8314.4621/28.9647
    Hv = 1.0

    u_tp = np.linspace(0, 81.0, 100, endpoint=True)
    result = np.zeros((100, 2))
    counter = 0
## Nozzle exit conditions
    beta = get_beta(M, g, npr)
    pr_exit = get_shock_jump(M, g, npr)
    te_exit = 300.0/(1 + (g-1)*0.5*M**2)
    pe_exit = 101325.0/pr_exit
    re_exit = pe_exit/(R*te_exit)
    ae_exit = np.sqrt(g*pe_exit/re_exit)
    
    M0 = M
    p0 = pe_exit
    t0 = te_exit
    r0 = re_exit
    a0 = ae_exit

    # Values after the incident shock wave
    M1 = F(M, beta, g)
    beta1 = beta
    theta1 = arctan(S(M, beta, g))
    pr1 = G(M, beta, g)
    ar1 = A(M, beta, g)
    rr1 = H(M, beta, g)

    # moving to the reflected shock which is considered in the moving frame
    # step 1 convert the flow values behind the incident shock to moving frame

    for U_tp in u_tp:
        data = (M0, p0, r0, a0, beta, U_tp)
        I0 = [2, 2, np.pi/4, 2, 2, np.pi/8, 0.9, 2, np.pi/2, 2, 2, np.pi/8, 2, 4, np.pi/8, np.pi/20, np.pi/4, np.pi/2, 2, 0.8]
        ztpmoving = fsolve(tp_threeshock, I0, args=data)

        theta2 = ztpmoving[14]
        theta3 = ztpmoving[15]
    
        beta2 = ztpmoving[16]
        beta3 = ztpmoving[17]

        M2 = ztpmoving[18]
        M3 = ztpmoving[19]

        p30 = ztpmoving[7]
        r30 = ztpmoving[9]
        a30 = ztpmoving[10]

    # Solving the expansion fan region
        nu_M2, mu_M2, MD, nu_MD, mu_MD = expansionregion(*(M2, g, theta3))

    # subsonic pocket calculation
        Mbar = subsonicpocket(*(M0, g, M3, r30, a30, a0, r0, p0))


    # calculation of mach stem height
        zLinBendor= machstem(*(Hv, mu_MD, beta1, theta1, beta2, mu_M2, theta3, g, Mbar))
        result[counter, 0] = zLinBendor[1]
        result[counter, 1] = U_tp/a0
        counter = counter + 1 
        print(counter, zLinBendor[1], U_tp/a0)

    np.savetxt('utp_50_794.csv', result, fmt='%.18e', delimiter=',', newline='\n', header='Hm, Utp', footer='', comments='# ', encoding=None)

    pass