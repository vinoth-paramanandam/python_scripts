import numpy as np
from numpy import tan, cos, sin, arcsin, arccos, arctan, pi
from scipy.optimize import fsolve

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

    vx1 = Md*ad*cos(phi-theta) - Utp*cos(psi)
    vy1 = Md*ad*sin(phi-theta) - Utp*sin(psi)
    
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


if __name__ == "__main__":
    g = 1.4
    M = 5.0
    npr = 40.0
    R = 287.14

    U_tp = 0.0
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

    ## Solution to three shock theory
    

    # Values after the incident shock wave
    M1 = F(M, beta, g)
    beta1 = beta
    theta1 = arctan(S(M, beta, g))
    pr1 = G(M, beta, g)
    ar1 = A(M, beta, g)
    rr1 = H(M, beta, g)

    M1star = PSI(M1, ar1*a0, U_tp, -theta1, pi-beta1)
    M0star = PSI(M0, a0, U_tp, 0, pi-beta1)
    
    # Assume a beta2star and beta3star
    beta2star = Bu(M1, ar1*a0, U_tp, -theta1, pi - beta1, beta2)
    beta3star = Bd(M0, a0, U_tp, 0, pi-beta1, beta3)

    # moving to the reflected shock which is considered in the moving frame
    # step 1 convert the flow values behind the incident shock to moving frame

    # data = (M0, p0, r0, a0, beta, U_tp)
    # I0 = [2, 2, np.pi/4, 2, 2, np.pi/8, 1.2, 2, np.pi/2, 2, 2, np.pi/8, 2, 4, np.pi/8, np.pi/8, np.pi/4, np.pi/2, 2, 0.8]
    # ztpmoving = fsolve(tp_threeshock, I0, args=data)
    # print(ztpmoving[18], ztpmoving[19])
    pass