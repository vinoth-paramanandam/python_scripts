import numpy as np

if __name__ == "__main__":

    # gamma = 1.4
    # M = 5.0

    # betainit = 30.0
    # betafinal = 40.0

    # beta = np.linspace(betainit, betafinal, num=21, endpoint=True)
    # pratio = (1.0 - gamma + 2.0*gamma*M*M*np.sin(beta*np.pi/180)**2)/(gamma + 1)
    # popj = (1 + (gamma-1)*0.5*M**2)**(gamma/(gamma-1.0))

    # for i in range(21):
    #     print (i, beta[i], pratio[i], popj/pratio[i], popj)

    gamma = 1.4
    M = 5.0

    npr = 41.0

    p0pj = (1 + 0.5*(gamma-1.0)*M*M)**(gamma/(gamma-1.0))
    p0 = 25*101325.0

    print (p0/p0pj, p0/npr)

    pass