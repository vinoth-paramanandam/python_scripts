import numpy as np
from numpy import sin, cos, tan, arcsin, arctan, sqrt

def get_shock_strength(M, g, npr):
    pope = (1.0 + 0.5*(g-1)*M**2)**(g/(g-1))
    return pope/npr
def get_temperature(M, g, npr):
    dummy = np.ones(npr.shape[0])
    return (1.0 + 0.5*(g-1)*M**2)*dummy

if __name__ == "__main__":
    
    # add the constants here for the main file
    g = 1.4
    pamb = 101325.0
    tstag = 300.0
    Rair = 8314.0/28.9647
    ramb = pamb/(Rair*tstag)

    # choose the mach number for the flow
    M = 2.44

    # choose the range of the nozzle pressure ratio
    npr = np.arange(5.0, 6.1, 0.1)

    # find the shock jump associated with the npr
    shock_jump = get_shock_strength(M, g, npr)

    pexit = pamb/shock_jump
    texit = 300.0/get_temperature(M, g, npr)
    
    # creating a new numpy array to store the relevant data and save it in csv file
    result = np.empty((npr.shape[0], 7))
    result[:, 0] = npr
    result[:, 1] = pexit
    result[:, 2] = texit
    result[:, 3] = pexit/(Rair*texit)
    result[:, 4] = g*Rair*texit
    result[:, 5] = pamb/(g*pexit)
    result[:, 6] = ramb/(pexit/(Rair*texit))
    # saving in csv file
    np.savetxt('nozzle_exit_2_44.csv', result,  delimiter=',', fmt='%15.5f')
    print('Done')
    pass