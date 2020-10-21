# import numpy as np
# from scipy.optimize import fsolve

# def basicvalues(uval, vval):
#     g = 1.2 # Gamma value for the air
#     R = 320.0
#     T = 3000.0
#     anot = g*R*T
#     Vval = np.sqrt(uval*uval + vval*vval)
#     thetaval = np.arctan(vval/uval)
#     aval = np.sqrt(anot - (g-1)*0.5*Vval*Vval)
#     alphaval = np.arcsin(aval/Vval)

#     return Vval, thetaval, aval, alphaval

# def charavalues(thetaval, alphaval, uval, vval, aval, yval, plus):
#     if(plus):
#         lambdaval = np.tan(thetaval + alphaval)
#         qval = uval*uval + aval*aval
#     else:
#         lambdaval = np.tan(thetaval - alphaval)
#         qval = uval*uval - aval*aval

#     rval = 2.0*uval*vval - qval*lambdaval
#     sval = aval*aval*vval/yval

#     return lambdaval, qval, rval, sval

# def interior_pt(upt2, upt1, vpt2, vpt1, ypt2, ypt1, xpt2, xpt1):
#     # calculate the co-efficents for the predictor
#     # uplus and uminus are set as point 2 and 1 respectively
#     uplus = upt2
#     uminus = upt1

#     vplus = vpt2
#     vminus = vpt1

#     yplus = ypt2
#     yminus = ypt1

#     xplus = xpt2
#     xminus = xpt1

#     Vplus, thetaplus, aplus, alphaplus = basicvalues(uplus, vplus)
#     Vminus, thetaminus, aminus, alphaminus = basicvalues(uminus, vminus)

#     Lamdaplus, Q_plus, R_Plus, S_plus = charavalues(thetaplus, alphaplus, uplus, vplus, aplus, yplus, True)
#     Lamdaminus, Q_minus, R_minus, S_minus = charavalues(thetaminus, alphaminus, uminus, vminus, aminus, yminus, False)

#     a11 = 1
#     a12 = -Lamdaplus
#     a21 = 1
#     a22 = -Lamdaminus

#     b11 = yplus - Lamdaplus*xplus
#     b22 = yminus - Lamdaminus*xminus

#     A = np.array([[a11, a12],[a21, a22]])
#     B = np.array([b11, b22])

#     xnew = np.linalg.solve(A, B)
    
#     return xnew


# """
# Prandtl-Meyer expansion relation to calculate the flow PM function for given mach number
# """
# def pm(*data):
#     M, g = data

#     a = np.sqrt((g+1)/(g-1))
#     b = np.arctan(np.sqrt((M**2 -1)/a**2))
#     c = np.arctan(np.sqrt(M**2-1))

#     return a*b - c

# """
# Function to calcualate the mach angle mu
# """
# def ma(*data):
#     M= data
#     return np.arcsin(1.0/M)

# """
# A function to calculate the interior point values using newton rhapson method
# """
# def interior(I0, *data):
#     x1, y1, x2, y2, slope1, mu1, slope2, mu2, Mcost, Msint = data

#     x3 = I0[0]
#     y3 = I0[1]
#     # Mcost = I0[2]
#     # Msint = I0[3]

#     theta3 = np.arctan2(Msint, Mcost)
#     mu3 = np.arcsin(np.cos(theta3)/Mcost)
#     lamda3plus = np.tan(theta3 + mu3)
#     lamda3minus = np.tan(theta3 - mu3)

#     lamdaplus = (np.tan(slope2 + mu2) + lamda3plus)*0.5
#     lamdaminus = (np.tan(slope1 - mu1) + lamda3minus)*0.5
 
#     f = np.zeros((2),dtype=float)

#     f[0] = y3 - y2 - lamdaplus*(x3 - x2)
#     f[1] = y3 - y1 - lamdaminus*(x3 - x1)

#     return f

# """
# New interior point calculation based of the Gerard Boberg implementation
# """
# def new_interior_pt(*data1):

#     g = 1.2    # No units
#     R = 320.0  # J/(Kg-K)
#     T = 3000.0 # IN kelvins

#     x1, y1, slope1, Mach1, x2, y2, slope2, Mach2 = data1
    
#     mu1 = np.arcsin(1./Mach1)
#     mu2 = np.arcsin(1./Mach2)

#     #order of magnitude of dy
#     dy = abs(y2 - y1)
#     theta12 = 0.5*(slope1 + slope2)
#     M12 = 0.5*(Mach1 + Mach2)

#     data = (x1, y1, x2, y2, slope1, mu1, slope2, mu2, M12*np.cos(theta12), M12*np.sin(theta12))
#     I0 = ([x1+dy, 0.5*(y1+y2)])
#     xfunc = fsolve(interior, I0, args=data, xtol=1.e-9)
    
#     return xfunc
#    # Calculate the Prandtl-Meyer angle 


# if __name__ == "__main__":
#     # Program to carry out the method of chatacteristics
#     g = 1.2
#     R = 320.0
#     T = 3000.0
#     # Program will be now used to test the interior point calculation

#     x_p1 = 0.131460
#     y_p1 = 0.040118

#     x_p2 = 0.135683
#     y_p2 = 0.037123

#     u_p1 = 2473.4
#     v_p1 = 812.8

#     u_p2 = 2502.8
#     v_p2 = 737.6

#     x_new = interior_pt(u_p2, u_p1, v_p2, v_p1, y_p2, y_p1, x_p2, x_p1)
#     print(x_new)

#     slope1 = np.arctan2(v_p1, u_p1)
#     slope2 = np.arctan2(v_p2, u_p2)

#     anot = g*R*T
#     a1 = np.sqrt(anot - (g-1)*0.5*(u_p1**2 + v_p1**2))
#     a2 = np.sqrt(anot - (g-1)*0.5*(u_p2**2 + v_p2**2))

#     Mach1 = np.sqrt(u_p1**2 + v_p1**2)/(a1)
#     Mach2 = np.sqrt(u_p2**2 + v_p2**2)/(a2)


#     data1 = (x_p1, y_p1, slope1, Mach1, x_p2, y_p2, slope2, Mach2)
#     Xnew = new_interior_pt(*data1)
#     print(Xnew)
#     pass

import numpy as np
from scipy.optimize import fsolve

def interior_eqns(I0, *data):
    xminus, yminus, uminus, vminus, xplus, yplus, uplus, vplus = data

    x3 = I0[0]
    y3 = I0[1]
    u3 = I0[2]
    v3 = I0[3]

    f = np.zeros((4),dtype=float)

    Vplus = np.sqrt(uplus**2 + vplus**2)
    Vminus = np.sqrt(uminus**2 + vminus**2)

    # V3 = 0.5*(Vplus + Vminus)

    anot = g*R*T
    aplus = np.sqrt(anot - (g-1.0)*0.5*(Vplus**2))
    aminus = np.sqrt(anot - (g-1.0)*0.5*(Vminus**2))

    # a3 = 0.5*(aplus + aminus)

    muplus = np.arcsin(aplus/Vplus)
    muminus = np.arcsin(aminus/Vminus)
    # mu3 = np.arcsin(a3/V3)

    thetaplus = np.arctan2(vplus, uplus)
    thetaminus = np.arctan2(vminus, uminus)
    # theta3 = np.arctan((vplus + vminus)/(uplus + uminus))

    lamdaplus = np.tan(thetaplus + muplus)
    lamdaminus = np.tan(thetaminus - muminus)

    qplus = uplus**2 - aplus**2
    rplus = 2.0*uplus*vplus - lamdaplus*qplus

    qminus = uminus**2 - aminus**2
    rminus = 2.0*uminus*vminus - lamdaminus*qminus

    tplus = qplus*uplus + rplus*vplus
    tminus = qminus*uminus + rminus*vminus

    # print (lamdaplus)

    f[0] = y3 - yplus - lamdaplus*(x3 - xplus)
    f[1] = y3 - yminus - lamdaminus*(x3 - xminus)
    f[2] = qplus*u3 + rplus*v3 - tplus
    f[2] = qminus*u3 + rminus*v3 - tminus

    return f

def interior(*data):
    
    xp1, yp1, up1, vp1, xp2, yp2, up2, vp2 = data

    #initial conditions for solving the point 
    xinit = 0.5*(xp1 + xp2)
    yinit = 0.5*(yp2 + yp1)
    uinit = ((up1 + up2))*0.5
    vinit = ((vp1 + vp2))*0.5

    I0 = ([xinit, yinit, uinit, vinit])
    retVal = fsolve(interior_eqns, I0, args=data)
    return retVal

if __name__ == "__main__":
    global g, R, T, anot

    g = 1.2
    R = 320.0
    T = 3000.0 

    #Specify the initial points
    xp1 = 0.131460
    yp1 = 0.040118

    xp2 = 0.135683
    yp2 = 0.037123

    up1 = 2473.4
    vp1 = 812.8
    
    up2 = 2502.8
    vp2 = 737.6

    errorVal = 1.0
    counter = 0

    while errorVal > 1.0e-9:
        if counter == 0:
            # mu1 = np.arcsin(1./M1)
            # mu2 = np.arcsin(1./M2)

            data = (xp1, yp1, up1, vp1, xp2, yp2, up2, vp2)
            intVal = interior(*data)
            xnew = intVal[0]
            ynew = intVal[1]
            unew = intVal[2]
            vnew = intVal[3]
            errorVal = 1.0e-2
            counter =+1 
            pass
        else:
            xplus = (xp2 + xnew)*0.5
            yplus = (yp2 + ynew)*0.5
            uplus = (up2 + unew)*0.5
            vplus = (vp2 + vnew)*0.5

            xminus = (xp1 + xnew)*0.5
            yminus = (yp1 + ynew)*0.5
            uminus = (up1 + unew)*0.5
            vminus = (vp1 + vnew)*0.5

            data = (xminus, yminus, uminus, vminus, xplus, yplus, uplus, vplus)
            intVal = interior(*data)
            xerr = abs(xnew - intVal[0])
            yerr = abs(ynew - intVal[1])
            uerr = abs(unew - intVal[2])
            verr = abs(vnew - intVal[3])

            xnew = intVal[0]
            ynew = intVal[1]
            unew = intVal[2]
            vnew = intVal[3]

            errorVal = max(xerr, yerr, uerr, verr)
            print (counter, errorVal, xerr, yerr, uerr, verr)
            counter = counter + 1
            
            # errorVal = errorVal/10
            pass
        pass

    print (intVal, counter)

    pass