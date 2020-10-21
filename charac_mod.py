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

def lin_interior_eqns(*data):
    xminus, yminus, uminus, vminus, xplus, yplus, uplus, vplus = data

    Vplus = np.sqrt(uplus**2 + vplus**2)
    Vminus = np.sqrt(uminus**2 + vminus**2)

    anot = g*R*T
    aplus = np.sqrt(anot - (g-1.0)*0.5*(Vplus**2))
    aminus = np.sqrt(anot - (g-1.0)*0.5*(Vminus**2))

    muplus = np.arcsin(aplus/Vplus)
    muminus = np.arcsin(aminus/Vminus)

    thetaplus = np.arctan2(vplus, uplus)
    thetaminus = np.arctan2(vminus, uminus)

    lamdaplus = np.tan(thetaplus + muplus)
    lamdaminus = np.tan(thetaminus - muminus)

    qplus = uplus**2 - aplus**2
    rplus = 2.0*uplus*vplus - lamdaplus*qplus

    qminus = uminus**2 - aminus**2
    rminus = 2.0*uminus*vminus - lamdaminus*qminus

    tplus = qplus*uplus + rplus*vplus
    tminus = qminus*uminus + rminus*vminus
    
    A = np.array([[-lamdaplus, 1.0, 0.0, 0.0], [-lamdaminus, 1.0, 0.0, 0.0], [0.0, 0.0, qplus, rplus], [0.0, 0.0, qminus, rminus]])
    b = np.array([yplus - lamdaplus*xplus, yminus - lamdaminus*xminus, tplus, tminus])

    x = np.linalg.solve(A, b)

    return x

def lin_axis_eqns(*data):
    xminus, yminus, uminus, vminus  = data

    Vminus = np.sqrt(uminus**2 + vminus**2)

    anot = g*R*T
    aminus = np.sqrt(anot - (g-1.0)*0.5*(Vminus**2))

    muminus = np.arcsin(aminus/Vminus)

    thetaminus = np.arctan2(vminus, uminus)

    lamdaminus = np.tan(thetaminus - muminus)

    qminus = uminus**2 - aminus**2
    rminus = 2.0*uminus*vminus - lamdaminus*qminus

    tminus = qminus*uminus + rminus*vminus
    
    x = np.empty(2)
    x[0] = ((yminus - lamdaminus*xminus)/lamdaminus)*-1.0
    x[1] = (tminus)/qminus

    return x

def lin_freepr_eqns(*data):
    xminus, yminus, uminus, vminus, xplus, yplus, uplus, vplus, p4, t0, p0 = data

    Vfinal = np.sqrt(((2*g*R*t0)/(g-1.0))*(1.0 - (p4/p0)**((g-1)/g)))
    Vplus = np.sqrt(uplus**2 + vplus**2)
    # Vminus = np.sqrt(uminus**2 + vminus**2)

    anot = g*R*T
    aplus = np.sqrt(anot - (g-1.0)*0.5*(Vplus**2))
    # aminus = np.sqrt(anot - (g-1.0)*0.5*(Vminus**2))

    muplus = np.arcsin(aplus/Vplus)
    # muminus = np.arcsin(aminus/Vminus)

    thetaplus = np.arctan2(vplus, uplus)
    thetaminus = np.arctan2(vminus, uminus)

    lamdaplus = np.tan(thetaplus + muplus)
    lamdaminus = np.tan(thetaminus) # since minus characteristic is from the streamline

    qplus = uplus**2 - aplus**2
    rplus = 2.0*uplus*vplus - lamdaplus*qplus

    # qminus = uminus**2 - aminus**2
    # rminus = 2.0*uminus*vminus - lamdaminus*qminus

    tplus = qplus*uplus + rplus*vplus
    # tminus = qminus*uminus + rminus*vminus
    
    A = np.array([[-lamdaplus, 1.0], [-lamdaminus, 1.0]])
    b = np.array([yplus - lamdaplus*xplus, yminus - lamdaminus*xminus])

    x = np.linalg.solve(A, b)

    u4 = (qplus*tplus - rplus*np.sqrt(Vfinal**2*(qplus**2 + rplus**2) - tplus**2))/(qplus**2 + rplus**2)
    v4 = np.sqrt(Vfinal**2 - u4**2)

    new = np.append(x, [u4, v4])
    
    return new

def interior_pt(*data):
    xminus, yminus, uminus, vminus, xplus, yplus, uplus, vplus = data

    counter = 0
    errVal = 1.0
    while errVal > 1.0e-9:
        
        # Calculate the initial condition for the newton rhapson method
        # xinit = 0.5*(xplus + xminus)
        # yinit = 0.5*(yplus + yminus)
        # uinit = 0.5*(uplus + uminus)
        # vinit = 0.5*(vplus + vminus)

        data = (xminus, yminus, uminus, vminus, xplus, yplus, uplus, vplus)
        # I0 = ([xinit, yinit, uinit, vinit])
        # newVal = fsolve(interior_eqns, I0, args=data)
        newVal = lin_interior_eqns(*data)

        xnew = newVal[0]
        ynew = newVal[1]
        unew = newVal[2]
        vnew = newVal[3]

        if counter == 0:
            errVal = 1.0
            #saving the data into a new variable to calculate the error
            oldVal = newVal
        else:
            err = abs(oldVal - newVal)
            oldVal = newVal
            errVal = max(err[0], err[1], err[2], err[3])
        
        xplus = 0.5*(xplus + xnew)
        yplus = 0.5*(yplus + ynew)
        uplus = 0.5*(uplus + unew)
        vplus = 0.5*(vplus + vnew)

        xminus = 0.5*(xminus + xnew)
        yminus = 0.5*(yminus + ynew)
        uminus = 0.5*(uminus + unew)
        vminus = 0.5*(vminus + vnew)

        counter = counter + 1 

        if counter > 25:
            break

    return newVal

def free_pressure_pt(*data):
    xminus, yminus, uminus, vminus, xplus, yplus, uplus, vplus, p4, t0, p0 = data
    
    counter = 0
    errVal = 1.0
    while errVal > 1.0e-9:
        
        # Calculate the initial condition for the newton rhapson method
        # xinit = 0.5*(xplus + xminus)
        # yinit = 0.5*(yplus + yminus)
        # uinit = 0.5*(uplus + uminus)
        # vinit = 0.5*(vplus + vminus)

        data = (xminus, yminus, uminus, vminus, xplus, yplus, uplus, vplus, p4, t0, p0)
        # I0 = ([xinit, yinit, uinit, vinit])
        # newVal = fsolve(interior_eqns, I0, args=data)
        newVal = lin_freepr_eqns(*data)
        
        xnew = newVal[0]
        ynew = newVal[1]
        unew = newVal[2]
        vnew = newVal[3]

        if counter == 0:
            errVal = 1.0
            #saving the data into a new variable to calculate the error
            oldVal = newVal
        else:
            err = abs(oldVal - newVal)
            oldVal = newVal
            errVal = max(err[0], err[1], err[2], err[3])
        
        xplus = 0.5*(xplus + xnew)
        yplus = 0.5*(yplus + ynew)
        uplus = 0.5*(uplus + unew)
        vplus = 0.5*(vplus + vnew)

        xminus = 0.5*(xminus + xnew)
        yminus = 0.5*(yminus + ynew)
        uminus = 0.5*(uminus + unew)
        vminus = 0.5*(vminus + vnew)

        counter = counter + 1

        if counter > 25:
            break

    return newVal

def axis_pt(*data):
    xminus, yminus, uminus, vminus = data

    counter = 0
    errVal = 1.0
    while errVal > 1.0e-9:
    
        data = (xminus, yminus, uminus, vminus)
        # I0 = ([xinit, yinit, uinit, vinit])
        # # newVal = fsolve(interior_eqns, I0, args=data)
        newVal = lin_axis_eqns(*data)

        xnew = newVal[0]
        ynew = 0.0
        unew = newVal[1]
        vnew = 0.0

        if counter == 0:
            errVal = 1.0
            #saving the data into a new variable to calculate the error
            oldVal = newVal
        else:
            err = abs(oldVal - newVal)
            oldVal = newVal
            errVal = max(err[0], err[1])

        xminus = 0.5*(xminus + xnew)
        yminus = 0.5*(yminus + ynew)
        uminus = 0.5*(uminus + unew)
        vminus = 0.5*(vminus + vnew)

        counter = counter + 1 

        if counter > 25:
            break

    return newVal


if __name__ == "__main__":
    global g, R, T

    g = 1.2
    R = 320.0
    T = 3000.0

    # Specify the initial points
    xp1 = 0.131460
    yp1 = 0.040118

    xp2 = 0.135683
    yp2 = 0.037123

    up1 = 2473.4
    vp1 = 812.8
    
    up2 = 2502.8
    vp2 = 737.6

    data = (xp1, yp1, up1, vp1, xp2, yp2, up2, vp2)
    newValue = interior_pt(*data)
    print (newValue)

    xp1 = 79.625e-3
    yp1 = 1.290e-3

    up1 = 2306.1
    vp1 = 35.7

    data= (xp1, yp1, up1, vp1)
    newValue = axis_pt(*data)
    print(newValue)

    xp2 = 0.34226
    yp2 = 0.12312

    up2 = 2455.3
    vp2 = 619.7

    xp1 = 0.32979
    yp1 = 0.12351

    up1 = 2422.5
    vp1 = 663.5

    p4 = 60000.0
    p0 = 7.0e6
    t0 = 3000.0
    data = (xp1, yp1, up1, vp1, xp2, yp2, up2, vp2, p4, t0, p0)
    newValue = free_pressure_pt(*data)
    print(newValue)
    pass
