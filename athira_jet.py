import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

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
    # print ( Vfinal**2*(qplus**2 + rplus**2) - tplus**2, tplus, uplus, vplus, qplus, rplus )
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
        
        data = (xminus, yminus, uminus, vminus, xplus, yplus, uplus, vplus)
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
    
        data = (xminus, yminus, uminus, vminus, xplus, yplus, uplus, vplus, p4, t0, p0)
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

    newValue =np.zeros(4)
    newValue[0] = newVal[0]
    newValue[1] = 0.0
    newValue[2] = newVal[1]
    newValue[3] = 0.0
    return newValue

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

"""
Prandtl-Meyer expansion relation to calculate the flow PM function for given mach number
"""
def pm(*data):
    M, g = data

    a = np.sqrt((g+1)/(g-1))
    b = np.arctan(np.sqrt((M**2 -1)/a**2))
    c = np.arctan(np.sqrt(M**2-1))

    return a*b - c

if __name__ == "__main__":
    global g, R, T

    g = 1.4
    R = 287.14

    # exit plane conditions at the nozzle
    Me = 1.073
    pj = 560474.171 # 200kpa

    # ambient conditions at the nozzle exit/ back pressure conditions
    pa = 128762.121 # 100kpa
    
    # assumptions of t0 and the other total  and exit jet conditions 
    t0 = 688.4226
    T = t0 # just for the sake of using the old code here
    tj = t0/(1 + (g-1.0)*0.5*Me**2)
    p0 = 1156525.398#pj*(1 + (g-1.0)*0.5*Me**2)**(g/(g-1.0))
    ue = Me*np.sqrt(g*R*tj)
    ve = 0.0
    
    data_frp = np.array((pa, t0, p0))
    # Exit plane height
    xe = 0.0
    ye = 0.012
    
    # Total number of characteristics at the exit of the nozzle lip
    waves = 4
    
    # calculation for the first point corresponding to the first wave
    # data = (xe, ye, ue, ve)
    # newValue = axis_pt(*data)

    # mach number behind the expansion fan at the inlet based on the p_jet and p_ambient
    Mf = np.sqrt((2.0/(g-1))*((p0/pa)**((g-1)/g)  - 1.0))
    delnu = (pm(*(Mf,g)) - pm(*(Me,g)))/waves
    

# creating a array for storing the point values
    points = np.zeros((3*waves, waves+1, 4), dtype=float)

# Initial conditions for the various expansion waves 
    M = Me
    for i in range(waves):
        if i > 0:
            nu_M = pm(*(M, g)) + delnu
            M = fsolve(PMfunction, Me*2, args=(nu_M, g))
            theta = delnu*i
            tjet = t0/(1 + (g-1.0)*0.5*M**2)
            ajet = np.sqrt(g*R*tjet)
            u = M*ajet*np.cos(theta)
            v = M*ajet*np.sin(theta)
            points[i, 0] = np.array((xe, ye, u, v))
        else:
            points[i, 0] = np.array((xe, ye, ue, ve))
            
    # Counter will indicate the wave number
    counter = 0

    # # Calculation of the wave interaction points
    # while(True):
        
    #     if counter < waves:
    #         pts = counter - 1
    #         # print(counter, pts)
    #         cnt = 0
    #         while(pts >= 0):
    #             # iterate over the points 
    #             # this will be the interior points
    #             # print(counter, cnt, pts)
    #             data = np.append(points[counter, cnt], points[counter-1, cnt+1])
    #             intValue = interior_pt(*data)
    #             points[counter, cnt+1] = intValue
                 
    #             print(counter, cnt+1)
    #             cnt += 1
    #             pts = pts - 1
                
    #         if pts < 0:
    #             # do the axis point
    #             # if counter == 0:
    #             #     # points[counter, 0] = np.array((xe, ye, ue, ve))
    #             #     data = points[counter, counter]#(xe, ye,ue, ve)

    #             #     axsValue = axis_pt(*data)
    #             #     points[counter, counter+1] = axsValue
    #             #     # points[waves, 0] = axsValue
    #             # else:
    #             #     points[counter, pts+1]
    #             #     pass
    #             data = points[counter, counter]
    #             axsValue = axis_pt(*data)
    #             points[counter, counter+1] = axsValue
    #             points[waves+counter, 0] = axsValue
    #             pass 
    #     elif counter >= waves and counter < 2*waves:
    #         print(counter)
    #     else:
    #         break
    #     # anew = np.sqrt(g*R*T - (g-1.0)*0.5*(points[counter, counter+1][2]**2))
    #     # print(counter, points[counter, counter+1], points[counter, counter+1][2]/anew)
    #     counter += 1
    # pass

    for i in range(waves):
        for j in range(waves):
            print('Solving for point ', i, j)
            if i == j:
                data = points[i, j]
                axsValue = axis_pt(*data)
                points[i, j+1] = axsValue
                break
            else:
                data = np.append(points[i, j], points[i-1,j+1])
                intValue = interior_pt(*data)
                points[i, j+1] = intValue
                pass

    for i in range(waves, 2*waves):
        for j in range(waves+1):
            if i - waves + j < waves:
                print(i, j , '-----', j + (i - waves), i-waves+1)
                points[i, j] = points[j + (i - waves), i-waves+1]
            else:
                if j == waves:
                    if i == waves:
                        data_val = np.append(points[waves-1, 0], points[i, j-1])
                        data = np.append(data_val, data_frp)
                        frpValue = free_pressure_pt(*data)
                        points[i, j] = frpValue
                        print(i, j, 'axis')
                        print(points[i, j])
                    else:
                        data_val = np.append(points[i-1, j], points[i, j-1])
                        data = np.append(data_val, data_frp)
                        frpValue = free_pressure_pt(*data)
                        points[i, j] = frpValue
                        print(i, j, 'axis')
                        print(points[i,j])
                else:
                    data = np.append(points[i-1, j+1], points[i, j-1])
                    intValue = interior_pt(*data)
                    points[i, j] = intValue
                    
                # break
            # points[i, j] = points[j, j+1]
            # print(i, j, j, i-2)
    # for i in range(waves):
    #     for j in range(waves):
    #         if (i == j):
    #             # points[waves + j, i] = points[i, j+1]
    #             print(i, j, waves+j, i)
    #             break
    #         else:
    #             # points[waves + j, i] = points[i, j]
    #             print(i, j, waves+j, i)
    #             pass
    # for i in range(waves):
    #     for j in range(1, waves+1):
    #         if i == waves-1:
    #             points[waves+i, j] = points[, j+1]
    #             print(i, j)
    #         else:
    #             points[waves + i, j] = points[i, j+1]
    #             break
