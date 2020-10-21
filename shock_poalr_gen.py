import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

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

if __name__ == "__main__":
    global g, M

    M = 5.0
    g = 1.4
    NPR = 70

    beta = get_beta(M, g, NPR)*180.0/np.pi
    pr = get_shock_jump(M, g, NPR)
    pj = 101325.0/pr
    print(pj)
    tj = 300.0/(1 + 0.5*(g-1)*M**2)
    print(beta)
    # beta_init = np.arcsin(1/M)
    # beta_final = np.pi/2
 
    # beta = np.linspace(beta_init, beta_final, num=1000, endpoint=True)
    # theta = get_theta(M, g, beta)
    # pr = get_oblique_pr(M, g, beta)

    # beta1 = get_beta(M, g, NPR)
    # M1 = get_M(M, g, beta1)
    # thetaadd = get_theta(M, g, beta1)
    # pradd = get_oblique_pr(M, g, beta1)
    
    # beta_init = np.arcsin(1./M)
    # beta_final = np.pi/2
    # beta = np.linspace(beta_init, beta_final, num=1000, endpoint=True)

    # theta1 = get_theta(M1, g, beta)
    # pr1 = get_oblique_pr(M1, g, beta1)

#     fig, ax = plt.subplots(figsize=(6,5))

#     _ = ax.plot(theta*180.0/np.pi, pr, linestyle='-', color='black', linewidth=0.75)
   
#     _ = ax.plot(-theta*180.0/np.pi, pr, linestyle='-', color='black', linewidth=0.75)

# # ax.text(4.4, 100, 'Over Expanded Regime',
# #          {'style':'italic', 'color': 'k', 'fontsize': 8, 'ha': 'center', 'va': 'center'})

# # ax.text(4.5, 70, 'RR',
# #          {'style':'italic', 'color': 'k', 'fontsize': 8, 'ha': 'center', 'va': 'center'})

# # ax.text(4.5, 20, 'MR',
# #          {'style':'italic', 'color': 'k', 'fontsize': 8, 'ha': 'center', 'va': 'center'})

# # ax.text(4.75, 45, 'DD',
# #          {'style':'italic', 'color': 'k', 'fontsize': 8, 'ha': 'center', 'va': 'center'})

# # ax.text(3, 250, 'Under Expanded Regime',
# #          {'style':'italic', 'color': 'k', 'fontsize': 8, 'ha': 'center', 'va': 'center'})


#     ax.set_xlabel(r'$\theta$')
#     ax.set_ylabel(r'\bf{$\frac{p_j}{p_i}$}', rotation=0, fontsize=10)
#     # ax.legend(('Incident Polar', 'Reflected Polar - RR solution', 'Reflected Polar - MR solution', 'Reflected Polar - Detachment solution', 'Reflected Polar - von-Neumann solution'), fontsize='x-small')

# # Move left y-axis to centre, passing through (0,0)
#     ax.spines['left'].set_position('center')

# # Eliminate upper and right axes
#     ax.spines['right'].set_color('none')
#     ax.spines['top'].set_color('none')

# # Show ticks in the left and lower axes only
#     ax.xaxis.set_ticks_position('bottom')
#     ax.yaxis.set_ticks_position('left')

#     ax.set_ylim(1, 18)
#     ax.set_xlim(-35, 35)
#     ax.set_aspect(3.025)
# # ax.grid('on', linestyle='-.')


#     fig.tight_layout()

# # fig.show()
#     fig.savefig('shock_polar.pdf')  
    pass