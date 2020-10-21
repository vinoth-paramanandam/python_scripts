import numpy as np
import matplotlib.pyplot as plt

reg = np.genfromtxt('regimes.csv', delimiter =',')

m = reg[:, 0]
sc = reg[:, 1]
so = reg[:, 2]
de = reg[:, 3]
vn = reg[:, 4]
tc = reg[:, 5]

plt.rc('font', family='serif', serif='Times')
plt.rc('text', usetex=False)
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('axes', labelsize=12)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

ax.set_xlabel('Mach Number')
ax.set_ylabel('Nozzle Pressure Ratio')

# ax.text(3, 50, 'Perfectly expanded Jet',
#          {'style':'italic', 'color': 'k', 'fontsize': 10, 'ha': 'center', 'va': 'center',
#           'bbox': dict(boxstyle="round", fc="w", ec="k", pad=0.5)})

ax.text(3.75, 40, 'OE',
         {'style':'italic', 'color': 'k', 'fontsize': 10, 'ha': 'center', 'va': 'center',
          'bbox': dict(boxstyle="round", fc="w", ec="k", pad=0.25)})

ax.text(3.74, 65, 'RR',
         {'style':'italic', 'color': 'k', 'fontsize': 10, 'ha': 'center', 'va': 'center',
          'bbox': dict(boxstyle="round", fc="w", ec="k", pad=0.25)})

ax.text(3.8, 12, 'MR',
         {'style':'italic', 'color': 'k', 'fontsize': 10, 'ha': 'center', 'va': 'center',
          'bbox': dict(boxstyle="round", fc="w", ec="k", pad=0.25)})

ax.text(3, 100, 'UE',
         {'style':'italic', 'color': 'k', 'fontsize': 10, 'ha': 'center', 'va': 'center',
          'bbox': dict(boxstyle="round", fc="w", ec="k", pad=0.25)})

# ax.legend(('k.' ,'k--', 'k-.', 'k'),('Second Critical NPR', 'Detachment NPR', 'vonNuemann NPR', 'Third Critical NPR'))
# ax.text(1.78, 80, 'OE - Over Expanded Regime',
#          {'color': 'k', 'fontsize': 10, 'ha': 'center', 'va': 'center'})

# ax.text(1.8, 70, 'UE - Under Expanded Regime',
#          {'color': 'k', 'fontsize': 10, 'ha': 'center', 'va': 'center'})

# ax.text(1.79, 60, 'RR - Regular Reflection Regime',
#          {'color': 'k', 'fontsize': 10, 'ha': 'center', 'va': 'center'})

# ax.text(1.8, 50, 'MR - Mach Reflection Regime',
#          {'color': 'k', 'fontsize': 10, 'ha': 'center', 'va': 'center'})


ax.plot(m, sc, 'ko', markevery=10)
#ax.plot(m, so, 'k')
ax.plot(m, de, 'k--')
ax.plot(m, vn, 'k-.')
ax.plot(m, tc, 'k')

ax.legend(('Second Critical NPR', 'Detachment NPR', 'vonNeumann NPR', 'Perfectly expanded NPR'), )

ax.set_xlim(1.1, 4)
ax.set_ylim(1, 150)

plt.show()
fig.savefig('myimage.svg', format='svg', dpi=1200)