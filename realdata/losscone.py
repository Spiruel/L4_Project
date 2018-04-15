import numpy as np
import matplotlib.pyplot as plt
import json
import re

from fig_housestyle import font
plt.rcParams["font.family"] = "Andale Modern"

#pre-impact orbits
met_a = [0.416, 0.602, 0.534, 0.671, 0.505, 0.541, 0.417, 0.395, 0.435, 1.175, 0.7644, 0.8130, 0.45, 0.490, 0.571, 0.405, 0.369, 0.386, 0.478, 0.581, 0.503]
met_a = [1/i for i in met_a] #because values given as 1/a in paper
met_e = [0.6711, 0.417, 0.4732, 0.41, 0.55, 0.47, 0.670, 0.690, 0.63, 0.245, 0.31206, 0.22, 0.8, 0.518, 0.431, 0.6023, 0.647, 0.824, 0.526, 0.571, 0.69]
met_i = [10.482, 12.0, 12.27, 4.9, 2.0, 32.2, 11.41, 3.2, 0.0, 9.07, 2.5422, 25.5, 0.26, 28.07, 9.6, 0.832, 2.0, 2.38, 5.51, 4.98, 14.7] 

phas = np.loadtxt('phas.csv', delimiter=',', skiprows=1, usecols=(1,2,3))
asteroids = np.loadtxt('asteroids.csv', delimiter=',', skiprows=1, usecols=(1,2,3))
comets = np.loadtxt('comets.csv', delimiter=',', skiprows=1, usecols=(1,2,3))

plt.plot(asteroids[:,0],asteroids[:,1],'b.', ms=2, label='Asteroids')
plt.plot(comets[:,0],comets[:,1],'r.', ms=3, label='Comets')

plt.plot(phas[:,0],phas[:,1],'y.', ms=3, label='PHAs')
plt.plot(met_a, met_e, 'go', ms=5, markeredgecolor=None, label='Known Impacts')

#plt.legend(frameon=1, loc=4, fancybox=True, framealpha=0.5)
#plt.plot(a[pha].astype(np.float),e[pha].astype(np.float),'b.', ms=3)

axis = plt.gca()
axis.set_xscale('log'); axis.set_yscale('log')
plt.ylim(10**(-2),10**(0.5)); plt.xlim(10**(-1),10**(2))
xs = np.logspace(-1, 2, num=1000)
plt.plot(xs, (1 - (0.985/xs)), 'k--')
plt.plot(xs, ((1.01/xs)-1), 'k--')

#from matplotlib.ticker import AutoMinorLocator
#minorLocator1 = AutoMinorLocator()
#minorLocator2 = AutoMinorLocator()
#plt.gca().xaxis.set_minor_locator(minorLocator1)
#plt.gca().yaxis.set_minor_locator(minorLocator2)

#plt.tick_params(which='both')
#plt.tick_params(which='major', length=7)
#plt.tick_params(which='minor', length=4)

plt.xlabel('a / AU', **font); plt.ylabel('e', **font)
plt.savefig('losscone.pdf', dpi=300, bbox_inches='tight')
plt.savefig('losscone.png', dpi=300, bbox_inches='tight')
