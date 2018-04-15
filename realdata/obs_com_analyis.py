import numpy as np
import json
from fig_housestyle import font

import matplotlib.pyplot as plt

#comet_els = np.genfromtxt('CometEls.txt', dtype='str', delimiter=[12,6,3,8,10,10,10,10,10,10,6,5,59,9])
#a,b = comet_els.shape
#for i in range(a):
#	for j in range(b):
#		comet_els[i][j] = comet_els[i][j].strip()

#fig = plt.figure(figsize=(9,9))

json_data=open('cometels.json').read()
comets = json.loads(json_data)
print 'loaded com!'
json_data=open('mpcorb_extended.json').read()
asteroids = json.loads(json_data)
print 'loaded ast!'

for a in asteroids[:50000]:
	e = a['e']
	if e == 1.:
		continue
	a = a['Perihelion_dist'] / (1-e)
	plt.plot(a,e,'b.',ms=1)
for c in comets:
	e = c['e']
	if e == 1.:
		continue
	a = c['Perihelion_dist'] / (1-e)
	plt.plot(a,e,'r.')

plt.xlim(0,8); plt.ylim(0,1)
#for a in comets:
#	qs = np.append(qs, i['Perihelion_dist'])

sma = np.linspace(1.8,8,200)
#plt.plot(sma, 1 - (4.9501/sma), 'k-')
#plt.plot(sma, (5.4588/sma) - 1, 'k-')

a_j = 5.2044
a_m = 1.524

plt.plot(sma, 1-((a_j/(4*sma)) * (3 - (a_j/sma))**2),  'k--', lw=2) #tj = 3
sma = np.linspace(2.5,8,1000)
plt.plot(sma, 1-((a_j/(4*sma)) * (2 - (a_j/sma))**2),  'k--', lw=2) #tj = 2

rot1 = 1-((a_j/(4*np.array([3.5,4.5]))) * (2 - (a_j/np.array([3.5,4.5])))**2)
rot2 = 1-((a_j/(4*np.array([6.1,7.1]))) * (3 - (a_j/np.array([6.1,7.1])))**2)
rot1 = np.degrees(np.arctan2(rot1[1]-rot1[0], 1))

rot2 = 1-((a_j/(4*np.array([6.2,7.4]))) * (3 - (a_j/np.array([6.2,7.4])))**2)
rot2 = np.degrees(np.arctan2(rot2[1]-rot2[0], 1))

plt.text(3.5,0.92, r'$T_J = 2$', rotation=rot1*4, size=13)
plt.text(6.2,0.085, r'$T_J = 3$', rotation=rot2*4, size=13)

plt.axvline(a_j, ymax=0.8, alpha=.5, color='k', lw=2)
plt.text(a_j-0.15,0.82, r'$a_J$', size=13)

plt.axvline(a_m, ymax=0.8, alpha=.5, color='k', lw=2)
plt.text(a_m-0.15,0.82, r'$a_M$', size=13)

#plt.axvline(2.5, ymax=0.9, lw=2, ls='--')
#plt.text(2.5,0.92, r'3:1', size=13)
#plt.axvline(2.82, ymax=0.9, lw=2, ls='--')
#plt.text(2.82,0.92, r'5:2', size=13)
#plt.axvline(2.95, ymax=0.9, lw=2, ls='--')
#plt.text(2.95,0.92, r'7:3', size=13)
#plt.axvline(3.27, ymax=0.9, lw=2, ls='--')
#plt.text(3.27,0.92, r'2:1', size=13)

from matplotlib.ticker import AutoMinorLocator
minorLocator1 = AutoMinorLocator()
minorLocator2 = AutoMinorLocator()
plt.gca().xaxis.set_minor_locator(minorLocator1)
plt.gca().yaxis.set_minor_locator(minorLocator2)

plt.tick_params(which='both')
plt.tick_params(which='major', length=7)
plt.tick_params(which='minor', length=4)

plt.xlabel('a / AU', **font); plt.ylabel('e', **font)
plt.savefig('a_e_tisserand.pdf')

