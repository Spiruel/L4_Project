import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import sklearn as sk

from sklearn import cross_validation
from sklearn import datasets
from sklearn import svm

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import confusion_matrix
from sklearn.ensemble import RandomForestClassifier
from sklearn import linear_model

class Plotter():
	def __init__(self):
		self.archs = 33 #number of archives (number starting from 1)
		self.t_end = 10. #length of sim in MYr
		self.plotdir = 'figures/'
		self.fig_ext = '.pdf'
	
		self.Noutputs = np.NaN

		self.vel = [np.NaN for i in range(self.archs)]
		self.sma = [np.NaN for i in range(self.archs)]
		self.ecc = [np.NaN for i in range(self.archs)]
		self.inc = [np.NaN for i in range(self.archs)]
		self.x = [np.NaN for i in range(self.archs)]
		self.y = [np.NaN for i in range(self.archs)]
		self.z = [np.NaN for i in range(self.archs)]
		self.pomega = [np.NaN for i in range(self.archs)]
		self.omega = [np.NaN for i in range(self.archs)]
		self.P = [np.NaN for i in range(self.archs)]
		self.error = [np.NaN for i in range(self.archs)]

		self.load_data()
		#self.load_data2()

		self.snaps = len(self.sma[0]) #initialise no. of snaps after file load
		self.ps_count = np.array([len(self.sma[arch][0])-8 for arch in range(self.archs)]).sum() #total test particle count across all simulations

	def load_data(self):
		for arch in range(self.archs):
			#arch = 0
			directory = '/media/spiruel/7A65-DA8A/good_sims/archive-' + str(arch)
			print 'loading data from: ' + directory
			if os.path.isfile(directory + '/sma' + str(arch) + '.npy') == True:
				self.x[arch] = np.load(directory + '/x' + str(arch) + '.npy')
				self.y[arch] = np.load(directory + '/y' + str(arch) + '.npy')
				self.z[arch] = np.load(directory + '/z' + str(arch) + '.npy')

				AU_m = 149597870700
				vel_conv = AU_m*2.*np.pi/(365*24*60*60) #NEED TO CONVERT FROM SIM UNITS
				self.vel[arch] = np.load(directory + '/vel' + str(arch) + '.npy')*vel_conv

				self.sma[arch] = np.load(directory + '/sma' + str(arch) + '.npy')
				self.ecc[arch] = np.load(directory + '/ecc' + str(arch) + '.npy')
				self.inc[arch] = np.pi - np.load(directory + '/inc' + str(arch) + '.npy') #180deg-inc CORRECTION
				
				#self.pomega[arch] = np.load(directory + '/pomega' + str(arch) + '.npy')
				self.omega[arch] = np.load(directory + '/omega' + str(arch) + '.npy')
				self.P[arch] = np.load(directory + '/P' + str(arch) + '.npy')

				self.error[arch] = np.load(directory + '/error' + str(arch) + '.npy')
			else:
				print 'COULD NOT LOAD DATA FROM: ' + directory




plotter = Plotter()

sma = np.hstack(plotter.sma)
ecc = np.hstack(plotter.ecc)
inc = np.hstack(plotter.inc)
omega = np.hstack(plotter.omega)
#Omega = np.hstack(plotter.sma)
#q = sma*(1-ecc)

ps = sma.shape[1]
selection = np.zeros((ps), dtype=bool)
for p in range(ps):
    if sma[0,p] > 0 and ecc[:,p][sma[:,p] > 0][-1] > 0.8:
        if (sma[:,p] < 5.2044).sum() < 10000:
            print p
            selection[p] = 1
            continue
    selection[p] = 0
        
p = np.where(selection)[0]

#print np.var(sma[::1,p])
plt.plot(sma[:,p])
plt.show()


tp_sel = np.random.randint(0,plotter.snaps,len(p))
print sma[tp_sel,selection]

df = pd.DataFrame({'a':sma[tp_sel,selection], 'e':ecc[tp_sel,selection], 'i':inc[tp_sel,selection], 
                   'w':omega[tp_sel,selection]})

print df

