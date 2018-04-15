import h5py
import numpy as np

import os
from fig_housestyle import font
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

def moving_average(x, n, type='simple'):
	"""
	compute an n period moving average.

	type is 'simple' | 'exponential'

	"""
	x = np.asarray(x)
	if type == 'simple':
		weights = np.ones(n)
	else:
		weights = np.exp(np.linspace(-1., 0., n))

	weights /= weights.sum()

	a = np.convolve(x, weights, mode='full')[:len(x)]
	a[:n] = a[n]
	return a

def for_all_methods(decorator):
    def decorate(cls):
        for attr in cls.__dict__: # there's propably a better way to do this
            if callable(getattr(cls, attr)):
                setattr(cls, attr, decorator(getattr(cls, attr)))
        return cls
    return decorate


#@for_all_methods(mpltex.acs_decorator)
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

	def stats(self):
		bound_ps = np.array([(self.ecc[arch][0,:] > 1).sum() for arch in range(self.archs)]).sum()
		end_bound_ps = np.array([(self.ecc[arch][self.snaps-1,:] < 1).sum() for arch in range(self.archs)]).sum()
		print '''
		Total number of particles: '''+str(self.ps_count)+'''
		Number of initially bound particles: '''+str(self.ps_count-bound_ps),str(format(float(self.ps_count-bound_ps)*100/self.ps_count, '.2f'))+'''%
		Number of bound particles at end of simulation: '''+str(end_bound_ps)
		return self.ps_count-bound_ps

	def retrograde(self):
		#check retrogrades
		retro = np.zeros((self.snaps))
		prog = np.zeros((self.snaps))

		for t in range(self.snaps):
			for arch in range(self.archs):
				for p in range(8,len(self.sma[arch][0])):
					# 90 degrees in rads
					if self.sma[arch][t,p] > 0:
						if self.inc[arch][t,p] > 90 * (np.pi/180) and self.inc[arch][t,p] < 270 * (np.pi/180):
							retro[t] += 1
						elif self.inc[arch][t,p] < 90 * (np.pi/180) or self.inc[arch][t,p] > 270 * (np.pi/180):
							prog[t] += 1

		plt.plot(prog/prog.max(), label='prog')
		plt.plot(retro/retro.max(), label='retro')
		plt.plot((prog/prog.max())/(retro/retro.max()))
		plt.legend()
		plt.show()
		print "Mean number retrogrades particles", np.mean(retro)

	def getaway(self):
		#check how many particles leave the solar system
		out, tot, vel_tot, vel_out = 0,0,0,0	
		for i in range(self.archs):
			tot += len(self.sma[i][0,10:])
			vel_tot += len(self.vel[i])
			for particle in range(len(self.sma[i][self.snaps-1,:])):
				if self.sma[i][self.snaps-1,particle] < 0 and self.ecc[i][self.snaps-1,particle] > 1:
					out += 1
					if particle in self.vel[i]:
						vel_out +=1
	
		print "Total number of meteoroids that leave the solar system: ", out, "/", tot
		#print "Total number of cometary meteoroids that leave the solar system:", vel_out, "/", vel_tot
	

	def leaverate(self):
		#plot of rate at which meteoroids leave the solar system
		time = np.linspace(0,self.t_end,self.snaps)
		leave = np.zeros((self.snaps))
		for i in range(self.archs):
			for j in range(self.snaps):
				for particle in range(8,len(self.sma[i][0,:])):
					if self.sma[i][j,particle] < 0 and self.ecc[i][j,particle] > 1:
						leave[j] += 1
	
		fig, ax = plt.subplots()
		plt.plot(time,leave, lw=1, zorder=0)
	
		#fit the leaverate from the moment it is constant
		y=[]
		#SHOULD BE 4000 NOT 40
		coeffs, error = np.polyfit(time[400:],leave[400:], 1, cov=True)
		for x in time:
			y.append(coeffs[0]*x + coeffs[1])
	
		print "Intersect at x=0, y=", coeffs[1] 
		print 'fit is: ', coeffs[0], '*x + ', coeffs[1] 
		print 'errors are: ', np.sqrt(error[0][0]), '*x + ', np.sqrt(error[1][1]) 
		#ax.annotate('- - $\mathrm{N_{leave}=(88.19 \pm 0.05)t + (826.9 \pm 0.4)}$', xy=(1,1600), xytext=(1,1600), color='purple', fontsize=18)
	
		plt.plot(time, y,'--', c='purple', zorder=1, alpha=1) #prolongation of plot above
		#plt.ylim(0,1800)
		plt.xlabel("t / Myr", **font)
		plt.ylabel("N$_{leave}$", **font)
		#plt.xticks(fontsize=22)
		#plt.yticks(fontsize=22)
		plt.savefig(self.plotdir + 'leaverate' + self.fig_ext)
		plt.show()
	
		return leave

	def timelosscone(self):
		#time particles stay in loss cone before hitting Earth
		leave1 = np.zeros((self.snaps)) 				#amount of particles that leave loss cone for the first time
		first_leave = np.array([])
		for i in range(self.archs):
			for k in range(8,len(self.sma[i][0,8:])):
				for j in range(self.snaps):
					a = self.sma[i][j,k]
					e = self.ecc[i][j,k]
					if a > 1 and a < 100:
						if e < -(1.01/a) + 1:
							leave1[j] += 1
							first_leave = np.append(first_leave, j)
							break
											
					elif a < 1 and a > 0:
						if e < (0.985/a) - 1:
							leave1[j] += 1
							first_leave = np.append(first_leave, j)
							break
	
		self.cum_leave1 = np.cumsum(leave1)
		leave1[leave1==0] = np.NaN
		print "Median time a particle is in the loss cone before leaving it for the first time:", np.nanmedian(first_leave) * (self.t_end/self.snaps), "Myr"

	
		time = np.arange(self.snaps)
		leave2 = np.zeros((self.snaps)) #amount of particles that leave loss cone for the last time
		last_leave = np.append(first_leave, j)
		for i in range(self.archs):
			for k in range(8,len(self.sma[i][0,8:])):
				for j in time[::-1]: 
					a = self.sma[i][j,k]
					e = self.ecc[i][j,k]
					#only looks at bound particles
					if a > 1 and a < 100:
						if e >= -(1.01/a) + 1:
							leave2[j] += 1
							last_leave = np.append(last_leave, j)
							break
											
					elif a < 1 and a > 0:
						if e >= (0.985/a) - 1:
							leave2[j] += 1
							last_leave = np.append(last_leave, j)
							break
	
		self.cum_leave2 = np.cumsum(leave2)
		leave2[leave2==0] = np.NaN
		print "Median time a particle is in the loss cone before leaving it for the last time:", np.nanmedian(last_leave) * (self.t_end/self.snaps), "Myr"


	def cumulative_plot(self):
			self.timelosscone() #gets cum_leave1 and cum_leave2
	
			#cumulative plot
			fig, ax = plt.subplots(figsize=(14, 9))
			t = np.linspace(0,self.t_end,self.snaps)
			plt.plot(t, self.cum_leave1, color='purple', label='First leave')
			plt.plot(t[:-1], self.cum_leave2[:-1], color='blue', label='Last leave')
		
			lgd = plt.legend(scatterpoints=1, loc=4,ncol=1, markerscale=3, frameon=False)
			ltext = plt.gca().get_legend().get_texts()
			plt.setp(ltext[0], color = 'purple')	
			plt.setp(ltext[1], color = 'blue')	
			plt.setp(ltext)
		
			#plt.xticks(fontsize=22)
			#plt.yticks(fontsize=22)
			plt.xlabel("$t_{lc}$ / Myr", labelpad=14, **font)
			plt.ylabel("$N_{cum}$", **font)
			plt.tight_layout()
			plt.savefig(self.plotdir + 'leavecone' + self.fig_ext)
			plt.show()

	def aeplot(self):
		from matplotlib.colors import LogNorm
		fig, ax = plt.subplots(2,3, sharex=True, sharey=True)
		plt.subplots_adjust(wspace=.05, hspace=.05)

		#[0,self.snaps/5,2*self.snaps/5,3*self.snaps/5,4*self.snaps/5,self.snaps-1]
		for pos,i in enumerate([0,self.snaps/40,self.snaps/20, self.snaps/10, self.snaps/2, self.snaps-1]):
			for arch in range(self.archs):
				a = self.sma[arch] #only consider bound particles
				e = self.ecc[arch]
				inc = self.inc[arch]

				axis = ax.flat[pos]

				axis.plot((a[i,:8]),(e[i,:8]), 'kx', ms=4)
				sc = axis.scatter((a[i][8:]),(e[i][8:]),c=(inc[i][8:]*(180/np.pi)), s=2, edgecolors='none',cmap='viridis')
			
			axis.set_xscale('log'); axis.set_yscale('log')
			axis.text(0.7,1.25, 't = %s Myr' % float('%.2g' % (-i *float(self.t_end)/self.snaps)), zorder=0, fontsize=9, **font)
			axis.set_ylim(10**(-2),10**(0.5)); axis.set_xlim(10**(-1.5),10**(2.5))
			xs = np.logspace(-1, 2, num=1000)
			axis.plot(xs, (1 - (0.985/xs)), 'k--')
			axis.plot(xs, ((1.01/xs)-1), 'k--')
			

		fig.add_subplot(111, frameon=False)
		# hide tick and tick label of the big axes
		plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
		
		#cbaxes = fig.add_axes([0.8, 0.1, 0.03, 0.8]) 
		cb = fig.colorbar(mappable=sc, orientation="horizontal", ax=ax.ravel().tolist(), cmap='viridis')

		cb.set_label('i / $^\degree$', **font)

		plt.xlabel("a / AU", **font)
		plt.ylabel("e", **font)
		
		plt.savefig(self.plotdir + 'ae_migration' + self.fig_ext)
		plt.show()

	def param_timeplot(self):
		time = np.linspace(0,self.t_end, self.snaps)

		fig, ax = plt.subplots(2,1, sharey=True)

		mmrs = {'4:1':2.06,
			'3:1':2.5,
			'5:2':2.82,
			'7:3':2.95,
			'2:1':3.27
			}

		exiters = {'83':1,
			'16':1,
			'97':2,
			'63':4,
			'33':5,
			'79':6,
			'24':7,
			'59':8,
			'75':9,
			'61':10,
			'53':11,
			'22':12,
			'29':13
			}
		#num = [624, 586, 668, 223]
		
		color=iter(plt.cm.spectral_r(np.linspace(0,1,20)))

		num = np.arange(0,100,1)#np.random.randint(9,self.sma[0].shape[1]-8,5)
		for p in num:
			xs = 1/(self.sma[0][:,p])
			xs[self.sma[0][:,p] < 0] = 0
			if np.diff(xs).max() > 0.01:
				c = next(color)
				ax.flat[0].plot(time, xs,  '-', label=str(p), c=c)
				ax.flat[1].plot(time, xs,  '-', label=str(p), c=c)

				exit = np.where(xs < 10**(-4))[0]
				if len(exit) == 0:
					exit = -1
				else:
					exit = exit[0]

				exit_time = time[exit]

				if p == 59: exit_time = 4.02899 #bodge

				if xs[exit:].max() < 10**(-4) or p in [83,59] and p not in [83]:
					if (self.inc[0][0,p]) > 90 * (np.pi/180) and (self.inc[0][0,p]) < 270 * (np.pi/180) and exit_time > 0.05:
						print 'PROGRADE', str(exiters[str(p)])
					if exit_time > 2:
						ax.flat[0].plot(exit_time, 10**(-4), 'ko', ms=3)
						ax.flat[0].text(exit_time-0.15, 10**(-3.5), str(exiters[str(p)]), bbox={'facecolor':'white', 'edgecolor':'k', 'pad':2, 'alpha':0.5}, **font)
					ax.flat[0].plot(exit_time, 10**(-4), 'ko', ms=3)
					ax.flat[1].plot(exit_time, 10**(-4), 'ko', ms=3)
					if exit_time > 0.05 and exit_time < 2:
						if exiters[str(p)] == 1:
							y_displace = -.1
						elif exiters[str(p)] == 2:
							y_displace = .75
						else: 
							y_displace = .5
						
						ax.flat[1].text(exit_time-0.02, 10**(-4+y_displace), str(exiters[str(p)]), bbox={'facecolor':'white', 'edgecolor':'k', 'pad':2, 'alpha':0.5}, **font)

		#bodge
		ax.flat[1].plot(0.272075*2, 10**(-4), 'ko', ms=3)
		ax.flat[1].text((0.272075*2)-0.02, 10**(-3.5), '3', bbox={'facecolor':'white', 'edgecolor':'k', 'pad':2, 'alpha':0.5}, **font)

		for a in ax.flat:
			a.set_yscale('log')

			a.axhline(0, color='k', ls='--')
			a.axhline(10**(-4), color='k', lw=1, ls='--')
			a.axhline(0.0292, color='k', lw=1, ls='-.')
		
		plt.tick_params(which='both')
		plt.tick_params(which='major', length=7)
		plt.tick_params(which='minor', length=4)
		
		ax.flat[0].set_ylim(10**(-4.5),)
		ax.flat[1].set_xlim(0,2)
		#plt.legend()
		#plt.xlim(0,len(xs)*1e-4)#; plt.ylim(-0.04,0.01)

		#for m in mmrs.values():
		#	plt.axhline(1/m, color='k', ls=':')

		fig = ax.flat[0].get_figure()
		ax = fig.add_subplot(111, frameon=False)
		# hide tick and tick label of the big axes
		ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
		
		ax.set_xlabel('t / Myr', **font); ax.set_ylabel('energy / AU$^{-1}$', **font)

		plt.savefig(self.plotdir + 'ae_timeplot' + self.fig_ext)
		plt.show()

	def sma_histogram(self):
		f, axes = plt.subplots(6, figsize=(14,9), sharex=True, sharey=True)
		
		for pos,i in enumerate(range(0,self.snaps,int(0.2*self.snaps))):
			ax = axes[pos]

			a = np.array([])
			for arch in range(self.archs):
				a = np.append(a, self.sma[arch][i,:])

			y,x = np.histogram(a, bins=np.arange(0,11,.05)) #only considering +ve
			ax.plot(x[:-1],y,'k-')
			#ax.set_ylim(-10,90)
			ax.set_xlim(0,10)
			ax.text(6,0.4*max(y),'t = ' + str(i* float(self.t_end*1e6/self.snaps)) + ' yrs')

		plt.xlabel('a / AU',)
		plt.text(-0.9,280,'N')#plt.ylabel('N', fontsize=16)
		# Fine-tune figure; make subplots close to each other and hide x ticks for
		# all but bottom plot.
		f.subplots_adjust(hspace=0)
		plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
		plt.savefig(self.plotdir + 'sma_histogram' + self.fig_ext)
		plt.show()

	def error_plot(self):
		xs = np.linspace(0,self.t_end,self.snaps)
		#plt.plot(xs, 1e-10*np.sqrt(xs),"--",color="black", lw=2) brouwers law not apply?
		plt.gca().set_yscale("log")

		#error for each individual arch
		#for i in range(self.archs):
		#	plt.plot(xs, self.error[i], marker=".",color="#555555",alpha=0.05)
		#	plt.plot(xs, moving_average(self.error[i], 40), color="black", alpha=0.4)
		#	plt.show()

		#mean of all errors in one plot
		plt.plot(xs, np.array(self.error).mean(axis=0), marker=".",color="#555555",alpha=0.05)
		plt.plot(xs, moving_average(np.array(self.error).mean(axis=0), 40), color="black")

		plt.ylabel('relative energy error', **font)
		plt.xlabel('t / Myr', **font)
		plt.savefig(self.plotdir + 'error' + self.fig_ext)
		plt.show()

	def dens_plot(self):
		fig = plt.figure()
		axis = plt.gca()

		sma = np.array([])
		ecc = np.array([])
		for arch in range(self.archs):
			mask = self.sma[arch] > 0
			sma = np.append(sma, self.sma[arch][mask])
			ecc = np.append(ecc, self.ecc[arch][mask])

		#axis.plot(sma.flatten()[:200000], ecc.flatten()[:200000], 'r.')
		nbins = 200
		plt.hist2d(np.log(sma.flatten()), np.log(ecc.flatten()), bins=nbins, cmap='viridis', range=[[-1, 2], [-2, 0]], normed=1)
		
		axis.set_ylim((-2),(0.0)); axis.set_xlim((-1),(2))

		plt.colorbar()
		#xs = np.logspace(-1, 2, num=1000)
		#axis.plot(xs, (1 - (0.985/xs)), 'k--')
		#axis.plot(xs, ((1.01/xs)-1), 'k--')
		axis.set_xlabel('log(a) / AU', **font); axis.set_ylabel('log(e)', **font)
		plt.savefig(self.plotdir + 'timeweighted_ae_plot' + self.fig_ext)
		plt.show()

	def a_timeplot(self):
		fig = plt.figure()
		axis = plt.gca()

		sma = np.array([])
		for arch in range(self.archs):
			mask = self.sma[arch] > 0
			sma = np.append(sma, self.sma[arch][mask])

		#axis.plot(sma.flatten()[:200000], ecc.flatten()[:200000], 'r.')
		#nbins = 200
		times = np.linspace(0,self.t_end,self.snaps)
		repeat = len(sma.flatten())/float(self.snaps)
		times = np.repeat(times,repeat)

		nbins=1000
		print len(sma), len(times)	
		plt.hist2d(sma.flatten()[:3520176], times, bins=nbins, cmap='viridis', range=[[0, 4], [0,10]], normed=1)
		#
		#axis.set_ylim((-2),(0.0)); axis.set_xlim((-1),(2))

		#plt.colorbar()
		#xs = np.logspace(-1, 2, num=1000)
		#axis.plot(xs, (1 - (0.985/xs)), 'k--')
		#axis.plot(xs, ((1.01/xs)-1), 'k--')
		#axis.set_xlabel('log(a) / AU'); axis.set_ylabel('log(e)')
		#plt.savefig(self.plotdir + 'timeweighted_ae_plot' + self.fig_ext)
		plt.show()

	def close_approach(self, target=4):
		fig = plt.figure()
		ax = plt.subplot(111)
		ax.set_xlabel("time / Myr", **font)
		ax.set_ylabel("distance / AU", **font)
	
		rand = np.random.randint(0,100)
		while self.sma[0][0,rand] < 0:
			rand = np.random.randint(0,100) #don't get an initial hyperbolic particle

		distance = np.sqrt(np.square(self.x[0][:,rand]-self.x[0][:,target])+np.square(self.y[0][:,rand]-self.y[0][:,target]))
		times = np.linspace(0, self.t_end, len(distance))
		plt.plot(times, distance)
		plt.plot(times, moving_average(distance,50), lw=2)
		closeencountertime = times[np.argmin(distance)]
		print("Minimum distance (%f AU) occurred at time: %f years." % (np.min(distance),closeencountertime))

		#for vel in self.vel:
		#	print (vel[-1,:] * (30 / vel[-1,3])).max()
		plt.show()


	def exotic(self):
		
		leave = self.leaverate()
		#check how many particles are exotic
		exo = np.zeros((self.snaps))
		#amount of exotic per unit time

		par = [np.zeros((self.snaps)) for i in range(self.archs)] #amount of time being exotic per particle
		aphe, peri = 0,0  #amount of exotic particles left and right
		curr_snap = 10000 #choose a random snapshot poin for analysis

		returns = np.array([])
		for arch in range(self.archs):
			indiv_return = np.zeros((self.snaps, len(self.sma[arch][0])))
			for snap in range(self.snaps):
				for particle in range(8,len(self.sma[arch][0])):
					a = self.sma[arch][snap, particle]
					e = self.ecc[arch][snap, particle]

					if a > 1 and a < 100:
						if e < -(1.01/a) + 1:
							exo[snap] += 1
							par[arch][particle] += 1
							if snap==curr_snap:
								aphe += 1
							indiv_return[snap,particle] = 1
					
					elif a < 1 and a > 0:
						if e < (0.985/a) - 1:
							exo[snap] += 1
							par[arch][particle] += 1
							if snap==curr_snap:
								peri += 1
							indiv_return[snap,particle] = 1
					else:
						#'reject unbound particle!'
						continue
			
			for particle in range(8,len(self.sma[arch][0])):
				returns = np.append(returns, (np.diff(indiv_return[:,particle]) == 1).sum())
			# count number of times particle goes back to loss cone (from 0 to 1)

		print 'min, max, mean, median number of returns to loss cone:', np.min(returns), np.max(returns), np.mean(returns), np.median(returns)

		#Calculate mean time particles are exotic
		tot_exotic = 0
		mean = []
		for arch in range(self.archs):
			for particle in range(8,len(self.sma[arch][0])):
				if par[arch][particle] > 0:
					tot_exotic += 1
					mean.append(par[arch][particle])

		print "Total particles that are exotic at some point: ", tot_exotic
		print "Mean time particles are exotic: ", np.mean(mean)/1000, " Myr"
	
		print "At time: - " + str(format(float(curr_snap)/(self.snaps/10),'.3f')) + " Myr --------------------------------------"
		print "Amount of exotic particles in total:", exo[curr_snap]
		print "Amount of exotic particles with Earth as aphelion:", aphe
		print "Amount of exotic particles with Earth as perihelion:", peri

		
		#plot time vs. ratio exotic/total
		exo = np.array(exo)

		tot = np.array([self.ps_count - leave[i] for i in range(self.snaps)])	#total amount of particles in sol sys per time
		exo_ratio = (tot-exo) / self.ps_count											#exotic particles ratio

		fig, ax = plt.subplots()
		time = np.linspace(0,self.t_end,self.snaps)
	
		#fit
		y=[]
		#SHOULD BE 2000 NOT 200
		coeffs, error = np.polyfit(time[10000:],exo_ratio[10000:], 1, cov=True)
		for x in time:
			y.append(coeffs[0]*x + coeffs[1])
	
		print 'fit is: ', coeffs[0], '*x + ', coeffs[1]
		print 'errors are: ', np.sqrt(error[0][0]), '*x + ', np.sqrt(error[1][1])
		#ax.annotate('- - $\mathrm{N_{out}/N_{tot}=}$(-2.7E-5 $\pm$ 1.6E-5)t + (0.25 $\pm$ 1.0E-4)', xy=(4,0.26), xytext=(4,0.26), color='purple', fontsize=18)
	
		plt.plot(time, y,'--', c='purple', zorder=2, alpha=1)
		#plt.plot(time, tot, color='y', zorder=1, lw=1)
		#plt.plot(time, tot-exo, color='r', zorder=1, lw=1)
		plt.plot(time, exo_ratio, zorder=1, lw=1, alpha=0.2)
		plt.plot(time, moving_average(exo_ratio,100), zorder=1, lw=2)
		#plt.xlim(0,10)
		#plt.ylim(0,210)
		#plt.xticks(fontsize=22)
		#plt.yticks(fontsize=22)
		plt.xlabel("t / Myr")
		plt.ylabel("N$_{haz}$/N$_{tot}$")

		minorLocator1 = AutoMinorLocator()
		minorLocator2 = AutoMinorLocator()
		axis = plt.gca()
		axis.xaxis.set_minor_locator(minorLocator1)
		axis.yaxis.set_minor_locator(minorLocator2)

		plt.tick_params(which='both')
		plt.tick_params(which='major', length=7)
		plt.tick_params(which='minor', length=4)


		plt.savefig(self.plotdir + 'exotic_ratio' + self.fig_ext)
		plt.show()

	def get_tisserand(self):
		tiss = [np.NaN for i in range(self.archs)]
		a_j = 5.2044 #j sma in AU

		for arch in range(self.archs):
			tiss_param = (a_j/np.abs(self.sma[arch])) + 2*np.sqrt( (1-self.ecc[arch]**2)*(np.abs(self.sma[arch])/a_j) ) * np.cos(self.inc[arch])

			tiss[arch] = tiss_param
			
		return tiss
		#plt.hist(self.tiss[0])
		#plt.show()

	def infoloss_plot(self):
		t = np.linspace(0, self.t_end, self.snaps)

		fig, ax = plt.subplots(4,2, sharey=True, sharex=True)
		plt.subplots_adjust(wspace=0.08, hspace=0.1)

		planets = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
		for p in range(8):
			axis = ax.flat[p]
			init = np.sqrt(self.x[0][:,p]**2 + self.y[0][:,p]**2 + self.z[0][:,p]**2)
			maxy = 0

			ys = [np.nan for arch in range(1,self.archs)]
			maxy = [np.nan for arch in range(1,self.archs)]

			for arch in range(1,self.archs):
				y = moving_average( np.abs( np.sqrt(self.x[arch][:,p]**2 + self.y[arch][:,p]**2 + self.z[arch][:,p]**2) - init ), 500)
				
				ys[arch-1] = y
				maxy[arch-1] = y.max()
		
			for i in range(len(ys)):
				axis.plot(t,ys[i]/maxy[i], '-', alpha=0.6)
			axis.text(0.5, 0.7, planets[p], bbox={'facecolor':'white', 'edgecolor':'none', 'pad':2, 'alpha':0.5}, **font)

			minorLocator1 = AutoMinorLocator()
			minorLocator2 = AutoMinorLocator()
			axis.xaxis.set_minor_locator(minorLocator1)
			axis.yaxis.set_minor_locator(minorLocator2)

			plt.tick_params(which='both')
			plt.tick_params(which='major', length=7)
			plt.tick_params(which='minor', length=4)

		plt.ylim(0,1.1)

		fig.add_subplot(111, frameon=False)
		# hide tick and tick label of the big axes
		plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
		plt.grid(False)
		plt.xlabel("time / Myr", **font)
		plt.ylabel("normalised separation distance", **font)

		plt.savefig(self.plotdir + 'info_loss' + self.fig_ext)
		plt.show()

	def migration_path(self):
		import matplotlib

		for arch in range(self.archs):
			ns = np.random.random_integers(9,100,5)
			for n in ns:
			    	#plt.plot(self.sma[arch][:,n],self.ecc[arch][:,n],'k:', zorder=-999)
				plt.scatter((self.sma[arch][:,n]),(self.ecc[arch][:,n]),c=np.log10(self.inc[arch][:,n]), s=10, vmin=0, vmax=1, edgecolors='none')

		plt.gca().set_xscale('log'); plt.gca().set_yscale('log')

		#plt.plot((self.sma[arch][:8,i]),(self.ecc[arch][:8,i]), 'kx', ms=1)
		plt.ylim(10**(-2),10**(1)); plt.xlim(10**(-2),10**(2))
		cb = plt.colorbar(norm=matplotlib.colors.LogNorm())
		cb.set_label('$\log_{10}(i)$ $/deg$')

		xs = np.logspace(-1, 2, num=1000)
		plt.plot(xs, (1 - (0.985/xs)), 'k--')
		plt.plot(xs, ((1.01/xs)-1), 'k--')
	
		#plt.plot(xs, (1 - (1.3814/xs)), 'k--')
		#plt.plot(xs, ((1.6660/xs)-1), 'k--')
		plt.plot(xs, (1 - (4.9501/xs)), 'k--')
		plt.plot(xs, ((5.4588/xs)-1), 'k--')


		plt.savefig(self.plotdir + 'migration' + self.fig_ext)
		plt.show()

	def dynam_life(self):
		#init_bound = frq.sum()#float(self.stats())

		dls = np.array([]) #dynamical lifetimes

		for arch in range(self.archs):
			ts,ps = self.sma[arch].shape
			for p in range(ps):
				for t in range(ts): #for every time for every particle
					sma = self.sma[arch][t,p]
					if sma < 0 and t < 1: #ignore bound particles
						break
					elif sma < 0:
						dls = np.append(dls, t-1) #choose last timestep as dynamical lifetime
						break
				#else:
					#dls = np.append(dls, ts)

		norm1 = float(self.t_end)/self.snaps

		frq, edges = np.histogram(dls*norm1, bins=20)
		plt.bar(edges[:-1], frq/float(frq.sum()), width=np.diff(edges), ec="k", align="edge", alpha=0.75, color='gray', edgecolor='None')
		plt.axvline(np.median(dls)*norm1,color='k',ls='--')
		print 'MEDIAN DYN LIFETIME:', np.median(dls)*norm1
		plt.xlabel('dynamical lifetime / Myr'); plt.ylabel('$N/N_{bound}$')

		plt.savefig(self.plotdir + 'dynamical_lifetime' + self.fig_ext)
		plt.show()

		norm2 = self.t_end*(200./self.snaps)
		frq, edges = np.histogram(dls[dls<200]*norm2, bins=10)
		plt.axvline(np.median(dls[dls<200])*norm2,color='k',ls='--')
		plt.bar(edges[:-1], frq/float(frq.sum()), width=np.diff(edges), ec="k", align="edge", alpha=0.75, color='gray', edgecolor='None')
		plt.show()

	def period_plot(self):
		t = np.linspace(0, self.t_end, self.snaps)
		period = (self.P[0][:,8:])/(2*np.pi)
		
		plt.plot(t, period/12., '-', label='jupiter')
		plt.plot(t, period/165., '--', label='neptune')
		plt.show()

	def comettypes_plot(self):
		times = np.linspace(0, self.t_end, self.snaps)

		tiss = self.get_tisserand()

		jfc = np.zeros((self.snaps))
		hfc = np.zeros((self.snaps))
		nic = np.zeros((self.snaps))
		ec = np.zeros((self.snaps))

		htc = np.zeros((self.snaps))

		for arch in range(self.archs):
			for t in range(self.snaps):
				hfc[t] = np.logical_and(tiss[arch][t,8:] < 2, self.P[arch][t,8:]/(2*np.pi) < 200).sum()
				jfc[t] = np.logical_and(tiss[arch][t,8:] < 3, np.logical_and(tiss[arch][t,8:] > 2, self.P[arch][t,8:]/(2*np.pi) < 20)).sum()

				nic[t] = (tiss[arch][t,8:] < 2).sum()
				ec[t] = (tiss[arch][t,8:] > 2).sum()

				htc[t] = np.logical_and(self.P[arch][t,8:]/(2*np.pi) > 20, self.P[arch][t,8:]/(2*np.pi) < 200).sum()

		plt.plot(times, moving_average(hfc, 200), label='hfc')
		plt.plot(times, moving_average(jfc, 200), label='jfc')
		plt.legend()
		plt.show()

		plt.plot(times, moving_average(ec, 200), label='ec')
		plt.plot(times, moving_average(nic, 200), label='nic')
		plt.legend()
		plt.show()

		plt.plot(times, moving_average(htc, 200), label='halley type')
		plt.legend()
		plt.show()

	def resonance_plot(self):
		smoothed_smas = np.array([])
		for arch in range(self.archs):
			for p in range(self.sma[arch].shape[-1]):
				sma = self.sma[arch][:,p]
				smoothed_smas = np.append(smoothed_smas, moving_average(sma, 10**4))

		plt.hist(smoothed_smas, range=(0,32), bins=320, facecolor='green', edgecolor='green', normed=True, histtype='stepfilled')

		a_j = 5.2044 #j sma in AU
		a_mmr = lambda s,r: a_j*(float(s)/float(r))**(2./3.)
		
		for s in range(1,30,2):
			plt.text(a_mmr(s,1)-0.61,0.8,'1:'+str(s))
			plt.plot(a_mmr(s,1),0.8-0.05,'ko', ms=2)


			plt.text(a_mmr(s,2)-0.61,0.7,'2:'+str(s))
			plt.plot(a_mmr(s,2),0.7-0.05,'ko', ms=2)

			#plt.text(a_mmr(s,3)-0.61,0.5,'3:'+str(s))
			#plt.plot(a_mmr(s,3),0.5,'ko')

			#interior resonances
			plt.text(a_mmr(1,s)-0.61,0.6,str(s)+':1')
			plt.plot(a_mmr(1,s),0.6-0.05,'ko', ms=2)
			plt.text(a_mmr(3,s)-0.61,0.5,str(s)+':3')
			plt.plot(a_mmr(3,s),0.5-0.05,'ko', ms=2)
		

		minorLocator1 = AutoMinorLocator()
		minorLocator2 = AutoMinorLocator()
		axis=plt.gca()
		axis.xaxis.set_minor_locator(minorLocator1)
		#axis.yaxis.set_minor_locator(minorLocator2)

		plt.tick_params(which='both')
		plt.tick_params(which='major', length=7)
		plt.tick_params(which='minor', length=4)

		plt.xlim(0,32)
		plt.xlabel('<a> / AU'); plt.ylabel('Probability')
		plt.savefig(self.plotdir + 'resonance' + self.fig_ext)
		plt.show()

		print format(np.logical_and(smoothed_smas > (a_j - 1), smoothed_smas < (a_j + 1)).sum() * (100 / float(len(smoothed_smas))), '.2f') + '% cumulative time in 1:1'

		print format(np.logical_and(smoothed_smas > (a_mmr(7,1) - 1), smoothed_smas < (a_mmr(7,1) + 1)).sum() * (100 / float(len(smoothed_smas))), '.2f') + '% cumulative time in 1:7'

	def comet_class(self):
		a_j = 5.2044 #j sma in AU

		jfc = np.array([])
		hfc = np.array([])
		centaurs = np.array([])
		tnos = np.array([])
		encke = np.array([])
		hyper = np.array([])
	
		encke_dyn_life = np.array([])
		
		for arch in range(self.archs):
			for p in range(self.sma[arch].shape[1]):
				tiss = np.empty((self.snaps))
				tiss[:] = np.nan
				for t in range(self.snaps):
					sma = self.sma[arch][t,p]
					ecc = self.ecc[arch][t,p]
					inc = self.inc[arch][t,p]

					if sma > 0 and ecc < 1:
						param = (a_j/np.abs(sma)) + 2*np.sqrt( (1-ecc**2)*(np.abs(sma)/a_j) ) * np.cos(inc)
						tiss[t] = param

				P = self.P[arch][:,p]/(2*np.pi)
				a = self.sma[arch][:,p]
				dyn_life = (a > 0).sum()

				jfc = np.append(jfc, np.logical_and(np.logical_and(tiss < 3, tiss > 2), P  < 20).sum()/float(dyn_life) )
				hfc = np.append(hfc, np.logical_and(tiss < 2, P  < 200).sum()/float(dyn_life) )
				centaurs = np.append(centaurs, np.logical_and(np.logical_and((1-self.ecc[arch][:,p])*a > 5.2, a < 30.1), P > 0).sum()/float(dyn_life))
				tnos = np.append(tnos, np.logical_and(a > 30.1, P > 0).sum()/float(dyn_life))
				encke = np.append(encke, np.logical_and(P>0, np.logical_and(tiss > 3, a < a_j)).sum()/float(dyn_life))
				hyper = np.append(hyper, (a< 0).sum()/float(self.snaps))

				if encke[-1] > 0:
					encke_dyn_life = np.append(encke_dyn_life, dyn_life)
			
		for name, com_type in zip(['jfc', 'hfc', 'centaurs', 'tnos', 'encke'], [jfc, hfc, centaurs, tnos, encke]):
			print name + ': ' + format( com_type[com_type>0].mean()*100, '.2f') + '%'
		print 'encke_dyn_life:', np.median(encke_dyn_life)*(float(self.t_end)/float(self.snaps))

	def radiant_plot(self):
		from radiant_plot import deg2hr, gal2equ, ecl2equ, conv2ra
		from scipy.interpolate import interp1d

		# ------------------------------------
		# SET UP -- basic landmarks/boundaries
		# ------------------------------------

		# Galactic disk
		gal_lon = np.linspace(0,360,4000)
		gal_lat = np.zeros(gal_lon.shape)
		galplane = conv2ra(gal2equ, zip(gal_lon, gal_lat))
		galplane_plus10 = conv2ra(gal2equ, zip(gal_lon, +10 * np.ones(gal_lon.shape)))
		galplane_minus10 = conv2ra(gal2equ, zip(gal_lon, -10 * np.ones(gal_lon.shape)))

		# Ecliptic plane, band of sky
		ecl_lon = np.linspace(0,360,4000)
		ecl_lat = np.zeros(ecl_lon.shape)
		ecliptic = conv2ra(ecl2equ, zip(ecl_lon, ecl_lat))
		ecliptic_plus35 = conv2ra(ecl2equ, zip(ecl_lon, +5 * np.ones(ecl_lon.shape)))
		ecliptic_minus35 = conv2ra(ecl2equ, zip(ecl_lon, -5 * np.ones(ecl_lon.shape)))

		# Plot galactic plane, ecliptic plane bands
		plt.plot(deg2hr * galplane[:,0], galplane[:,1], 'k-', lw=5, alpha=1, zorder=5)
		plt.plot(deg2hr * galplane_plus10[:,0], galplane_plus10[:,1], 'k--', alpha=1, zorder=5)
		plt.plot(deg2hr * galplane_minus10[:,0], galplane_minus10[:,1], 'k--', alpha=1, zorder=5)

		plt.plot(deg2hr * ecliptic[:,0], ecliptic[:,1], 'g-', lw=2, alpha=1, zorder=5)
		plt.plot(deg2hr * ecliptic_plus35[:,0], ecliptic_plus35[:,1], 'g--', alpha=1, zorder=5)
		plt.plot(deg2hr * ecliptic_minus35[:,0], ecliptic_minus35[:,1], 'g--', alpha=1, zorder=5)

		# Fill between bands
		# Must interpolate because x-vector (RA) is different between bounding curves
		ra_vec = np.linspace(0,360,4000)
		# Note: interpolation fails at edges.  No big deal
		f_gal_top = interp1d(galplane_plus10[:,0], galplane_plus10[:,1], bounds_error=False)
		f_gal_bot = interp1d(galplane_minus10[:,0], galplane_minus10[:,1], bounds_error=False)
		plt.fill_between(deg2hr*ra_vec, f_gal_bot(ra_vec), f_gal_top(ra_vec), color='gray', alpha=0.2, zorder=5)

		f_ecl_top = interp1d(ecliptic_plus35[:,0], ecliptic_plus35[:,1], bounds_error=False)
		f_ecl_bot = interp1d(ecliptic_minus35[:,0], ecliptic_minus35[:,1], bounds_error=False)
		plt.fill_between(deg2hr*ra_vec, f_ecl_bot(ra_vec), f_ecl_top(ra_vec), color='green', alpha=0.2, zorder=5)

		# Mark galactic centre
		plt.scatter(deg2hr*gal2equ([0,0])[0], gal2equ([0,0])[1], marker='o',c='k', s=200, alpha=1, zorder=5)

		#plt.scatter(deg2hr*ecl2equ([0,0])[0], ecl2equ([0,0])[1], marker='o',c='g', s=200, alpha=1, zorder=5)

		AU_m = 149597870700
		R = 4E8/AU_m


		earth_x = 0.351458957848
		earth_y=0.927984574948
		earth_z=-0.000164653819752

		for arch in range(self.archs):

			############################

			jfcs = np.zeros((self.sma[arch].shape[1]), dtype=bool)
			hfcs = np.zeros((self.sma[arch].shape[1]), dtype=bool)
			centaurs = np.zeros((self.sma[arch].shape[1]), dtype=bool)
			tnos = np.zeros((self.sma[arch].shape[1]), dtype=bool)
			enckes = np.zeros((self.sma[arch].shape[1]), dtype=bool)
			hypers = np.zeros((self.sma[arch].shape[1]), dtype=bool)

			for p in range(self.sma[arch].shape[1]):
				tiss = np.empty((self.snaps))
				tiss[:] = np.nan
				for t in range(self.snaps):
					sma = self.sma[arch][t,p]
					ecc = self.ecc[arch][t,p]
					inc = self.inc[arch][t,p]

					a_j = 5.2044
					if sma > 0 and ecc < 1:
						param = (a_j/np.abs(sma)) + 2*np.sqrt( (1-ecc**2)*(np.abs(sma)/a_j) ) * np.cos(inc)
						tiss[t] = param

				P = self.P[arch][:,p]/(2*np.pi)
				a = self.sma[arch][:,p]
				dyn_life = (a > 0).sum()

				jfc = np.logical_and(np.logical_and(tiss < 3, tiss > 2), P  < 20).sum()
				hfc = np.logical_and(tiss < 2, P  < 200).sum()
				centaur = np.logical_and(np.logical_and((1-self.ecc[arch][:,p])*a > 5.2, a < 30.1), P > 0).sum()
				tno = np.logical_and(a > 30.1, P > 0).sum()
				encke = np.logical_and(P>0, np.logical_and(tiss > 3, a < a_j)).sum()
				hyper = (a < 0).sum()

				max_class = max(jfc, hfc, centaur, tno, encke, hyper)
				if max_class == jfc:
					#print 'jfc'
					jfcs[p] = 1
				elif max_class == hfc:
					#print 'hfc'
					hfcs[p] = 1
				elif max_class == centaur:
					#print 'centaurs'
					centaurs[p] = 1
				elif max_class == tno:
					#print 'tnos'
					tnos[p] = 1
				elif max_class == encke:
					#print 'encke'
					enckes[p] = 1
				elif max_class == hyper:
					#print 'hyper'
					hypers[p] = 1

			##################################

			mask = self.sma[arch][10000,:] > 0
			x = (self.x[arch][0,:] - earth_x) * AU_m
			y = (self.y[arch][0,:] - earth_y) * AU_m
			z = (self.z[arch][0,:] - earth_z) * AU_m

			r = np.sqrt(x**2 + y**2 + z**2)
			theta = np.arcsin(z/r)*(180/np.pi)#- (np.arccos(z/r) - (np.pi/2))
			phi = np.arctan2(y,x)*(180/np.pi) + 180 #longitude

			#theta = np.array([ np.arccos(2*np.random.random() - 1) for i in range(2000)])*(180/np.pi) - 90
			#phi = np.array([ np.random.random() * 2 * np.pi for i in range(2000)])*(180/np.pi)

			print 'theta',theta.min(), theta.max()
			print 'phi',phi.min(), phi.max()

			radec = conv2ra(ecl2equ, zip(phi, theta))
	
			#from scipy.stats import gaussian_kde	
			#xy = np.vstack([radec[:,0],radec[:,1]])
			#z = gaussian_kde(xy)(xy)

			#plt.scatter(radec[:,0]*deg2hr, radec[:,1], c=z, edgecolor='none')

			plt.plot(radec[:,0][jfcs]*deg2hr, radec[:,1][jfcs], '^', ms=2, zorder=99, markeredgewidth=0)
			plt.plot(radec[:,0][hfcs]*deg2hr, radec[:,1][hfcs], 'x', ms=2, zorder=99, markeredgewidth=0)
			plt.plot(radec[:,0][centaurs]*deg2hr, radec[:,1][centaurs], 'o', ms=2, zorder=99, markeredgewidth=0)
			plt.plot(radec[:,0][tnos]*deg2hr, radec[:,1][tnos], 's', ms=2, zorder=99, markeredgewidth=0)
			plt.plot(radec[:,0][enckes]*deg2hr, radec[:,1][enckes], 'D', ms=2, zorder=99, markeredgewidth=0)
			#plt.plot(radec[:,0][hypers]*deg2hr, radec[:,1][hypers], '*', ms=2, zorder=99)
	
			import seaborn as sns
			sns.kdeplot(radec[:,0]*deg2hr, radec[:,1], shade=True, cmap="binary")

		plt.xlim([24,0])
		plt.xticks([24,21,18,15,12,9,6,3,0],[str(i)+'$^h$' for i in [24,21,18,15,12,9,6,3,0]])
		plt.ylim([-90,90])
		plt.yticks(np.linspace(-80,80,9),[str(i)+'$^\degree$' for i in np.linspace(-80,80,9)])

		plt.xlabel('RA', **font); plt.ylabel('Dec', **font)
		plt.savefig(self.plotdir + 'radiants' + self.fig_ext)
		plt.show()

plotter = Plotter()

plotter.stats()
#plotter.retrograde()
#plotter.getaway()
#plotter.leaverate()
#plotter.exotic()
#plotter.timelosscone()
#plotter.cumulative_plot()
#plotter.aeplot()
#plotter.param_timeplot()
#plotter.sma_histogram()
#plotter.error_plot()
#plotter.dens_plot()
#plotter.a_timeplot()
#plotter.close_approach()

#plotter.get_tisserand()
#plotter.infoloss_plot()

#plotter.migration_path()
#plotter.dynam_life()
#plotter.period_plot()
#plotter.comettypes_plot()
#plotter.resonance_plot()

#plotter.comet_class()
plotter.radiant_plot()

#q = (1-plotter.ecc[0])*plotter.sma[0]
#plt.plot(plotter.omega[0][:,8:], q[:,8:], '.')
#plt.show()
