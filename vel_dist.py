import numpy as np
import matplotlib.pyplot as plt
#plt.style.use('mphys')
import scipy.stats as st #for maxwell distribution
from scipy import interpolate
from fig_housestyle import font

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

def normalise(v):
    norm=np.linalg.norm(v, ord=1)
    if norm==0:
        norm=np.finfo(v.dtype).eps
    return v/norm

def comet_vel_dist(size,verbose=False):
	#params = np.loadtxt('realdata/comets.csv', delimiter=',', skiprows=1, usecols=[1,2,3])
	#a,e,i = params[:,0], params[:,1], params[:,2]
	#i *= np.pi/180
	#q = (1-e)*a
	#v = np.sqrt( 887.21 * (3-2*np.cos(i)*(q*(1+e))**(0.5) - ((1-e)/q)) )
	#v_i = np.sqrt(124.99 + 887.21*(3-2*np.cos(i)*(q*(1+e))**(0.5) - ((1-e)/q)))
	#v_i = [value for value in v_i if not np.isnan(value)]

	data = np.loadtxt('jeffers_vel_dist.csv', delimiter=',')
	x,y = data[:,0], data[:,1]
	norm_y = normalise(y)
	xs = np.random.choice(x, p=norm_y, size=size)
	if verbose:
		plt.plot(x,norm_y, 'k-', lw=2)
		print '@@@@', (x*norm_y).sum()
		f = interpolate.interp1d(x, norm_y)
		#plt.plot(xs, f(xs),'bo',label='randomly selected velocities from distribution')
		#plt.plot(x, norm_y, 'ro')

		from matplotlib.ticker import AutoMinorLocator
		minorLocator1 = AutoMinorLocator()
		minorLocator2 = AutoMinorLocator()
		plt.gca().xaxis.set_minor_locator(minorLocator1)
		plt.gca().yaxis.set_minor_locator(minorLocator2)

		plt.tick_params(which='both')
		plt.tick_params(which='major', length=7)
		plt.tick_params(which='minor', length=4)

		plt.xlabel('velocity / km$\;$s$^{{-1}}$', **font); plt.ylabel('relative probability', **font)
		plt.savefig('vel_dist.pdf')
		plt.show()
	return xs*1e3
	
#makes a maxwell distribution that fits the velocity distribution
def make_maxwell(loc,scale,correction,xlim,ylim):
	x = np.linspace(xlim, ylim, 4000)
	y = st.maxwell.pdf(x, loc=loc, scale=scale) * correction	
	return x, y


#plot maxwell fits of velocity distribution from Brown et al. (2004)
def plotdistribution(x1,y1,x_ast,y_ast,x_comet,y_comet):
	fig, ax = plt.subplots(figsize=(14, 9))
	ax.plot(x1,np.log10(y1),c="k",lw=3, label='Brown et al. 2004')

	ax.plot(x_ast, np.log10(y_ast), lw=1.5, c="purple", ls="--", label='Maxwell distribution of ast velocity')
	ax.plot(x_comet, np.log10(y_comet), lw=1.5, c="purple", ls="--", label='Maxwell distribution of comet velocity')	
	
	plt.xlim(0,80)
	plt.ylim(-4,2)
	plt.xlabel("$v$ / $\mathrm{km}\;\mathrm{s}^{-1}$", fontsize=30) #velocity
	plt.ylabel("$\mathrm{log_{10}(N)\: / per\:\mathrm{km}\;\mathrm{s}^{-1}}$", fontsize=30) #Relative Number
	plt.xticks(fontsize=22)
	plt.yticks(fontsize=22)
	plt.savefig('maxwell.pdf')
	plt.show()
	

#returns the normalisation of the two fits by dividing the area under both
def integral(x_ast,y_ast,x_comet,y_comet):
	a_ast = np.trapz(y_ast, x=x_ast) #integrate the fit to get the area
	a_comet = np.trapz(y_comet, x=x_comet)
	return a_comet/a_ast


#returns the initial velocities
def calculate_velocities(loc,scale,num,min_vel,max_vel):
	vel = np.zeros(num)
	i=0
	while i < num :
		vel[i] = st.maxwell.rvs(loc=loc, scale=scale)
		if vel[i] > min_vel and vel[i] < max_vel: 
			i+=1

	return vel


#plot histogram of initial velocities
def plotvelocities(vel_ast,vel_comet):
	fig, ax = plt.subplots()
	plt.hist(vel_ast, bins=20, color='purple')
	plt.hist(vel_comet, bins=15, color='deeppink')
	
	plt.xlim(0,90)
	plt.gca().set_yscale("log")
	plt.xlabel("$v$ / $\mathrm{km}\;\mathrm{s}^{-1}$", fontsize=16)
	plt.ylabel("Number of objects", fontsize=16)
	#plt.rc('text', usetex=True)
	#plt.rc('font', family='serif')
	plt.savefig('velocities.png')
	plt.show()


def make_velocities(number):
	veldis_parameters = np.loadtxt('veldis_parameters.csv', delimiter=',')

	x = veldis_parameters[:,0]
	y = veldis_parameters[:,1]

	#make maxwell fit for distribution (Brown et al. 2005)
	#parameters ast
	loc_ast = -30. #right start
	scale_ast = 19. #inclination	
	correction_ast = 115. #altitude
	xlim_ast = -1.
	ylim_ast = 80.
		
	#parameters comet
	loc_comet = 32.
	scale_comet = 14.
	correction_comet = 0.28
	xlim_comet = 32.
	ylim_comet = 90.
	
	
	x_ast, y_ast = make_maxwell(loc_ast,scale_ast, correction_ast, xlim_ast, ylim_ast)
	x_comet, y_comet = make_maxwell(loc_comet,scale_comet, correction_comet, xlim_comet, ylim_comet)


	#plot the maxwell fits for the velocity distribution
	plotdistribution(x,y,x_ast, y_ast,x_comet, y_comet)

	#get the normalisation of the two fits
	norm = integral(x_ast,y_ast,x_comet,y_comet)
	
	#get the velocities
	ast_number = int(np.round(number/(1+norm),0))
	comet_number = int(np.round(norm*ast_number,0))
	print "Number of asts: ", ast_number - 20 #minus the extra ones
	print "Number of comets: ", comet_number
	print "------------------------------------------"
	
	vel_ast, vel_ast = np.array([]), np.array([])
	if ast_number > 0: 
		vel_ast = calculate_velocities(loc_ast,scale_ast,ast_number,17.,50.)	
	if comet_number > 0: 
		vel_comet = calculate_velocities(loc_comet,scale_ast,comet_number,50.,72.)
	vel = np.append(vel_ast, vel_comet)
	
	plotvelocities(vel_ast,vel_comet) #histogram of velocities to check distribution	
	
	#difference in v_esc at 400.000km and at 6461km (R_earth + 90km) is 11102 m/s - 1411 m/s = 9.691 km/s
	difference = 9.691
	vel = vel - difference

	return vel

if __name__ == '__main__':
	number = 10000
	#vel = make_velocities(number + 20)
	comet_vel_dist(100,verbose=1)
