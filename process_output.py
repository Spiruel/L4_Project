import numpy as np
import h5py #to read hdf5 files

# Run the process
def run_archive(arch):
	directory = '/media/spiruel/7A65-DA8A/good_sims/archive-' + str(arch)

	print 'Archive-', arch
	f = h5py.File(directory + '/out.hdf5', 'r')   # 'r' means that the file is open in read-only mode

	#define variables
	data = f

	steps = len(data) #amount of time steps
	number = len( np.array(data['0/x']) ) #amount of particles (solar system + meteoroids)

	x = [np.NaN for k in range(steps)]
	y = [np.NaN for k in range(steps)]
	z = [np.NaN for k in range(steps)]
	vx = [np.NaN for k in range(steps)]
	vy = [np.NaN for k in range(steps)]
	vz = [np.NaN for k in range(steps)]
	sma = [np.NaN for k in range(steps)]
	ecc = [np.NaN for k in range(steps)]
	inc = [np.NaN for k in range(steps)]
	d = [np.NaN for k in range(steps)]
	P = [np.NaN for k in range(steps)]
	n = [np.NaN for k in range(steps)]
	Omega = [np.NaN for k in range(steps)]
	omega = [np.NaN for k in range(steps)]
	pomega = [np.NaN for k in range(steps)]
	true_anom = [np.NaN for k in range(steps)]
	M = [np.NaN for k in range(steps)]
	h = [np.NaN for k in range(steps)]

	v = [np.NaN for k in range(steps)]

	error = np.array([])
	
	'''#Attributes
	d 	(float) radial distance from reference
	v 	(float) velocity relative to central object's velocity
	h	(float) specific angular momentum
	P 	(float) orbital period (negative if hyperbolic)
	n 	(float) mean motion (negative if hyperbolic)
	a 	(float) semimajor axis
	e 	(float) eccentricity
	inc 	(float) inclination
	Omega 	(float) longitude of ascending node
	omega 	(float) argument of pericenter
	pomega 	(float) longitude of pericenter
	f 	(float) true anomaly
	M 	(float) mean anomaly
	l 	(float) mean longitude = Omega + omega + M
	theta 	(float) true longitude = Omega + omega + f
	T 	(float) time of pericenter passage'''

	step = 0 #year
	for name in data:
		print name
		step = int(name)
		x[step] = np.array(data[name + '/x'])
		y[step] = np.array(data[name + '/y'])
		z[step] = np.array(data[name + '/z'])

		vx[step] = np.array(data[name + '/vx'])
		vy[step] = np.array(data[name + '/vy'])
		vz[step] = np.array(data[name + '/vz'])

		sma[step] = np.array(data[name + '/a'])
		ecc[step] = np.array(data[name + '/e'])
		inc[step] = np.array(data[name + '/inc'])

		d[step] = np.array(data[name + '/d'])
		P[step] = np.array(data[name + '/P'])
		n[step] = np.array(data[name + '/n'])

		Omega[step] = np.array(data[name + '/Omega'])
		omega[step] = np.array(data[name + '/omega'])
		pomega[step] = np.array(data[name + '/pomega'])

		true_anom[step] = np.array(data[name + '/f'])
		M[step] = np.array(data[name + '/M'])
		h[step] = np.array(data[name + '/h'])

		v[step] = np.array(data[name + '/v'])

		error = np.append(error, data[name + '/error'])

	#save data

	np.save(directory + '/sma' + str(arch) + '.npy', sma)
	np.save(directory + '/ecc' + str(arch) + '.npy', ecc)
	np.save(directory + '/inc' + str(arch) + '.npy', inc)
	np.save(directory + '/x' + str(arch) + '.npy', x)
	np.save(directory + '/y' + str(arch) + '.npy', y)
	np.save(directory + '/z' + str(arch) + '.npy', z)
	np.save(directory + '/vx' + str(arch) + '.npy', vx)
	np.save(directory + '/vy' + str(arch) + '.npy', vy)
	np.save(directory + '/vz' + str(arch) + '.npy', vz)
	np.save(directory + '/vel' + str(arch) + '.npy', v)
	np.save(directory + '/d' + str(arch) + '.npy', d)
	np.save(directory + '/P' + str(arch) + '.npy', P)
	np.save(directory + '/n' + str(arch) + '.npy', n)
	np.save(directory + '/Omega' + str(arch) + '.npy', Omega)
	np.save(directory + '/omega' + str(arch) + '.npy', omega)
	np.save(directory + '/pomega' + str(arch) + '.npy', pomega)
	np.save(directory + '/f' + str(arch) + '.npy', f)
	np.save(directory + '/M' + str(arch) + '.npy', M)
	np.save(directory + '/h' + str(arch) + '.npy', h)
	np.save(directory + '/error' + str(arch) + '.npy', error)

		
	f.close() #close file




# ---------------------------- MAIN ------------------------------------


from joblib import Parallel, delayed
import multiprocessing

num_cores = multiprocessing.cpu_count() - 2 # number of cores to use for parallelization
print 'Using %i cores for parallelization' %(num_cores)

#Parallel(n_jobs = num_cores)(delayed(run_archive) (arch) for arch in range(21,39))

[run_archive(arch) for arch in range(35,38)]
