import numpy as np
import rebound

import sys
import h5py

def write_set_to_file(filename, snap, sim):
	part_num = sim.N
	particles = sim.particles
	orbits = sim.calculate_orbits()

	print 'now saving data...',

	f = h5py.File(filename,'a')
	f.create_group(str(snap))

	f.create_dataset(str(snap)+'/x', data=np.array([particles[j].x for j in range(part_num-1)]))
	f.create_dataset(str(snap)+'/y', data=np.array([particles[j].y for j in range(part_num-1)]))
	f.create_dataset(str(snap)+'/z', data=np.array([particles[j].z for j in range(part_num-1)]))

	f.create_dataset(str(snap)+'/vx', data=np.array([particles[j].vx for j in range(part_num-1)]))
	f.create_dataset(str(snap)+'/vy', data=np.array([particles[j].vy for j in range(part_num-1)]))
	f.create_dataset(str(snap)+'/vz', data=np.array([particles[j].vz for j in range(part_num-1)]))

	f.create_dataset(str(snap)+'/ax', data=np.array([particles[j].ax for j in range(part_num-1)]))
	f.create_dataset(str(snap)+'/ay', data=np.array([particles[j].ay for j in range(part_num-1)]))
	f.create_dataset(str(snap)+'/az', data=np.array([particles[j].az for j in range(part_num-1)]))

	f.create_dataset(str(snap)+'/d', data=np.array([orbits[j].d for j in range(part_num-1)]))
	f.create_dataset(str(snap)+'/v', data=np.array([orbits[j].v for j in range(part_num-1)]))
	f.create_dataset(str(snap)+'/h', data=np.array([orbits[j].h for j in range(part_num-1)]))
	f.create_dataset(str(snap)+'/P', data=np.array([orbits[j].P for j in range(part_num-1)]))
	f.create_dataset(str(snap)+'/n', data=np.array([orbits[j].n for j in range(part_num-1)]))

	f.create_dataset(str(snap)+'/a', data=np.array([orbits[j].a for j in range(part_num-1)]))
	f.create_dataset(str(snap)+'/e', data=np.array([orbits[j].e for j in range(part_num-1)]))
	f.create_dataset(str(snap)+'/inc', data=np.array([orbits[j].inc for j in range(part_num-1)]))

	f.create_dataset(str(snap)+'/Omega', data=np.array([orbits[j].Omega for j in range(part_num-1)]))
	f.create_dataset(str(snap)+'/omega', data=np.array([orbits[j].omega for j in range(part_num-1)]))
	f.create_dataset(str(snap)+'/pomega', data=np.array([orbits[j].pomega for j in range(part_num-1)]))

	f.create_dataset(str(snap)+'/f', data=np.array([orbits[j].f for j in range(part_num-1)]))
	f.create_dataset(str(snap)+'/M', data=np.array([orbits[j].M for j in range(part_num-1)]))
	f.create_dataset(str(snap)+'/theta', data=np.array([orbits[j].theta for j in range(part_num-1)]))
	f.create_dataset(str(snap)+'/T', data=np.array([orbits[j].T for j in range(part_num-1)]))

	f.create_dataset(str(snap)+'/error', data=abs(((sim.calculate_energy()) - E0)/E0))
	f.create_dataset(str(snap)+'/time', data=sim.t)
	
	print 'saved snap,', snap

	f.close()


directory = sys.argv[1].split('/')[0]+'/'
if 'output' not in directory:
	print 'Not found'
	raise SystemError

sim = rebound.Simulation.from_file(directory+"/rewound_data.bin")

hf = h5py.File(directory+"out.hdf5", "r")
print hf.keys()
print hf.get('0/error')[0]
Noutputs = 20000
t_end = 1e7
raise SystemError

year = 2.*np.pi
print 'Integrating to', t_end, 'yrs (' + str(t_end*year) +') with stepsize of', str(sim.dt)
times = np.linspace(0.,t_end*year, Noutputs+1)

E0 = sim.calculate_energy()#+gr_potential(sim)

for i,time in enumerate(times):
    if time <= sim.t:
	print 'done', time
	pass
    #sim.integrate(time)
    print 'reached time,', time
    if time != 0: 
        print'{:.2f}%'.format((sim.t/float(t_end*year))*100)

    #write_set_to_file(directory+'out.hdf5', i, sim)

    if i%10 == 0:
	print 'checkpoint...'
	#save_file = directory+"rewound_data.bin"
	#sim.save(save_file)
