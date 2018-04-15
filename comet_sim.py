import numpy as np
import rebound

from os import mkdir
from time import strftime
import h5py

from vel_dist import comet_vel_dist

def write_set_to_file(filename, snap, sim):
	part_num = sim.N
	particles = sim.particles
	orbits = sim.calculate_orbits()

	print 'now saving data...',

	f = h5py.File(output+filename,'a')
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

sim = rebound.Simulation.from_file("checkpoint_mercurius.bin")
sim.move_to_com()

#reverse direction of solar system
for p in sim.particles:
	p.vx *= -1
	p.vy *= -1
	p.vz *= -1

sim.status()

earth = sim.particles[3]
#define units
year = 2.*np.pi
AU_m = 149597870700
R = 4E8/AU_m
vel_conv = AU_m*2.*np.pi/(365*24*60*60) #convert from sim units to m/s

#6 day timestep
sim.dt = (6/365.)*year
   
output = 'output_'+(strftime("%d%m%Y%H%m%s"))+'/'
mkdir(output)
print 'OUTPUT:', output

def make_particles(N, velocities):
	for i in range(N):
		launch = False
		while not launch:
			#pick random point on surface of a sphere the size of earth-moon system
			theta = np.arccos(2*np.random.random() - 1)
			phi = np.random.random() * 2 * np.pi
			x = R * np.sin(theta) * np.cos(phi)
			y = R * np.sin(theta) * np.sin(phi)
			z = R * np.cos(theta)

			dir_vec = np.array([x,y,z])
		
			rand = velocities[i]
			vx, vy, vz = rand * (dir_vec/np.sqrt(dir_vec.dot(dir_vec)))

			res_x, res_y, res_z = earth.x+x, earth.y+y, earth.z+z #resultant coordinates of particle relative to solar system in AU
			res_vx, res_vy, res_vz = (earth.vx*vel_conv) + vx, (earth.vy*vel_conv) + vy, (earth.vz*vel_conv) + vz #result velocities of particles relative to solar system in m/s
		
			length_velocity = np.sqrt((res_vx)**2 + (res_vy)**2 + (res_vz)**2)

			#ensure meteoroid vel is less than escape vel and that it does not move towards the earth, but rather moves increasingly far from it
			if length_velocity < 71e3 and np.sqrt((x*AU_m)**2 + (y*AU_m)**2 + (z*AU_m)**2) - np.sqrt(((x*AU_m)+res_vx)**2 + ((y*AU_m)+res_vy)**2 + ((z*AU_m)+res_vz)**2) < 0:
				launch = True

		sim.add(m=0,x=res_x, y=res_y, z=res_z, vx=res_vx/vel_conv, vy=res_vy/vel_conv, vz=res_vz/vel_conv)

number = 100
t = 0
launch_time = 11.8*year
velocities = comet_vel_dist(number)

count = 0

#sim.integrate( np.random.uniform(0,launch_time) )

while t <= launch_time:
	#sim.integrate(t)

	#launch new particles
	if t < launch_time:
		#choose random dt close to 0.1 years
		dt = np.round(np.random.normal(loc = 0.2, scale = 0.1), decimals = 2) * year	

		#choose random amount of particles to launch with a normal distribution
		launch_now = int(np.round((np.random.normal(loc = number / ((launch_time) / dt), scale = 2)), 0))
		#make particles
		if launch_now>0:
			print count, launch_now
			count += launch_now
			make_particles(launch_now, velocities)

		if t + dt > launch_time:
			dt = launch_time - t #to get exactly the same time snaps
		
	t += dt
	
print 'reached time', t/year, 'and launched', len(sim.particles)-9, 'particles'

Noutputs = 20000
t_end = 1e7
print 'Integrating to', t_end, 'yrs (' + str(t_end*year) +') with stepsize of', str(sim.dt)
times = np.linspace(0.,t_end*year, Noutputs+1)

E0 = sim.calculate_energy()#+gr_potential(sim)

for i,time in enumerate(times):
    sim.integrate(time)
    print 'reached time,', time
    if time != 0: 
        print'{:.2f}%'.format((sim.t/float(t_end*year))*100)

    write_set_to_file('out.hdf5', i, sim)

    if i%10 == 0:
	print 'checkpoint...'
	save_file = output+"rewound_data.bin"
	sim.save(save_file)
