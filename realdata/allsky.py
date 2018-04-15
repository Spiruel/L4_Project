import matplotlib.pyplot as plt
import numpy as np #casue it rocks

# get random values to plot
#needs to be in radians from (-pi,pi) & (-pi/2, pi/2)
a = np.arccos(2*np.random.rand(200)-1) #np.random.rand(1000)*2*np.pi - np.pi
d = 2*np.random.rand(200)*2*np.pi

a1 = np.random.rand(100)*2*np.pi - np.pi
d1 = np.random.rand(100)*np.pi - np.pi/2
#todo(make values astropy.cooridnates and use wrap_at(180*u.degree). Allows for more variable data.

num_pts = 200
indices = np.arange(0, num_pts, dtype=float) + 0.5
phi = np.arccos(1 - 2*indices/num_pts)
theta = np.pi * (1 + 5**0.5) * indices

import mpl_toolkits.mplot3d
fig = plt.figure()
x, y, z = np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)
x1, y1, z1 = np.cos(a) * np.sin(d), np.sin(a) * np.sin(d), np.cos(d)
ax = fig.add_subplot(121, projection='3d')
ax.scatter(x, y, z, color='r')
ax.scatter(x1, y1, z1, color='y')
plt.axis('off')
ax.set_aspect('equal')
#plt.show()

#make figure
ax = fig.add_subplot(122, projection='lambert' ) # or mollweide, aitoff, or lambert. [example](http://matplotlib.org/examples/pylab_examples/geo_demo.html) 
ax.scatter(a,d, marker='*', color='y', s=20)
ax.scatter(phi,theta, marker='*', color='r', s=20)
ax.grid(True)
#ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']) #use if you want to change to a better version of RA.
for label in ax.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)
plt.xlabel(r'$\theta$', fontsize=16)
plt.ylabel(r'$\phi$', fontsize=16)
plt.tight_layout()
plt.show()
fig.savefig('allsky.pdf')
