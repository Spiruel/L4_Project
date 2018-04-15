import numpy as np

from math import atan2, asin, cos, sin, radians, degrees
from scipy.interpolate import interp1d

# Quick conv
deg2hr = 24./360.

import math
PI=math.pi
HALFPI = PI/2.0
D2R = PI/180.0
R2D = 1.0/D2R

#
# Copied from https://gist.github.com/barentsen/2367839
# alternately can use astropy tools?
#
def gal2equ(ga):
	"""
	Convert Galactic to Equatorial coordinates (J2000.0)
	(use at own risk)

	Input: [l,b] in decimal degrees
	Returns: [ra,dec] in decimal degrees

	Source:
	- Book: "Practical astronomy with your calculator" (Peter Duffett-Smith)
	- Wikipedia "Galactic coordinates"

	Tests (examples given on the Wikipedia page):
	>>> gal2equ([0.0, 0.0]).round(3)
	array([ 266.405,  -28.936])
	>>> gal2equ([359.9443056, -0.0461944444]).round(3)
	array([ 266.417,  -29.008])
	"""
	l = radians(ga[0])
	b = radians(ga[1])

	# North galactic pole (J2000) -- according to Wikipedia
	pole_ra = radians(192.859508)
	pole_dec = radians(27.128336)
	posangle = radians(122.932-90.0)

	ra = atan2( (cos(b)*cos(l-posangle)), (sin(b)*cos(pole_dec) - cos(b)*sin(pole_dec)*sin(l-posangle)) ) + pole_ra
	dec = asin( cos(b)*cos(pole_dec)*sin(l-posangle) + sin(b)*sin(pole_dec) )

	return np.array([np.mod(degrees(ra),360), degrees(dec)])


def ecl2equ(ecl):
	"""
	Convert ecliptic to equatorial coordinates (use at own risk)
	Formulae from https://en.wikibooks.org/wiki/General_Astronomy/Coordinate_Systems

	Input: [lon,lat] in decimal degrees
	Returns: [ra,dec] in decimal degrees
	"""
	lon = radians(ecl[0])
	lat = radians(ecl[1])

	# Earth's axial tilt w.r.t. ecliptic plane
	# http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
	eps = radians(23.44)

	dec = asin( sin(lat)*cos(eps) + cos(lat)*sin(eps)*sin(lon) )
	#ra = asin( (sin(dec)*cos(eps) - sin(lat)) / (cos(dec)*sin(eps)) )
	ra = lon  # this is a gross approximation.
	return np.array([np.mod(degrees(ra),360), degrees(dec)])

# Convert lists of pairs, some coordinate system y
# to list of pairs in equatorial coordinates (ra/dec)
# using function f2ra (which maps y -> equatorial)
def conv2ra(f2ra, coords):
	arr = np.array(map(lambda x : f2ra(x), coords))
	ra = arr[:,0]
	dec = arr[:,1]
	ra_sort = np.argsort(ra)
	return np.array(zip(ra[ra_sort], dec[ra_sort]))
