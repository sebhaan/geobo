"""
Calcuation of forward models of gravity and magnetic sensitivity as well as drill core transfromation.

The forward models transform the localized measurement of a remote sensor grid into a 3D representation 
of geophysical properties of a region, here  gravity and magnetic forward models: 
The gravity forward model is defined by using Li's tractable approximation for a 3-D field 
of constant density prisms ( Li and Oldenbur, 3D-inversion of gravity data, 1998}) 
and can be determined analytically. The induced magnetic field calculation uses Li's tractable approximation 
for a 3-D field of prisms of constant magnetic susceptibility, 
which depends on the magnetic mineral content below the surface and is measured 
by the response of their magnetic dipoles induced by the Earth's magnetic field. 
The joint GP inversion takes into account a covariance that exists between density and magnetic susceptibility.

Copyright 2020 Sebastian Haan

This file is part of GeoBO.

GeoBO is free software made available under the AGPL License. 
For details see the LICENSE file.

@author: Sebastian Haan
"""

import numpy as np
from .config_loader import *



def A_sens(magneticField, locations, Edges, func):
	"""
	calculate gravity and magnetic forward model matrix

	PARAMETER

	:param magneticField: ambient magnetic field vector as array with length 3 (x,y,z)
	:param locations: X,Y,Z sensor position with shape (N Sensors, 3)
	:param Edges: Cube edge positions
	:param func: 'grav' for gravity, 'magn' for magnetic

	RETURNS

	Forward model
	eZ
	"""
	Edges = np.asarray(Edges)
	xEdges = Edges[0]
	yEdges = Edges[1]
	zEdges = Edges[2]
	bx = magneticField[0]
	by = magneticField[1]
	bz = magneticField[2]
	nprism = xNcube* yNcube * zNcube
	nedge = (xNcube+1)* (yNcube+1) * (zNcube+1)
	sens = np.zeros((1 * xNcube * yNcube, nprism))
	result_ez = np.zeros((1 * xNcube * yNcube, nedge))

	# For each sensor location
	for n in range(1 * xNcube * yNcube):
		x0 = xEdges - locations[n, 0]
		y0 = yEdges - locations[n, 1]
		z0 = zEdges - locations[n, 2] # z-axis already inverted

		# Lazy edge padding for both grav and mag
		aLongWay = 1e6 # Metres, chosen as in Obsedian
		x0[0] -= aLongWay
		y0[0] -= aLongWay
		x0[-1] += aLongWay
		y0[-1] += aLongWay

		#Precompute eZ for each point
		if func == 'grav':
			eZ = grav_func(x0, y0, z0) # shape (Ncube+1,Ncube+1,Ncube+1)
		elif func == 'magn':
			eZ = magn_func(x0, y0, z0, bx, by, bz)
		else:
			print('function not supported')
		result_ez[n,:] = eZ.reshape((xNcube+1)*(yNcube+1)*(zNcube+1))

		# Compute the sensitivities	
		idx = 0
		for i in range(xEdges.shape[0]-1):
		    for j in range(yEdges.shape[1]-1):
		        for k in range(zEdges.shape[2]-1):
		        	sens[n, idx] = -((eZ[i + 1, j + 1, k + 1] - eZ[i + 1, j + 1, k] - eZ[i + 1, j, k + 1] + eZ[i + 1, j, k]) 
		        		- (eZ[i, j + 1, k + 1] - eZ[i, j + 1, k] - eZ[i, j, k + 1] + eZ[i, j, k]))
		        	idx += 1

	if func == 'grav':
		sens = c_MILLIGALS_UNITS * sens / fcor_grav
	if func == 'magn':
		sens = sens / fcor_mag

	return sens, result_ez


def grav_func(x, y, z):
	"""
	Compute the vertical component of gravity
	Computes the sensitivity for a particular point in the gravity

	PARAMETER

	:param x: x coordinate of the position.
    :param y: y coordinate of the position.
    :param z: z coordinate of the position.
	"""
	eps = 1e-9 
	r = np.sqrt(x**2 + y**2 + z**2)
	func =  x * np.log(y + r) + y * np.log(x + r) - z * np.arctan((x * y) / (z * r + eps))
	return func


def magn_func(x,y,z,bx,by,bz):
	"""
	Compute magnetic forward model
	Calculating the magnetic sensitivity at a particular position relative to the origin.

	PARAMETER

    :param x: x coordinate of the position.
    :param y: y coordinate of the position.
    :param z: z coordinate of the position.
    :param bx: The magnetic field in x-direction at this position
    :param by: The magnetic field in y-direction at this position
    :param bz: The magnetic field in z-direction at this position
	"""
	r = np.sqrt(x**2 + y**2 + z**2)
	#Compute the normalisation factor for the magnetic field
	normB = np.sqrt(bx * bx + by * by + bz * bz)
	func = 1./ normB * ((2. * by * bz * np.log(x + r)) + (2. * bz * bx * np.log(y + r)) + (2. * by * bx * np.log(z + r))
            + (bz * bz - by * by) * np.arctan((x * z) / (y * r)) + (bz * bz - bx * bx) * np.arctan((y * z) / (x * r))) 
	res = -func
	return res


def A_drill(loc, voxelpos):
	"""
	Transform the voxel cube into filter matrix for drill hole with sensitivity 1

	PARAMETER
	
	:param loc: x,y,z drillcore coordinates in shape (Ndrillcore, 3)
	:param voxelpos: x,y,z coordinates in shape (3, Nvoxel)
	"""
	x = voxelpos[0].flatten() #- 1
	y = voxelpos[1].flatten() #- 1
	z = voxelpos[2].flatten() #- 1
	sens = np.zeros((loc.shape[1], xNcube * yNcube * zNcube))
	for i in range(loc.shape[1]):
		coord = loc[:,i]#.astype(int)
		sel = np.where((x == coord[0]) & (y == coord[1]) & (z == coord[2]))
		sens[i, sel] = 1 
	return sens
