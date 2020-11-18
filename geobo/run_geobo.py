"""
This script creates reconstructed Cubes with mean subtracted properties of density, magnetic susceptibility, 
and  drill core properties plus their predicted variance cubes (see inversion.py for more details)

See settings.yaml for specifications.

Copyright 2020 Sebastian Haan

This file is part of GeoBO.

GeoBO is free software made available under the AGPL License. 
For details see the LICENSE file.

@author: Sebastian Haan
"""
import os
import sys
import numpy as np
import pandas as pd
import rasterio
from scipy import interpolate
from scipy.ndimage.interpolation import zoom
from scipy.optimize import minimize, least_squares, shgo
from scipy.stats import norm, pareto
from scipy import interpolate
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm


def read_surveydata(plot = True):
	"""
	Reading gravity data and magnetic data, including cropping and downsampling to Cube size

	PARAMETER
	:param plot: boolean, if True plots survey data and downsampled

	RETURN
	gravity data
	magnetic data
	sensor locations
	"""
	if FNAME_gravsurvey is not None:
		gravimg = rasterio.open(os.path.join(inpath,FNAME_gravsurvey))
		grav = gravimg.read(1)
	else:
		grav = None
	if FNAME_magsurvey is not None:
		magimg = rasterio.open(os.path.join(inpath, FNAME_magsurvey))
		mag = magimg.read(1)
	else:
		mag = None
	# downsample to Cube voxelsize
	zoomfac = xNcube *1./grav.shape[1]
	grav2 = zoom(grav, zoomfac)
	assert grav2.shape == (yNcube, xNcube)
	zoomfac = xNcube *1./mag.shape[1]
	mag2 = zoom(mag, zoomfac)
	assert mag2.shape == (yNcube, xNcube)
	# Create survey sensor coordinates:
	x_s = np.linspace(0.5, xNcube - 0.5, xNcube)  * xvoxsize
	y_s = np.linspace(0.5, yNcube - 0.5, yNcube)  * yvoxsize
	z_s = zmax + zoff
	x_sensor, y_sensor, z_sensor = np.meshgrid(x_s,y_s,z_s)
	locations = np.asarray([x_sensor.flatten(), y_sensor.flatten(), z_sensor.flatten()]).T
	if not os.path.exists(outpath):
		os.makedirs(outpath)
	if plot:
		extent=[xmin,xmax,ymin,ymax]
		plt.imshow(grav, aspect = 'equal', cmap = 'viridis', extent=extent, origin='lower')
		plt.colorbar()
		plt.savefig(os.path.join(outpath,'gravfield.png'))
		plt.clf()
		plt.imshow(mag, aspect = 'equal', cmap = 'viridis', extent=extent, origin='lower')
		plt.colorbar()
		plt.savefig(os.path.join(outpath, 'magfield.png'))
		plt.clf()
		plt.imshow(grav2, aspect = 'equal', cmap = 'viridis', extent=extent, origin='lower')
		plt.colorbar()
		plt.savefig(os.path.join(outpath, 'gravfield_downsampled.png'))
		plt.clf()
		plt.imshow(mag2, aspect = 'equal', cmap = 'viridis', extent=extent, origin='lower')
		plt.colorbar()
		plt.savefig(os.path.join(outpath,'magfield_downsampled.png'))
		plt.close("all")
	return grav2.flatten(), mag2.flatten(), locations


def read_drilldata(features):
	"""
	Reading drill data

	PARAMETER
	:param features: List of drill features of interest (see settings.yaml)

	RETURN
	Drill data
	x,y,z coordinates of drilldata
	x,y,z min,max arrays of drilldata
	"""
	if FNAME_drilldata is not None:
		drill = pd.read_csv(os.path.join(inpath,FNAME_drilldata))
		drill = drill[(drill.x >= xmin) & (drill.x <= xmax) & (drill.y >= ymin) & (drill.y <= ymax) & (drill.z <= zmax) & (drill.z >= zmin)].copy()
	else: 
		drill = None
	# Select data only within extent
	# Convert to local coordinate system with origin X,Y = 0,0
	drill['x'] = drill['x'] - xmin
	drill['y'] = drill['y'] - ymin
	# Set up voxel coordinate grid
	xdrill = drill['x'].values
	ydrill = drill['y'].values
	zdrill = drill['z'].values 
	try:
		drillfirst = drill.groupby('SiteID').first()
		drilllast = drill.groupby('SiteID').last()
		xdrillminmax = np.asarray([drillfirst.x.values, drilllast.x.values]).T
		ydrillminmax = np.asarray([drillfirst.y.values, drilllast.y.values]).T
		zdrillminmax = np.asarray([drillfirst.z.values, drilllast.z.values]).T
	except:
		xdrillminmax = 0.#np.asarray([xdrill.min(), xdrill.min()])
		ydrillminmax = 0.#np.asarray([ydrill.min(), ydrill.min()])
		zdrillminmax = 0.#np.asarray([zdrill.min(), zdrill.min()])
	coord = np.vstack([xdrill, ydrill, zdrill]).T
	drilldata = []
	for feature in features:
		data = drill[feature]
		drilldata.append(align_drill(coord, data))
	return np.asarray(drilldata), coord, (xdrillminmax, ydrillminmax, zdrillminmax)


def align_drill(coord, data):
	"""
	Convert drill-core data in Model Cube shape with iteration over all voxels

	PARAMETER
	param coord: xyz Coordinates of drillcore, shape (N_drill, 3)
	param data: 1 dim data for drillcore, shape (N_drill)

	RETURN
	drilldata in voxel shape
	"""
	dx = xvoxsize
	dy = yvoxsize
	dz = zvoxsize
	data = np.asarray(data)
	res = np.zeros_like(xxx)
	for ix in range(xxx.shape[0]):
		for iy in range(xxx.shape[1]):
			for iz in range(xxx.shape[2]):
				sel = np.where(((xxx[ix,iy,iz] - dx) <= coord[:, 0]) & (coord[:, 0] < (xxx[ix,iy,iz] + dx))
					& ((yyy[ix,iy,iz] - dy) <= coord[:, 1]) & (coord[:, 1] < (yyy[ix,iy,iz] + dy))
					& ((zzz[ix,iy,iz] - dz) <= coord[:, 2]) & (coord[:, 2] < (zzz[ix,iy,iz] + dz)))
				if np.size(data[sel]) > 0:
					#print(sel)
					m = np.nanmean(data[sel])
					if np.isfinite(m):
						res[ix,iy,iz] = m
	return res



"""
Below: defineition of acquisition function functions for BO

The key of BO is the acquisition function, which typically has to balance between 
a) exploration, i.e., querying points that maximise the information gain and minimize the uncertainty of a model, 
b) exploitation, i.e. querying points that maximise the reward 
(e.g. concentrating search in the vicinity locations with high value such as minerals), and 
c) minimize the number of samples given an expensive cost function for any new measurement. 

Exploration-exploitation and cost trade-off parameters can be set in settings.yaml
"""

def futility_vertical(params, costs = None):
	""" 
	Utility/Aquisition function for bayesian optimisation assuming vertical drillcores

	PARAMETER
	params: parameters for drillhole (x_drill,y_drill) 
	param costs: cube of costs with same shape as reconstructed cube

	RETURN
	Output of utility function (scalar)
	"""
	if costs is None:
		costs = drill_rec * 0.
	params = np.asarray(params)
	xmaxvox = drill_rec.shape[0] - 1 
	ymaxvox = drill_rec.shape[1] - 1
	if np.isfinite(params).all():
		xd = int(np.round(params[0]))
		yd = int(np.round(params[1]))
		if (xd > 0) & (xd < xmaxvox) & (yd > 0) & (yd < ymaxvox):
			func = np.sum(drill_rec[xd, yd, :]) + kappa * np.sqrt(np.sum(drill_var[xd, yd, :])) - beta * np.sum(costs[xd,yd, :])
		else:
			func = -np.inf
	else:
		func = -np.inf
	return -func # (negative if optimiser use minimize)


def futility_drill(params, costs = None): 
	"""
	Calculates utility function for proposed non-vertical drillcore, 
	which takes into account the azimuth and dip in addition to location and drill-core length

	PARAMETER
	param params: [x0, y0, azimuth, dip] x0 and y0 in meters; azimuth and dip in degree 
	param length_newdrill: length of drillcore in meters
	param costs: cube of costs with same shape as reconstructed cube

	RETURN
	Output of utility function (scalar)
	"""
	if costs is None:
		costs = drill_rec * 0.
	length_newdrill = zLcube
	x0, y0, azimuth, dip = params
	nstep = int(2 * length_newdrill/np.min([xvoxsize, yvoxsize, zvoxsize]))
	rladder = np.linspace(0,length_newdrill, nstep)
	x0 = rladder * 0 + x0
	y0 = rladder * 0 + y0
	z0 = rladder * 0 + zmax
	azimuth = rladder * 0 + azimuth
	dip = rladder * 0 + dip
	try:
		xdrillnew, ydrillnew, zdrillnew = spherical2cartes(x0, y0, z0, azimuth * np.pi/180., (180 - dip)* np.pi/180., rladder)
		xnew = (xdrillnew/xvoxsize).astype(int)
		ynew = (ydrillnew/yvoxsize).astype(int)
		znew = (-zdrillnew/zvoxsize).astype(int)
		funct = np.sum(drill_rec[xnew,ynew,znew])  + kappa * np.sqrt(np.sum(drill_var[xnew,ynew,znew])) - beta * np.sum(costs[xnew,ynew,znew])
	except:
		funct = 0.
	return -funct

"""

Chose in settings.yaml sepcifications and whether vertical or non-vertical drillcores

The optimisation of the acquisitioan function is performed with a global optimiser: SHGO with sobol sampling
This has the advantage of proposing and sorting multiple local maxima

"""

def bayesopt_vert(drillcoord = None):
	"""
	New Drillcore proposal based on mean, uncertainty, and costs as defined in acquisitian function futility_vertical

	The optimisation of the acquisitioan function is performed with a global optimiser: SHGO with sobol sampling
	This has the advantage of proposing and sorting multiple local maxima

	PARAMETER
	param drillcoord: exisiting x, y, x drill coordinates

	RETURN
	Saves list of proposal of new drillcore coordinates in file newdrill_vertical_proposals.csv

	"""
	print('Calculating propoals list of new vertical drillholes...')
	# Define boundaries for optimisation:
	bp0 = np.asarray([yNcube/2, xNcube/2])
	blb = np.asarray([1, yNcube-1])
	bub = np.asarray([1, xNcube-1])
	#bopt_res = minimize(futility_vertical, bp0, bounds = (blb, bub),  method='SLSQP', options={'maxiter': 500}) #tol =1e-6, method='SLSQP'
	bopt_res = shgo(futility_vertical, bounds = ((1,yNcube - 1), (1,xNcube- 1)), n=20, iters=20, sampling_method='sobol' ) #tol =1e-6, method='SLSQP'
	if not bopt_res.success:
		print('WARNING: ' + bopt_res.message) 
	print('')
	print('New vertical Drillcore Proposal:')
	print('_______________________')
	print('EASTING [meters]: ', np.round(bopt_res.x[1] * xvoxsize) + xmin + 0.5 * xvoxsize)
	print('NORTHING [meters]: ', np.round(bopt_res.x[0] * yvoxsize) + ymin + 0.5 * yvoxsize)
	print(' ')

	# Write new drillcore proposals as csv file
	print('Saving all new drillhole proposals in file: newdrill_proposals_vertical.csv')
	newdf_values = bopt_res.xl 
	newdf_head = ['NORTHING','EASTING']
	df_newdrill = pd.DataFrame(np.round(newdf_values,2), columns = newdf_head)
	df_newdrill['EASTING'] = np.round(df_newdrill['EASTING']) * xvoxsize + xmin + 0.5 * xvoxsize
	df_newdrill['NORTHING'] = np.round(df_newdrill['NORTHING']) * yvoxsize + ymin + 0.5 * yvoxsize
	df_newdrill['BO_GAIN'] = -np.round(bopt_res.funl,4)
	df_newdrill.to_csv(os.path.join(outpath,'newdrill_proposals_vertical.csv'), index=False)

	# Create image of proposed drillcore positions
	densimg = density_rec.mean(axis = 2)
	magsusimg = magsus_rec.mean(axis = 2)
	drillimg = drill_rec.mean(axis = 2)
	drillvarimg = drill_var.mean(axis = 2)
	#extent=[xmin + xvoxsize, xmax - xvoxsize,ymin + yvoxsize,ymax- yvoxsize]
	extent=[xmin, xmax , ymin , ymax]
	plt.clf()
	plt.imshow(drillimg, aspect = 'equal', cmap = 'viridis', extent=extent, origin='lower')
	plt.xlabel('EASTING')
	plt.ylabel('NORTHING')
	if drillcoord is not None:
		xdrill, ydrill, zdrill = drillcoord.T
		plt.scatter(xdrill+xmin, ydrill+ymin, color='k')
	plt.scatter(df_newdrill['EASTING'].values, df_newdrill['NORTHING'].values, color='white')
	plt.scatter(df_newdrill.loc[0,'EASTING'], df_newdrill.loc[0,'NORTHING'], color='red')
	plt.title('Proposed Vertical Drillcores')
	plt.tight_layout()
	plt.savefig(os.path.join(outpath, 'newdrill_vertical_proposals.png'))
	plt.close("all")


def bayesopt_nonvert(drillcoord = None):
	"""
	New Drillcore proposal based on mean, uncertainty, and costs as defined in acquisitian function futility_vertical

	The optimisation of the acquisitioan function is performed with a global optimiser: SHGO with sobol sampling
	This has the advantage of proposing and sorting multiple local maxima

	PARAMETER
	param drillcoord: exisiting x, y, x drill coordinates

	RETURN
	Saves list of proposal of new drillcore coordinates in file newdrill_vertical_proposals.csv
	"""

	bnds = ((yvoxsize,yLcube - yvoxsize), (xvoxsize,xLcube- xvoxsize), (0,360),(30,90))
	bopt_res = shgo(futility_drill, bnds, n=10, iters=500, sampling_method='sobol')
	print('')
	print('New non-vertical Drillcore Proposal:')
	print('_______________________')
	print('EASTING [meters]: ', np.round(bopt_res.x[1] + ymin))
	print('NORTHING [meters]: ', np.round(bopt_res.x[0] + xmin))
	print('Azimuth Angle [degree]: ', np.round(bopt_res.x[2],1))
	print('Dip Angle [degree]: ', np.round(bopt_res.x[3],1))
	print(' ')

	# Write new drillcore proposals as csv file
	print('Saving all new drillhole proposals in file: newdrill_proposals_non-vertical.csv')
	newdf_values = bopt_res.xl
	newdf_head = ['NORTHING', 'EASTING', 'AZIMUTH', 'DIP']
	df_newdrill = pd.DataFrame(np.round(newdf_values,2), columns = newdf_head)
	df_newdrill['EASTING'] = np.round(df_newdrill['EASTING'] + xmin,1)
	df_newdrill['NORTHING'] = np.round(df_newdrill['NORTHING'] + ymin,1)
	df_newdrill['BO_GAIN'] = -np.round(bopt_res.funl,4)
	df_newdrill.to_csv(os.path.join(outpath,'newdrill_proposals_non-vertical.csv'), index=False)

	# Create image of proposed drillcore positions
	densimg = density_rec.mean(axis = 2)
	magsusimg = magsus_rec.mean(axis = 2)
	drillimg = drill_rec.mean(axis = 2)
	drillvarimg = drill_var.mean(axis = 2)
	#extent=[xmin + xvoxsize, xmax - xvoxsize,ymin + yvoxsize,ymax- yvoxsize]
	extent=[xmin, xmax, ymin, ymax]
	plt.clf()
	plt.imshow(drillimg, aspect = 'equal', cmap = 'viridis', extent=extent, origin='lower')
	plt.xlabel('EASTING')
	plt.ylabel('NORTHING')
	if drillcoord is not None:
		xdrill, ydrill, zdrill = drillcoord.T
		plt.scatter(xdrill+xmin, ydrill+ymin, color='k')
	plt.scatter(df_newdrill['EASTING'].values, df_newdrill['NORTHING'].values, color='white')
	plt.scatter(df_newdrill.loc[0,'EASTING'], df_newdrill.loc[0,'NORTHING'], color='red')
	plt.title('Proposed Drillcores')
	plt.tight_layout()
	plt.savefig(outpath + 'newdrill_proposals.png')
	plt.close('all')


def create_costcube():
	"""
	User function to create costcube for drilling. Need to be same cube shape as reconstructed cube.
	By default set to zero costs.

	RETURN
	costcube
	"""
	# Here costs are set to zero, change below
	costcube = np.zeros((xNcube, yNcube, zNcube))
	return costcube


# Run main functions:

if len(sys.argv) < 2:
	fname_settings = input("Please provide settings file name (including path):  ")
	sys.argv.append(fname_settings)

if len(sys.argv) == 2:
	from .config_loader import *  # loads settings
	from .utils import *
	from . import cubeshow as cs 
	from . import inversion
	from . import simcube
	
	if gen_simulation:
		# Create first simulated datacube, see settings.yaml
		simcube.create_simdata(modelname)  

	# Initiate inversion  class
	inv = inversion.Inversion()

	# Create new cube geometry
	voxelpos = inv.create_cubegeometry()
	xxx, yyy, zzz = voxelpos
	xxx = inv.xxx = xxx.reshape(xNcube, yNcube, zNcube)
	yyy = inv.yyy = yyy.reshape(xNcube, yNcube, zNcube)
	zzz = inv.zzz = zzz.reshape(xNcube, yNcube, zNcube)

	# Read in survey data
	gravfield, magfield, sensor_locations = read_surveydata()

	# # Read in existing drillcore data
	drilldata, drillcoord, drillminmax = read_drilldata(drill_features)
	drilldata0 = drilldata[ifeature]
	drillfield = drilldata0[drilldata0 != 0]
	xdrillminmax, ydrillminmax, zdrillminmax = drillminmax

	# Joint Inversion and reconmstyrcuting of cube with geophysical properties:
	density_rec, magsus_rec, drill_rec, density_var, magsus_var, drill_var = inv.cubing(gravfield, magfield, drillfield, sensor_locations, drilldata0)

	### Create VTK cubes or reconstructed cubes:  
	origin = (voxelpos[0].min(), voxelpos[1].min(), voxelpos[2].min())
	voxelsize = (xvoxsize, yvoxsize,zvoxsize)
	cs.create_vtkcube(density_rec, origin, voxelsize, fname = outpath + 'cube_density.vtk')
	cs.create_vtkcube(magsus_rec, origin, voxelsize, fname = outpath + 'cube_magsus.vtk')
	cs.create_vtkcube(drill_rec, origin, voxelsize, fname = outpath + 'cube_drill.vtk')
	cs.create_vtkcube(density_var, origin, voxelsize, fname = outpath + 'cube_density_variance.vtk')
	cs.create_vtkcube(magsus_var, origin, voxelsize, fname = outpath + 'cube_magsus_variance.vtk')
	cs.create_vtkcube(drill_var, origin, voxelsize, fname = outpath + 'cube_drill_variance.vtk')


	# Create plots of 2D maps of vertically integrated cube properties:
	if plot_vertical:	
		densimg = density_rec.mean(axis = 2)
		magsusimg = magsus_rec.mean(axis = 2)
		drillimg = drill_rec.mean(axis = 2)
		extent=[xmin + xvoxsize, xmax - xvoxsize,ymin + yvoxsize,ymax- yvoxsize]
		plt.clf()
		plt.imshow(densimg, aspect = 'equal', cmap = 'viridis', extent=extent, origin='lower')
		plt.colorbar()
		plt.savefig(os.path.join(outpath,'dens_rec2D_loc2.png'))
		plt.clf()
		plt.imshow(magsusimg, aspect = 'equal', cmap = 'viridis', extent=extent, origin='lower')
		plt.colorbar()
		plt.savefig(os.path.join(outpath,'magsus_rec2D_loc2.png'))
		plt.clf()
		plt.imshow(drillimg, aspect = 'equal', cmap = 'viridis', extent=extent, origin='lower')
		plt.colorbar()
		plt.savefig(os.path.join(outpath, 'drill_rec2D_loc2.png'))
		plt.close('all')


	# Create 3D Cube plot if needed (same data as in VTK cube)
	if plot3d:
			plt.clf()
			#cs.skplot3(density_rec, Nsize = (xNcube, yNcube, zNcube), drill = (yyydrill/xvoxsize,xxxdrill/yvoxsize,-zzzdrill/zvoxsize), sensor = (x_sensor/xvoxsize, y_sensor/yvoxsize, x_sensor * 0.), show = False, path_out = outpath, filename = 'density-drill-mesh.png')
			cs.skplot3(density_rec, Nsize = (yNcube, xNcube, zNcube), drill = (ydrillminmax/xvoxsize, xdrillminmax/yvoxsize,-zdrillminmax/zvoxsize), sensor = (sensor_locations[1]/xvoxsize,sensor_locations[0]/yvoxsize, sensor_locations[2] * 0. + zmax), show = False, path_out = outpath, filename = 'density-mesh3D.png')
			plt.clf()
			cs.skplot3(magsus_rec, Nsize = (yNcube, xNcube, zNcube), drill = (ydrillminmax/xvoxsize, xdrillminmax/yvoxsize,-zdrillminmax/zvoxsize), sensor = (sensor_locations[1]/xvoxsize,sensor_locations[0]/yvoxsize, sensor_locations[2] * 0. + zmax), show = False, path_out = outpath, filename = 'magsus-mesh3D.png')
			plt.clf()
			cs.skplot3(drill_rec, Nsize = (yNcube, xNcube, zNcube), drill = (ydrillminmax/xvoxsize, xdrillminmax/yvoxsize,-zdrillminmax/zvoxsize), sensor = (sensor_locations[1]/xvoxsize,sensor_locations[0]/yvoxsize, sensor_locations[2] * 0. + zmax), show = False, path_out = outpath, filename = 'drill-mesh3D.png')
			plt.close('all')


	# Default: no costs, change in function create_costcube to specify costs 
	costcube = create_costcube

	# Propose new drill-core, chose in settings.yaml sepcifications and whether vertical or non-vertical drillcores
	if bayesopt_vertical:
		bayesopt_vert(drillcoord)

	if bayesopt_nonvertical:
		bayesopt_nonvert(drillcoord)






