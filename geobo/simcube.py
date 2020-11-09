""" 
Script for generating simulated cubes and sensors

Set of multiple distinct synthetic geophysical models can be created, 
including two-dipping body ("cylinders") and  layered models. 
For each model  a 3D voxel cube with geological structures is generated 
given by their density and magnetic susceptibility 
as well as the corresponding 2D gravity and magnetic remote sensor measurements. 

Other custom models can be included by adding a new model in function create_syncube()

Copyright 2020 Sebastian Haan

This file is part of GeoBO.

GeoBO is free software made available under the AGPL License. 
For details see the LICENSE file.

@author: Sebastian Haan
"""
import os
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import random
import rasterio
from .config_loader import *
from . import inversion
from .sensormodel import * 
from . import cubeshow as cs


def create_syncube(modelname, voxelpos):
	"""Creates synthetic cube for density and magnetic susceptibility
	
	Generates two output files, one vtk cube and one csv file

	PARAMETER

	param modelname: String, options: "layers_2", "layers_3", "cylinders" 
	param voxelpos: voxel positions (x, y, z)

	RETURN

	density cube
	susceptibility cube 
	"""
	print("Creating simulated cube data ...")
	xxx, yyy, zzz = voxelpos
	x3 = xxx.reshape(yNcube, xNcube, zNcube)
	y3 = yyy.reshape(yNcube, xNcube, zNcube)
	z3 = zzz.reshape(yNcube, xNcube, zNcube)
	if modelname == 'layers_2':
		zshift = zLcube/8. * 1. / (1 + np.exp(2.*(-y3 + zLcube/2)))
		layer1 = 4. * (1./(1+np.exp(-2*(-z3 - zLcube * 0.3 + zshift) )) - 1./(1+np.exp(-2*(-z3 - zLcube * 0.325 + zshift))))
		cut1 = np.percentile(layer1,90)
		layer1[layer1 < cut1] = 0.
		layer1[layer1 >= cut1] = layer1.max()
		layer2 = 8. * (1./(1+np.exp(-2*(-z3 - zLcube * 0.25 + zshift) )) - 1./(1+np.exp(-2*(-z3 - zLcube * 0.275 + zshift))))
		cut2 = np.percentile(layer2,90)
		layer2[layer2 < cut2] = 0.
		layer2[layer2 >= cut2] = layer2.max()
		density = 0.5 + layer1 + layer2 #
		# assume simple correlation between magnetism and density:
		magsus = gp_coeff[1] * density
	if modelname == 'layers_3':
		zshift = zLcube/8. * 1. / (1 + np.exp(2.*(-y3 + yLcube/2.)))
		layer3 = 6. * (1./(1+np.exp(-2*(-z3 - zLcube * 0.35 + zshift) )) - 1./(1+np.exp(-2*(-z3 - zLcube * 0.375 + zshift))))
		cut3 = np.percentile(layer3,90)
		layer3[layer3 < cut3] = 0.
		layer3[layer3 >= cut3] = layer3.max()
		layer1 = 4. * (1./(1+np.exp(-2*(-z3 - zLcube * 0.3 + zshift) )) - 1./(1+np.exp(-2*(-z3 - zLcube * 0.325 + zshift))))
		cut1 = np.percentile(layer1,90)
		layer1[layer1 < cut1] = 0.
		layer1[layer1 >= cut1] = layer1.max()
		layer2 = 8. * (1./(1+np.exp(-2*(-z3 - zLcube * 0.25 + zshift) )) - 1./(1+np.exp(-2*(-z3 - zLcube * 0.275 + zshift))))
		cut2 = np.percentile(layer2,90)
		layer2[layer2 < cut2] = 0.
		layer2[layer2 >= cut2] = layer2.max()
		density = 0.5 + layer1 + layer2 + layer3
		magsus = gp_coeff[1] * density
	if modelname == 'cylinders':
		rad = yLcube/18.
		rc1 = ((y3-yLcube/1.3 - rad)**2 ) + ((z3 + zLcube/(4)  - rad)**2)
		rc2 = ((y3-yLcube/4. - rad)**2 ) + ((z3 + zLcube/(4)  - rad)**2)
		density = x3 * 0. + 0.1
		density[rc2 <= rad**2] = 1.
		density[rc1 <= rad**2] = 1.
		density[(x3<xLcube/5.) | (x3>xLcube * 4./5.)] = 0.1
		#magnetic = density
		magsus = gp_coeff[1] * density

	# Create simulated VTK cube
	origin = (voxelpos[0].min(), voxelpos[1].min(), voxelpos[2].min())
	voxelsize = (xvoxsize, yvoxsize,zvoxsize)
	cs.create_vtkcube(density, origin, voxelsize, fname = os.path.join(inpath,'simcube_' + modelname + '.vtk'))

	# Save simulated data as csv file:
	newdf_head = ['x','y','z', 'DENSITY', 'MAGSUS']
	data = np.asarray([x3.flatten(), y3.flatten(), z3.flatten(), density.flatten(), magsus.flatten()])
	df= pd.DataFrame(data.T, columns = newdf_head)
	df.to_csv(os.path.join(inpath,'simcube_' + modelname + '.csv'), index = False)

	# Create simulated drill data from 4 drillcores:
	#select four random x,y coordinates
	dfdrill = df.copy()
	xdrill = [random.randint(2,xNcube -2) for p in range(2)] 
	ydrill = [random.randrange(2,yNcube -2) for p in range(2)]
	xdrill = np.asarray(xdrill) * xvoxsize + 0.5 * xvoxsize
	ydrill = np.asarray(ydrill) * yvoxsize+ 0.5 * yvoxsize
	dfdrill = dfdrill.loc[(dfdrill['x'].isin(xdrill)) & (dfdrill['y'].isin(ydrill))]
	dfdrill['SiteID'] = 'SiteID_' + dfdrill['x'].astype(str) + dfdrill['y'].astype(str)
	dfdrill.to_csv(os.path.join(inpath,'simdrill_' + modelname + '.csv'), index = False)

	return density, magsus


def create_synsurvey(modelname, density, magsus):
	"""
	Generates synthetic gravity and magntics survey data. Sensors are positiones on top of cube data

	PARAMETER

	param density: density cube data
	param magsus: magnetic susceptibiity cube data

	RETURN

	gravity 2D array
	magnetic 2D array
	"""
	print("Creating simulated sensor data...")
	# Create voxel edges
	xedge = np.linspace(0, xNcube, xNcube + 1) * xvoxsize
	yedge = np.linspace(0, yNcube, yNcube + 1) * yvoxsize
	zedge = np.linspace(0, -zNcube, zNcube + 1) * zvoxsize + zmax
	xEdges, yEdges, zEdges = np.meshgrid(xedge, yedge, zedge)
	Edges = np.asarray([xEdges, yEdges, -zEdges])
	# Create sensor positions
	xnew = np.arange(xvoxsize/2., xLcube + xvoxsize/2., xvoxsize)
	ynew = np.arange(yvoxsize/2., yLcube + yvoxsize/2., yvoxsize)
	xx ,yy = np.meshgrid(xnew,ynew)
	zz = xx * 0. + zoff
	sensor_locations = np.asarray([xx.flatten(), yy.flatten(), zz.flatten()]).T
	# Calculate sensitivities
	gravsens, _ = A_sens(magneticField * 0., sensor_locations, Edges, 'grav')
	magsens, _ = A_sens(magneticField, sensor_locations, Edges, 'magn')
	gravfield = np.dot(gravsens, density.flatten()) # shape Nsensor, equivalent to (gravsens * properties).sum(axis = 1)
	magfield =  np.dot(magsens, magsus.flatten())
	grav2D = gravfield.reshape(xx.shape)
	magn2D = magfield.reshape(xx.shape)
	# Write csv file
	newdf_head = ['X','Y', 'GRAVITY', 'MAGNETIC']
	data = np.asarray([xx.flatten(), yy.flatten(), gravfield, magfield])
	df= pd.DataFrame(data.T, columns = newdf_head)
	df.to_csv(os.path.join(inpath,'simsurveydata_' + modelname + '.csv'), index = False)

	return grav2D, magn2D


def create_simdata(modelname = "cylinders", plot = True):
	""" 
	Generates two simulated 3D cubes (density and magnetic susceptibility) 
	plus their corresponding gravity and magnetics 2D sensor data above surface.

	Sensor data is calculated via forward models as specified in sensormodel.py.
	Other settings such as cube geometry, earth's magnetic field, and sensor height 
	above ground are specified in settings. Sensor x,y positions are by default centered 
	as grid on top of voxel centers.

	PARAMETER
	param modelname: String ["layers_2", "layers_3", "cylinders"], defaults to "cylinders" 
	param plot: boolean, if True, 2D plots of sensor and simulated data are created

	RETURN
	The simulated data is saved in output directory (settings) as csv files (one for cube and one for sensor data) 
	The Two cubes are aslo saved in addition as VTK format files.
	Gravity and magentic survey data as tif file.
	Plots of sensor and vertically integrated cube data are saved in output directory.
	"""
	# Check that input path exists
	os.makedirs(inpath, exist_ok=True)

	inv = inversion.Inversion()

	voxelpos = inv.create_cubegeometry()

	# Create cubes
	density, magsus = create_syncube(modelname, voxelpos)

	# Create sensor data
	grav2D, magn2D = create_synsurvey(modelname, density, magsus)

	# Create tif survey files:
	with rasterio.open(os.path.join(inpath,'gravity_simdata_' + modelname +'.tif'), 'w',driver = 'GTiff', width = grav2D.shape[1], height=grav2D.shape[0], count=1, dtype='float32') as dst: 
		dst.write(grav2D.astype(rasterio.float32), 1) 
	with rasterio.open(os.path.join(inpath,'magnetic_simdata_' + modelname +'.tif'), 'w',driver = 'GTiff', width = magn2D.shape[1], height=magn2D.shape[0], count=1, dtype='float32') as dst: 
		dst.write(magn2D.astype(rasterio.float32), 1) 

	# Plot sensor data and integrated cube data along vertical axis
	if plot:
		print("Creating plots of simulated data ...")
		extent = np.asarray([0, xLcube, 0, yLcube])
		fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(10, 8))
		axs[0,0].imshow(grav2D, extent = extent)
		axs[0,0].set_title('Gravity Measurements')
		axs[0,0].grid(True)
		axs[0,1].imshow(magn2D, extent = extent)
		axs[0,1].set_title('Magnetic Meauerements')
		axs[0,1].grid(True)
		axs[1,0].imshow(np.sum(density, axis =2), extent = extent)
		axs[1,0].set_title('Vertical Sum Density')
		axs[1,0].grid(True)
		axs[1,1].imshow(np.sum(magsus, axis =2), extent = extent)
		axs[1,1].set_title('Vertical Sum Magnetic Susceptibility')
		axs[1,1].grid(True)
		plt.tight_layout()
		plt.savefig(os.path.join(inpath,'figure_simdata_' + modelname +'.png'), dpi = 300)
		plt.close()
