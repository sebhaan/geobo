"""
Utility functions for coordinate conversion and reading ion drillcore and survey data

Copyright 2020 Sebastian Haan

This file is part of GeoBO.

GeoBO is free software made available under the AGPL License. 
For details see the LICENSE file.

@author: Sebastian Haan
"""

import numpy as np
import pandas as pd
import rasterio
import os
from scipy.ndimage.interpolation import zoom
from .config_loader import *

def spherical2cartes(x0, y0, z0, phi, theta, r):
	""" Conversion from spherical coordinates to cartesian

	PARAMETER
	param x0, y0,z0: coordinates of origin
	param phi: azimuthal angle
	param theta: polar angle
	param r: radial length

	RETURN
	x ,y, z coordinates
	"""
	x = x0 + r * np.sin(theta) * np.cos(phi)
	y = y0 + r * np.sin(theta) * np.sin(phi)
	z = z0 + r * np.cos(theta)
	return x, y, z


def cartes2spherial(x0, y0, z0, x1, y1, z1):
	""" Conversion from cartesian coordinates to spehrical

	PARAMETER
	param x0, y0,z0: coordinates of origin
	param x1, y1,z1: coordinates of end point

	RETURN
	radius, polar angle, azimuthal angle
	"""
	r = np.sqrt((x0-x1)**2 + (y0-y1)**2 + (z0-z1)**2)
	theta = np.arccos((z1-z0)/r)
	phi = np.arctan2((y1-y0),(x1-x0))
	return r, theta, phi


def align_drill2(coord, data):
	"""Convert drill-core data in Model Cube shape with iteratinon over all voxels

	PARAMETER

	param coord: xyz Coordinates of drillcore, shape (N_drill, 3)
	param data: 1 dim data for drillcore, shape (N_drill)

	RETURN

	drillcore voxel data cube
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
				if np.size(sel) > 0:
					#print(sel)
					m = np.nanmean(data[sel])
					if np.isfinite(m):
						res[ix,iy,iz] = m
	return res


def normalize(x):
	""" Normalise x
	:param x: input 1D array to normalise
	
	Return
	normalized array
	"""
	if abs(x.max() - x.min()) > 0:
		norm = (x- x.min()) / (x.max() - x.min())  
	else:
		norm = x / x.max()
	return norm


def readcsv_drill(fname, prop_names, pos_names, cartesian = True):
	""" Read csv file with drill-core data and ectract data for density, magsus, and mineral content
	Drill-core data is converted in casrtesian xyz coordinates if only xy position, length and azimuth and dip angle available

	PARAMETER

	:param fname: Path and filename of csv file
	:param prop_names: string array of property names in header of file, e.g ['density', 'magsus', 'mineral']
	:param pos_names: string array of position names in header ['x', 'y', 'z'] 
		or ['East', 'North', 'Elev', DepthFrom', 'DepthTo', 'Azimuth', 'Dip'] with Easting and Northings in meter and Azimuth and Drip in degrees
	:param cartesian: if 'True' (Default), the drill positions are provided in xyz fromat in meters, 
		if 'False' the drill positions are cualcuated from drill depth, azimuth and dip


	RETURN

	returns: pandas array with positions and drill properties
	"""

	names = prop_names + pos_names
	drill = pd.read_csv(fname_drilldata, names = names)
	if cartesian:
		xdrill, ydrill, zdrill = drill[pos_names[0]].values, drill[pos_names[1]].values, drill[pos_names[2]].values
	else:
		xdrill,ydrill,zdrill = spherical2cartes(drill.East.values, drill.North.values, drill.Elev.values, 
			drill.Azimuth.values * np.pi/180., (90 - drill.Dip.values)* np.pi/180., 0.5*(drill.DepthFrom.values + drill.DepthTo.values))
	data = pd.DataFrame()
	data['x'] = xdrill
	data['y'] = ydrill
	data['z'] = zdrill
	data['density'] = drill[prop_names[0]].values
	data['magsus'] = drill[prop_names[1]].values
	data['mineral'] = drill[prop_names[2]].values
	return data


def readraster_survey(fname, pixres = None, clipext = None):
	""" Read geotif rasterfile with one band per property

	PARAMETER
	:param fname: Path and filename of csv file
	:param pixres: desired pixel resolution in meters (assuming same resolution for x and y) 
	:param clipext: provide boundary box [xmin,ymin,xmax,ymax] to crop data to extent


	RETURN
	returns: numpy array with survey data and x,y pixelcoordinates
	"""

	img = rasterio.open(fname)
	bounds = img.bounds #format: [xmin, ymin, xmax, ymax]
	img_shape = img.shape
	# Native pixel resolutiojn:
	ypixres, xpixres = img.res
	# if ypixsize == xpixsize:
	# 	print("WARNING, Input raster pixelsize (resolution) is not the same in vertical and horizontal direction")
	# Read in data into numpy array:
	array = img.read(1)
	# Regrid array to desired 
	if pixres is not None:
		assert ypixres == xpixres, " Raster pixelsize (resolution) is not the same in vertical and horizontal direction"
		array  = zoom(array, ypixsize / pixres)
		xpixres = ypixres = pixres
	# Define pixel coordinates of original image
	xx, yy = np.meshgrid(np.linspace(bounds[0] + 0.5*xpixres, bounds[2] - 0.5*xpixres, array.shape[1]), 
				np.linspace(bounds[1] + 0.5*ypixres, bounds[3] - 0.5*ypixres, array.shape[0]))
	if clipext is not None:
		clipext = np.asarray(clipext)
		# check if clip extent is within boudning box';
		if ((bounds[0]<= clipext[0] <= bounds[2]) &  
			(bounds[1]<= clipext[1] <= bounds[3]) &
			(bounds[0]<= clipext[2] <= bounds[2]) &
			(bounds[1]<= clipext[3] <= bounds[3])):
			# Clip array:
			array = array[((xx - 0.5*xpixres) >= clipext[0]) & ((xx +0.5*xpixres) < clipext[2])
			 & ((yy - 0.5*ypixres) >= clipext[1]) & ((yy + 0.5*ypixres) >= clipext[3])]
			xx2 = xx[((xx - 0.5*xpixres) >= clipext[0]) & ((xx +0.5*xpixres) < clipext[2])
			 & ((yy - 0.5*ypixres) >= clipext[1]) & ((yy + 0.5*ypixres) >= clipext[3])]
			yy2 = yy[((xx - 0.5*xpixres) >= clipext[0]) & ((xx +0.5*xpixres) < clipext[2])
			 & ((yy - 0.5*ypixres) >= clipext[1]) & ((yy + 0.5*ypixres) >= clipext[3])]
			array = array.reshape( np.round((yy2.max()-yy2.min())/ypixres + 1).astype(int), np.round((xx2.max()-xx.min())/xpixres + 1).astype(int))
			xx2 = xx2.reshape( np.round((yy2.max()-yy2.min())/ypixres + 1).astype(int), np.round((xx2.max()-xx.min())/xpixres + 1).astype(int))
			yy2 = yy2.reshape( np.round((yy2.max()-yy2.min())/ypixres + 1).astype(int), np.round((xx2.max()-xx.min())/xpixres + 1).astype(int))
		else:
			print('WARNING: Clip extent exceeds image boundary!... No clipping is applied')
	else: 
		xx2, yy2 = xx, yy
	return array, np.asarray([xx2, yy2])
