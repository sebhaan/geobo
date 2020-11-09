"""
Framework for plotting 3D volumetric data (3D cube) as voxel grid.

Includes export of datacube as vtk file

Copyright 2020 Sebastian Haan

This file is part of GeoBO.

GeoBO is free software made available under the AGPL License. 
For details see the LICENSE file.

@author: Sebastian Haan
"""

import os 
import numpy as np
import os
import sys
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LightSource
import pyvista as pv # helper module for the Visualization Toolkit (VTK)
#import plotly.plotly as py
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage import measure



def skplot2(density, Nsize, drill = [], sensor = [], savefig = True, show = True, path_out = '', filename = 'density-drill-mesh.png'):
    """
    Makes volumetric 3D mesh plot at different density levels
    and plots vertical drillholes as well
    :params drill: x and y drillpositions
    """
    if savefig:
        if not os.path.exists(path_out):
            raise ValueError("Path not found: " + path_out)
        else:
            outfile = os.path.join(path_out, filename)
    Nsize = np.asarray(Nsize)
    if len(drill) > 0:
        drill = np.asarray(drill)
        xdrill = drill[0]
        ydrill = drill[1]
    if len(sensor) > 0:
        sensor = np.asarray(sensor)
        x_sensor = sensor[0]
        y_sensor = sensor[1]
        z_sensor = sensor[2]
    colorscheme = cm.viridis
    fstd = density.std()
    fmean = density.mean()
    fmax = density.max()  
    if fmean < 0:
        sfaces = np.linspace(0+0.25*fstd, fmax*0.99, 5)
    else:
        sfaces = np.linspace(fmean+0.25*fstd, fmax*0.99, 5)
    cfaces = colorscheme(sfaces/sfaces.max())
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    # plot density
    sfact = sfaces.min()/sfaces.max()
    cfaces[:,3] = 0.2*sfaces/sfaces.max() # introduced in new library
    for i in range(len(sfaces)):
        if sfaces[i] > 0.:
            verts, faces, normals, values = measure.marching_cubes_lewiner(density, sfaces[i])
            mesh = Poly3DCollection(verts[faces])
            mesh.set_facecolor(cfaces[i])
            #mesh.set_alpha(0.2*sfaces[i]/sfaces.max())
            mesh.set_edgecolor(cfaces[i])
            ax.add_collection3d(mesh)
    # plot drillhole
    if len(drill) > 0:
        for i in range(len(xdrill)):
            ax.plot([xdrill[i],xdrill[i]],[ydrill[i],ydrill[i]],[0,density.shape[2]-1], 'k')
    # plot grav and mag sensors
    if len(sensor) > 0:
        ax.scatter(x_sensor.flatten(), y_sensor.flatten(), z_sensor.flatten(), c='r', marker='o')

    ax.set_xlabel("Y")
    ax.set_ylabel("X")
    ax.set_zlabel("-Z")

    ax.set_xlim(1, Nsize[0])  # a = 6 (times two for 2nd ellipsoid)
    ax.set_ylim(1, Nsize[1])  # b = 10
    ax.set_zlim(0, Nsize[2]-1)  # c = 16
    ax.invert_zaxis()
    ax.set_title("Drillholes: " + str(len(xdrill)))
    m = cm.ScalarMappable(cmap = colorscheme)
    m.set_array(sfaces)
    #cbar = plt.colorbar(m)
    plt.tight_layout()
    if savefig:
        plt.savefig(outfile)
    if show:
        plt.show()


def skplot3(density, Nsize, drill = [], sensor = [], savefig = True, show = True, path_out = '', filename = 'density-drill-mesh.png'):
    """
    Makes columetric 3D mesh plot at different density levels using MArching cubes levinger algorithm
    and plots non-vertical drillholes
    :params drill: x, y, z drillpositions
    """
    if savefig:
        if not os.path.exists(path_out):
            raise ValueError("Path not found: " + path_out)
        else:
            outfile = os.path.join(path_out,filename)
    Nsize = np.asarray(Nsize)
    if len(drill) > 0:
        drill = np.asarray(drill)
        xdrill = drill[0]
        ydrill = drill[1]
        zdrill = drill[2]
    if len(sensor) > 0:
        sensor = np.asarray(sensor)
        x_sensor = sensor[0]
        y_sensor = sensor[1]
        z_sensor = sensor[2]
    colorscheme = cm.viridis
    fstd = density.std()
    fmean = density.mean()
    fmax = np.percentile(density, 99)
    #fmax = density.max() 
    fmin = np.percentile(density, 1)
    if fmean < 0:
        sfaces = np.linspace(0+0.25*fstd, fmax, 5)
    else:
        sfaces = np.linspace(fmin+0.25*fstd, fmax, 5)
    cfaces = colorscheme((sfaces - fmin)/(sfaces - fmin).max())
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    # plot density
    sfact = sfaces.min()/sfaces.max()
    cfaces[:,3] = 0.2*sfaces/sfaces.max() # introduced in new library
    for i in range(len(sfaces)):
        if sfaces[i] > 0.:
            verts, faces, normals, values = measure.marching_cubes_lewiner(density, sfaces[i])
            mesh = Poly3DCollection(verts[faces])
            mesh.set_facecolor(cfaces[i])
            #mesh.set_alpha(0.2*sfaces[i]/sfaces.max())
            mesh.set_edgecolor(cfaces[i])
            ax.add_collection3d(mesh)
    # plot drillhole
    if len(drill) > 0:
        #ax.scatter(xdrill, ydrill, zdrill, c = 'k', marker = 'o')
        for i in range(len(xdrill)):
            ax.plot([xdrill[i, 0],xdrill[i, 1]],[ydrill[i, 0],ydrill[i, 1]],[zdrill[i, 0],zdrill[i, 1]], 'k')
    # plot grav and mag sensors
    if len(sensor) > 0:
        ax.scatter(x_sensor.flatten(), y_sensor.flatten(), z_sensor.flatten(), c='r', marker='o')

    ax.set_xlabel("Y")
    ax.set_ylabel("X")
    ax.set_zlabel("-Z")

    ax.set_xlim(1, Nsize[0])  # a = 6 (times two for 2nd ellipsoid)
    ax.set_ylim(1, Nsize[1])  # b = 10
    ax.set_zlim(0, Nsize[2]-1)  # c = 16
    ax.invert_zaxis()
    ax.set_title("Drillholes: " + str(len(xdrill)))
    m = cm.ScalarMappable(cmap = colorscheme)
    m.set_array(sfaces)
    #cbar = plt.colorbar(m)
    plt.tight_layout()
    if savefig:
        plt.savefig(outfile)
    if show:
        plt.show()


def create_vtkcube(density, origin, voxelsize, fname):
    """
    Export Cube as VTK file (can be used in e.g. ParaView)
    and create a range of 3D cube plots with pyvista
    :param density: 3D cube in shape (xdim, ydim, zdim)
    :param origin: origin cooridnates of cube
    :param voxelsize: voxel sizes in (xsize, ysize, zsize)
    :param fname: path + filename for files
    """
    grid = pv.UniformGrid()
    grid.dimensions = np.array(density.shape) + 1
    grid.origin = origin
    grid.spacing = voxelsize
    grid.cell_arrays["values"] = density.flatten(order="F")
    grid.save(fname)
