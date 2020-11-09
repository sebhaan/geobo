"""
Kernel library for Gaussian Processes including sparse kernels and cross-covariance terms

The choice for an appropriate covariance function is important, 
as the GP's output directly depends on it. These parameters of the covariance function are 
referred to as the hyperparameters of the GP, which can be either given by a fixed covariance scale 
and noise, or learned from data by optimising the marginal likelihood. To handle the computational problem 
of inverting a large covariance matrix, sparse covariance function are included here as well.
To take fully into account cross-covariances between multiple model parameters 
(e.g., rock density and magnetic susceptibility), we construct cross-covariance terms (here labeled with edning "2")
between all kernel pairs. One important requirement for constructing cross-covariance terms is that they must be defined 
to be both positive semi-definite and informative; for an overview how to construct such as matrix in detail see Melkumyan 2011.

Copyright 2020 Sebastian Haan

This file is part of GeoBO.

GeoBO is free software made available under the AGPL License. 
For details see the LICENSE file.

@author: Sebastian Haan
"""
from scipy import reshape, sqrt, identity
import numpy as np


def calcGridPoints3D(Lpix, pixscale):
    """
    returns grid points for distance matrix calculation.
    :param Lpix: number of pixels in each dimension as array (xLpix, yLpix, zLpix)
    :param pixscale: pixelscale in each dimension as array (xpixscale, ypixscale, zpixscale)
    """
    Lpix = np.asarray(Lpix)
    pixscale = np.asarray(pixscale)
    xLpix, yLpix, zLpix = Lpix[0], Lpix[1], Lpix[2]
    xpixscale, ypixscale, zpixscale = pixscale[0], pixscale[1], pixscale[2]
    xrange = np.arange(1, xLpix+1) * xpixscale
    yrange = np.arange(1, yLpix+1) * ypixscale
    zrange = np.arange(1, zLpix+1) * zpixscale
    _xg, _yg, _zg = np.meshgrid(xrange, yrange, zrange)
    xr, yr, zr = _xg.ravel(), _yg.ravel(), _zg.ravel()
    return np.asarray([xr, yr, zr]).T


def calcDistanceMatrix(nDimPoints, 
                       distFunc=lambda deltaPoint: np.sum(deltaPoint[d]**2 for d in range(len(deltaPoint)))):
    """ Returns the matrix of squared distances from one coordinate to any other
    :param nDimPoints: list of n-dim tuples
    :param distFunc: calculates the distance based on the differences
    """
    nDimPoints = np.array(nDimPoints)
    dim = len(nDimPoints[0])
    delta = [None]*dim
    for d in range(dim):
        data = nDimPoints[:,d]
        delta[d] = data - np.reshape(data,(len(data),1)) # computes all possible combinations

    dist = distFunc(delta)
    #dist = dist + np.identity(len(data))*dist.max() # eliminate self matching
    # returns  squared distance:
    return dist 


def calc_square_distances2D(Lpix, pixscale):
    """
    Initialize (squared) distance matrix for stationary kernel.
    """
    Lpix = np.asarray(Lpix)
    pixscale = np.asarray(pixscale)
    xLpix, yLpix = Lpix[0], Lpix[1]
    xpixscale, ypixscale = pixscale[0], pixscale[1]
    xrange = (np.arange(0, xLpix) - xLpix/2.0) * xpixscale
    yrange = (np.arange(0, yLpix) - xLpix/2.0) * ypixscale
    _xg, _yg = np.meshgrid(xrange, yrange)
    xr, yr = _xg.ravel(), _yg.ravel()
    Dx = xr[:, np.newaxis] - xr[np.newaxis,:]
    Dy = yr[:, np.newaxis] - yr[np.newaxis,:]
    return Dx**2 + Dy**2


def gpkernel(D2, gamma):
    """2-D round  RBF kernel, with length scale = standard deviation of
    the PSF of a Gaussian process scene drawn from this kernel.
    Squared Exponential kernel
    :param D2: pairwise square distances
    :param gamma: kernel length scale
    """
    return np.exp(-0.5 * D2/gamma**2)

def gpkernel2(D2, gammas):
    """exp squared x qxp squared kernel
    the PSF of a Gaussian process scene drawn from this kernel.
    Squared Exponential kernel
    :param D2: pairwise square distances
    :param gamma: kernel length scales (gamma1, gamma2)
    """
    l1 = gammas[0]
    l2 =gammas[1]
    return np.sqrt(2.*l1 * l2/(l1**2 + l2**2)) * np.exp(- D2/(l1**2 + l2**2))

def gpkernel_sparse(D2, gamma):
    """2-D round sparse RBF kernel, defined in Melkumyan and Ramos, 2009
    lengthscale is roughly equivlanet to 4 times the lengthcale of squared exponential
    Same as in gpcubesolve_r1.py
    :param D2: pairwise square distances
    :param gamma: kernel length scale
    """
    D2 = np.sqrt(D2)
    #gamma = 4 * gamma
    res = np.zeros_like(D2)
    res[D2 < gamma] = (2 + np.cos(2*np.pi * D2[D2 < gamma] /gamma))/3.*(1-D2[D2 < gamma] /gamma) + 1/(2.*np.pi) * np.sin(2*np.pi*D2[D2 < gamma] /gamma)
    # Remove floating errors
    res[res<0.] = 0.
    return res

def gpkernel_sparse2(D2, gammas):
    """ sparse matrix x sparse matrix
    :param D2: pairwise square distances
    :param gamma: kernel length scales (gamma1, gamma2)
    """
    D2 = np.sqrt(D2)
    l1 = gammas[0]
    l2 = gammas[1]
    # offset for avoiding non-zero terms
    if l1 == l2:
        l2 += 1e-3 * l2
    lmean = np.mean([l1, l2])
    lmin = np.min([l1,l2])
    lmax = np.max([l1,l2])
    #if D2 >= lmean:
    res = np.zeros_like(D2)
    #if D2 <= abs(l2 - l1) / 2.:
    res[D2 <= abs(l2 - l1) / 2.] = 2./(3*np.sqrt(l1*l2)) * (lmin + 1/np.pi*lmax**3/(lmax**2 - lmin**2) * np.sin(np.pi*lmin/lmax*np.cos(2*np.pi*D2[D2 <= abs(l2 - l1) / 2.]/lmax)))
    #if (D2 >= abs(l2 - l1)/2.) & (D2<=(l1 + l2)/2.):
    res[(D2 >= abs(l2 - l1)/2.) & (D2<=(l1 + l2)/2.)] = 2./(3*np.sqrt(l1*l2)) * (lmean - D2[(D2 >= abs(l2 - l1)/2.) & (D2<=(l1 + l2)/2.)] + l1**3*np.sin(np.pi*(l2-2.*D2[(D2 >= abs(l2 - l1)/2.) & (D2<=(l1 + l2)/2.)])/l1)/(2*np.pi*(l1**2-l2**2)) - l2**3*np.sin(np.pi*(l1-2.*D2[(D2 >= abs(l2 - l1)/2.) & (D2<=(l1 + l2)/2.)])/l2)/(2*np.pi*(l1**2-l2**2)))
    # Remove floating errors
    res[res<0.] = 0.
    return res

def gpkernel_matern32(D2, gamma):
    ''' Matern3/2 kernel
    :param D2: pairwise square distances
    :param gamma: kernel length scale
    '''
    nu = np.sqrt(3) * np.sqrt(D2) / gamma    
    return (1 + nu) * np.exp(-nu)

def gpkernel_matern32_2(D2, gammas):
    ''' Matern3/2 x Matern3/2 kernel
    :param D2: pairwise square distances
    :param gammas: kernel length scales (gamma1, gamma2)
    '''
    l1 = gammas[0]
    l2 = gammas[1]
    norm = 2*np.sqrt(l1 * l2)/(l1**2 - l2**2)
    return  norm * (l1 * np.exp(- np.sqrt(3*D2)/l1) - l2 * np.exp(- np.sqrt(3*D2)/l2))

def create_cov(D2, gplength, crossweights = [1,1,1], fkernel = 'sparse'):
    """
    Compute cross-correlation covariance matrix
    :param D2: distance matrix
    :param gplengths: scalar or array with lengthscales for each of the properties (up to three)
    :param crossweights: array with cross-correlation coefficents correspdoning to the output property (w1,w2,w3)
    :param fkernel: kernel function, which can be either 'sparse' (Default), 'exp' (squared exponential) or 'matern32'

    Cross weights (w1,w2,w3) defined as 
    w1: density - drill correltion
    w2: magentic - drill correlation 
    w3: density - magnetic correlation

    Return: cross-correlation matrix
    """
    # first calculate kernel
    params= np.asarray(gplength)
    if params[1] == params[0]:
        params[1] = 1.01 * params[0]
    if params[2] == params[0]:
        params[1] = 1.02 * params[0]
    if params[2] == params[1]:
        params[2] = 1.01 * params[1]
    w1, w2, w3 = np.asarray(crossweights)
    #kcov = np.asarray([gpkernel(D2, params[0]), gpkernel2(D2, params[0:2]), gpkernel2(D2, params[0:2]), gpkernel(D2, params[1])]).T
    if fkernel == 'matern32':
        kcov1 = np.vstack([gpkernel_matern32(D2, params[0]), w3 * gpkernel_matern32_2(D2, params[[0,1]]), w1 * gpkernel_matern32_2(D2, params[[0,2]])])
        kcov2 = np.vstack([w3 * gpkernel_matern32_2(D2, params[[1,0]]), gpkernel_matern32(D2, params[1]), w2 * gpkernel_matern32_2(D2, params[[1,2]])])
        kcov3 = np.vstack([w1* gpkernel_matern32_2(D2, params[[2,0]]), w2 * gpkernel_matern32_2(D2, params[[2,1]]), gpkernel_matern32(D2, params[2])])
    if fkernel == 'sparse':
        kcov1 = np.vstack([gpkernel_sparse(D2, params[0]), w3 * gpkernel_sparse2(D2, params[[0,1]]), w1 * gpkernel_sparse2(D2, params[[0,2]])])
        kcov2 = np.vstack([w3 * gpkernel_sparse2(D2, params[[1,0]]), gpkernel_sparse(D2, params[1]), w2 * gpkernel_sparse2(D2, params[[1,2]])])
        kcov3 = np.vstack([w1 * gpkernel_sparse2(D2, params[[2,0]]), w2 * gpkernel_sparse2(D2, params[[2,1]]), gpkernel_sparse(D2, params[2])])
    if fkernel == 'exp':
        kcov1 = np.vstack([gpkernel(D2, params[0]), w3 * gpkernel2(D2, params[[0,1]]), w1 * gpkernel2(D2, params[[0,2]])])
        kcov2 = np.vstack([w3 * gpkernel2(D2, params[[1,0]]), gpkernel(D2, params[1]), w2 * gpkernel2(D2, params[[1,2]])])
        kcov3 = np.vstack([w1 * gpkernel2(D2, params[[2,0]]), w2 * gpkernel2(D2, params[[2,1]]), gpkernel(D2, params[2])])
    return np.hstack([kcov1, kcov2, kcov3])

