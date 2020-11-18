#GeoBO: A Python package for Multi-Objective Bayesian Optimisation and Joint Inversion in Geosciences
#
#Copyright 2020 Sebastian Haan
#
#This file is part of GeoBO.
#
#GeoBO is free software made available under the AGPL License. 
#For details see the LICENSE file.
#
#@author: Sebastian Haan

"""
GeoBO is build upon a probabilistic framework using Gaussian Process (GP) priors to jointly solve multi-linear forward models. 
This software generates multi-output 3D cubes of geophysical properties (e.g. density, magnetic susceptibility, mineral concentrations) and their uncertainties from 2D survey data (e.g. magnetics and gravity) and any pre-existing drillcore measurements. 
The reconstructed 3D model is then used to query the next most promising measurement location given an expensive cost function (e.g. for drillcores). 
A ranked list of new measurements is proposed based on user-defined objectives as defined in the acquisitian function which typically aims to optimize exploration (reducing global model uncertainty) and exploitation (focusing on highly promising regions) while minimizing costs. 

## Functionality

GeoBO's probablistic framework includes all steps from  prior selection, data fusion and inversion, to sensor optimisation and real world model output. The main functionalities of GeoBO are summarised in the following:
 - Joint probabilistic inversion tool by solving simultanously multi-linear forward models (e.g. gravity, magnetics) using cross-variances between geophysical properties (cross-variance terms can be specified by user)
 - Output 1: Generation of cubes and computation of complete posterior distribution for all geophysical properties (described by their mean and variance value at each location (cubecell aka voxel). 
 - Output 2: Generation of ranked proposal list for new most promising drillcores based on global optimisation of acquisitian function
 - Templates for acquisitian function to use in Bayesian Optimisation
 - Flexible parameter settings for exploration-exploitation trade-off and inclusion of local 3D cost function in acquisitian function 


Other features are:
 - Generation of simulated geophysical data with a choice of three different models
 - Package includes geological survey/drillcore sample as well as synthetic data and functions for synthetic data generation
 - Generation of 2D/3D visualisation plots of reconstructed cubes and survey data
 - 3D Cube export in VTK format (for subsequent analysis, e.g., in Python or with ParaView)
 - Options to include any pre-existing drillcore data 
 - Library of Gaussian Process (GP) kernels including sparse GP kernels
 - Flexible settings for any cube geometry and resolution
 - (Optional) Optimization of GP hyperparameters and cross-correlation coeffcients via computation of marginal GP likelihood

Example outputs can be found in the directory `examples/results/`.

## Installation And Requirements


### Installation

To install GeoBO locally using setuptools: 

python setup.py build
python setup.py install

or using pip:

pip3 install geobo


### Requirements
- python >=3.6
- numpy
- matplotlib
- scikit_image
- scipy
- rasterio
- pandas
- pyvista
- skimage
- PyYAML


## Usage and Settings

1) Change the main settings such as filenames and parameters in the settings file settings.yaml. These settings specify:

2) Then run geobo
python run_main.py settings.yaml
"""
__version__ = "0.1.1"
__author__ = "Sebastian Haan"