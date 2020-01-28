GEOBO
==========================================

Python package for Multi-Objective Bayesian Optimisation and Joint Inversion in Geosciences

``GEOBO`` is build upon a probabilistic framework using Gaussian Process (GP) priors to jointly solve multi-linear forward models. This software generates multi-output 3D cubes of geophysical properties (e.g. density, magnetic susceptibility, mineral concentrations) and their uncertainties from 2D survey data (e.g. magnetics and gravity) and any pre-existing drillcore measurements. The reconstructed 3D model is then used to query the next most promising measurement location given an expensive cost function (e.g. for drillcores). A ranked list of new measurements is proposed based on user-defined objectives as defined in the acquisitian function which typically aims to optimize exploration (reducing global model uncertainty) and exploitation (focusing on highly promising regions) while minimizing costs.

The probablistic framework includes all steps from  prior selection, data fusion and inversion, to sensor optimisation and real world model output.

Core features are:

 - Joint probabilistic inversion tool by solving simultanously multi-linear forward models (e.g. gravity, magnetics) using cross-variances between geophysical properties (cross-variance terms can be specified by user)
 - Computation of complete posterior distribution for all geophysical properties (described by their mean and variance value at each location/voxel) 
 - Generation of ranked proposal list for new most promising drillcores based on global optimisation of acquisitian function
 - Templates for acquisitian function to use in Bayesian Optimisation
 - Flexible parameter settings for exploration-exploitation trade-off and inclusion of local 3D cost function in acquisitian function 


Other features are:
 - Package includes geological survey/drillcore sample as well as synthetic data and functions for synthetic data generation
 - Generation of 2D/3D visualisation plots of reconstructed cubes and survey data
 - 3D Cube export in VTK format (forsubsequent analysis, e.g., with ParaView)
 - Options to include any pre-existing drillcore data 
 - Library of Gaussian Process (GP) kernels including sparse GP kernels
 - Flexible settings for any cube geometry and resolution
 - Optional output: marginal GP likelihood for subsequent optimization of GP inversion process (not automated here)


INSTALLATION AND REQUIREMENTS
-----------------------------

Requires Python 3, see requirements.txt for additional Python libraries

Test installation by running script with included synthetic data.



GETTING STARTED 
---------------

1) Change the main settings such as filenames and parameters in settings.yaml

2) python run_geobo.py 


The main functions for acquisitian function  can be found in run_geobo.py
Visualisation functions and VTK export are defined in cubeshow.py
Inversion functions are defined in inversion.py 


