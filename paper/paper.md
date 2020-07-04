---
title: 'GeoBO: Python package for Multi-Objective Bayesian Optimisation and Joint Inversion in Geosciences'
tags:
  - Python
  - Geoscience
  - Bayesian Optimisation
  - Inversion
  - Gaussian Processes
authors:
  - name: Sebastian Haan
    orcid: 0000-0002-5994-5637
    affiliation: "1"
affiliations:
  - name: Sydney Informatics Hub, The University of Sydney, Australia
    index: 1
date: 2 July 2020
bibliography: paper.bib
---
<!-- pandoc -V geometry:margin=1in -V fontsize:11pt --filter pandoc-citeproc -o paper.pdf paper.md -->

# Summary

A critical decision process in geophysical data collection is how to efficiently combine a variety of sensor types and where to collect new observations, in particular if measurements are very costly or resources are limited. This may include drill-core placements or the layout of a geophysical survey. Bayesian Optimisation (BO) solves this decision making problem by finding global optimal solutions given multiple, often distinct measurements and model uncertainties. One of the major challenges for applying BO in geoscience still lies in constructing a 3D model from sparse, typically 2D sensor data that can be computationally efficiently evaluated and that transforms the combined data and uncertainties from multiple sensor measurements coherently into a posterior function approximating the true model; this is also know as an inversion problem.

``GeoBO`` is a Python package that solves both, Multi-Objective Optimisation and Joint Inversion, within one probabilistic framework: from prior selection and input, data fusion and inversion, to the final sensor optimisation and real world model output. Fig.1 shows a graphical overview model of ``GeoBO``; the two main process steps are summarized in the following. 

First, ``GeoBO`` jointly solves multi-linear forward models of 2D survey data (e.g. magnetics and gravity) to 3D-geophysical properties using Gaussian Process kernels [see e.g. @Rasmussen:2006; @Melkumyan:2009], without the requirement of prior geological knowledge. In addition, sparse drillhole measurements can be taken into account. One of the key features is to consider correlations  between geophysical properties by solving simultaneously multiple forward models [@Melkumyan:2011; @Reid:2013]. Thus, this approach provides a better constrained combined solution of multiple, distinct sensor types, rather than taking their individual solutions. By using non-parametric Gaussian Process priors, the solution of the inverse problem provides a complete posterior distribution for all geophysical properties which are described by their mean and variance value at each location (voxel). While this inversion method is limited to linear forward models, the reconstructed posterior distribution of the geophysical properties can be further refined using site-specific knowledge or more complex inversion methods for sampling non-linear forward models, such as hIPPYlib [@Hippylib], Obsidian [@Obsidian], or GemPy [@GemPy].

In a second step, ``GeoBO`` allows the user to define objectives in an acquisition function [@Brochu:2010] which typically has to balance between a) exploration, i.e., querying points that maximise the information gain and minimize the uncertainty of a geophysical site model, b) exploitation, i.e. querying points that maximise the reward (e.g. concentrating search in the vicinity locations with high value such as minerals), and c) minimize the number of samples given a cost function for any new measurement (e.g. costly drill-cores). The reconstructed 3D model output of the joint inversion model is then used to query for the next most promising point based on the aquisitian function, which guides the search for a user-defined optimum.

``GeoBO`` allows applications to specify priors as additional options, such as a) the typical correlation lengthscale of the geological structure via GP hyperparameters, b) the gain/cost parameters in the BO acquisition function, and c) the correlation coefficients between different geophysical properties. This package includes forward models for gravity and magnetics surveys [@Li:1998], drillcores, and test-sets of geological models and data. 

While it has been tested on solving the problem of allocating iteratively new expensive drill-core samples [@Haan:2020], the same method can be applied to a wide range of optimisation and sensor fusion problems. 


![Overview of the probabilistic framework for multi-objective optimisation and joint inversion; the process stages are: 1) Joint inversion via Gaussian Processes based on multiple sensor data fusion and forward model,  2) Posterior model generation of 3D multi-geophysical properties. 3) Maximisation of acquisition function to allocate optimal new sensor location and type. 4)  New acquired data is combined with existing data; process repeats until maximum number of iterations is achieved. For more details see @Haan:2020](graphmodel2.png)

## Installation, Dependencies and Usage
``GeoBO`` can be installed using setuptools or pip and requires the following python packages: numpy, matplotlib, scikit_image, scipy, rasterio, pandas, pyvista, skimage, PyYAML(for details see requirements.txt). The installation can be tested by running the script with included synthetic data and default settings.

To use ``GeoBO``, change the main settings such as filenames and parameters in settings.yaml (see description there), and run 
```sh
cd geobo/
python main.py settings.yaml
```

Core features of ``GeoBO``:

 - Joint probabilistic inversion tool by solving simultaneously multi-linear forward models (e.g. gravity, magnetics) using cross-variances between geophysical properties (cross-variance terms can be specified by use in settings.yaml)
 - Computation of complete posterior distribution for all geophysical properties (described by their mean and variance value at each location/voxel) 
 - Generation of ranked proposal list for new most promising drillcores (optional: non-vertical drillcores) based on global optimisation of acquisition function (see settings.yaml)
 - Templates for acquisition function to use in Bayesian Optimisation (see 'futility' functions in run_geobo.py)
 - Flexible parameter settings (settings.yaml) for exploration-exploitation trade-off and inclusion of local 3D cost function in acquisition function (see function 'create_costcube' in run_geobo.py)

Other features are:

 - Package includes geological survey/drillcore sample as well as synthetic data and functions for synthetic data generation
 - Generation of 2D/3D visualisation plots of reconstructed cubes and survey data (see plot settings in settings.yaml)
 - 3D Cube export in VTK format (for subsequent analysis, e.g., with ParaView)
 - Options to include any pre-existing drillcore data (see settings.yaml for feature names as well as format of sample drilldata)
 - Library of Gaussian Process (GP) kernels including sparse GP kernels (kernels.py)
 - Flexible settings for any cube geometry and resolution (cube length/width/height and voxel resolutions)
 - Optional  marginal GP likelihood for  optimization of GP hyperparameters and inversion process



# Acknowledgements
Sebastian Haan would like to acknowledge Dietmar Muller, Fabio Ramos, and Ben Mather for their very valuable contributions to this project. This research was supported by the Sydney Informatics Hub, a Core Research Facility of the University of Sydney.


# References
