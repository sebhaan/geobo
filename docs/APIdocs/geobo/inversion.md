::: {role="main"}
::: {#section-intro .section}
Script for running inversion and reconstructing 3D cubes from 2D sensor
data using Gaussina processes

Author: Sebastian Haan

Expand source code

    """
    Script for running inversion and reconstructing 3D cubes from 2D sensor data using Gaussina processes

    Author: Sebastian Haan
    """

    import numpy as np
    import sys
    from scipy.linalg import pinv, solve, cholesky, solve_triangular
    from config_loader import *
    import kernels as kernel #for Gaussian Process prior
    import sensormodel

    class Inversion:
            """
            Class for Inversion and reconstructon of 3D cubes from 2D sensor data.

            To solve the inverse problem for linear systems, we deploy a Bayesian framework with Gaussian process priors 
            over the physical properties of interest. Gaussian processes (GP) are a flexible, probabilistic approach 
            using kernel machines with non-parametric priors (see e.g. Rasmussen 2006).  An important advantage of 
            the Bayesian method is that it generates a predictive distribution with a mean and variance, 
            which are prerequisites for Bayesian optimisation with respect to, e.g., information gain from a new measurement. 
            Another advantage is that the GP marginal likelihood function is well defined by the values of their hyper-parameters, 
            which allows it to optimise them exactly. This reasoning about functions under uncertainty and their well-tuned 
            interpolation character allows GPs to work extremely well for sparse data as well as for big data 
            (using sparse covariance matrix, see e.g. Melkumyan 2000).

            To take fully into account cross-covariances between multiple model parameters 
            (e.g., rock density and magnetic susceptibility), we construct sparse cross-covariance terms 
            between all kernel pairs following Melkumyan 2011. The full covariance matrix is then constructed 
            by setting sparse covariance functions on all diagonal block elements and sparse-sparse cross-covariance functions 
            on all non-diagonal block elements. Cross-covariance amplitude terms are given by the linear correlation coefficients 
            between the corresponding geophysical properties.

            Call function cubing() for running inversion.
            """
            def __init__(self):
                    # Set Lengthscale and amplitude for GP
                    self.gp_length = gp_lengthscale * np.asarray([xvoxsize, xvoxsize, xvoxsize]) # 2*voxelsize seems to have max logl
                    self.gp_sigma = np.asarray(gp_err)
                    self.coeffm = np.asarray(gp_coeff) # coefficients for GP kernel mix


            def create_cubegeometry(self):
                    """
                    Create Cube geometry and coordinates

                    See settings.yaml for cube specifications
                    """
                    # Create voxel Edge coordinates:
                    xedge = np.linspace(0, xNcube, xNcube + 1) * xvoxsize
                    yedge = np.linspace(0, yNcube, yNcube + 1) * yvoxsize
                    zedge = np.linspace(0, -zNcube, zNcube + 1) * zvoxsize + zmax
                    #zedge = np.linspace(0.,Ncube/zfactor, Ncube/zfactor+1) * voxsize
                    xEdges, yEdges, zEdges = np.meshgrid(xedge, yedge, zedge)
                    self.Edges = np.asarray([xEdges, yEdges, -zEdges])
                    # Create voxel coordinates
                    xnew = np.arange(xvoxsize/2., xLcube + xvoxsize/2., xvoxsize)
                    ynew = np.arange(yvoxsize/2., yLcube + yvoxsize/2., yvoxsize)
                    znew = zmax - np.arange(zvoxsize/2., zLcube + zvoxsize/2., zvoxsize)
                    xx ,yy = np.meshgrid(xnew,ynew)
                    self.xxx, self.yyy, self.zzz = np.meshgrid(xnew, ynew, znew)
                    self.voxelpos = np.vstack([self.xxx.flatten(), self.yyy.flatten(), self.zzz.flatten()])
                    return self.voxelpos


            def predict3(self, calclogl = False):
                    """
                    Calculate mean, variance, and likelihood likelihood for GP with 3x3 kernel block matrix 
                
                PARAMETERS

                :param calclogl: if True calculate marginal GP loglikelihood, if False logl return is set to 0

                RETURN: 
                mean
                covariance
                log-likelihood
                    """
                    # Calculate kernel 3x3 block matrix:
                    self.datastd = np.mean([np.nanstd(self.gravfield), np.nanstd(self.magfield), np.nanstd(self.drillfield)])
                    self.kcov = kernel.create_cov(self.D2, self.gp_length, self.coeffm, fkernel = kernelfunc)
                    #Ak = np.dot(self.Asens3, kcov) # shape (3*Nsensor, 3*Ncube^3)
                    yerr = np.hstack((self.gravfield * 0. + self.gp_sigma[0], self.magfield * 0. + self.gp_sigma[1], self.drillfield *0. + self.gp_sigma[2]))
                    #try:
                    AkA = np.dot(self.Asens3, np.dot(self.kcov, self.Asens3.T)) + np.diag(yerr**2)
                    #AkA = np.dot(Ak, self.Asens3.T) + np.diag(yerr**2) # shape(2*Nsensor, 2*Nsensor)
                    # Cholesky decomposition
                    try:
                            AkA_chol = cholesky(AkA, lower=True)
                    except:
                            print("Cholesky decompostion failed, AkA matrix i likely not positive semitive.")
                            print("Change GP parameter settings")
                            sys.exit(1)
                    usolve = solve_triangular(AkA_chol, self.Fs3, lower=True) #shape(2*Nsensor)
                    # Calculate Likelihood if option is set
                    if calclogl:
                            log_det_AkA = np.log(np.diag(AkA_chol)**2).sum()
                            n_log_2pi = xNcube * yNcube * zNcube * np.log(2 * np.pi)
                            logl = -0.5 * (np.dot(usolve, usolve) + log_det_AkA + n_log_2pi)
                    else:
                            logl = 0.
                    # predicted mean
                    self._V = solve_triangular(AkA_chol, np.dot(self.Asens3, self.kcov), lower=True) #shape(2*Nsensor, 2* Ncube^3)
                    mu = np.dot(self._V.T, usolve) #* _fstd + _fmean
                    # covariance
                    covar = self.kcov - np.dot(self._V.T, self._V)
                    #var = np.diagonal(covar) #* _fstd**2
                    #except:
                    #       print("Warning: GP Inversion failed!")
                    #       mu, covar, logl = np.zeros(3*xNcube*yNcube*zNcube), kcov * 1e9, np.inf 
                    return mu, covar, logl

            
            def cubing(self, gravfield, magfield, drillfield, sensor_locations, drilldata0):
                    """
                    Joint Inversion and Cubing of sensor data

                    PARAMETERS

                    :param gravfield: 2D gravitational survey data
                    :param magfield: 2D magnetic survey data
                    :param drillfield: drilldata

                    RETURN

                    density_rec: cube with mean density 
                    magsus_rec: cube with mean magnetic susceptibility
                    drill_rec: cube with mean drill-core property  
                    density_var: varicance cube with density 
                    magsus_var: variance cube with magnetic susceptibility  
                    drill_var: varianbce cube with drill-core property  


                    """
                    self.gravfield = gravfield
                    self.magfield = magfield
                    self.drillfield = drillfield
                    self.sensor_locations = sensor_locations
                    self.drilldata0 = drilldata0
                    # Normalize data
                    grav_mean, grav_std = self.gravfield.mean(), self.gravfield.std()
                    gravfield_norm = (self.gravfield - grav_mean) / grav_std
                    magn_mean, magn_std = self.magfield.mean(), self.magfield.std()
                    magfield_norm = (self.magfield - magn_mean) / magn_std
                    drill_mean, drill_std = self.drillfield.mean(), self.drillfield.std()
                    drillfield_norm = (self.drillfield - drill_mean) / drill_std
                    # Create kernel distance matrix
                    self.points3D = kernel.calcGridPoints3D((xNcube, yNcube, zNcube), (xvoxsize, yvoxsize, zvoxsize))
                    self.D2 = kernel.calcDistanceMatrix(self.points3D)
                    # Combine Data and Forward models
                    voxelpos_drill = np.vstack([self.xxx[self.drilldata0 != 0], self.yyy[self.drilldata0 != 0], self.zzz[self.drilldata0 != 0]])
                    # Combine data:
                    self.Fs3 = np.hstack((gravfield_norm, magfield_norm, drillfield_norm))
                    # Calculate Sensitivity matrix
                    Asens_grav, _ = A_sens(magneticField * 0., self.sensor_locations, self.Edges, 'grav')
                    Asens_mag, _ = A_sens(magneticField, self.sensor_locations, self.Edges, 'magn')
                    Asens_drill = A_drill(voxelpos_drill, self.voxelpos)
                    # Calculate total transformation matrix
                    Asens_g= np.vstack([Asens_grav, np.zeros_like(Asens_mag), np.zeros_like(Asens_drill)])
                    Asens_m= np.vstack([np.zeros_like(Asens_grav), Asens_mag, np.zeros_like(Asens_drill)])
                    Asens_d= np.vstack([np.zeros_like(Asens_grav), np.zeros_like(Asens_mag), Asens_drill])
                    self.Asens3 = np.hstack([Asens_g, Asens_m, Asens_d])
                    # Run actual GP inversion
                    self.mu_rec, self.cov_rec, self.logl = self.predict3(calclogl = True)
                    # Reshape results and multiply with stadnard deviatio 
                    results_rec = self.mu_rec.reshape(3, yNcube, xNcube, zNcube)
                    results_var = np.diag(self.cov_rec).reshape(3, yNcube, xNcube, zNcube)
                    #Cut off outer voxel rows
                    #results_rec = results_rec[:, 1:yNcube-1, 1:xNcube-1, :]
                    #results_var = results_var[:, 1:yNcube-1, 1:xNcube-1, :]
                    density_rec = results_rec[0] * grav_std  # Model represents deviation from the mean
                    density_var = results_var[0] * grav_std**2
                    magsus_rec = results_rec[1] * magn_std # Model represents deviation from the mean
                    magsus_var = results_var[1] * magn_std**2
                    drill_rec = results_rec[2] * drill_std  # Model represents deviation from the mean
                    drill_var = results_var[2] * drill_std**2
                    return density_rec, magsus_rec, drill_rec, density_var, magsus_var, drill_var
:::

::: {.section}
:::

::: {.section}
:::

::: {.section}
:::

::: {.section}
Classes {#header-classes .section-title}
-------

` class Inversion`{.flex .name .class}

:   ::: {.section .desc}
    Class for Inversion and reconstructon of 3D cubes from 2D sensor
    data.

    To solve the inverse problem for linear systems, we deploy a
    Bayesian framework with Gaussian process priors over the physical
    properties of interest. Gaussian processes (GP) are a flexible,
    probabilistic approach using kernel machines with non-parametric
    priors (see e.g. Rasmussen 2006). An important advantage of the
    Bayesian method is that it generates a predictive distribution with
    a mean and variance, which are prerequisites for Bayesian
    optimisation with respect to, e.g., information gain from a new
    measurement. Another advantage is that the GP marginal likelihood
    function is well defined by the values of their hyper-parameters,
    which allows it to optimise them exactly. This reasoning about
    functions under uncertainty and their well-tuned interpolation
    character allows GPs to work extremely well for sparse data as well
    as for big data (using sparse covariance matrix, see e.g. Melkumyan
    2000).

    To take fully into account cross-covariances between multiple model
    parameters (e.g., rock density and magnetic susceptibility), we
    construct sparse cross-covariance terms between all kernel pairs
    following Melkumyan 2011. The full covariance matrix is then
    constructed by setting sparse covariance functions on all diagonal
    block elements and sparse-sparse cross-covariance functions on all
    non-diagonal block elements. Cross-covariance amplitude terms are
    given by the linear correlation coefficients between the
    corresponding geophysical properties.

    Call function cubing() for running inversion.
    :::

    Expand source code

        class Inversion:
                """
                Class for Inversion and reconstructon of 3D cubes from 2D sensor data.

                To solve the inverse problem for linear systems, we deploy a Bayesian framework with Gaussian process priors 
                over the physical properties of interest. Gaussian processes (GP) are a flexible, probabilistic approach 
                using kernel machines with non-parametric priors (see e.g. Rasmussen 2006).  An important advantage of 
                the Bayesian method is that it generates a predictive distribution with a mean and variance, 
                which are prerequisites for Bayesian optimisation with respect to, e.g., information gain from a new measurement. 
                Another advantage is that the GP marginal likelihood function is well defined by the values of their hyper-parameters, 
                which allows it to optimise them exactly. This reasoning about functions under uncertainty and their well-tuned 
                interpolation character allows GPs to work extremely well for sparse data as well as for big data 
                (using sparse covariance matrix, see e.g. Melkumyan 2000).

                To take fully into account cross-covariances between multiple model parameters 
                (e.g., rock density and magnetic susceptibility), we construct sparse cross-covariance terms 
                between all kernel pairs following Melkumyan 2011. The full covariance matrix is then constructed 
                by setting sparse covariance functions on all diagonal block elements and sparse-sparse cross-covariance functions 
                on all non-diagonal block elements. Cross-covariance amplitude terms are given by the linear correlation coefficients 
                between the corresponding geophysical properties.

                Call function cubing() for running inversion.
                """
                def __init__(self):
                        # Set Lengthscale and amplitude for GP
                        self.gp_length = gp_lengthscale * np.asarray([xvoxsize, xvoxsize, xvoxsize]) # 2*voxelsize seems to have max logl
                        self.gp_sigma = np.asarray(gp_err)
                        self.coeffm = np.asarray(gp_coeff) # coefficients for GP kernel mix


                def create_cubegeometry(self):
                        """
                        Create Cube geometry and coordinates

                        See settings.yaml for cube specifications
                        """
                        # Create voxel Edge coordinates:
                        xedge = np.linspace(0, xNcube, xNcube + 1) * xvoxsize
                        yedge = np.linspace(0, yNcube, yNcube + 1) * yvoxsize
                        zedge = np.linspace(0, -zNcube, zNcube + 1) * zvoxsize + zmax
                        #zedge = np.linspace(0.,Ncube/zfactor, Ncube/zfactor+1) * voxsize
                        xEdges, yEdges, zEdges = np.meshgrid(xedge, yedge, zedge)
                        self.Edges = np.asarray([xEdges, yEdges, -zEdges])
                        # Create voxel coordinates
                        xnew = np.arange(xvoxsize/2., xLcube + xvoxsize/2., xvoxsize)
                        ynew = np.arange(yvoxsize/2., yLcube + yvoxsize/2., yvoxsize)
                        znew = zmax - np.arange(zvoxsize/2., zLcube + zvoxsize/2., zvoxsize)
                        xx ,yy = np.meshgrid(xnew,ynew)
                        self.xxx, self.yyy, self.zzz = np.meshgrid(xnew, ynew, znew)
                        self.voxelpos = np.vstack([self.xxx.flatten(), self.yyy.flatten(), self.zzz.flatten()])
                        return self.voxelpos


                def predict3(self, calclogl = False):
                        """
                        Calculate mean, variance, and likelihood likelihood for GP with 3x3 kernel block matrix 
                    
                    PARAMETERS

                    :param calclogl: if True calculate marginal GP loglikelihood, if False logl return is set to 0

                    RETURN: 
                    mean
                    covariance
                    log-likelihood
                        """
                        # Calculate kernel 3x3 block matrix:
                        self.datastd = np.mean([np.nanstd(self.gravfield), np.nanstd(self.magfield), np.nanstd(self.drillfield)])
                        self.kcov = kernel.create_cov(self.D2, self.gp_length, self.coeffm, fkernel = kernelfunc)
                        #Ak = np.dot(self.Asens3, kcov) # shape (3*Nsensor, 3*Ncube^3)
                        yerr = np.hstack((self.gravfield * 0. + self.gp_sigma[0], self.magfield * 0. + self.gp_sigma[1], self.drillfield *0. + self.gp_sigma[2]))
                        #try:
                        AkA = np.dot(self.Asens3, np.dot(self.kcov, self.Asens3.T)) + np.diag(yerr**2)
                        #AkA = np.dot(Ak, self.Asens3.T) + np.diag(yerr**2) # shape(2*Nsensor, 2*Nsensor)
                        # Cholesky decomposition
                        try:
                                AkA_chol = cholesky(AkA, lower=True)
                        except:
                                print("Cholesky decompostion failed, AkA matrix i likely not positive semitive.")
                                print("Change GP parameter settings")
                                sys.exit(1)
                        usolve = solve_triangular(AkA_chol, self.Fs3, lower=True) #shape(2*Nsensor)
                        # Calculate Likelihood if option is set
                        if calclogl:
                                log_det_AkA = np.log(np.diag(AkA_chol)**2).sum()
                                n_log_2pi = xNcube * yNcube * zNcube * np.log(2 * np.pi)
                                logl = -0.5 * (np.dot(usolve, usolve) + log_det_AkA + n_log_2pi)
                        else:
                                logl = 0.
                        # predicted mean
                        self._V = solve_triangular(AkA_chol, np.dot(self.Asens3, self.kcov), lower=True) #shape(2*Nsensor, 2* Ncube^3)
                        mu = np.dot(self._V.T, usolve) #* _fstd + _fmean
                        # covariance
                        covar = self.kcov - np.dot(self._V.T, self._V)
                        #var = np.diagonal(covar) #* _fstd**2
                        #except:
                        #       print("Warning: GP Inversion failed!")
                        #       mu, covar, logl = np.zeros(3*xNcube*yNcube*zNcube), kcov * 1e9, np.inf 
                        return mu, covar, logl

                
                def cubing(self, gravfield, magfield, drillfield, sensor_locations, drilldata0):
                        """
                        Joint Inversion and Cubing of sensor data

                        PARAMETERS

                        :param gravfield: 2D gravitational survey data
                        :param magfield: 2D magnetic survey data
                        :param drillfield: drilldata

                        RETURN

                        density_rec: cube with mean density 
                        magsus_rec: cube with mean magnetic susceptibility
                        drill_rec: cube with mean drill-core property  
                        density_var: varicance cube with density 
                        magsus_var: variance cube with magnetic susceptibility  
                        drill_var: varianbce cube with drill-core property  


                        """
                        self.gravfield = gravfield
                        self.magfield = magfield
                        self.drillfield = drillfield
                        self.sensor_locations = sensor_locations
                        self.drilldata0 = drilldata0
                        # Normalize data
                        grav_mean, grav_std = self.gravfield.mean(), self.gravfield.std()
                        gravfield_norm = (self.gravfield - grav_mean) / grav_std
                        magn_mean, magn_std = self.magfield.mean(), self.magfield.std()
                        magfield_norm = (self.magfield - magn_mean) / magn_std
                        drill_mean, drill_std = self.drillfield.mean(), self.drillfield.std()
                        drillfield_norm = (self.drillfield - drill_mean) / drill_std
                        # Create kernel distance matrix
                        self.points3D = kernel.calcGridPoints3D((xNcube, yNcube, zNcube), (xvoxsize, yvoxsize, zvoxsize))
                        self.D2 = kernel.calcDistanceMatrix(self.points3D)
                        # Combine Data and Forward models
                        voxelpos_drill = np.vstack([self.xxx[self.drilldata0 != 0], self.yyy[self.drilldata0 != 0], self.zzz[self.drilldata0 != 0]])
                        # Combine data:
                        self.Fs3 = np.hstack((gravfield_norm, magfield_norm, drillfield_norm))
                        # Calculate Sensitivity matrix
                        Asens_grav, _ = A_sens(magneticField * 0., self.sensor_locations, self.Edges, 'grav')
                        Asens_mag, _ = A_sens(magneticField, self.sensor_locations, self.Edges, 'magn')
                        Asens_drill = A_drill(voxelpos_drill, self.voxelpos)
                        # Calculate total transformation matrix
                        Asens_g= np.vstack([Asens_grav, np.zeros_like(Asens_mag), np.zeros_like(Asens_drill)])
                        Asens_m= np.vstack([np.zeros_like(Asens_grav), Asens_mag, np.zeros_like(Asens_drill)])
                        Asens_d= np.vstack([np.zeros_like(Asens_grav), np.zeros_like(Asens_mag), Asens_drill])
                        self.Asens3 = np.hstack([Asens_g, Asens_m, Asens_d])
                        # Run actual GP inversion
                        self.mu_rec, self.cov_rec, self.logl = self.predict3(calclogl = True)
                        # Reshape results and multiply with stadnard deviatio 
                        results_rec = self.mu_rec.reshape(3, yNcube, xNcube, zNcube)
                        results_var = np.diag(self.cov_rec).reshape(3, yNcube, xNcube, zNcube)
                        #Cut off outer voxel rows
                        #results_rec = results_rec[:, 1:yNcube-1, 1:xNcube-1, :]
                        #results_var = results_var[:, 1:yNcube-1, 1:xNcube-1, :]
                        density_rec = results_rec[0] * grav_std  # Model represents deviation from the mean
                        density_var = results_var[0] * grav_std**2
                        magsus_rec = results_rec[1] * magn_std # Model represents deviation from the mean
                        magsus_var = results_var[1] * magn_std**2
                        drill_rec = results_rec[2] * drill_std  # Model represents deviation from the mean
                        drill_var = results_var[2] * drill_std**2
                        return density_rec, magsus_rec, drill_rec, density_var, magsus_var, drill_var

    ### Methods

    ` def create_cubegeometry(self)`{.name .flex}

    :   ::: {.section .desc}
        Create Cube geometry and coordinates

        See settings.yaml for cube specifications
        :::

        Expand source code

            def create_cubegeometry(self):
                    """
                    Create Cube geometry and coordinates

                    See settings.yaml for cube specifications
                    """
                    # Create voxel Edge coordinates:
                    xedge = np.linspace(0, xNcube, xNcube + 1) * xvoxsize
                    yedge = np.linspace(0, yNcube, yNcube + 1) * yvoxsize
                    zedge = np.linspace(0, -zNcube, zNcube + 1) * zvoxsize + zmax
                    #zedge = np.linspace(0.,Ncube/zfactor, Ncube/zfactor+1) * voxsize
                    xEdges, yEdges, zEdges = np.meshgrid(xedge, yedge, zedge)
                    self.Edges = np.asarray([xEdges, yEdges, -zEdges])
                    # Create voxel coordinates
                    xnew = np.arange(xvoxsize/2., xLcube + xvoxsize/2., xvoxsize)
                    ynew = np.arange(yvoxsize/2., yLcube + yvoxsize/2., yvoxsize)
                    znew = zmax - np.arange(zvoxsize/2., zLcube + zvoxsize/2., zvoxsize)
                    xx ,yy = np.meshgrid(xnew,ynew)
                    self.xxx, self.yyy, self.zzz = np.meshgrid(xnew, ynew, znew)
                    self.voxelpos = np.vstack([self.xxx.flatten(), self.yyy.flatten(), self.zzz.flatten()])
                    return self.voxelpos

    ` def cubing(self, gravfield, magfield, drillfield, sensor_locations, drilldata0)`{.name .flex}

    :   ::: {.section .desc}
        Joint Inversion and Cubing of sensor data

        PARAMETERS

        :param gravfield: 2D gravitational survey data :param magfield:
        2D magnetic survey data :param drillfield: drilldata

        RETURN

        density\_rec: cube with mean density magsus\_rec: cube with mean
        magnetic susceptibility drill\_rec: cube with mean drill-core
        property\
        density\_var: varicance cube with density magsus\_var: variance
        cube with magnetic susceptibility\
        drill\_var: varianbce cube with drill-core property
        :::

        Expand source code

            def cubing(self, gravfield, magfield, drillfield, sensor_locations, drilldata0):
                    """
                    Joint Inversion and Cubing of sensor data

                    PARAMETERS

                    :param gravfield: 2D gravitational survey data
                    :param magfield: 2D magnetic survey data
                    :param drillfield: drilldata

                    RETURN

                    density_rec: cube with mean density 
                    magsus_rec: cube with mean magnetic susceptibility
                    drill_rec: cube with mean drill-core property  
                    density_var: varicance cube with density 
                    magsus_var: variance cube with magnetic susceptibility  
                    drill_var: varianbce cube with drill-core property  


                    """
                    self.gravfield = gravfield
                    self.magfield = magfield
                    self.drillfield = drillfield
                    self.sensor_locations = sensor_locations
                    self.drilldata0 = drilldata0
                    # Normalize data
                    grav_mean, grav_std = self.gravfield.mean(), self.gravfield.std()
                    gravfield_norm = (self.gravfield - grav_mean) / grav_std
                    magn_mean, magn_std = self.magfield.mean(), self.magfield.std()
                    magfield_norm = (self.magfield - magn_mean) / magn_std
                    drill_mean, drill_std = self.drillfield.mean(), self.drillfield.std()
                    drillfield_norm = (self.drillfield - drill_mean) / drill_std
                    # Create kernel distance matrix
                    self.points3D = kernel.calcGridPoints3D((xNcube, yNcube, zNcube), (xvoxsize, yvoxsize, zvoxsize))
                    self.D2 = kernel.calcDistanceMatrix(self.points3D)
                    # Combine Data and Forward models
                    voxelpos_drill = np.vstack([self.xxx[self.drilldata0 != 0], self.yyy[self.drilldata0 != 0], self.zzz[self.drilldata0 != 0]])
                    # Combine data:
                    self.Fs3 = np.hstack((gravfield_norm, magfield_norm, drillfield_norm))
                    # Calculate Sensitivity matrix
                    Asens_grav, _ = A_sens(magneticField * 0., self.sensor_locations, self.Edges, 'grav')
                    Asens_mag, _ = A_sens(magneticField, self.sensor_locations, self.Edges, 'magn')
                    Asens_drill = A_drill(voxelpos_drill, self.voxelpos)
                    # Calculate total transformation matrix
                    Asens_g= np.vstack([Asens_grav, np.zeros_like(Asens_mag), np.zeros_like(Asens_drill)])
                    Asens_m= np.vstack([np.zeros_like(Asens_grav), Asens_mag, np.zeros_like(Asens_drill)])
                    Asens_d= np.vstack([np.zeros_like(Asens_grav), np.zeros_like(Asens_mag), Asens_drill])
                    self.Asens3 = np.hstack([Asens_g, Asens_m, Asens_d])
                    # Run actual GP inversion
                    self.mu_rec, self.cov_rec, self.logl = self.predict3(calclogl = True)
                    # Reshape results and multiply with stadnard deviatio 
                    results_rec = self.mu_rec.reshape(3, yNcube, xNcube, zNcube)
                    results_var = np.diag(self.cov_rec).reshape(3, yNcube, xNcube, zNcube)
                    #Cut off outer voxel rows
                    #results_rec = results_rec[:, 1:yNcube-1, 1:xNcube-1, :]
                    #results_var = results_var[:, 1:yNcube-1, 1:xNcube-1, :]
                    density_rec = results_rec[0] * grav_std  # Model represents deviation from the mean
                    density_var = results_var[0] * grav_std**2
                    magsus_rec = results_rec[1] * magn_std # Model represents deviation from the mean
                    magsus_var = results_var[1] * magn_std**2
                    drill_rec = results_rec[2] * drill_std  # Model represents deviation from the mean
                    drill_var = results_var[2] * drill_std**2
                    return density_rec, magsus_rec, drill_rec, density_var, magsus_var, drill_var

    ` def predict3(self, calclogl=False)`{.name .flex}

    :   ::: {.section .desc}
        Calculate mean, variance, and likelihood likelihood for GP with
        3x3 kernel block matrix

        PARAMETERS

        :param calclogl: if True calculate marginal GP loglikelihood, if
        False logl return is set to 0

        RETURN: mean covariance log-likelihood
        :::

        Expand source code

            def predict3(self, calclogl = False):
                    """
                    Calculate mean, variance, and likelihood likelihood for GP with 3x3 kernel block matrix 
                
                PARAMETERS

                :param calclogl: if True calculate marginal GP loglikelihood, if False logl return is set to 0

                RETURN: 
                mean
                covariance
                log-likelihood
                    """
                    # Calculate kernel 3x3 block matrix:
                    self.datastd = np.mean([np.nanstd(self.gravfield), np.nanstd(self.magfield), np.nanstd(self.drillfield)])
                    self.kcov = kernel.create_cov(self.D2, self.gp_length, self.coeffm, fkernel = kernelfunc)
                    #Ak = np.dot(self.Asens3, kcov) # shape (3*Nsensor, 3*Ncube^3)
                    yerr = np.hstack((self.gravfield * 0. + self.gp_sigma[0], self.magfield * 0. + self.gp_sigma[1], self.drillfield *0. + self.gp_sigma[2]))
                    #try:
                    AkA = np.dot(self.Asens3, np.dot(self.kcov, self.Asens3.T)) + np.diag(yerr**2)
                    #AkA = np.dot(Ak, self.Asens3.T) + np.diag(yerr**2) # shape(2*Nsensor, 2*Nsensor)
                    # Cholesky decomposition
                    try:
                            AkA_chol = cholesky(AkA, lower=True)
                    except:
                            print("Cholesky decompostion failed, AkA matrix i likely not positive semitive.")
                            print("Change GP parameter settings")
                            sys.exit(1)
                    usolve = solve_triangular(AkA_chol, self.Fs3, lower=True) #shape(2*Nsensor)
                    # Calculate Likelihood if option is set
                    if calclogl:
                            log_det_AkA = np.log(np.diag(AkA_chol)**2).sum()
                            n_log_2pi = xNcube * yNcube * zNcube * np.log(2 * np.pi)
                            logl = -0.5 * (np.dot(usolve, usolve) + log_det_AkA + n_log_2pi)
                    else:
                            logl = 0.
                    # predicted mean
                    self._V = solve_triangular(AkA_chol, np.dot(self.Asens3, self.kcov), lower=True) #shape(2*Nsensor, 2* Ncube^3)
                    mu = np.dot(self._V.T, usolve) #* _fstd + _fmean
                    # covariance
                    covar = self.kcov - np.dot(self._V.T, self._V)
                    #var = np.diagonal(covar) #* _fstd**2
                    #except:
                    #       print("Warning: GP Inversion failed!")
                    #       mu, covar, logl = np.zeros(3*xNcube*yNcube*zNcube), kcov * 1e9, np.inf 
                    return mu, covar, logl
:::

Index
=====

::: {.toc}
:::

-   ### Super-module

    -   `geobo`

-   ### [Classes](#header-classes) {#classes}

    -   #### `Inversion`

        -   `create_cubegeometry`
        -   `cubing`
        -   `predict3`
:::

Generated by [pdoc 0.7.4](https://pdoc3.github.io/pdoc).
