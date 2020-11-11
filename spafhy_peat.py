# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 16:18:57 2016

@author: slauniai & khaahti

"""

import numpy as np
import pandas as pd
from canopygrid import CanopyGrid
from soilprofile import SoilGrid

eps = np.finfo(float).eps  # machine epsilon

""" ************** SpaFHy v1.0 ************************************

Simple spatial hydrology and catchment water balance model.

CONSISTS OF THREE CLASSES, defined in separate modules:
    CanopyGrid - vegetation and snowpack water storages and flows
    BucketGrid - topsoil bucket model (root zone / topsoil water storage)
    Topmodel - integration to catchment scale using Topmodel -concept
HELPER FUNCTIONS:
    spafhy_parameters - parameter definition file
    spafhy_io - utility functions for data input & output

MAIN PROGRAM:
    spafhy_driver is main program, call it as
    outargs = spathy_driver(spathyparamfile, args)

    spathyparamfile - path to parameter file, default is 'spathy_default.ini'
    soil type dependent parameters are in 'soilparam.ini'

NEEDS 2D gis rasters in ascii-grid format

CanopyGrid & BucketGrid can be initialized from gis-data or set to be spatially constant

ToDo:
    CanopyGrid:
        -include topographic shading to radiation received at canopy top
        -radiation-based snowmelt coefficient
        -add simple GPP-model; 2-step Farquhar or LUE-based approach
    BucketGrid:
        -make soil hydrologic properties more realistic e.g. using pedotransfer functions
        -kasvupaikkatyyppi (multi-NFI) --> soil properties
        -add soil frost model, simplest would be Stefan equation with coefficients modified based on snow insulation
          --> we need snow density algorithm: SWE <-----> depth
    Topmodel:
        -think of definging 'relative m & to grids' (soil-type & elevation-dependent?) and calibrate 'catchment averages'
        -topmodel gives 'saturated zone storage deficit in [m]'. This can be converted to gwl proxy (?) if:
        local water retention characteristics are known & hydrostatic equilibrium assumes.
        Look which water retention model was analytically integrable (Campbell, brooks-corey?)

    Graphics and analysis of results:
        -make ready functions


(C) Samuli Launiainen 10/2016-->

VERSION 05.10.2018 / equations correspond to GMDD paper

"""



"""
******************************************************************************
            ----------- SpaFHy model class --------
******************************************************************************
"""


class SpaFHy():
    """
    SpaFHy model class
    """
    def __init__(self, pgen, pcpy, psoil):

        self.dt = pgen['dt']  # s

        """--- initialize CanopyGrid ---"""
        self.cpy = CanopyGrid(pcpy, pcpy['state'])

        """--- initialize SoilGrid ---"""
        self.soil = SoilGrid(psoil)

    def run_timestep(self, forc):
        """
        Runs SpaFHy for one timestep starting from current state
        Args:
            forc - dictionary or pd.DataFrame containing forcing values for the timestep
            ncf - netCDF -file handle, for outputs
            flx - returns flux and state grids to caller as dict
            ave_flx - returns averaged fluxes and states to caller as dict
        Returns:
            dict of results from canopy and soil models
        """
        doy = forc['doy'].values
        ta = forc['air_temperature'].values
        vpd = forc['vapor_pressure_deficit'].values + eps
        rg = forc['global_radiation'].values
        par = forc['par'].values + eps
        prec = forc['precipitation'].values
        co2 = forc['CO2'].values
        u = forc['wind_speed'].values + eps

        # run CanopyGrid
        canopy_results = self.cpy.run_timestep(
                doy, self.dt, ta, prec, rg, par, vpd, U=u, CO2=co2,
                beta=self.soil.Ree, Rew=self.soil.Rew, P=101300.0)

        # run Soilprofile water balance
        soil_results = self.soil.watbal(
                dt=self.dt,
                rr=1e-3*canopy_results['potential_infiltration'],
                tr=1e-3*canopy_results['transpiration'],
                evap=1e-3*canopy_results['forestfloor_evaporation'])

        return canopy_results, soil_results
