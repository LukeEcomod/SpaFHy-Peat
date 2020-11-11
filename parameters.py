# -*- coding: utf-8 -*-
"""
PARAMETERS

@author: slauniai & khaahti
"""
import pathlib

def parameters(folder=''):

    pgen = {'description': 'testcase',  # description written in result file
            'start_date': '1981-01-01',
            'end_date': '1991-12-31',
            'spinup_end': '1982-01-01',
            'dt': 86400.0,
            'spatial_cpy': True,  # if False uses parameters from cpy['state']
            # else needs cf.dat, hc.dat, LAI_decid.dat, LAI_spruce.dat, LAI_pine.dat, (cmask.dat)
            'spatial_soil': True,  # if False uses soil_id, ditch_depth, ditch_spacing from psp
            # else needs soil_id.dat, ditch_depth.dat, ditch_spacing.dat
            'spatial_forcing': True,  # if False uses forcing from forcing file with pgen['forcing_id'] and cpy['loc']
            # else needs Ncoord.dat, Ecoord.dat, forcing_id.dat
            'stand_development': False,  # if True stand characteristics change annually accoording to input,
            # give input (cf.dat, hc.dat, LAI_decid.dat, LAI_spruce.dat, LAI_pine.dat) for each year in columns
            'gis_folder': str(pathlib.Path(folder+r'/parameters')),
            'forcing_file': str(pathlib.Path(folder+r'/forcing/Weather_id_[forcing_id].csv')),
            'forcing_id': 0,  # used if spatial_forcing == False
            'ncf_file': folder + r'.nc',
            'results_folder': r'results/',
            'save_interval': 366, # interval for writing results to file (decreases need for memory during computation)
            'variables':[ # list of output variables (rows can be commented out if not all variables are of interest)
                    ['parameters_lai_conif', 'leaf area index of conifers [m2 m-2]'],
                    ['parameters_lai_decid_max', 'leaf area index of decidious trees [m2 m-2]'],
                    # ['parameters_hc', 'canopy height [m]'],
                    ['parameters_cf', 'canopy closure [-]'],
                    # ['parameters_soil_id', 'soil class index'],
                    ['parameters_ditch_depth', 'ditch depth [m]'],
                    ['parameters_ditch_spacing', 'ditch spacing [m]'],
                    ['parameters_lat', 'latitude [deg]'],
                    ['parameters_lon', 'longitude [deg]'],
                    ['forcing_air_temperature', 'air temperature [degC]'],
                    ['forcing_precipitation', 'precipitation [mm d-1]'],
                    ['forcing_vapor_pressure_deficit', 'vapor pressure deficit [kPa]'],
                    ['forcing_global_radiation', 'global radiation [Wm-2]'],
                    ['forcing_CO2', 'CO2 mixing ratio [ppm]'],
                    # ['forcing_wind_speed','wind speed [m s-1]'],
                    # ['forcing_snow_depth', 'snow depth [cm]'],
                    # ['soil_pond_storage', 'pond storage [m]'],
                    ['soil_ground_water_level', 'ground water level [m]'],
                    # ['soil_infiltration', 'infiltration [mm d-1]'],
                    # ['soil_surface_runoff', 'surface runoff [mm d-1]'],
                    ['soil_evaporation', 'evaporation from soil surface [mm d-1]'],
                    # ['soil_drainage', 'subsurface drainage [mm d-1]'],
                    # ['soil_moisture_top', 'volumetric water content of moss layer [m3 m-3]'],
                    # ['soil_rootzone_moisture', 'volumetric water content of rootzone [m3 m-3]'],
                    # ['soil_water_closure', 'soil water balance error [mm d-1]'],
                    # ['soil_transpiration_limitation', 'transpiration limitation [-]'],
                    # ['canopy_interception', 'canopy interception [mm d-1]'],
                    ['canopy_evaporation', 'evaporation from interception storage [mm d-1]'],
                    ['canopy_transpiration','transpiration [mm d-1]'],
                    # ['canopy_stomatal_conductance','stomatal conductance [m s-1]'],
                    # ['canopy_gs_raw','stomatal conductance [m s-1]'],
                    # ['canopy_throughfall', 'throughfall to moss or snow [mm d-1]'],
                    ['canopy_snow_water_equivalent', 'snow water equivalent [mm]'],
                    # ['canopy_water_closure', 'canopy water balance error [mm d-1]'],
                    # ['canopy_phenostate', 'canopy phenological state [-]'],
                    ['canopy_leaf_area_index', 'canopy leaf area index [m2 m-2]'],
                    ['canopy_degree_day_sum', 'sum of degree days [degC]'],
                    ]
             }

    # canopygrid
    pcpy = {'flow' : {  # flow field
                     'zmeas': 2.0,
                     'zground': 0.5,
                     'zo_ground': 0.01
                     },
            'interc': {  # interception
                        'wmax': 1.5,  # storage capacity for rain (mm/LAI)
                        'wmaxsnow': 4.5,  # storage capacity for snow (mm/LAI)
                        'c_snow': 1.0,  #correctioon for snow fall (-)
                        'Tmin': 0.0,  # temperature below which all is snow [degC]
                        'Tmax': 2.0,  # temperature above which all is water [degC]- Koivusalo & Kokkonen 2002
                        },
            'snow': {  # degree-day snow model
                    'kmelt': 2.5,  # melt coefficient in open (mm/d)
                    'kfreeze': 0.5,  # freezing coefficient (mm/d)
                    'r': 0.05  # maximum fraction of liquid in snow (-)
                    },
            'physpara': {
                        # canopy conductance
                        'amax': 10.0, # maximum photosynthetic rate (umolm-2(leaf)s-1)
                        'g1_conif': 2.1, # stomatal parameter, conifers
                        'g1_decid': 3.5, # stomatal parameter, deciduous
                        'q50': 50.0, # light response parameter (Wm-2)
                        'kp': 0.6, # light attenuation parameter (-)
                        # soil evaporation
                        'gsoil': 1e-2 # soil surface conductance if soil is fully wet (m/s)
                        },
            'phenopara': {
                        # seasonal cycle of physiology: smax [degC], tau[d], xo[degC],fmin[-](residual photocapasity)
                        'smax': 18.5, # degC
                        'tau': 13.0, # days
                        'xo': -4.0, # degC
                        'fmin': 0.05, # minimum photosynthetic capacity in winter (-)
                        # deciduos phenology
                        'lai_decid_min': 0.1, # minimum relative LAI (-)
                        'ddo': 45.0, # degree-days for bud-burst (5degC threshold)
                        'ddur': 23.0, # duration of leaf development (days)
                        'sdl': 9.0, # daylength for senescence start (h)
                        'sdur': 30.0, # duration of leaf senescence (days),
                         },
            'state': {  # following properties are used if spatial_cpy == False
                       'lai_conif': 3.5, # conifer 1-sided LAI (m2 m-2)
                       'lai_decid_max': 0.5, # maximum annual deciduous 1-sided LAI (m2 m-2)
                       'hc': 16.0, # canopy height (m)
                       'cf': 0.6, # canopy closure fraction (-)
                       #initial state of canopy storage [mm] and snow water equivalent [mm]
                       'w': 0.0, # canopy storage mm
                       'swe': 0.0, # snow water equivalent mm
                       },
            'loc': {  # following coordinates used if spatial_forcing == False
                    'lat': 61.4,  # decimal degrees
                    'lon': 23.4
                    }
            }

    # soil profile
    psp = {
            # soil profile, following properties are used if spatial_soil = False
            'soil_id': 2.0,
            # drainage parameters, following properties are used if spatial_soil = False
            'ditch_depth': 1.0,  # ditch depth [m]
            'ditch_spacing': 45.0,  # ditch spacing [m]
            'ditch_width': 1.0,  # ditch width [m]
            # organic (moss) layer
            'org_depth': 0.04, # depth of organic top layer (m)
            'org_poros': 0.9, # porosity (-)
            'org_fc': 0.3, # field capacity (-)
            'org_rw': 0.24, # critical vol. moisture content (-) for decreasing phase in Ef
            'pond_storage_max': 0.05,  # maximum pond depth [m]
            # initial states
            'ground_water_level': -0.2,  # groundwater depth [m]
            'org_sat': 1.0, # organic top layer saturation ratio (-)
            'pond_storage': 0.0  # initial pond depth at surface [m]
            }

    return pgen, pcpy, psp

def peat_soilprofiles():
    """
    Properties of typical peat profiles
    """
    peatp = {
        'sphagnum':{
            'soil_id': 1.0,
            'z': [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1., -1.5, -2.0],
            'pF': {  # vanGenuchten water retention parameters
                    'ThetaS': [0.945, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918],
                    'ThetaR': [0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098],
                    'alpha': [0.338, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072],
                    'n': [1.402, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371]},
            'saturated_conductivity': [30*8.99E-05, 20*2.98E-05, 10*9.86E-06, 3.27E-06, 1.08E-06, 3.58E-07, 1.19E-07, 1.16E-07, 1.16E-07, 1.16E-07, 1.16E-07, 1.16E-07],
                },
        'carex': {
            'soil_id': 2.0,
            'z': [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1., -1.5, -2.0],
            'pF': {  # vanGenuchten water retention parameters
                    'ThetaS': [0.943, 0.874, 0.874, 0.874, 0.874, 0.874, 0.874, 0.874, 0.874, 0.874, 0.874, 0.874],
                    'ThetaR': [0.002, 0.198, 0.198, 0.198, 0.198, 0.198, 0.198, 0.198, 0.198, 0.198, 0.198, 0.198],
                    'alpha': [0.202, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030],
                    'n': [1.349, 1.491, 1.491, 1.491, 1.491, 1.491, 1.491, 1.491, 1.491, 1.491, 1.491, 1.491]},
            'saturated_conductivity': [30*4.97E-05, 20*3.21E-05, 10*2.07E-05, 1.34E-05, 8.63E-06, 5.57E-06, 3.60E-06, 2.32E-06, 1.50E-06, 9.68E-07, 2.61E-07, 1.16E-07],
                },
            }
    return peatp
