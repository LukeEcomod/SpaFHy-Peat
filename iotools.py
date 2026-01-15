# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 10:34:37 2016

@author: slauniai & khaahti

"""
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from soilprofile import gwl_Wsto, gwl_Ksat, nan_function

eps = np.finfo(float).eps  # machine epsilon
workdir = os.getcwd()

def read_soil_gisdata(fpath, plotgrids=False):
    """
    reads gis-data grids and returns numpy 2d-arrays
    Args:
        fpath - relative path to data folder (str)
        plotgrids - True plots
    Returns:
        gis - dict of gis-data rasters
            cmask
            soil_id
            ditch_depth
            ditch_spacing
    """
    fpath = os.path.join(workdir, fpath)

    # soil classification
    soilclass, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'soil_id.dat'))

    # ditch depth and spacing
    ditch_depth, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'ditch_depth.dat'))
    ditch_spacing, _, _, _, _ = read_AsciiGrid(os.path.join(fpath,'ditch_spacing.dat'))

    # catchment mask cmask[i,j] == 1, np.NaN outside
    if os.path.isfile(os.path.join(fpath, 'cmask.dat')):
        cmask, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'cmask.dat'))
    else:
        cmask = np.ones(np.shape(soilclass))

    # dict of all rasters
    gis = {'cmask': cmask,
           'soilclass': soilclass,
           'ditch_depth': ditch_depth,
           'ditch_spacing': ditch_spacing
           }

    for key in gis.keys():
        gis[key] *= cmask

    if plotgrids is True:
        plt.figure()
        plt.subplot(311); plt.imshow(soilclass); plt.colorbar(); plt.title('soiltype')
        plt.subplot(312); plt.imshow(ditch_depth); plt.colorbar();
        plt.title('ditch depth (m)')
        plt.subplot(313); plt.imshow(ditch_spacing); plt.colorbar();
        plt.title('ditch spacing (m)')

    return gis

def read_cpy_gisdata(fpath, plotgrids=False):
    """
    reads gis-data grids and returns numpy 2d-arrays
    Args:
        fpath - relative path to data folder (str)
        plotgrids - True plots
    Returns:
        gis - dict of gis-data rasters
            cmask
            LAI_pine, LAI_spruce - pine and spruce LAI (m2m-2)
            LAI_conif - conifer total annual max LAI (m2m-2)
            LAI_dedid - deciduous annual max LAI (m2m-2)
            cf - canopy closure (-)
            hc - mean stand height (m)

    """
    fpath = os.path.join(workdir, fpath)

    # tree height [m]
    hc, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'hc.dat'))

    # canopy closure [-]
    cf, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'cf.dat'))

    # leaf area indices
    try:
        LAI_pine, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'LAI_pine.dat'))
        LAI_spruce, _, _, _, _ = read_AsciiGrid(os.path.join(fpath,'LAI_spruce.dat'))
        LAI_conif = LAI_pine + LAI_spruce
    except:
        LAI_conif, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'LAI_conif.dat'))
    LAI_decid, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'LAI_decid.dat'))

    # catchment mask cmask[i,j] == 1, np.NaN outside
    if os.path.isfile(os.path.join(fpath, 'cmask.dat')):
        cmask, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'cmask.dat'))
        # if cmask.shape != hc.shape: DO SOMETHING !!
    else:
        cmask = np.ones(np.shape(hc))

    # dict of all rasters
    gis = {'cmask': cmask,
           'LAI_conif': LAI_conif,
           'LAI_decid': LAI_decid, 'hc': hc, 'cf': cf}

    for key in gis.keys():
        gis[key] *= cmask

    if plotgrids is True:

        plt.figure()
        plt.subplot(221); plt.imshow(LAI_pine+LAI_spruce); plt.colorbar();
        plt.title('LAI conif (m2/m2)')
        plt.subplot(222); plt.imshow(LAI_decid); plt.colorbar();
        plt.title('LAI decid (m2/m2)')
        plt.subplot(223); plt.imshow(hc); plt.colorbar(); plt.title('hc (m)')
        plt.subplot(224); plt.imshow(cf); plt.colorbar(); plt.title('cf (-)')

    return gis

def read_forcing_gisdata(fpath):
    """
    reads gis-data grids and returns numpy 2d-arrays
    Args:
        fpath - relative path to data folder (str)
    Returns:
        gis - dict of gis-data rasters
            cmask
            lat
            lon
            forcing_id

    """
    fpath = os.path.join(workdir, fpath)

    # latitude and longitude
    lat, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'latitude.dat'))
    lon, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'longitude.dat'))

    forcing_id, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'forcing_id.dat'))

    # catchment mask cmask[i,j] == 1, np.NaN outside
    if os.path.isfile(os.path.join(fpath, 'cmask.dat')):
        cmask, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'cmask.dat'))
    else:
        cmask = np.ones(np.shape(lat))

    # dict of all rasters
    gis = {'cmask': cmask, 'lat': lat, 'lon': lon, 'forcing_id': forcing_id}

    for key in gis.keys():
        gis[key] *= cmask

    return gis

def preprocess_soildata(psp, peatp, gisdata, spatial=True):
    """
    creates input dictionary for initializing SoilGrid
    Args:
        soil parameters
        soiltype parameters
        gisdata
            cmask
            soilclass
        spatial
    """
    # create dict for initializing soil profile.
    # copy pbu into sdata and make each value np.array(np.shape(cmask))
    data = psp.copy()
    data.update((x, y * gisdata['cmask']) for x, y in data.items())

    data.update({'soiltype': np.empty(np.shape(gisdata['cmask']),dtype=object),
                 'depth_id': np.empty(np.shape(gisdata['cmask']),dtype=int)})

    if spatial == False:
        data['soilclass'] = psp['soil_id'] * gisdata['cmask']
    else:
        data['soilclass'] = gisdata['soilclass']
        data['ditch_depth'] = gisdata['ditch_depth']
        data['ditch_spacing'] = gisdata['ditch_spacing']

    data['gwl_to_Ksat'] = np.full(
            len(peatp)*len(np.unique(np.round(data['ditch_depth'],2))),
            nan_function, dtype=object)

    soil_ids = []
    for key, value in peatp.items():
        soil_ids.append(value['soil_id'])

    if set(soil_ids) >= set(np.unique(data['soilclass']).tolist()):
        # no problems
        pass
    else:
        raise ValueError("Soil id in inputs not specified in parameters.py")

    i = 0
    for key, value in peatp.items():
        c = value['soil_id']
        ix = np.where(data['soilclass'] == c)
        data['soiltype'][ix] = key
        # interpolation function between wsto and gwl
        value.update(gwl_Wsto(value['z'], value['pF']))
        # interpolation function between root_wsto and gwl
        value.update(gwl_Wsto(value['z'][:2], {key: value['pF'][key][:2] for key in value['pF'].keys()}, root=True))

        for depth in np.unique(np.round(data['ditch_depth'][ix],2)):
            data['gwl_to_Ksat'][i] = gwl_Ksat(value['z'],
                    value['saturated_conductivity'], depth)
            ixx = np.where((np.round(data['ditch_depth'],2) == depth) &
                           (data['soiltype'] == key))
            data['depth_id'][ixx] = i
            i=i+1

    data['gwl_to_Ksat'] = data['gwl_to_Ksat'][:i]

    data['wtso_to_gwl'] = {soiltype: peatp[soiltype]['to_gwl'] for soiltype in peatp.keys()}
    data['gwl_to_wsto'] = {soiltype: peatp[soiltype]['to_wsto'] for soiltype in peatp.keys()}
    data['gwl_to_rootmoist'] = {soiltype: peatp[soiltype]['to_rootmoist'] for soiltype in peatp.keys()}

    return data

def preprocess_cpydata(pcpy, gisdata, spatial=True):
    """
    creates input dictionary for initializing CanopyGrid
    Args:
        canopy parameters
        gisdata
            cmask
            LAI_pine, LAI_spruce - pine and spruce LAI (m2m-2)
            LAI_conif - conifer total annual max LAI (m2m-2)
            LAI_dedid - deciduous annual max LAI (m2m-2)
            cf - canopy closure (-)
            hc - mean stand height (m)
            (lat, lon)
        spatial
    """
    # inputs for CanopyGrid initialization: update pcpy using spatial data
    cstate = pcpy['state'].copy()

    if spatial:
        cstate['lai_conif'] = gisdata['LAI_conif']
        cstate['lai_decid_max'] = gisdata['LAI_decid']
        cstate['cf'] = gisdata['cf']
        cstate['hc'] = gisdata['hc']
        for key in ['w', 'swe']:
            cstate[key] *= gisdata['cmask']
        if {'lat','lon'}.issubset(gisdata.keys()):
            pcpy['loc']['lat'] = gisdata['lat']
            pcpy['loc']['lon'] = gisdata['lon']
    else:
        for key in cstate.keys():
            cstate[key] *= gisdata['cmask']

    pcpy['state'] = cstate

    return pcpy

def read_FMI_weather(start_date, end_date, sourcefile, CO2=400.0, U=2.0):
    """
    reads FMI interpolated daily weather data from file
    """

    sourcefile = os.path.join(sourcefile)

    # import forcing data
    try:
        fmi = pd.read_csv(sourcefile, sep=';', header='infer',
                          usecols=['OmaTunniste', 'Kunta', 'aika','vuosi','kk','paiva',
                          'longitude','latitude', 't_mean', 't_max', 't_min', 'rainfall',
                          'radiation', 'hpa', 'lamposumma_v', 'rainfall_v'],
                          parse_dates=['aika'],encoding="ISO-8859-1")

        fmi['aika'] = pd.to_datetime({'year': fmi['vuosi'],
                                    'month': fmi['kk'],
                                    'day': fmi['paiva']})

        fmi = fmi.rename(columns={'aika': 'date',
                                  'OmaTunniste': 'ID',
                                  't_mean': 'air_temperature',
                                  'rainfall': 'precipitation',
                                  'radiation': 'global_radiation',
                                  'hpa': 'h2o'})

        time = fmi['date']
    except:
        try:
            fmi = pd.read_csv(sourcefile, sep=';', header='infer',
                              usecols=['x','y','date','temp_avg','prec',
                              'wind_speed_avg','global_rad','vapour_press'],
                              parse_dates=['date'],encoding="ISO-8859-1")

            fmi = fmi.rename(columns={'temp_avg': 'air_temperature',
                                      'prec': 'precipitation',
                                      'global_rad': 'global_radiation',
                                      'vapour_press': 'h2o',
                                      'wind_speed_avg':'wind_speed'})

            time = pd.to_datetime(fmi['date'], format='%Y-%m-%d')
        except:
            raise ValueError('Problem reading forcing data')

    fmi.index = time
    # get desired period
    fmi = fmi[(fmi.index >= start_date) & (fmi.index <= end_date)]

    fmi['h2o'] = 1e-1*fmi['h2o']  # hPa-->kPa
    fmi['global_radiation'] = 1e3 / 86400.0*fmi['global_radiation']  # kJ/m2/d-1 to Wm-2

    # saturated vapor pressure
    esa = 0.6112*np.exp(
            (17.67*fmi['air_temperature']) / (fmi['air_temperature'] + 273.16 - 29.66))  # kPa
    vpd = esa - fmi['h2o']  # kPa
    vpd[vpd < 0] = 0.0
    rh = 100.0*fmi['h2o'] / esa
    rh[rh < 0] = 0.0
    rh[rh > 100] = 100.0

    fmi['RH'] = rh
    fmi['esa'] = esa
    fmi['vapor_pressure_deficit'] = vpd

    fmi['doy'] = fmi.index.dayofyear
    # replace nan's in prec with 0.0
    fmi['precipitation'] = fmi['precipitation'].fillna(0.0)

    fmi['par'] = 0.45*fmi['global_radiation']
    fmi.loc[fmi['vapor_pressure_deficit'] < 0.0, 'vapor_pressure_deficit'] = 0.0

    # add CO2 and wind speed concentration to dataframe
    # print('CO2 set constant: ' + str(CO2) + ' ppm')
    fmi['CO2'] = float(CO2)
    if 'wind_speed' not in fmi:
        fmi['wind_speed'] = float(U)

    fmi['wind_speed'] = fmi['wind_speed'].fillna(U)

    dates = pd.date_range(start_date, end_date).tolist()

    if len(dates) != len(fmi):
        print(str(len(dates) - len(fmi)) + ' days missing from forcing file, interpolated')
    forcing = pd.DataFrame(index=dates, columns=[])
    forcing = forcing.merge(fmi, how='outer', left_index=True, right_index=True)
    forcing = forcing.fillna(method='ffill')

    return forcing

def initialize_netcdf(pgen, cmask, filepath, filename, description):
    """
    netCDF4 format output file initialization

    Args:
        variables (list): list of variables to be saved in netCDF4
        cmask
        filepath: path for saving results
        filename: filename
        description: description
    """
    from netCDF4 import Dataset, date2num
    from datetime import datetime

    # dimensions
    date_dimension = None
    i_dimension, j_dimension = np.shape(cmask)

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    ff = os.path.join(filepath, filename)

    # create dataset and dimensions
    ncf = Dataset(ff, 'w')
    ncf.description = 'SpaFHy results : ' + description
    ncf.history = 'created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    ncf.source = 'modified SpaFHy v.1.0'

    ncf.createDimension('date', date_dimension)
    ncf.createDimension('i', i_dimension)
    ncf.createDimension('j', j_dimension)

    date = ncf.createVariable('date', 'f8', ('date',))
    date.units = 'days since 0001-01-01 00:00:00.0'
    date.calendar = 'standard'
    tvec = pd.date_range(pgen['spinup_end'], pgen['end_date']).tolist()[1:]
    date[:] = date2num(tvec, units=date.units, calendar=date.calendar)

    for var in pgen['variables']:

        var_name = var[0]
        var_unit = var[1]

        if (var_name.split('_')[0] == 'forcing' and
            pgen['spatial_forcing'] == False):
            var_dim = ('date')
        elif var_name.split('_')[0] == 'parameters':
            var_dim = ('i', 'j')
        else:
            var_dim = ('date','i', 'j')

        variable = ncf.createVariable(
                var_name, 'f4', var_dim)

        variable.units = var_unit

    return ncf, ff

def write_ncf(results, ncf, steps=None):
    """
    Writes model simultaion results in netCDF4-file

    Args:
        index (int): model loop index
        results (dict): calculation results from group
        ncf (object): netCDF4-file handle
    """

    keys = results.keys()
    variables = ncf.variables.keys()

    for key in keys:

        if key in variables and key != 'date':
            if len(ncf[key].shape) > 2:
                if steps==None:
                    ncf[key][:,:,:] = results[key]
                else:
                    ncf[key][steps[0]:steps[1],:,:] = results[key][0:steps[1]-steps[0],:,:]
            elif len(ncf[key].shape) > 1:
                ncf[key][:,:] = results[key]
            else:
                if steps==None:
                    ncf[key][:] = results[key]
                else:
                    ncf[key][steps[0]:steps[1]] = results[key][0:steps[1]-steps[0]]

def read_AsciiGrid(fname, setnans=True):
    """
    reads AsciiGrid format in fixed format as below:
        ncols         750
        nrows         375
        xllcorner     350000
        yllcorner     6696000
        cellsize      16
        NODATA_value  -9999
        -9999 -9999 -9999 -9999 -9999
        -9999 4.694741 5.537514 4.551162
        -9999 4.759177 5.588773 4.767114
    IN:
        fname - filename (incl. path)
    OUT:
        data - 2D numpy array
        info - 6 first lines as list of strings
        (xloc,yloc) - lower left corner coordinates (tuple)
        cellsize - cellsize (in meters?)
        nodata - value of nodata in 'data'
    Samuli Launiainen Luke 7.9.2016
    """
    import numpy as np

    fid = open(fname, 'r')
    info = fid.readlines()[0:6]
    fid.close()

    # print info
    # conversion to float is needed for non-integers read from file...
    xloc = float(info[2].split(' ')[-1])
    yloc = float(info[3].split(' ')[-1])
    cellsize = float(info[4].split(' ')[-1])
    nodata = float(info[5].split(' ')[-1])

    # read rest to 2D numpy array
    data = np.loadtxt(fname, skiprows=6)

    if setnans is True:
        data[data == nodata] = np.nan
        nodata = np.nan

    data = np.array(data, ndmin=2)

    return data, info, (xloc, yloc), cellsize, nodata

def read_results(outputfile):
    """
    Opens simulation results netcdf4 dataset in xarray
    Args:
        outputfile (str): outputfilename
    Returns:
        results (xarray): simulation results from given outputfile
    """

    import xarray as xr

    result = xr.open_dataset(outputfile)
    result.coords['i'] = -np.arange(0,result.dims['i'])
    result.coords['j'] = np.arange(0,result.dims['j'])

    return result
