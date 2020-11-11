# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 13:52:46 2019

@author: khaaahti
"""

import time
import numpy as np
import pandas as pd
from spafhy_peat import SpaFHy
from iotools import read_FMI_weather, initialize_netcdf, write_ncf
import matplotlib.pyplot as plt
from copy import deepcopy as copy

eps = np.finfo(float).eps

def driver(create_ncf=False, output=True, folder=''):
    """
    Model driver: sets up model, runs it and saves results to file (create_ncf==True)
    or return dictionary of results.
    """

    """ set up model """
    running_time = time.time()

    # load and process parameters parameter
    pgen, pcpy_all, psoil, cmask = preprocess_parameters(folder)

    # if stand development is enabled pcpy['state'] includes values for all years
    pcpy = copy(pcpy_all)
    if pgen['stand_development']:
        if psoil['soil_id'].shape[0] != 1:
            raise ValueError("When stand development enabled, other than stand charactristics should be given in one column")
        for key in ['lai_conif', 'lai_decid_max', 'hc', 'cf']:
            if pcpy_all['state'][key].shape[0] == 1:
                raise ValueError("When stand development enabled, give stand characteristics for each year in columns")
            pcpy['state'][key] = np.array(pcpy_all['state'][key][:,0], ndmin=2)

    # initialize SpaFHy
    spa = SpaFHy(pgen, pcpy, psoil)

    # read forcing data
    forcing = preprocess_forcing(pgen)

    Nsteps = len(forcing['date'])
    Nspin = (pd.to_datetime(pgen['spinup_end']) - pd.to_datetime(pgen['start_date'])).days + 1

    # results dictionary to accumulate simulation results
    # for save_interval at a time
    if create_ncf:
        save_interval = min(pgen['save_interval'], Nsteps - Nspin)
        results = _create_results(pgen, cmask, save_interval)
    else:
        save_interval = Nsteps - Nspin
        results = _create_results(pgen, cmask, Nsteps - Nspin)

    Nsaveresults = list(np.arange(Nspin, Nsteps, save_interval)[1:] - 1)

    # save parameters to results
    results = _append_results('parameters', pcpy['state'], results)
    results = _append_results('parameters', pcpy['loc'], results)
    results = _append_results('parameters', psoil, results)

    if create_ncf:
        ncf, outputfile = initialize_netcdf(
                pgen=pgen,
                cmask=cmask,
                filepath=pgen['results_folder'],
                filename=pgen['ncf_file'],
                description=pgen['description'])

    print('*** Running model ***')

    interval = 0
    Nsaved = Nspin - 1
    Nyear = 0

    for k in range(0, Nsteps):

        if pgen['stand_development']:
            if forcing['date.dayofyear'].values[k] == 1 and k > 10:
                Nyear += 1
                for key in ['lai_conif', 'lai_decid_max', 'hc', 'cf']:
                    pcpy['state'][key] = np.array(pcpy_all['state'][key][:,Nyear], ndmin=2)
                # update canopy state:
                spa.cpy.update_state(pcpy['state'])
                print('Stand characteristics updated ' + forcing['date'].dt.strftime("%d-%m-%Y").values[k])

        canopy_results, soil_results = spa.run_timestep(forcing.isel(date=k))

        if k >= Nspin:  # save results after spinup done
            results = _append_results('canopy', canopy_results, results, k - Nsaved - 1)
            results = _append_results('soil', soil_results, results, k - Nsaved - 1)

            if k in Nsaveresults and create_ncf:
                interval += 1
                print('*** Writing results to netCDF4-file, subset %.0f/%.0f ***' % (interval, len(Nsaveresults)+1))
                # save forcing to results
                results = _append_results('forcing', forcing[dict(
                        date=slice(Nsaved + 1, k + 1))], results)
                write_ncf(results=results, ncf=ncf, steps=[Nsaved + 1 - Nspin, k + 1 - Nspin])
                Nsaved = k

    if create_ncf:
        interval += 1
        print('*** Writing results to netCDF4-file, subset %.0f/%.0f ***' % (interval, len(Nsaveresults)+1))
        # save forcing to results
        results = _append_results('forcing', forcing[dict(
                date=slice(Nsaved + 1, k + 1))], results)
        write_ncf(results=results, ncf=ncf, steps=[Nsaved + 1 - Nspin, k + 1 - Nspin])
        ncf.close()
        print('--- Running time %.2f seconds ---' % (time.time() - running_time))
        print('--- Results are in file: ' + outputfile + ' ---')
        if output:
            return outputfile
    else:
        print('--- Running time %.2f seconds ---' % (time.time() - running_time))
        if output:
            return results

def preprocess_parameters(folder=''):
    """
    Reading gisdata if applicable and preprocesses parameters
    """

    from iotools import read_soil_gisdata, read_cpy_gisdata, read_forcing_gisdata
    from iotools import preprocess_soildata, preprocess_cpydata
    from parameters import peat_soilprofiles, parameters

    pgen, pcpy, psp= parameters(folder)
    peatp = peat_soilprofiles()
    gisdata = {}

    if pgen['spatial_soil']:
        gisdata.update(read_soil_gisdata(pgen['gis_folder']))
    if pgen['spatial_cpy']:
        gisdata.update(read_cpy_gisdata(pgen['gis_folder']))
    if pgen['spatial_forcing']:
        gisdata.update(read_forcing_gisdata(pgen['gis_folder']))
        pgen.update({'forcing_id': gisdata['forcing_id']})
    if (pgen['spatial_cpy'] == False and
        pgen['spatial_soil'] == False and
        pgen['spatial_forcing'] == False):
        gisdata = {'cmask': np.ones((1,1))}

    soildata = preprocess_soildata(psp, peatp, gisdata, pgen['spatial_soil'])

    cpydata = preprocess_cpydata(pcpy, gisdata, pgen['spatial_cpy'])

    return pgen, cpydata, soildata, gisdata['cmask']

def preprocess_forcing(pgen):
    """
    Reads forcing file(s) based on indices in pgen['forcing_id']
    Creates xarray dataset of forcing data
    """

    import xarray as xr

    indices = np.array(pgen['forcing_id'],ndmin=2)

    variables = ['doy',
                 'air_temperature',
                 'vapor_pressure_deficit',
                 'global_radiation',
                 'par',
                 'precipitation',
                 'CO2',
                 'wind_speed']

    dims = ['date','i','j']
    dates = pd.date_range(pgen['start_date'], pgen['end_date']).tolist()
    empty_array = np.ones((len(dates),) + np.shape(indices)) * np.nan

    ddict = {var: (dims, empty_array.copy()) for var in variables}

    for index in np.unique(indices):
        if np.isnan(index):
            break
        fp = pgen['forcing_file'].replace('[forcing_id]',str(int(index)))
        df = read_FMI_weather(pgen['start_date'],
                              pgen['end_date'],
                              sourcefile=fp)
        ix = np.where(indices==index)
        for var in variables:
            ddict[var][1][:,ix[0],ix[1]] = np.matmul(
                    df[var].values.reshape(len(dates),1),np.ones((1,len(ix[0]))))

    ds = xr.Dataset(ddict, coords={'date': dates})

    return ds

def _create_results(pgen, cmask, Nsteps):
    """
    Creates results dictionary to accumulate simulation results
    """
    i, j = np.shape(cmask)

    results = {}

    for var in pgen['variables']:

        var_shape = []
        var_name = var[0]

        if var_name.split('_')[0] != 'parameters':
            var_shape.append(Nsteps)
        if (var_name.split('_')[0] != 'forcing' or
            pgen['spatial_forcing'] == True):
            var_shape.append(i)
            var_shape.append(j)

        results[var_name] = np.full(var_shape, np.NAN)

    return results

def _append_results(group, step_results, results, step=None):
    """
    Adds results from each simulation steps to temporary results dictionary
    """

    results_keys = results.keys()

    for key in results_keys:
        var = key.split('_',1)[-1]

        if var in step_results and key.split('_',1)[0] == group:
            if group == 'forcing':
                res = step_results[var].values
            else:
                res = step_results[var]
            if step==None:
                results[key] = res
            else:
                results[key][step] = res

    return results

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', help='parameter folder', type=str)

    args = parser.parse_args()

    outputfile = driver(create_ncf=True, folder=args.folder)

    print(outputfile)