# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 12:58:54 2018

@author: khaahti
"""
import numpy as np
eps = np.finfo(float).eps
from scipy.interpolate import interp1d

class Soil(object):
    """
    Soil profile model
    Simulates moss/organic layer with interception and evaporation,
    soil water storage, drainage to ditches and pond storage on top of soil.
    """
    def __init__(self, spara):
        """
        Initializes Soil:
        Args:
            spara (dict):
                # organic (moss) layer
                'org_depth': depth of organic top layer [m]
                'org_poros': porosity [m3 m-3]
                'org_fc': field capacity [m3 m-3]
                'org_rw': critical vol. moisture content [m3 m-3] for decreasing phase in Ef
                'pond_storage_max': maximum pond depth [m]
                # layerwise soil para (lists)
                'z': position of soil layer bottom below soil surface <0 [m]
                'pF' (dict):  # vanGenuchten water retention parameters
                    'ThetaS': porosity [m3 m-3]
                    'ThetaR': residual water content [m3 m-3]
                    'alpha': air entry suction [cm-1]
                    'n': pore size distribution [-]
                'saturated_conductivity': saturated soil hydraulic conductivity [m s-1]
                # drainage parameters
                'ditch_depth':ditch depth [m]
                'ditch_spacing': ditch spacing [m]
                # initial states
                'ground_water_level': groundwater depth [m]
                'org_sat': organic top layer saturation ratio [-]
                'pond_storage': initial pond depth at surface [m]
        """
        # top layer is an interception storage, which capacity is depends on its depth and field capacity
        self.poros_top = spara['org_poros']  # porosity, m3 m-3
        self.fc_top = spara['org_fc']  # field capacity m3 m-3
        self.rw_top = spara['org_rw']  # ree parameter m3 m-3
        self.Wsto_top_max = spara['org_fc'] * spara['org_depth']  # maximum storage m

        # pond storage [m]
        self.h_pond_max = spara['pond_storage_max']

        # drainage parameters [m]
        self.ditch_depth = spara['ditch_depth']
        self.ditch_spacing = spara['ditch_spacing']

        # interpolated functions [m]
        # soil column ground water dpeth vs. water storage
        self.wsto_to_gwl, self.gwl_to_wsto = gwl_Wsto(spara['z'], spara['pF'])
        # Saturated hydraulooc conductivity above ditch depth
        self.gwl_to_Ksat = gwl_Ksat(spara['z'], spara['saturated_conductivity'], spara['ditch_depth'])

        # maximum water storage [m] of soil column (i.e. gwl=0)
        self.Wsto_max = self.gwl_to_wsto(0.0)

        # initialize state
        # ground water level [m]
        self.gwl = spara['ground_water_level']
        # water storage in soil profile [m]
        self.Wsto = self.gwl_to_wsto(self.gwl)
        # toplayer water storage [m] 
        self.Wsto_top = self.Wsto_top_max * spara['org_sat']
        # relative conductance for evaporation
        self.Wliq_top = self.poros_top * self.Wsto_top / self.Wsto_top_max
        self.Ree = np.maximum(0.0, np.minimum(
                0.98*self.Wliq_top / self.rw_top, 1.0)) # relative evaporation rate (-)
        # pond storage
        self.h_pond = spara['pond_storage']

    def watbal(self, dt=86400, pot_inf=0.0, tr=0.0, evap=0.0):

        r""" Solves water balance assuming hydrostatic equilibrium in soil profile.

        Args:
            dt (float): solution timestep [s]
            pot_inf (float): potential infiltration during dt [m]
            tr (float): transpiration from root zone during dt [m]
            evap (float): evaporation from top layer during dt [m]
        Returns:
            results (dict)

        """

        initial_state = self.Wsto + self.Wsto_top + self.h_pond + pot_inf

        # pond starage is added to potential infiltration
        pot_inf += self.h_pond
        self.h_pond = 0.0

        # top layer interception: asymptotic approach of saturation
        interc = np.maximum(0.0, (self.Wsto_top_max - self.Wsto_top))\
                    * (1.0 - np.exp(-(pot_inf / self.Wsto_top_max)))
        # potential infiltration to soil profile
        pot_inf -= interc  
        # update top layer water storage
        self.Wsto_top += interc
        # evaporation from top layer limited by top layer water storage
        evap = np.minimum(evap, self.Wsto_top)
        # update top layer water storage
        self.Wsto_top -= evap

        # drainage to ditches
        # saturated hydraulic conductivity of soil profile above ditch depth [m s-1]
        self.Ksat = self.gwl_to_Ksat(self.gwl)
        # saturated soil depth above ditch depth [m]
        Hdr = np.minimum(np.maximum(0, self.gwl + self.ditch_depth), self.ditch_depth)
        # drainage [m]: Hooghoudt (1940) equation (neglecting drainage from below the ditches)
        drain = 4 * self.Ksat * Hdr**2 / (self.ditch_spacing**2) * dt

        """ soil column water balance """
        # potential net flow to soil profile during dt [m]
        Qin = pot_inf - drain - tr
        # airvolume available in soil profile after previous time step [m]
        Airvol = np.maximum(0.0, self.Wsto_max - self.Wsto)

        # actual net flow to soil profile during dt limited by airvolume in column [m]
        Qin = np.minimum(Qin, Airvol)
        self.Wsto += Qin

        # actual infiltration [m]
        infiltration = Qin + drain + tr

        # inflow excess [m], water not fitting in soil due to saturation
        exfil = pot_inf - infiltration
        # first fill top layer
        # water that can fit in top layer
        to_top_layer = np.minimum(exfil, self.Wsto_top_max - self.Wsto_top)
        self.Wsto_top += to_top_layer
        # then pond storage
        to_pond = np.minimum(exfil - to_top_layer, self.h_pond_max - self.h_pond)
        self.h_pond += to_pond
        # and route remaining to surface runoff
        surface_runoff = exfil - to_top_layer - to_pond

        """ update state """
        # soil profile ground water level [m]
        self.gwl = self.wsto_to_gwl(self.Wsto)

        # organic top layer; maximum that can be hold is Fc
        self.Wliq_top = self.fc_top * self.Wsto_top / self.Wsto_top_max
        self.Ree = np.maximum(0.0, np.minimum(0.98*self.Wliq_top / self.rw_top, 1.0))

        # mass balance error [m]
        mbe = (initial_state  - self.Wsto_top - self.Wsto - self.h_pond -
               drain - tr - surface_runoff - evap)

        results = {
                'pond_storage': self.h_pond,  # [m]
                'ground_water_level': self.gwl,  # [m]
                'infiltration': infiltration * 1e3,  # [mm]
                'surface_runoff': surface_runoff * 1e3,  # [mm]
                'evaporation': evap * 1e3,  # [mm]
                'drainage': drain * 1e3,  # [mm]
                'moisture_top': self.Wliq_top,  # [m3 m-3]
                'water_closure':  mbe,  # [mm]
                }

        return results

def gwl_Wsto(z, pF, grid_step=0.01):
    r""" Forms interpolated function for soil column ground water dpeth, <0 [m], as a
    function of water storage [m] and vice versa

    Args:
        z (list): position of soil layer bottom below soil surface <0 [m]
        pF (dict of lists):
            'ThetaS' saturated water content [m3 m-3]
            'ThetaR' residual water content [m3 m-3]
            'alpha' air entry suction [cm-1]
            'n' pore size distribution [-]
    Returns:
        WstoToGwl: interpolated function for gwl(Wsto)
        GwlToWsto: interpolated function for Wsto(gwl)
    """

    # finer grid (J-P's code)
    z_fine= np.arange(0, min(z), -grid_step)
    dz_fine = z_fine*0.0 + grid_step
    z_mid_fine = dz_fine / 2 - np.cumsum(dz_fine)

    # pF parameter arrays for finer grid
    ix = np.zeros(len(z_fine))
    for depth in z:
        # below makes sure floating point precision doesnt mess with the ix
        ix += np.where((z_fine < depth) & ~np.isclose(z_fine, depth, atol=1e-9), 1, 0)
    pF_fine={}
    for key in pF.keys():
        pp = []
        for i in range(len(z_fine)):
            pp.append(pF[key][int(ix[i])])
        pF_fine.update({key: np.array(pp)})

    # fine resolution array of gwl values for which corresponding water storage will be determined
    gwl = np.arange(0.0, -5, -grid_step)
    gwl[-1] = -150

    # solve water storage corresponding to gwls
    Wsto = [sum(h_to_cellmoist(pF_fine, g - z_mid_fine, dz_fine) * dz_fine) for g in gwl]

    # interpolate functions
    WstoToGwl = interp1d(np.array(Wsto), np.array(gwl), fill_value='extrapolate')
    GwlToWsto = interp1d(np.array(gwl), np.array(Wsto), fill_value='extrapolate')

    return WstoToGwl, GwlToWsto

def h_to_cellmoist(pF, h, dz):
    r""" Cell moisture based on vanGenuchten-Mualem soil water retention model.
    Partly saturated cells calculated as thickness weigthed average of
    saturated and unsaturated parts.

    ! cell thickness should be small, because no interpolation is done here.

    Args:
        pF (dict):
            'ThetaS' saturated water content [m3 m-3]
            'ThetaR' residual water content [m3 m-3]
            'alpha' air entry suction [cm-1]
            'n' pore size distribution [-]
        h (array): pressure head [m]
        dz (array): soil conpartment thichness, node in center [m]
    Returns:
        theta (array): volumetric water content of cell [m3 m-3]
    """

    # water retention parameters
    Ts = np.array(pF['ThetaS'])
    Tr = np.array(pF['ThetaR'])
    alfa = np.array(pF['alpha'])
    n = np.array(pF['n'])
    m = 1.0 - np.divide(1.0, n)

    # moisture based on cell center head
    x = np.minimum(h, 0)
    theta = Tr + (Ts - Tr) / (1 + abs(alfa * 100 * x)**n)**m

    # correct moisture of partly saturated cells
    ix = np.where(abs(h) < dz/2)
    if len(Ts) == 1:
        ixx = 0
    else:
        ixx = ix
    # moisture of unsaturated part
    x[ix] = -(dz[ix]/2 - h[ix]) / 2
    theta[ix] = Tr[ixx] + (Ts[ixx] - Tr[ixx]) / (1 + abs(alfa[ixx] * 100 * x[ix])**n[ixx])**m[ixx]
    # total moisture as weighted average
    theta[ix] = (theta[ix] * (dz[ix]/2 - h[ix]) + Ts[ixx] * (dz[ix]/2 + h[ix])) / (dz[ix])

    return theta

def gwl_Ksat(z, Ksat, DitchDepth, grid_step=0.01):
    r""" Forms interpolated function for hydraulic conductivity of saturated layer 
    above ditch bottom vs gwl
    """

    # gwl from soil surface gwl = 0 to gwl = -150 (finer resolution until ditch bottom)
    gwl = np.arange(0.0, - DitchDepth - 0.1, - grid_step)
    gwl[-1] = -150

    # solve water storage corresponding to gwls
    Ka = [Ksat_layer(z, Ksat, g, DitchDepth) for g in gwl]

    # interpolate functions
    GwlToKsat = interp1d(np.array(gwl), np.array(Ka), fill_value='extrapolate')

    return GwlToKsat

def Ksat_layer(z, Ksat, gwl, DitchDepth):
    r""" Calculates hydraulic conductivity of saturated layer 
    above ditch bottom.

    Args:
       dz (array): soil layer thichnesses, node in center [m]
       Ksat (array): horizontal saturated hydr. cond. of soil layers [m s-1]
       gwl (float): ground water level below soil surface, <0 [m]
       DitchDepth (float): depth of drainage ditch bottom, >0 [m]

    Returns:
       Ka (array): hydraulic conductivity of saturated layer above ditch bottom [m s-1]
    """

    z = np.array(z)
    dz = abs(z)
    dz[1:] = z[:-1] - z[1:]

    # each layers thickness above ditch bottom [m] (zero when layer completely below ditch bottom)
    dz_dd = np.minimum(dz, np.maximum(DitchDepth + z + dz, 0))
    # each layers saturated thickness above ditch bottom [m]
    dz_sat = np.minimum(np.maximum(gwl - np.maximum(-DitchDepth, z), 0), dz_dd)
    if sum(dz_sat) > 0:
        # transmissivity of layers  [m2 s-1]
        Trans = Ksat * dz_sat
        # effective hydraulic conductivity ms-1
        Ka = sum(Trans) / sum(dz_sat) 
    else:
        Ka = 0.0

    return Ka
