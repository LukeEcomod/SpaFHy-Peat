# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 12:58:54 2018

@author: khaahti
"""
import numpy as np
eps = np.finfo(float).eps
from scipy.interpolate import interp1d

apply_vectorized = np.vectorize(lambda f, x: f(x))

class SoilGrid(object):
    """
    Soil profile model for gridded uses (inputs as arrays)
    Simulates moss/organic layer with interception and evaporation,
    soil water storage, drainage to ditches and pond storage on top of soil.
    """
    def __init__(self, spara):
        """
        Initializes SoilProfile:
        Args:
            spara (dict):
                # scipy interpolation functions describing soil behavior
                'wtso_to_gwl'
                'gwl_to_wsto'
                'gwl_to_Ksat'
                # organic (moss) layer
                'org_depth': depth of organic top layer (m)
                'org_poros': porosity (-)
                'org_fc': field capacity (-)
                'org_rw': critical vol. moisture content (-) for decreasing phase in Ef
                'pond_storage_max': maximum pond depth [m]
                # drainage parameters
                'ditch_depth':ditch depth [m]
                'ditch_spacing': ditch spacing [m]
                'ditch_width': ditch width [m]
                # initial states
                'ground_water_level': groundwater depth [m]
                'org_sat': organic top layer saturation ratio (-)
                'pond_storage': initial pond depth at surface [m]
        """
        # top layer is interception storage, which capacity is depends on its depth [m]
        # and field capacity
        self.dz_top = spara['org_depth']  # depth, m3 m-3
        self.poros_top = spara['org_poros']  # porosity, m3 m-3
        self.fc_top = spara['org_fc']  # field capacity m3 m-3
        self.rw_top = spara['org_rw']  # ree parameter m3 m-3
        self.Wsto_top_max = self.fc_top * self.dz_top  # maximum storage m
        self.soiltype = spara['soiltype']

        # pond storage
        self.h_pond_max = spara['pond_storage_max']

        # interpolated functions for soil column ground water dpeth vs. water storage
        self.wsto_to_gwl = spara['wtso_to_gwl']
        self.gwl_to_wsto = spara['gwl_to_wsto']
        self.gwl_to_Ksat = spara['gwl_to_Ksat']

        self.Wsto_max = np.full_like(self.h_pond_max,0.0)
        for key, value in self.gwl_to_wsto.items():
            self.Wsto_max[self.soiltype == key] = value(0.0)

        self.gwl_to_rootmoist = spara['gwl_to_rootmoist']

        # drainage parameters
        self.ditch_depth = spara['ditch_depth']
        self.ditch_spacing = spara['ditch_spacing']
        self.depth_id = spara['depth_id']
        self.Ksat = np.full_like(self.h_pond_max, 0.0)

        # initialize state
        # soil profile
        self.gwl = spara['ground_water_level']
        self.Wsto = np.full_like(self.gwl, 0.0)
        self.rootmoist = np.full_like(self.gwl, 0.0)
        self.root_fc0 = np.full_like(self.gwl, 0.0)
        self.root_fc1 = np.full_like(self.gwl, 0.0)
        self.root_wp = np.full_like(self.gwl, 0.0)
        for key, value in self.gwl_to_wsto.items():
            self.Wsto[self.soiltype == key] = value(self.gwl[self.soiltype == key])
        for key, value in self.gwl_to_rootmoist.items():
            self.rootmoist[self.soiltype == key] = value(self.gwl[self.soiltype == key])
            self.root_fc0[self.soiltype == key] = value(-0.7 - 0.1)
            self.root_fc1[self.soiltype == key] = value(-1.2 - 0.1)
            self.root_wp[self.soiltype == key] = value(-150.0 - 0.1)

        self.Rew = 1.0

        # toplayer storage and relative conductance for evaporation
        self.Wsto_top = self.Wsto_top_max * spara['org_sat']

        self.Wliq_top = self.poros_top * self.Wsto_top / self.Wsto_top_max
        self.Ree = np.maximum(0.0, np.minimum(
                0.98*self.Wliq_top / self.rw_top, 1.0)) # relative evaporation rate (-)
        # pond storage
        self.h_pond = spara['pond_storage']

    def watbal(self, dt=1, rr=0.0, tr=0.0, evap=0.0):

        r""" Solves soil water storage in column assuming hydrostatic equilibrium.

        Args:
            dt (float): solution timestep [s]
            rr (float/array): potential infiltration [m]
            tr (float/array): transpiration from root zone [m]
            evap (float/array): evaporation from top layer [m]
        Returns:
            results (dict)

        """

        state0 = self.Wsto + self.Wsto_top + self.h_pond + rr
        rr += self.h_pond
        self.h_pond = 0.0

        # top layer interception & water balance
        interc = np.maximum(0.0, (self.Wsto_top_max - self.Wsto_top))\
                    * (1.0 - np.exp(-(rr / self.Wsto_top_max)))
        rr -= interc  # to soil profile
        self.Wsto_top += interc
        evap = np.minimum(evap, self.Wsto_top)
        self.Wsto_top -= evap

        # drainage [m]
        for i in range(len(self.gwl_to_Ksat)):
            self.Ksat[self.depth_id == i] = self.gwl_to_Ksat[i](self.gwl[self.depth_id == i])  # [m s-1]
        Hdr = np.minimum(np.maximum(0, self.gwl + self.ditch_depth), self.ditch_depth)
        drain = 4 * self.Ksat * Hdr**2 / (self.ditch_spacing**2) * dt  # [m]

        """ soil column water balance """
        # net flow to soil profile during dt
        Qin = rr - drain - tr
        # airvolume available in soil profile after previous time step
        Airvol = np.maximum(0.0, self.Wsto_max - self.Wsto)

        Qin = np.minimum(Qin, Airvol)
        self.Wsto += Qin

        infiltration = Qin + drain + tr

        # inflow excess
        exfil = rr - infiltration
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
        # soil profile
        # ground water depth corresponding to Wsto
        for key, value in self.wsto_to_gwl.items():
            self.gwl[self.soiltype == key] = value(self.Wsto[self.soiltype == key])

        # organic top layer; maximum that can be hold is Fc
        self.Wliq_top = self.fc_top * self.Wsto_top / self.Wsto_top_max
        self.Ree = np.maximum(0.0, np.minimum(0.98*self.Wliq_top / self.rw_top, 1.0))

        for key, value in self.gwl_to_rootmoist.items():
            self.rootmoist[self.soiltype == key] = value(self.gwl[self.soiltype == key])

        # Koivusalo et al. 2008 HESS without wet side limit
        self.Rew = np.where(self.rootmoist > self.root_fc1,
                            np.minimum(1.0, 0.5*(1 + (self.rootmoist - self.root_fc1)/(self.root_fc0 - self.root_fc1))),
                            np.maximum(0.0, 0.5*(self.rootmoist - self.root_wp)/(self.root_fc1 - self.root_wp))
                            )

        # mass balance error [m]
        mbe = (state0  - self.Wsto_top - self.Wsto - self.h_pond -
               drain - tr - surface_runoff - evap)

        results = {
                'pond_storage': self.h_pond,  # [m]
                'ground_water_level': self.gwl,  # [m]
                'infiltration': infiltration * 1e3,  # [mm d-1]
                'surface_runoff': surface_runoff * 1e3,  # [mm d-1]
                'evaporation': evap * 1e3,  # [mm d-1]
                'drainage': drain * 1e3,  # [mm d-1]
                'moisture_top': self.Wliq_top,  # [m3 m-3]
                'water_closure':  mbe,  # [mm d-1]
                'transpiration_limitation': self.Rew,  # [-]
                'rootzone_moisture': self.rootmoist,  # [m3 m-3]
                }

        return results

def gwl_Wsto(z, pF, root=False):
    r""" Forms interpolated function for soil column ground water dpeth, < 0 [m], as a
    function of water storage [m] and vice versa

    Args:
        pF (dict of arrays):
            'ThetaS' saturated water content [m\ :sup:`3` m\ :sup:`-3`\ ]
            'ThetaR' residual water content [m\ :sup:`3` m\ :sup:`-3`\ ]
            'alpha' air entry suction [cm\ :sup:`-1`]
            'n' pore size distribution [-]
        dz (np.arrays): soil conpartment thichness, node in center [m]
    Returns:
        (dict):
            'to_gwl': interpolated function for gwl(Wsto)
            'to_wsto': interpolated function for Wsto(gwl)
    """

    z = np.array(z)
    dz = abs(z)
    dz[1:] = z[:-1] - z[1:]
    z_mid = dz / 2 - np.cumsum(dz)

    # --------- connection between gwl and water storage------------
    # gwl from ground surface gwl = 0 to gwl = -5
    gwl = np.arange(0.0, -5, -1e-3)
    gwl[-1] = -150
    # solve water storage corresponding to gwls
    Wsto = [sum(h_to_cellmoist(pF, g - z_mid, dz) * dz) for g in gwl]

    if root:
        Wsto = Wsto/sum(dz)

    # interpolate functions
    WstoToGwl = interp1d(np.array(Wsto), np.array(gwl), fill_value='extrapolate')
    GwlToWsto = interp1d(np.array(gwl), np.array(Wsto), fill_value='extrapolate')

    del gwl, Wsto

    if root:
        return {'to_rootmoist': GwlToWsto}
    else:
        return {'to_gwl': WstoToGwl, 'to_wsto': GwlToWsto}

def h_to_cellmoist(pF, h, dz):
    r""" Cell moisture based on vanGenuchten-Mualem soil water retention model.
    Partly saturated cells calculated as thickness weigthed average of
    saturated and unsaturated parts.

    Args:
        pF (dict):
            'ThetaS' (array): saturated water content [m\ :sup:`3` m\ :sup:`-3`\ ]
            'ThetaR' (array): residual water content [m\ :sup:`3` m\ :sup:`-3`\ ]
            'alpha' (array): air entry suction [cm\ :sup:`-1`]
            'n' (array): pore size distribution [-]
        h (array): pressure head [m]
        dz (array): soil conpartment thichness, node in center [m]
    Returns:
        theta (array): volumetric water content of cell [m\ :sup:`3` m\ :sup:`-3`\ ]

    Kersti Haahti, Luke 8/1/2018
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

def gwl_Ksat(z, Ksat, DitchDepth):
    r""" Forms interpolated function for drainage vs gwl
    """
    z = np.array(z)
    dz = abs(z)
    dz[1:] = z[:-1] - z[1:]

    # --------- connection between gwl and drainage------------
    # gwl from ground surface gwl = 0 to gwl = -5
    gwl = np.arange(0.0, -DitchDepth-0.1, -1e-2)
    gwl[-1] = -5.0

    # solve water storage corresponding to gwls
    Ka = [Ksat_layer(dz, Ksat, g, DitchDepth) for g in gwl]

    # interpolate functions
    GwlToKsat = interp1d(np.array(gwl), np.array(Ka), fill_value='extrapolate')

    return GwlToKsat

def Ksat_layer(dz, Ksat, gwl, DitchDepth):
    r""" Calculates drainage to ditch using Hooghoud's drainage equation,
    accounts for drainage from saturated layers above and below ditch bottom.

    Args:
       dz (array):  soil conpartment thichness, node in center [m]
       Ksat (array): horizontal saturated hydr. cond. [ms-1]
       gwl (float): ground water level below surface, <0 [m]
       DitchDepth (float): depth of drainage ditch bottom, >0 [m]

    Returns:
       Qz_drain (array): drainage from each soil layer [m3 m-3 s-1]

    Reference:
       Follows Koivusalo, Lauren et al. FEMMA -document. Ref: El-Sadek et al., 2001.
       J. Irrig.& Drainage Engineering.

    Samuli Launiainen, Metla 3.11.2014.; converted to Python 14.9.2016
    Kersti Haahti, 29.12.2017. Code checked, small corrections
    """
    z = dz / 2 - np.cumsum(dz)
    Ka = 0.0

    Hdr = min(max(0, gwl + DitchDepth), DitchDepth)  # depth of saturated layer above ditch bottom

    """ drainage from saturated layers above ditch base """
    # layers above ditch bottom where drainage is possible
    ix = np.intersect1d(np.where((z - dz / 2)- gwl < 0), np.where(z + dz / 2 > -DitchDepth))

    if Hdr > 0:
        # saturated layer thickness [m]
        dz_sat = np.minimum(np.maximum(gwl - (z - dz / 2), 0), dz)
        # transmissivity of layers  [m2 s-1]
        Trans = Ksat * dz_sat

        """ drainage from saturated layers above ditch base """
        # layers above ditch bottom where drainage is possible
        ix = np.intersect1d(np.where((z - dz / 2)- gwl < 0), np.where(z + dz / 2 > -DitchDepth))

        if ix.size > 0:
            dz_sat[ix[-1]] = dz_sat[ix[-1]] + (z[ix][-1] - dz[ix][-1] / 2 + DitchDepth)
            if abs(sum(dz_sat[ix]) - Hdr) > eps:
                print(sum(dz_sat[ix]), Hdr, DitchDepth, gwl)
            Trans[ix[-1]] = Ksat[ix[-1]] * dz_sat[ix[-1]]
            Ka = sum(Trans[ix]) / sum(dz_sat[ix])  # effective hydraulic conductivity ms-1

    return Ka

nan_function = interp1d(np.array([np.nan, np.nan]),
                        np.array([np.nan, np.nan]),
                        fill_value='extrapolate')