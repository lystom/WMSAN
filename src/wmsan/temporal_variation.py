import os.path
import cartopy
import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd


from netCDF4 import Dataset, date2num
from math import radians, log
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from obspy.geodetics.base import gps2dist_azimuth
from datetime import datetime
from scipy.interpolate import interp1d
from tqdm import tqdm
from calendar import monthrange
from pyproj import Geod

from wmsan.read_hs_p2l import read_hs, read_p2l

__author__ = "Reza D.D. Esfahani" # mod. by Lisa Tomasetto 07/2024
__copyright__ = "Copyright 2024, UGA"
__credits__ = ["Reza D.D. Esfahani"] 
__version__ = "0.1"
__maintainer__ = "Lisa Tomasetto"
__email__ = "lisa.tomasetto@univ-grenoble-alpes.fr"

def temporal_evolution(paths, dpt1, zlon, zlat, date_vec=[2020, [], [], []], extent=[-180, 180, -90, 90],parameters= [1/12, 1/2],  c_file = '../../data/C.nc', prefix = 'WW3-GLOB-30M', **kwargs):
    """Compute the temporal evolution of the seismic sources in a given region.

    Args:
        paths (list): A list containing the paths to the file with bathymetry data and the local path for WW3 data.
        dpt1 (xarray.DataArray): The bathymetry data.
        zlon (xarray.DataArray): The longitude values of the bathymetry data.
        zlat (xarray.DataArray): The latitude values of the bathymetry data.
        date_vec (list, optional): A list containing the year, month, day, and hour of the dates to compute the temporal evolution.
        extent (list, optional): A list containing the longitude and latitude extent of the region.
        parameters (list, optional): A list containing the minimum and maximum frequencies for integration (f1 and f2).
        c_file (str, optional): The path to the file with the amplification coefficient data.
        prefix (str, optional): The prefix for the WW3 data file name.

    Returns:
        time (list): A list of datetime objects representing the time of each computation.
        temporal_variation (ndarray): An array containing the temporal variation of force of the seismic sources.
    """
    
    file_bathy = paths[0]
    ww3_local_path = paths[1]

    # Constants
    radius = 6.371*1e6 # radius of the earth in meters
    lg10 = log(10) # log of 10
    res_mod = radians(0.5) # angular resolution of the model
    #
    f1 = parameters[0]
    f2 = parameters[1]
    

    ## Initialize variables
    if 'temporal_resolution' in kwargs:
        temporal_resolution = kwargs['temporal_resolution']
        if temporal_resolution == 'hourly':
            temporal_resol_hourly = True
        else:
            temporal_resol_hourly = False
            F_f1_concat = []
            
        if temporal_resolution == 'daily':
            temporal_resol_daily = True
        else:
            temporal_resol_daily = False

        if temporal_resolution == 'monthly':
            temporal_resol_monthly = True
             
        else:
            temporal_resol_monthly = False            
    else:
        # Default temporal resolution
        temporal_resol_hourly = True
        
    ## Adapt latitude and longitude to values in parameters file
    lon_min = extent[0]
    lon_max = extent[1]
    lat_min = extent[2]
    lat_max = extent[3]
    
    ## Work on the pacific ocean
    if lon_min > lon_max:
        ## work on the pacific ocean
        lon_min = ((360 + (lon_min % 360)) % 360)
        lon_max = ((360 + (lon_max % 360)) % 360)
        
    ## Open bathymetry
    dpt1 = dpt1.sel(latitude = slice(lat_min, lat_max), longitude = slice(lon_min, lon_max))
    zlat = zlat.sel(latitude = slice(lat_min, lat_max))
    zlon = zlon.sel(longitude = slice(lon_min, lon_max))
    ## Open Amplification Coefficient
    # check for refined bathymetry
    res_bathy = abs(zlon[1] - zlon[0])
    amplification_coeff = xr.open_dataarray(c_file)
    ## Work on the pacific ocean and amplification isn't already in the right coordinate system
    if (extent[0] > extent[1]) and max(amplification_coeff.longitude) < 180:
        amplification_coeff = amplification_coeff.assign_coords(longitude=((360 + (amplification_coeff.longitude % 360)) % 360))
        amplification_coeff = amplification_coeff.roll(longitude=int(len(amplification_coeff['longitude']) / 2),roll_coords=True)
    
    if res_bathy == 0.5:
        refined = False
    else:
        print("Refined bathymetry grid \n PLEASE RUN amplification_coefficients.ipynb before running this script")
        refined = True

    ## Surface Element
    msin = np.array([np.sin(np.pi/2 - np.radians(zlat))]).T
    ones = np.ones((1, len(zlon)))
    res_mod = radians(abs(zlat[1] - zlat[0]))
    dA = radius**2*res_mod**2*np.dot(msin,ones)
    
    ## Loop over dates
    YEAR = date_vec[0]
    MONTH = date_vec[1]
    DAY = date_vec[2]
    HOUR = date_vec[3]

    temporal_variation = []
    time = []

    for iyear in np.array([YEAR]):
        if isinstance(MONTH, int):
            MONTH = np.array([MONTH])
        elif not len(MONTH):
            MONTH = np.arange(1, 13)
        else:
            MONTH = np.array(MONTH)
        for imonth in MONTH:
            TOTAL_month = np.zeros(dpt1.shape)  # Initiate monthly source of Rayleigh wave matrix
            daymax = monthrange(iyear,imonth)[1]
            filename_p2l = '%s/%s_%d%02d_p2l.nc'%(ww3_local_path, prefix, iyear, imonth)
            print("File WW3 ", filename_p2l)
            try:
                day = np.array(DAY)
                if day[0] > day[-1]:
                    index = np.squeeze(np.argwhere(day==daymax))
                    if imonth == MONTH[0]:
                        day = day[:index+1]
                    elif imonth == MONTH[-1]:
                        index = np.squeeze(np.argwhere(day==monthrange(iyear, imonth-1)[1]))
                        day = day[index+1:]
                    else:
                        day = np.arange(1,daymax+1)
            except:
                try:
                    day = int(np.squeeze(day))
                except:
                    day= np.arange(1,(monthrange(iyear,imonth)[1])+1)
            for iday in day:
                if isinstance(HOUR, int):
                    HOUR = np.array([HOUR])
                elif not len(HOUR):
                    HOUR = np.arange(0,24,3)
                else:
                    HOUR = np.array(HOUR)
                for ih in HOUR:
                    
                    ## Open F_p3D 
                    (lati, longi, freq_ocean, p2l, unit1) = read_p2l(filename_p2l, [iyear, imonth, iday, ih], [extent[0], extent[1]], [lat_min, lat_max])
                    nf = len(freq_ocean)  # number of frequencies 
                    xfr = np.exp(np.log(freq_ocean[-1]/freq_ocean[0])/(nf-1))  # determines the xfr geometric progression factor
                    df = freq_ocean*0.5*(xfr-1/xfr)  # frequency interval in wave model times 2
                    freq_seismic = 2*freq_ocean  # ocean to seismic waves freq
                
                    ## Check units of the model, depends on version
                    if unit1 == 'log10(Pa2 m2 s+1E-12':
                        p2l = np.exp(lg10*p2l)  - (1e-12-1e-16)
                    elif unit1 == 'log10(m4s+0.01':
                        p2l = np.exp(lg10*p2l) - 0.009999
                    elif unit1 == 'log10(Pa2 m2 s+1E-12)':
                        p2l = np.exp(lg10*p2l)  - (1e-12-1e-16)
                        
                    ## Integral over a frequency band  
                    if f1 < f2:
                        index_freq = np.logical_and(f1 <= freq_seismic, freq_seismic <= f2)
                                        ## Single frequency
                    elif f1 == f2:
                        index_freq = np.squeeze(np.argmin(abs(freq-f1)))
                        print('unique frequency ', f1)
                        
                    ## Exception in parametrization of frequencies
                    else:
                        print('two frequencies with f2 < f1 were given')
                        
                    df = df[index_freq]
                    freq_seismic = freq_seismic[index_freq]
                    n_freq = len(df)
                    Fp = p2l.sel(frequency = freq_ocean[index_freq], latitude = slice(lat_min, lat_max), longitude = slice(lon_min, lon_max))
                    Fp = Fp.where(np.isfinite(Fp))
                    Fp = Fp.assign_coords({"freq": freq_seismic})
                    Fp = Fp.swap_dims({"frequency": "freq"})
                    Fp = Fp.drop_vars('frequency')
                    Fp = Fp.rename({"freq": "frequency"})
                    res_p2l = abs(Fp.latitude[1] - Fp.latitude[0])
                    if res_bathy != res_p2l:
                        # interpolate Fp
                        Fp = Fp.interp(longitude=zlon, latitude=zlat, method='linear')
                    ## Amplification coefficient
                    amplification_coeff = amplification_coeff.sel(frequency = freq_seismic, method='nearest')
                    amplification_coeff = amplification_coeff.sel(latitude = slice(lat_min, lat_max), longitude = slice(lon_min, lon_max))
                    amplification_coeff = amplification_coeff.reindex_like(Fp, method='nearest', tolerance=0.01)
                    F_f = Fp*amplification_coeff**2

                    ## Compute Equivalent Vertical Force
                    F = F_f
                    # Integrate over the frequency band
                    for ifq, fq in enumerate(freq_seismic.values):
                        df_fq = df[ifq]
                        F_fi = F_f[ifq, :, :]
                        F_fi *= dA
                        F[ifq, :, :] = F_fi*df_fq
                    # Force
                    F = 2*np.pi*np.sqrt(F.sum(axis = 0))

                    if temporal_resol_hourly:
                        time.append(datetime(iyear, imonth, iday, ih))
                        ## Concatenate temporal variation of sources + spatial averaging 
                        temporal_variation.append(np.nanmean(F))
                    else:
                        F_f1_concat.append(F)

                if temporal_resol_daily:
                    time.append(datetime(iyear, imonth, iday))
                    ## Concatenate temporal variation of  daily sources averaging + spatial averaging
                    temporal_variation.append(np.nanmean(F_f1_concat))
                    F_f1_concat = []

            if temporal_resol_monthly:
                # need to define the day! or chage datetime
                time.append(datetime(iyear, imonth, 15))
                ## Concatenate temporal variation of  monthly sources averaging + spatial averaging
                temporal_variation.append(np.nanmean(F_f1_concat))
                F_f1_concat = []

    return time, np.asarray(temporal_variation) 
