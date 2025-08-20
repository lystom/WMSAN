#!/usr/bin/env python3

# Preamble
#__author__ = "Lisa Tomasetto"
#__copyright__ = "Copyright 2024, UGA"
#__credits__ = ["Lisa Tomasetto"]
#__version__ = "2025.0.0"
#__maintainer__ = "Lisa Tomasetto"
#__email__ = "lisa.tomasetto@univ-grenoble-alpes.fr"
#__status__ = "Production"

"""This set of functions reads the netcdf files _hs.nc and _p2l.nc produced by WW3.

It contains four functions read_WWNC, read_WWNCf, read_hs and read_p2l:

- `read_WWNC (file_path, time_vect, lon1, lat1)`: 
    Read netcdf _hs.nc file and return a matrix with dimension lon x lat of significant wave height in meters.

- `read_WWNCf (file_path, time_vect, lon1, lat1)`: 
    Read netcdf _p2l.nc file and return latitude, longitude, frequency,
    p2l data which is the base 10 logarithm of power specral density of equivalent surface pressure and the units of p2l.

- `read_hs (file_path, time_vect, lon1, lat1)`:
    Read netcdf significant wave height (_hs.nc) file and
    return latitude, longitude, frequenc, hs data. Uses xarray.

- `read_p2l (file_path, time_vect, lon1, lat1)`: 
    Read netcdf _p2l.nc file and return latitude, longitude, frequency, p2l data the units of p2l. Uses xarray.
"""
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import argparse

from datetime import date
from datetime import datetime
from netCDF4 import Dataset
import netCDF4 as nc
import os, fnmatch, glob, shutil

def read_WWNC(file_path, time_vect, lon1, lat1):
    """Read netcdf _hs.nc file and return a matrix with dimension lon x lat of significant height of wind and swell waves in meters.
    
    Examples:
        >>> file_path = '../data/ftp.ifremer.fr/ifremer/ww3/HINDCAST/SISMO/GLOBAL05_2006_REF102040/data/WW3-GLOB-30M_200609_hs.nc'
        >>> year = 2006
        >>> month = 9
        >>> day = 4
        >>> hour = 12
        >>> lon1 = []
        >>> lat1 = []
        >>> time_vect = [year, month, day, hour]
        >>> A = read_WWNC(file_path, time_vect, [], [])
        >>> print(A)
    
    Args:
        file_path (str): path of the netcdf file
        time_vect (list): [year, month, day, hour]
        lon1 (list): [lon_min, lon_max]
        lat1 (list): [lat_min, lat_max]
    
    Returns:
        hs (numpy.ndarray): matrix with dimension lon x lat of significant height of wind and swell waves in meters 
    """

    ## open file
    f = Dataset(file_path, mode='r')

    #Coordinates variables
    lon = f.variables['longitude'][:]
    nx = len(lon)
    lat = f.variables['latitude'][:]
    ny = len(lat)

    #time variables
    # We assume that the date reference is 1 Jan 1990, this is normally written in the time attributes
    time0 = date.toordinal(date(1990, 1, 1)) + 366
    time = f.variables['time'][:] + time0
    nt = len(time)

    # Define indices for data subset
    if time_vect != []:  # if a date is specified select the closest time in data
        date1 = 366 + date.toordinal(date(time_vect[0], time_vect[1], time_vect[2])) + time_vect[3]/24
        kk = np.argmin(abs(time-date1))
        time = time[kk]
        KK = kk
        nk = 1
    else:
        KK = 0
        nk = nt

    if lon1 != []:  # if a longitude is specified select the closest longitude in data
        ii = np.argmin(abs(lon-lon1))
        lon = lon[ii]
        II = ii
        ni = 1
    else:
        II = 0
        ni = nx

    if lat1 != []:  # if a latitude is specified select the closest latitude in data
        jj = np.argmin(abs(lat-lat1))
        lat = lat[jj]
        JJ = jj
        nj = 1
    else:
        JJ = 0
        nj = ny

    # Extract data

    hs = f.variables['hs'][KK:KK+nk][JJ:JJ+nj][II:II+ni]
   
    f.close()
    hs = (np.squeeze(hs.filled(fill_value=np.nan))).T  # replace masked values in data by NaNs
    return hs

def read_WWNCf(file_path, time_vect, lon1, lat1):
    """Read netcdf _p2l.nc file and return latitude, longitude, frequenc, p2l data which is the base 10 logarithm of power specral density of equivalent surface pressure and the units of p2l.
    
    Examples:
        >>> file_path = '../data/ftp.ifremer.fr/ifremer/ww3/HINDCAST/SISMO/GLOBAL05_2006_REF102040/WW3-GLOB-30M_200609_p2l.nc'
        >>> year = 2006
        >>> month = 9
        >>> day = 4
        >>> hour = 12
        >>> lon1 = []
        >>> lat1 = []
        >>> time_vect = [year, month, day, hour]
        >>> (lat, lon, freq, p2l, unit1) = read_WWNCf(file_path, time_vect, [], [])
    
    Args:
        file_path (str): path of the netcdf file
        time_vect (list): [year, month, day, hour]
        lon1 (list): [lon_min, lon_max]
        lat1 (list): [lat_min, lat_max]
        
    Returns:
        lat (numpy.ndarray): latitude
        lon (numpy.ndarray): longitude
        freq (numpy.ndarray): frequency
        p2l (numpy.ndarray): p2l
        unit1 (str): unit of p2l   
    """

    ## open file
    f = Dataset(file_path, mode='r')

    #Coordinates variables
    lon = f.variables['longitude'][:]
    nx = len(lon)
    lat = f.variables['latitude'][:]
    ny = len(lat)

    #time variables
    # We assume that the date reference is 1 Jan 1990, this is normally written in the time attributes
    time0 = date.toordinal(date(1990, 1, 1)) + 366
    time = f.variables['time'][:] + time0
    nt = len(time)

    #Frequency variables
    freq = f.variables['f'][:]
    nf = len(freq)

    # Define indices for data subset
    if time_vect != []:  # if a date is specified select the closest time in data
        date1 = 366 + date.toordinal(date(time_vect[0], time_vect[1], time_vect[2])) + time_vect[3]/24
        kk = np.argmin(abs(time-date1))
        time = time[kk]
        KK = kk
        nk = 1
    else:
        KK = 0
        nk = nt

    if lon1 != []:  # if a longitude is specified select the closest longitude in data
        ii = np.argmin(abs(lon-lon1))
        lon = lon[ii]
        II = ii
        ni = 1
    else:
        II = 0
        ni = nx

    if lat1 != []:  # if a latitude is specified select the closest latitude in data
        jj = np.argmin(abs(lat-lat1))
        lat = lat[jj]
        JJ = jj
        nj = 1
    else:
        JJ = 0
        nj = ny

    LL = 0
    nl = nf

    # Extract data
    p2l = f.variables['p2l'][KK:KK+nk][LL:LL+nl][JJ:JJ+nj][II:II+ni]

 # Check units and convert to normal units in case of log scales
    unit1 = f.variables['p2l'].units
    scale = f.variables['p2l'].scale_factor
    p2l = (np.squeeze(p2l.filled(fill_value=np.nan))).T  # replace masked values in data by NaNs
    f.close()
    return lat, lon, freq, p2l, unit1, scale

def read_hs(file_path, time_vect, lon1 = (-180, 180), lat1 = (-90, 90)):
    """Read netcdf significant wave height (_hs.nc) file and return latitude, longitude, frequenc, hs data units are m. The output is an xarray of shape (lat, lon)
    
    Examples:
        >>> file_path = '../data/2006/WW3-GLOB-30M_200609_hs.nc'
        >>> year = 2006
        >>> month = 9
        >>> day = 4
        >>> hour = 12
        >>> lon1 = [-180, 180]
        >>> lat1 = [-90, 90]
        >>> time_vect = [year, month, day, hour]
        >>> hs = read_hs(file_path, time_vect, lon1, lat1)
        >>> fig = plt.figure(figsize=(9,6))
        >>> ax = plt.axes(projection=ccrs.Robinson())
        >>> ax.coastlines()
        >>> ax.gridlines()
        >>> hs.plot(ax=ax, transform=ccrs.PlateCarree(), cbar_kwargs={'shrink': 0.4})
        >>> plt.show()
    
    Args:
        file_path (str): path to _hs.nc file.
        time_vect (list): [year, month, day, hour] of the time step.
        lon1 (tuple, optional): (lon_min, lon_max) of the spatial extent.
        lat1 (tuple, optional): (lat_min, lat_max) of the spatial extent.
        
    Returns:
        hs (xarray.DataArray): array of shape (lat, lon) with significant height of wind and swell waves in meters
    """
    
    # open file
    try:
        ds = xr.open_dataset(file_path)
    except:
        print('Error opening file, please check that the file exists\n')
        print(file_path)
        exit()

    # datetime
    year = time_vect[0]
    month = time_vect[1]
    day = time_vect[2]
    hour = time_vect[3]
    timestep = datetime(year, month, day, hour)
    
    # spatial extent
    ## Check latitude within -90 90
    if np.abs(lat1[0]) > 90 or np.abs(lat1[1]) > 90:
        print("Latitude not correct, absolute value > 90")
        return 
    # check longitude
    if lon1[0] > lon1[1]:
        ## working near pacific
        print("pacific")
        ds = ds.assign_coords(longitude=((360 + (ds.longitude % 360)) % 360))
        ds = ds.roll(longitude=int(len(ds['longitude']) / 2),roll_coords=True)
        lon1[0] = ((360 + (lon1[0] % 360)) % 360)
        lon1[1] = ((360 + (lon1[1] % 360)) % 360)
    lon = ds.longitude[np.logical_and(ds.longitude >= lon1[0], ds.longitude <= lon1[1])]
    lat = ds.latitude[np.logical_and(ds.latitude >= lat1[0], ds.latitude <= lat1[1])]
    
    # extract data
    hs = ds.hs.sel(time= timestep, longitude = lon, latitude= lat, method = 'nearest')
    ds.close()
    return hs
    
def read_p2l(file_path, time_vect, lon1 = (-180, 180), lat1 = (-90, 90)):
    """Read netcdf _p2l.nc file and return latitude, longitude, frequency, p2l data which is the base 10 logarithm of power specral density of equivalent surface pressure and the units of p2l.
 
    Examples:
        >>> file_path = '../data/2006/WW3-GLOB-30M_200609_p2l.nc'
        >>> year = 2006
        >>> month = 9
        >>> day = 4
        >>> hour = 12
        >>> lon1 = []
        >>> lat1 = []
        >>> time_vect = [year, month, day, hour]
        >>> (lat, lon, freq, p2l, unit1) = read_p2l(file_path, time_vect, lon1, lat1)
    
    Args:
        file_path (str): path of the netcdf file _p2l.nc
        time_vect (list): [year, month, day, hour] of the time step.
        lon1 (tuple, optional): (lon_min, lon_max) of the spatial extent.
        lat1 (tuple, optional): (lat_min, lat_max) of the spatial extent.
    
    Returns:
        lat (numpy.ndarray): latitude
        lon (numpy.ndarray): longitude
        freq (numpy.ndarray): frequency
        p2l (numpy.ndarray): p2l
        unit1 (str): unit of p2l 
    """
    
    # open file
    try:
        ds = xr.open_dataset(file_path)
    except:
        print('Error opening file, please check that the file exists\n')
        print(file_path)
        exit()
    ds.assign_coords(f=("f", ds.f.data))
    ds = ds.rename({'f':'frequency'})
    # datetime
    year = time_vect[0]
    month = time_vect[1]
    day = time_vect[2]
    hour = time_vect[3]
    timestep = datetime(year, month, day, hour)

    # spatial extent
    if np.abs(lat1[0]) > 90 or np.abs(lat1[1]) > 90:
        print("Latitude not correct, absolute value > 90")
        return
    if lon1[0] > lon1[1]:
        print("pacific")
        ds = ds.assign_coords(longitude=((360 + (ds.longitude % 360)) % 360))
        ds = ds.roll(longitude=int(len(ds['longitude']) / 2),roll_coords=True)
        lon1[0] = ((360 + (lon1[0] % 360)) % 360)
        lon1[1] = ((360 + (lon1[1] % 360)) % 360)
        lon = ds.longitude[np.logical_and(ds.longitude >= lon1[0], ds.longitude <= lon1[1])]
    lon = ds.longitude[np.logical_and(ds.longitude >= lon1[0], ds.longitude <= lon1[1])]
    lat = ds.latitude[np.logical_and(ds.latitude >= lat1[0], ds.latitude <= lat1[1])]
    # frequency range
    freq = ds.frequency
    # extract data
    p2l = ds.p2l.sel(time= timestep, frequency = freq, longitude = lon, latitude= lat, method = 'nearest')

    # units
    unit1 = ds.p2l.units
    ds.close()
    return lat, lon, freq, p2l, unit1

def read_p2l_from_url(url, time_vect, lon1 = (-180, 180), lat1 = (-90, 90)):
    try:
        nc_ds.close()
        extract_ds.close()
    except:
        pass
    (lat_min, lat_max) = lat1
    (lon_min, lon_max) = lon1
    #load netcdf from url as netCDF4 dataset
    ncfile = nc.Dataset(url+'#mode=bytes')
    #load netCDF4 dataset as Xarray dataset
    nc_ds = xr.open_dataset(xr.backends.NetCDF4DataStore(ncfile))
    # extract from Jan 1st to Jan 5th included
    extract_ds=nc_ds.sel(time=datetime(time_vect[0], time_vect[1], time_vect[2], time_vect[3]))
    extract_ds = extract_ds.sel(latitude=slice(lat_min, lat_max), longitude=slice(lon_min, lon_max))
    extract_ds = extract_ds.rename({'f':'frequency'})

    return extract_ds.latitude, extract_ds.longitude, extract_ds.frequency, extract_ds.p2l, 'log10(Pa2 m2 s+1E-12)'

def read_hs_from_url(url, time_vect, lon1 = (-180, 180), lat1 = (-90, 90)):
    lat_min, lat_max = lat1
    lon_min, lon_max = lon1
    #load netcdf from url as netCDF4 dataset
    ncfile = nc.Dataset(url+'#mode=bytes')
    #load netCDF4 dataset as Xarray dataset
    nc_ds = xr.open_dataset(xr.backends.NetCDF4DataStore(ncfile))
    # extract from Jan 1st to Jan 5th included
    extract_ds=nc_ds.sel(time=datetime(time_vect[0], time_vect[1], time_vect[2], time_vect[3]))
    extract_ds = extract_ds.sel(latitude=slice(lat_min, lat_max), longitude=slice(lon_min, lon_max))
    extract_ds=extract_ds[['hs']]
    nc_ds.close()
    return extract_ds.hs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog= 'ReadHSP2L',
                                     description = 'Read WW3 hs.nc and p2l.nc files for python, given the year, month, day, hour')
    parser.add_argument('-y', '--year', help='year', type=int, default=2021)
    parser.add_argument('-m', '--month', help='month', type=int, default=7)
    parser.add_argument('-d', '--day', help='day', type=int, default=4)
    parser.add_argument('-hr', '--hour', help='hour multiple of 3', type=int, default=12)

    args = parser.parse_args()
    
    file_path_p2l = "/Volumes/LaCie/data/2021/WW3-GLOB-30M_202107_p2l.nc"
    file_path_hs = "/Volumes/LaCie/data/2021/WW3-GLOB-30M_202107_hs.nc"
    year = int(args.year)
    month = int(args.month)
    day = int(args.day)
    hour = int(args.hour)
    time_vect = [year, month, day, hour]
    lon_min, lon_max = -180, 180
    lat_min, lat_max = -90, 90
    
    ## read p2l.nc
    ### xarray

    
    url = 'https://data-ww3.ifremer.fr/PROJECT/SISMO/CCI_ERA5_CERSATSSMI_MERCAREA_T702_NOREF_05/2010/FIELD_NC/CCI_WW3-GLOB-30M_201001_p2l.nc'
    
    lat, lon, freq, p2l, unit1 = read_p2l_from_url(url, [2010, 1, 1, 6])
    
    plt.pcolormesh(lon, lat, p2l.sel(frequency=0.1, method='nearest'))
    plt.show()