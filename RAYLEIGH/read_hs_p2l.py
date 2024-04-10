#!/usr/bin/env python3

from datetime import date
from datetime import datetime
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
import xarray as xr
import cartopy.crs as ccrs

def read_WWNC(file_path, time_vect, lon1, lat1):
    """Read netcdf _hs.nc file and return a matrix with dimension lon x lat
    of significant height of wind and swell waves in meters
    example:
    file_path = '../data/ftp.ifremer.fr/ifremer/ww3/HINDCAST/SISMO/GLOBAL05_2006_REF102040/data/WW3-GLOB-30M_200609_hs.nc'
    year = 2006
    month = 9
    day = 4
    hour = 12
    lon1 = []
    lat1 = []
    time_vect = [year, month, day, hour]
    A = read_WWNC(file_path, time_vect, [], [])
    print(A) """

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
    scale = f.variables['hs'].scale_factor
    hs = f.variables['hs'][KK:KK+nk][JJ:JJ+nj][II:II+ni]
    #hs = hs * scale values are the same as in matlab without scaling Oo
    f.close()
    hs = (np.squeeze(hs.filled(fill_value=np.nan))).T  # replace masked values in data by NaNs
    #print(hs.shape)
    return hs

def read_WWNCf(file_path, time_vect, lon1, lat1):
    """Read netcdf _p2l.nc file and return latitude, longitude, frequenc, p2l data which is the base 10 logarithm
 of power specral density of equivalent surface pressure and the units of p2l
    example:
    file_path = '../data/ftp.ifremer.fr/ifremer/ww3/HINDCAST/SISMO/GLOBAL05_2006_REF102040/WW3-GLOB-30M_200609_p2l.nc'
    year = 2006
    month = 9
    day = 4
    hour = 12
    lon1 = []
    lat1 = []
    time_vect = [year, month, day, hour]
    (lat, lon, freq, p2l, unit1) = read_WWNCf(file_path, time_vect, [], [])"""

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
        #print(time_vect)
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
    """Read netcdf significant wave height (_hs.nc) file and return latitude, longitude, frequenc, hs data
    units are m
    the output is an xarray of shape (lat, lon)
    
    example:
    file_path = '../data/ftp.ifremer.fr/ifremer/ww3/HINDCAST/SISMO/GLOBAL05_2006_REF102040/WW3-GLOB-30M_200609_hs.nc'
    year = 2006
    month = 9
    day = 4
    hour = 12
    lon1 = [-180, 180]
    lat1 = [-90, 90]
    time_vect = [year, month, day, hour]
    hs = read_hs(file_path, time_vect, lon1, lat1)
    
    # plot
    fig = plt.figure(figsize=(9,6))
    ax = plt.axes(projection=ccrs.Robinson())
    ax.coastlines()
    ax.gridlines()
    hs.plot(ax=ax, transform=ccrs.PlateCarree(), cbar_kwargs={'shrink': 0.4})
    plt.show()"""
    
    # open file
    try:
        ds = xr.open_dataset(file_path)
    except:
        print('Error opening file, please check that the file exists\n')
        print(file_path)

    # datetime
    year = time_vect[0]
    month = time_vect[1]
    day = time_vect[2]
    hour = time_vect[3]
    timestep = datetime(year, month, day, hour)
    
    # spatial extent
    lon = ds.longitude[np.logical_and(ds.longitude >= lon1[0], ds.longitude <= lon1[1])]
    lat = ds.latitude[np.logical_and(ds.latitude >= lat1[0], ds.latitude <= lat1[1])]
    
    # extract data
    hs = ds.hs.sel(time= timestep, longitude = lon, latitude= lat, method = 'nearest')
    return hs
    
def read_p2l(file_path, time_vect, lon1 = [-180, 180], lat1 = [-90, 90]):
    """Read netcdf _p2l.nc file and return latitude, longitude, frequenc, p2l data which is the base 10 logarithm
 of power specral density of equivalent surface pressure and the units of p2l
    example:
    file_path = '../data/ftp.ifremer.fr/ifremer/ww3/HINDCAST/SISMO/GLOBAL05_2006_REF102040/WW3-GLOB-30M_200609_p2l.nc'
    year = 2006
    month = 9
    day = 4
    hour = 12
    lon1 = []
    lat1 = []
    time_vect = [year, month, day, hour]
    (lat, lon, freq, p2l, unit1) = read_WWNCf(file_path, time_vect, lon1, lat1)"""
    
    # open file
    try:
        ds = xr.open_dataset(file_path)
    except:
        print('Error opening file, please check that the file exists\n')
        print(file_path)
    ds.assign_coords(f=("f", ds.f.data))
    ds = ds.rename({'f':'frequency'})
    # datetime
    year = time_vect[0]
    month = time_vect[1]
    day = time_vect[2]
    hour = time_vect[3]
    timestep = datetime(year, month, day, hour)
    # spatial extent
    lon = ds.longitude[np.logical_and(ds.longitude >= lon1[0], ds.longitude <= lon1[1])]
    lat = ds.latitude[np.logical_and(ds.latitude >= lat1[0], ds.latitude <= lat1[1])]
    # frequency range
    freq = ds.frequency
    # extract data
    p2l = ds.p2l.sel(time= timestep, frequency = freq, longitude = lon, latitude= lat, method = 'nearest')
    # units
    unit1 = ds.p2l.units
    return lat, lon, freq, p2l, unit1
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog= 'ReadHSP2L',
                                     description = 'Read WW3 hs.nc and p2l.nc files for python, given the year, month, day, hour')
    parser.add_argument('-y', '--year', help='year', type=int, default=2021)
    parser.add_argument('-m', '--month', help='month', type=int, default=7)
    parser.add_argument('-d', '--day', help='day', type=int, default=4)
    parser.add_argument('-hr', '--hour', help='hour multiple of 3', type=int, default=12)

    args = parser.parse_args()
    
    file_path_p2l = "/Users/tomasetl/Documents/code/ocean_source/data/2021/WW3-GLOB-30M_202107_p2l.nc"
    file_path_hs = "/Users/tomasetl/Documents/code/ocean_source/data/2021/WW3-GLOB-30M_202107_hs.nc"
    year = int(args.year)
    month = int(args.month)
    day = int(args.day)
    hour = int(args.hour)
    time_vect = [year, month, day, hour]
    lon_min, lon_max = -8, 40
    lat_min, lat_max = 30, 48
    
    ## read p2l.nc
    ### xarray
    lat , lon, freq, p2l, unit1 = read_p2l(file_path_p2l, time_vect, [lon_min, lon_max], [lat_min, lat_max])
    print(p2l.shape)
    p2l_f = p2l.sel(frequency= 0.3, method = 'nearest')
    fig = plt.figure(figsize=(9,6))
    ax = plt.axes(projection=ccrs.Robinson())
    ax.coastlines()
    ax.gridlines()
    p2l_f.plot(ax=ax, transform=ccrs.PlateCarree(), cbar_kwargs={'shrink': 0.4})
    plt.show()
    print(unit1)
    
#     # ### not xarray
#     # (lat, lon, freq, p2l, unit1, scale) = read_WWNCf(file_path_p2l, time_vect, [lon_min, lon_max], [lat_min, lat_max])
    
#     # print(p2l.shape)
#     # print(unit1, scale)
#     ## read hs.nc
#     ### xarray
#     hs_xr = read_hs(file_path_hs, time_vect, [lon_min, lon_max], [lat_min, lat_max])
#     fig = plt.figure(figsize=(9,6))
#     ax = plt.axes(projection=ccrs.Robinson())
#     ax.coastlines()
#     ax.gridlines()
#     hs_xr.plot(ax=ax, transform=ccrs.PlateCarree(), cbar_kwargs={'shrink': 0.4})
#     plt.show()
    
#     ### not xarray
#     A = read_WWNC(file_path_hs, time_vect, [], [])
#     print(A.shape)
#     ind_lat = np.squeeze(np.argwhere((lat > lat_min)&(lat < lat_max))) 
#     ind_lon = np.squeeze(np.argwhere((lon > lon_min)&(lon < lon_max)))        

#     hs = A[ind_lon, :]
#     hs = hs [:, ind_lat]
#     hs = np.flipud(hs)
    
#     fig, ax = plt.subplots(figsize=(15, 5))
#     im = ax.imshow(np.flipud(hs_xr), extent=(lon_min, lon_max, lat_min, lat_max), rasterized=True)
#     fig.colorbar(im, ax=ax, label='Significant wave height (m)', shrink=1)
#     ax.set_title('%d-%02d-%02dT%02d:00:00.000000'%(year, month, day, hour))
#     plt.tight_layout()
#     plt.show()
#     #plt.savefig('hs_%d-%02d-%02dT%02d.png'%(year, month, day, hour), transparent=True)
