#!/usr/bin/env python3

# Preamble
#__author__ = "Lisa Tomasetto"
#__copyright__ = "Copyright 2024, UGA"
#__credits__ = ["Lisa Tomasetto"]
#__version__ = "0.1"
#__maintainer__ = "Lisa Tomasetto"
#__email__ = "lisa.tomasetto@univ-grenoble-alpes.fr"


"""This set of functions aims at modeling the ambient noise source in the secondary microseismic range for body waves.
Using Lucia Gualtieri's site effect computation and WAVEWATCHIII model.

It contains six functions:

- `download_ww3_local(YEAR, MONTH, ftp_path_to_files, ww3_local_path, prefix)`: download the WW3 files for the given month.

- `open_bathy(file_bathy, refined_bathymetry, extent)`: open the bathymetry file. Either default WW3 grid, ETOPOv2 or a custom grid.

- `subfctn_liquid_solid(p, mi, mt)`: compute the reflection and transmission coefficients between two media.

- `bathy(z, f, p, m)`: compute the amplification coefficient for P and S waves.

- `ampli(dpt1, f, rp, layers, theta)`: compute the amplification coefficient for P and S waves integrated over a range of takeoff angles.

- `loop_ww3_sources(paths, dpt1, zlon, zlat, wave_type, date_vec, extent, parameters, c_file, prefix, **kwargs)`: compute the proxy for the source force amplitude.
"""
##################################################################################

##Libraries
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
import cartopy
import os.path

from netCDF4 import Dataset, date2num
from math import radians, log
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from datetime import datetime
from calendar import monthrange
from numpy.lib.scimath import sqrt as csqrt

from wmsan.read_hs_p2l import read_p2l

plt.style.use("ggplot")
SMALL_SIZE = 18
MEDIUM_SIZE = 20
BIGGER_SIZE = 22

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
##################################################################################
############# DOWNLOAD WW3 FILES #################################################
##################################################################################

def download_ww3_local(YEAR, MONTH, ftp_path_to_files="ftp://ftp.ifremer.fr/ifremer/dataref/ww3/GLOBMULTI_ERA5_GLOBCUR_01/GLOB-30M/2020/FIELD_NC/", ww3_local_path= '../../data/ww3/', prefix = "WW3-GLOB-30M"):
    """Download WW3 files for a given year and month from the specified FTP path to a local directory.
    
    Args:
        YEAR (int): The year for which the files should be downloaded.
        MONTH (list): The month or months for which the files should be downloaded.
        ftp_path_to_files (str, optional): The FTP path to the WW3 files.
        ww3_local_path (str, optional): The local directory where the files should be saved.
        prefix (str, optional): The prefix for the WW3 files.
    
    """
    workdir = os.getcwd()
    # create directory if it does not exist
    try: 
        os.mkdir(ww3_local_path) 
    except OSError as error: 
        print(error)
    os.chdir(ww3_local_path)
    
    # download files
    if len(MONTH) == 0:
        print("-----------------------------------------------------------------\n")
        print("Downloading WW3 files for year %d\n"%(YEAR))
        print("-----------------------------------------------------------------\n")
        MONTH = np.arange(1, 13)
    for m in MONTH:
        print("Downloading can take some time...\n")
        file_p2l = ftp_path_to_files + "%s_%d%02d_p2l.nc"%(prefix, YEAR, m) # p2l file
        check_file_p2l = ww3_local_path + "%s_%d%02d_p2l.nc"%(prefix, YEAR, m) # p2l file
        
        if os.path.exists(check_file_p2l):
            print("-----------------------------------------------------------------\n")
            print(check_file_p2l + " already downloaded\n")
            print("-----------------------------------------------------------------\n")
        else:
            os.system("wget -nv -c %s"%(file_p2l))
            print("-----------------------------------------------------------------\n")
            print(file_p2l + " downloaded\n")
            print("-----------------------------------------------------------------\n")

    print( "WW3 files downloaded in %s"%(ww3_local_path))
    os.chdir(workdir)
    print("current directory : ", os.getcwd())

##################################################################################
################ OPEN BATHYMETRY FILE ############################################
##################################################################################

def open_bathy(file_bathy = '../../data/WW3-GLOB-30M_202002_p2l.nc', refined_bathymetry=False, extent=[-180, 180, -90, 90]):
    """Open bathymetry file and optionally refine bathymetry using ETOPOv2 dataset. 

    Args:
        file_bathy (str): Path to the bathymetry file.
        refined_bathymetry (bool, optional): Whether to use the refined ETOPOv2 dataset. Defaults to False.
        extent (list, optional): The geographical extent of the bathymetry data in the format [lon_min, lon_max, lat_min, lat_max].

    Returns:
        dpt1_mask (xarray.DataArray): Masked bathymetry data.
        zlon (xarray.DataArray): Longitude coordinates.
        zlat (xarray.DataArray): Latitude coordinates.
    """
    [lon_min, lon_max, lat_min, lat_max] = extent
    if np.abs(lat_min) > 90 or np.abs(lat_max) > 90:
        print("Latitude not correct, absolute value > 90")
        return
    
    if refined_bathymetry or file_bathy == '../../data/ETOPO_2022_v1_60s_N90W180_bed.nc':
        try:
            ds = xr.open_mfdataset(file_bathy, combine='by_coords')
            ds  = ds.rename({'lon':'longitude', 'lat': 'latitude'})
            if extent[0] > extent[1]:
                ## work on the pacific ocean
                ds = ds.assign_coords(longitude=((360 + (ds.longitude % 360)) % 360))
                ds = ds.roll(longitude=int(len(ds['longitude']) / 2),roll_coords=True)
                lon_min = ((360 + (lon_min % 360)) % 360)
                lon_max = ((360 + (lon_max % 360)) % 360)
            z = ds['z']
            z *= -1 # ETOPOv2 to Depth
            z = z.where(z>0, other=np.nan)
            dpt1 = z.sel(latitude = slice(lat_min, lat_max), longitude = slice(lon_min, lon_max))         
        except:
            try:
                # load refined bathymetry ETOPOv2
                file_bathy = '../../data/ETOPO_2022_v1_60s_N90W180_bed.nc'
                ds = xr.open_mfdataset(file_bathy, combine='by_coords')
                ds  = ds.rename({'lon':'longitude', 'lat': 'latitude'})
                if extent[0] > extent[1]:
                    ## work on the pacific ocean
                    ds = ds.assign_coords(longitude=((360 + (ds.longitude % 360)) % 360))
                    ds = ds.roll(longitude=int(len(ds['longitude']) / 2),roll_coords=True)
                    lon_min = ((360 + (lon_min % 360)) % 360)
                    lon_max = ((360 + (lon_max % 360)) % 360)
                z = ds['z']
                z *= -1 # ETOPOv2 to Depth
                z = z.where(z>0, other=np.nan)
                dpt1 = z.sel(latitude = slice(lat_min, lat_max), longitude = slice(lon_min, lon_max))    
            except:
                print("Refined bathymetry ETOPOv2 not found. \nYou can download it from:\n https://www.ngdc.noaa.gov/thredds/catalog/global/ETOPO2022/60s/60s_bed_elev_netcdf/catalog.html?dataset=globalDatasetScan/ETOPO2022/60s/60s_bed_elev_netcdf/ETOPO_2022_v1_60s_N90W180_bed.nc\nSave in ../data/")
                return None, None, None
    else:
        ds = xr.open_mfdataset(file_bathy, combine='by_coords')
        if lon_min > lon_max:
            ## work on the pacific ocean
            ds = ds.assign_coords(longitude=((360 + (ds.longitude % 360)) % 360))
            ds = ds.roll(longitude=int(len(ds['longitude']) / 2),roll_coords=True)
            lon_min = ((360 + (lon_min % 360)) % 360)
            lon_max = ((360 + (lon_max % 360)) % 360)
        dpt1 = ds['dpt'].squeeze(dim = 'time', drop=True)
        dpt1 = dpt1.sel(latitude = slice(lat_min, lat_max), longitude = slice(lon_min, lon_max))
    ## Mask nan values    
    dpt1_mask = dpt1.where(np.isfinite(dpt1))
    zlon = dpt1_mask.longitude
    zlat = dpt1_mask.latitude
    return dpt1_mask, zlon, zlat
##################################################################################
############################ SITE EFFECT #########################################
##################################################################################
def subfcn_liquid_solid(p, mi, mt):
    """Calculate the reflection and transmission coefficients for P and S waves between two media.
    Author: LI Lei, ll.ynyf@gmail.com modified by Pierre Boue 23/11/2020
    
    Args:
        p (float or array-like): The wave slowness.
        mi (tuple): The properties of the incident medium, consisting of the vertical P-wave velocity (vp1) and density (rho1).
        mt (tuple): The properties of the transmitted medium, consisting of the vertical P-wave velocity (vp2), shear wave velocity (vs2), and density (rho2).

    Returns:
        Rpp (float or array-like): The reflection coefficient for P waves.
        Tpp (float or array-like): The transmission coefficient for P waves.
        Tps (float or array-like): The transmission coefficient for S waves.
    """
    p2 = np.array(p)**2
    # vp, density of incident medium, vp, vs, density of transmitted medium
    (vp1, rho1, vp2, vs2, rho2) = (mi[0], mi[-1], mt[0], mt[1], mt[2])
    q1 = csqrt(1/vp1**2 - p2)  # vertical slowness
    q2p = csqrt(1/vp2**2 - p2)
    q2s = csqrt(1/vs2**2 - p2)

    a = rho2*q1*((1-2*p2*vs2**2)**2 + 4*vs2**4*p2*q2p*q2s)
    b = rho1*q2p
    D = a + b
    Rpp = (a-b)/D
    Tpp = 2*rho1*q1*(1-2*(p*vs2)**2)*np.reciprocal(D)
    Tps = 4*rho1*q1*q2p*p*(vs2**2)*np.reciprocal(D)
    return Rpp, Tpp, Tps

def bathy(z, f, p=[], m= [1500, 1000, 55400, 3200, 2500]):
    """Bathymetry secondary microseismic excitation coefficients (for P, S amplitude).
    Based on LI Lei, ll.ynyf@gmail.com modified by Pierre Boue 23/11/2020.
    
    Examples:
        >>> z = np.linspace(0, 25000, 5001)
        >>> f = 1/8
        >>> cP, cS = bathy(z, f, [], [])
        >>> x = f*z/1500
        >>> plt.figure()
        >>> plt.plot(x, abs(cP[0]), x, abs(cS[0]))
        >>> plt.show()
        
    Args:
        z (np.ndarray): thickness of water layer in meters, if z in m|km, v should be in m/s|km/s. All rhos must keep the same units.
        f (np.ndarray): seismic frequency in Hz
        p (np.ndarray, optional): slowness, if not specified return p values integral to 1/vp_crust
        m (list, optional): [vp_water, rho_water, vp_crust, vs_crust, rho_crust]

    Returns:
        cP (np.ndarray): excitation coefficients of P waves in shape of [[p,]f, z]
        cS (np.ndarray): excitation coefficients of S waves in shape of [[p,]f, z]
    
    """
    ##
    (rhow, rhoc) = (m[1], m[4])  # in kg/m³
    (vpw, vpc, vsc) = (m[0], m[2], m[3])  # in m/s
    if len(m) > 0:
        if (len(m) >= 1 and m[0] > 0):
            vpw = m[0]
        if (len(m) >= 2 and m[1] > 0):
            rhow = m[1]
        if (len(m) >= 3 and m[2] > 0):
            vpc = m[2]
        if (len(m) >= 4 and m[3] > 0):
            vsc = m[3]
        if (len(m) >= 5 and m[4] > 0):
            rhoc = m[4]

    elif np.max(z) < 50:
        (vpw, rhow, vpc, vsc, rhoc) = (1.5, 1, 5.54, 3.2, 2.5)
        print('[rho_w, rho_c] = [1, 2.5] g/cm³')
        print('[vpw, vpc, vsc] = [1.5, 5.54, 3.2] km/s')
    if len(p)==0:
        p = np.linspace(0, 0.995, 200)/vpc
        m = np.array([vpw, rhow, vpc, vsc, rhoc])
        a = np.arcsin(vpw*p)
        cP = np.zeros((np.size(f), np.size(z)), dtype=np.complex128)
        cS = np.zeros((np.size(f), np.size(z)), dtype=np.complex128)
        for i in range(np.size(f)):
            if np.size(f) == 1:
                c_P, c_S = bathy(z, f, p, m)
            else:
                c_P, c_S = bathy(z, f[i], p, m)
            cP[i, :] = cP[i, :] + np.trapz(abs(c_P)**2, x=a, axis=0)
            cS[i, :] = cS[i, :] + np.trapz(abs(c_S)**2, x=a, axis=0)
        (cP, cS) = (csqrt(cP), csqrt(cS))
        return cP, cS
    
    elif (p != np.linspace(0, 0.995, 200)/vpc).all():
        c_P = np.empty((np.size(p), np.size(f), np.size(z)), dtype=np.complex128)
        c_P[:] = np.nan
        c_S = np.empty((np.size(p), np.size(f), np.size(z)), dtype=np.complex128)
        c_S[:] = np.nan
        cP = np.zeros((np.size(f), np.size(z)), dtype=np.complex128)
        cS = np.zeros((np.size(f), np.size(z)), dtype=np.complex128)
        Rpp, Tpp, Tps = subfcn_liquid_solid(p, [vpw, rhow], [vpc, vsc, rhoc])
        qw = csqrt(1/vpw**2 - p**2)
        if np.size(f) == 1:
            phi = 4*np.pi*f*(z.flatten('F')).T
            for i in range(np.size(p)):
                C = 1*np.reciprocal(1 + Rpp[i]*np.exp(1j*phi*qw[i]))
                c_P[i, :, :] = np.dot(Tpp[i], C) 
                c_S[i, :, :] = np.dot(Tps[i], C)  
        else:
            for i in range(np.size(f)):
                c_Pi, c_Si = bathy(z, f[i], p, m)
                c_P[:, i, :] = c_Pi
                c_S[:, i, :] = c_Si

        for i in range((np.size(f))):
                a = np.arcsin(vpw*p)
                cP[i, :] = cP[i, :] + np.trapz(abs(c_P[:,i,:])**2, x=a, axis=0)
                cS[i, :] = cS[i, :] + np.trapz(abs(c_S[:, i, :])**2, x=a, axis=0)
        (cP, cS) = (np.squeeze(cP), np.squeeze(cS))
    Rpp, Tpp, Tps = subfcn_liquid_solid(p, [vpw, rhow], [vpc, vsc, rhoc])
    qw = csqrt(1/vpw**2 - p**2)
    if np.size(f) == 1:
        phi = 4*np.pi*f*(z.flatten('F')).T
    else:
        phi = 4*np.pi*np.dot(f.flatten('F'), (z.flatten('F')).T)
    cP = np.empty((np.size(p), np.size(f), np.size(z)), dtype=np.complex128)
    cP[:] = np.nan
    cS = cP.copy()
    for i in range(np.size(p)):
        C = 1*np.reciprocal(1 + Rpp[i]*np.exp(1j*phi*qw[i]))
        cP[i, :, :] = np.dot(Tpp[i], C)  
        cS[i, :, :] = np.dot(Tps[i], C)
    (cP, cS) = (np.squeeze(cP), np.squeeze(cS))
    return cP, cS

def ampli(dpt1, f, rp=[], layers=[1500, 1000, 5540, 3200, 2500], theta = radians(15.71)):
    """
    Compute amplification coefficient for P and S waves. 
    
    Args:
        dpt1 (ndarray): bathymetry grid in meters
        f (ndarray): frequency vector in Hz
        rp (ndarray, optional): ray parameter matrix to integrate over. 
        layers (list, optional): layers properties [Vp_w, rho_w, Vp_c, Vs_c, rho].
        theta (float, optional): limit angle of. Defaults to radians(15.71).
    
    Returns:
        cP (ndarray): Amplification coefficient for P waves.
        cS (ndarray): Amplification coefficient for S waves.
        bathy_ampli_P (ndarray): Amplification coefficient for P waves divided by theta.
        bathy_ampli_S (ndarray): Amplification coefficient for S waves divided by theta.
        
    """
    # compute amplification for Tp and rp following G14
    (cP, cS) = bathy(dpt1, f, rp, layers)
    if np.size(f)==1:
        cP = abs(cP.reshape(dpt1.shape, order='F').copy())
        cS = abs(cS.reshape(dpt1.shape, order='F').copy())
        bathy_ampli_P = cP/theta
        bathy_ampli_S = cS/theta
        return cP, cS, bathy_ampli_P, bathy_ampli_S
    else:
        new_cP = np.empty((np.size(f), dpt1.shape[0], dpt1.shape[1]))
        new_cS = np.empty((np.size(f), dpt1.shape[0], dpt1.shape[1]))
        for i in range(np.size(f)):
            new_cP[i, :] = abs(cP[i,:].reshape(dpt1.shape, order='F').copy())
            new_cS[i, :] = abs(cS[i,:].reshape(dpt1.shape, order='F').copy())
        bathy_ampli_P = new_cP/theta
        bathy_ampli_S = new_cS/theta
        return new_cP, new_cS, bathy_ampli_P, bathy_ampli_S

##################################################################################
################# LOOP WW3 SOURCES ###############################################
##################################################################################

def loop_ww3_sources(paths, dpt1, zlon, zlat, wave_type='P', date_vec=[2020, [], [], []], extent=[-180, 180, -90, 90],parameters= [1/12, 1/2], c_file = "../../data/cP.nc", prefix = "WW3-GLOB-30M", **kwargs):
    """Compute the Proxy for the Source Force on the seafloor for a given wave type (P or S), given a path to the ww3 p2l file, the bathymetry, the wave type, the date vector and the spatial extent.
    Saves in netcdf format the Proxy for the Source Force for each frequency if save argument True.
    Plots in PNG source maps of P/S waves at given intervals depending on plot variables.

    Args:
        paths (list): [file_bathy, ww3_local_path]: paths of additional files bathymetry, ww3 p2l file
        dpt1 (xarray.ndarray): bathymetry grid in m (depth) with dimensions lon x lat
        zlon (xarray.ndarray): longitude of bathymetry file (°)
        zlat (xarray.ndarray): latitude of bathymetry file (°)
        wave_type (str, optional): P or S waves.
        date_vec (list, optional): date vector [year, month, day, hour], with hour in [0, 3, 6, 9, 12, 15, 18, 21].
        extent (list, optional): spatial extent format [lon_min, lon_max, lat_min, lat_max].
        parameters (list, optional): parameters minimum frequency, maximum frequency.
        c_file (str, optional): path to amplification coefficient file.
        prefix (str, optional) : prefix of ww3 p2l file.
        plot (bool, optional): default: True
        plot_hourly (bool, optional): plot maps every 3-hours default: False
        plot_daily (bool, optional): plot maps every day, default: False
        plot_monthly (bool, optional): plot maps every month, default : True
        plot_yearly (bool, optional): plot map for the year average, default: False
        save (bool, optional): save 3-hourly matrix, default: False

    """
    
    ww3_local_path = paths[1]
    
    # Constants
    radius = 6.371*1e6 # radius of the earth in meters
    lg10 = log(10) # log of 10
    #
    f1 = parameters[0]
    f2 = parameters[1]
    ## Initialize variables
    if 'plot_type' in kwargs:
        plot_type = kwargs['plot_type']
        if plot_type == 'hourly':
            plot_hourly = True
        else:
            plot_hourly = False
        if plot_type == 'daily':
            plot_daily = True
            F_daily =  np.zeros(dpt1.shape)
        else:
            plot_daily = False
        if plot_type == 'monthly':
            plot_monthly = True
            F_monthly =  np.zeros(dpt1.shape)
        else:
            plot_monthly = False
        if plot_type == 'yearly':
            plot_yearly = True
            F_yearly =  np.zeros(dpt1.shape)
        else:
            plot_yearly = False
    
    if 'save' in kwargs:
        save = kwargs['save']
    else:
        save = False
        
    if 'vmin' in kwargs:
        vmin = kwargs['vmin']
    else:
        vmin = 0
        
    if 'vmax' in kwargs:
        vmax = kwargs['vmax']
    else:
        vmax = 1e10

    ## Adapt latitude and longitude to values in parameters file
    lon_min = extent[0]
    lon_max = extent[1]
    lat_min = extent[2]
    lat_max = extent[3]
    
    central_longitude = 0
    if np.abs(lat_min) > 90 or np.abs(lat_max) > 90:
        print("Latitude not correct, absolute value > 90")
        return
    
    if lon_min > lon_max:
        ## work on the pacific ocean
        lon_min = ((360 + (lon_min % 360)) % 360)
        lon_max = ((360 + (lon_max % 360)) % 360)
        central_longitude = 180
        
    dpt1 = dpt1.sel(latitude = slice(lat_min, lat_max), longitude = slice(lon_min, lon_max))
    zlat = zlat.sel(latitude = slice(lat_min, lat_max))
    zlon = zlon.sel(longitude = slice(lon_min, lon_max))
    
    ## Open Amplification Coefficient
    # check for refined bathymetry
    res_bathy = abs(zlon[1] - zlon[0])
    if res_bathy == 0.5:
        ds_ampli = xr.open_dataset('../../data/c%s.nc'%(wave_type)).astype('float64')

    else:
        print("Refined bathymetry grid \n PLEASE RUN amplification_coefficients.ipynb before running this script")
        ds_ampli = xr.open_dataset(c_file).astype('float64')

        
    if extent[0] > extent[1]:
        ds_ampli = ds_ampli.assign_coords(longitude=((360 + (ds_ampli.longitude % 360)) % 360))
        ds_ampli = ds_ampli.roll(longitude=int(len(ds_ampli['longitude']) / 2),roll_coords=True)
    amplification_coeff = ds_ampli['c%s'%wave_type]
    
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
    
    for iyear in np.array([YEAR]):
        if isinstance(MONTH, int):
            MONTH = np.array([MONTH])
        elif not len(MONTH):
            MONTH = np.arange(1, 13)
        else:
            MONTH = np.array(MONTH)
        for imonth in MONTH:
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
                        df = df[index_freq]
                        freq_seismic = freq_seismic[index_freq]
                        n_freq = len(df)
                        Fp = p2l.sel(frequency = freq_ocean[index_freq], latitude = slice(lat_min, lat_max), longitude = slice(lon_min, lon_max))
                        Fp = Fp.where(np.isfinite(Fp))
                        Fp = Fp.drop('frequency')
                        Fp.coords['frequency'] = freq_seismic
                        res_p2l = abs(Fp.longitude[1]-Fp.longitude[0])
                        if res_bathy != res_p2l:
                            # interpolate Fp
                            Fp = Fp.interp(longitude=zlon, latitude=zlat, method='linear')
                        amplification_coeff = amplification_coeff.sel(frequency = freq_seismic, method='nearest', tolerance=0.01)
                        amplification_coeff = amplification_coeff.sel(latitude = slice(lat_min, lat_max), longitude = slice(lon_min, lon_max))
                        amplification_coeff = amplification_coeff.reindex_like(Fp, method='nearest', tolerance=0.01)
                        F_f = Fp*amplification_coeff**2
                        F = F_f.copy()
                        for ifq, fq in enumerate(freq_seismic):
                            F_f[ifq, :, :] *= dA
                            F[ifq, :, :] = F_f[ifq, :, :]*df[ifq]
                        F = 2*np.pi*np.sqrt(F.sum(dim = 'frequency'))
                    ## Single frequency
                    elif f1 == f2:
                        index_freq = np.squeeze(np.argmin(abs(freq_seismic-f1)))
                        print('unique frequency ', f1)
                        return
                        
                    ## Exception in parametrization of frequencies
                    else:
                        print('two frequencies with f2 < f1 were given')
                        return

                    ## Save F to file
                    if save:
                        path_out = './F/'
                        if not os.path.exists(path_out):
                            print("make directory %s"%path_out)
                            os.makedirs(path_out)

                        # Create netCDF 3-hourly file
                        ncfile = Dataset(path_out+"F_%d%02d%02d%02d.nc"%(iyear, imonth, iday, ih), mode='w',format='NETCDF4_CLASSIC')
                        lat_dim = ncfile.createDimension('latitude', len(zlat))  # latitude axis
                        lon_dim = ncfile.createDimension('longitude', len(zlon))  # longitude axis
                        time_dim = ncfile.createDimension('time', daymax*8)  # unlimited axis (can be appended to).
                        freq_dim = ncfile.createDimension('frequency', n_freq)
                        ncfile.title='Proxy for the Source Force on %d-%02d-%02d-%02d'%(iyear, imonth, iday, ih)
                        ncfile.subtitle='Equivalent Force every 3 hours for the secondary microseismic peak'
                        lat = ncfile.createVariable('latitude', np.float64, ('latitude',))
                        lat.units = 'degrees_north'
                        lat.long_name = 'latitude'
                        lon = ncfile.createVariable('longitude', np.float64, ('longitude',))
                        lon.units = 'degrees_east'
                        lon.long_name = 'longitude'
                        time = ncfile.createVariable('time', np.float64, ('time',))
                        time.units = 'hours since 1990-01-01'
                        time.long_name = 'time'
                        freq_nc = ncfile.createVariable('frequency', np.float64, ('frequency',))
                        freq_nc.units = 'Hz'
                        freq_nc.long_name = 'frequency'
                        F_freq = ncfile.createVariable('F_f', np.float64, ('frequency','latitude', 'longitude'))
                        F_freq.units = 'N.s^{1/2}'
                        F_freq.long_name = 'Proxy for the Source Force spectrum'
                        lat[:] = zlat
                        lon[:] = zlon
                        freq_nc[:] = freq_seismic
                        F_freq[:] = F_f
                        time[:] = date2num(datetime(iyear, imonth, iday, ih), units='hours since 1990-01-01', calendar='standard')
                        ncfile.close()
                        
                    ## Plot F
                    if plot_hourly:
                        plt.close('all') 
                        F_plot = xr.DataArray(F, 
                                                coords={'latitude': zlat,'longitude': zlon}, 
                                                dims=["latitude", "longitude"],
                                                name = 'Frequency %.3f-%.3f Hz.%d-%02d-%02dT%02d\n %s waves.\n'%(f1, f2, iyear, imonth, iday, ih, wave_type))
                        fig = plt.figure(figsize=(9,6))
                        fig.suptitle('Frequency %.3f-%.3f Hz.%d-%02d-%02dT%02d\n %s waves.\n'%(f1, f2, iyear, imonth, iday, ih, wave_type))
                        ax = plt.axes(projection=ccrs.Robinson(central_longitude=central_longitude))
                        ax.coastlines()
                        gl = ax.gridlines()
                        gl.xformatter = LONGITUDE_FORMATTER
                        gl.yformatter = LATITUDE_FORMATTER
                        ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k', facecolor='linen')
                        if wave_type == 'P':
                            F_plot.plot(ax=ax, transform=ccrs.PlateCarree(),  cbar_kwargs={'label':'F (N)', 'orientation': 'horizontal'}, vmin=vmin, vmax=vmax)
                        else:
                            F_plot.plot(ax=ax, transform=ccrs.PlateCarree(),  cbar_kwargs={'label':'F (N)', 'orientation': 'horizontal'}, vmin=vmin, vmax=vmax)
                        plt.savefig('F_%s_%d%02d%02dT%02d.png'%(wave_type, iyear, imonth, iday, ih), dpi = 300, bbox_inches='tight')

                    ## Sum F
                    if plot_daily :
                        F_daily += F
                    if plot_monthly :
                        F_monthly += F
                    if plot_yearly :
                        F_yearly += F
                   
                if plot_daily:
                    plt.close('all') 
                    F_plot = xr.DataArray(F_daily, 
                        coords={'latitude': zlat,'longitude': zlon}, 
                        dims=["latitude", "longitude"],
                        name = 'Frequency %.3f-%.3f Hz.%d-%02d-%02d\n %s waves.\n'%(f1, f2, iyear, imonth, iday, wave_type))
                    fig = plt.figure(figsize=(9,6))
                    fig.suptitle('Frequency %.3f-%.3f Hz.%d-%02d-%02d.\n %s waves.\n'%(f1, f2, iyear, imonth, iday, wave_type))
                    ax = plt.axes(projection=ccrs.Robinson(central_longitude=central_longitude))
                    ax.coastlines()
                    gl = ax.gridlines()
                    gl.xformatter = LONGITUDE_FORMATTER
                    gl.yformatter = LATITUDE_FORMATTER
                    ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k', facecolor='linen')
                    if wave_type == 'P':
                        F_plot.plot(ax=ax, transform=ccrs.PlateCarree(),  cbar_kwargs={'label':'F (N)', 'orientation': 'horizontal'}, vmin=vmin, vmax=vmax),
                    else:
                        F_plot.plot(ax=ax, transform=ccrs.PlateCarree(),  cbar_kwargs={'label':'F (N)', 'orientation': 'horizontal'}, vmin=vmin, vmax=vmax)
                    plt.savefig('F_%s_%d%02d%02d.png'%(wave_type, iyear, imonth, iday), dpi = 300, bbox_inches='tight')
                    F_daily = np.zeros((dpt1.shape))
                    
            if plot_monthly :
                plt.close('all')
                F_plot = xr.DataArray(F_monthly, 
                    coords={'latitude': zlat,'longitude': zlon}, 
                    dims=["latitude", "longitude"],
                    name = 'Frequency %.3f-%.3f Hz.%d-%02d\n %s waves.\n'%(f1, f2, iyear, imonth, wave_type))
                fig = plt.figure(figsize=(9,6))
                fig.suptitle('Frequency %.3f-%.3f Hz.%d-%02d\n %s waves.\n'%(f1, f2, iyear, imonth, wave_type))
                ax = plt.axes(projection=ccrs.Robinson(central_longitude=central_longitude))
                ax.coastlines()
                gl = ax.gridlines()
                gl.xformatter = LONGITUDE_FORMATTER
                gl.yformatter = LATITUDE_FORMATTER
                ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k', facecolor='linen')
                if wave_type == 'P':
                    F_plot.plot(ax=ax, transform=ccrs.PlateCarree(),  cbar_kwargs={'label':'F (N)', 'orientation': 'horizontal'}, vmin=vmin, vmax=vmax)
                else:
                    F_plot.plot(ax=ax, transform=ccrs.PlateCarree(),  cbar_kwargs={'label':'F (N)', 'orientation': 'horizontal'}, vmin=vmin, vmax=vmax)
                plt.savefig('F_%s_%d%02d.png'%(wave_type, iyear, imonth), dpi = 300, bbox_inches='tight')
                F_monthly = np.zeros((dpt1.shape))
                    

        if plot_yearly:
            plt.close('all')
            F_plot = xr.DataArray(F_yearly,
                                coords={'latitude': zlat,'longitude': zlon}, 
                                dims=["latitude", "longitude"],
                                name = 'Frequency %.3f-%.3f Hz.%d\n %s waves.\n'%(f1, f2, iyear, wave_type))
            fig = plt.figure(figsize=(9,6))
            fig.suptitle('Frequency %.3f-%.3f Hz.%d\n %s waves.\n'%(f1, f2, iyear, wave_type))
            ax = plt.axes(projection=ccrs.Robinson(central_longitude=central_longitude))
            ax.coastlines()
            gl = ax.gridlines()
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k', facecolor='linen')
            if wave_type == 'P':
                F_plot.plot(ax=ax, transform=ccrs.PlateCarree(),  cbar_kwargs={'label':'F (N)', 'orientation': 'horizontal'}, vmin=vmin, vmax=vmax)
            else:
                F_plot.plot(ax=ax, transform=ccrs.PlateCarree(),  cbar_kwargs={'label':'F (N)', 'orientation': 'horizontal'}, vmin=vmin, vmax=vmax)
            plt.savefig('F_%s_%d.png'%(wave_type, iyear), dpi = 300, bbox_inches='tight')
            F_daily = np.zeros((dpt1.shape))
            F_yearly = np.zeros((dpt1.shape))
        plt.close('all')
    print('%s source maps done!'%wave_type)