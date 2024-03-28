#!/usr/bin/env python3

# Preamble
__author__ = "Lisa Tomasetto"
__copyright__ = "Copyright 2024, UGA"
__credits__ = ["Lisa Tomasetto"]
__version__ = "1.0"
__maintainer__ = "Lisa Tomasetto"
__email__ = "lisa.tomasetto@univ-grenoble-alpes.fr"
__status__ = "Production"

"""
This program aims at modeling the ambient noise source in the secondary microseismic range for Rayleigh waves.
Using Longuet-Higgins site effect and F.Ardhuin WW3 model.
"""
##################################################################################

##Libraries
import numpy as np
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import xarray as xr
import pandas as pd
import os.path

from netCDF4 import Dataset, date2num
from math import radians, log
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from obspy.geodetics.base import gps2dist_azimuth
from datetime import datetime
from scipy.interpolate import interp1d
from tqdm import tqdm
from calendar import monthrange
from pyproj import Geod

##Own librairies
from read_hs_p2l import read_hs, read_p2l
from readWW31 import read_dpt

## Set font size parameters to make readable figures
plt.style.use("ggplot")
SMALL_SIZE = 18
MEDIUM_SIZE = 22
BIGGER_SIZE = 24

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

plt.rcParams['xtick.direction'] = 'inout'
plt.rcParams['ytick.direction'] = 'inout'
plt.rcParams['font.family'] = "serif"

def site_effect(z, f, zlat, zlon, vs_crust=2.8, path='./longuet_higgins.txt'):
    """ Bathymetry secondary microseismic excitation coefficients (Rayleigh waves).
    Input:
    z : thickness of water layer
    f: seismic frequency in hertz
    vs_crust : shear waves velocity in the crust (sea bed)
      default value is 2.8 km/s
    path : the path to the Longuet Higgins file containing tabulated values of site effect coefficient for the 4th first modes.

    Output :
    Numerical value of the site effect in the shape (f.shape, z.shape)."""

    df = pd.read_csv('%s'%path, sep='\t', header =0, usecols=[0, 1, 2, 3, 4, 5, 6, 7], names = ['fh1', 'c1', 'fh2', 'c2', 'fh3', 'c3', 'fh4', 'c4'])
    fc1 = interp1d(df.fh1, df.c1, kind='nearest', bounds_error=False, fill_value=0)
    fc2 = interp1d(df.fh2, df.c2, kind='nearest', bounds_error=False, fill_value=0)
    fc3 = interp1d(df.fh3, df.c3, kind='nearest', bounds_error=False, fill_value=0)
    fc4 = interp1d(df.fh4, df.c4, kind='nearest', bounds_error=False, fill_value=0)
    try :
        n = len(f)
        x = z.shape[1]
        y = z.shape[0]
        #print(n, x, y)
        #print(z.shape)
        C = np.empty((n, y, x))
        for i, fq in enumerate(f):
            fh_v = 2*np.pi*fq*z/(vs_crust*1e3)
            C[i, :, :] = fc1(fh_v)**2 + fc2(fh_v)**2 + fc3(fh_v)**2 + fc4(fh_v)**2
    except:
        raise
    C = np.squeeze(C)
    ## C to xarray
    C = xr.DataArray(C, dims=('frequency','latitude', 'longitude'), coords={'frequency': f,'latitude': zlat, 'longitude': zlon})
    return C

def download_ww3_local(YEAR, MONTH, ftp_path_to_files="ftp://ftp.ifremer.fr/ifremer/dataref/ww3/GLOBMULTI_ERA5_GLOBCUR_01/GLOB-30M/2020/FIELD_NC/", ww3_local_path= './data/ww3/'):
    
    """
    Download WW3 files for a given year and month from the specified FTP path to a local directory.
    
    Input:
    - YEAR: the year for which the files should be downloaded
    - MONTH: the month or months for which the files should be downloaded
    - ftp_path_to_files: the FTP path to the WW3 files
    - ww3_local_path: the local directory where the files should be saved
    
    Output:
    None
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
        file_p2l = ftp_path_to_files + "LOPS_WW3-GLOB-30M_%d%02d_p2l.nc"%(YEAR, m) # p2l file
        check_file_p2l = ww3_local_path + "LOPS_WW3-GLOB-30M_%d%02d_p2l.nc"%(YEAR, m) # p2l file
        
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


def open_bathy(file_bathy, refined_bathymetry=False, extent=[-180, 180, -90, 90]):
    """
    Open bathymetry file and optionally refine bathymetry using ETOPOv2 dataset. 

    Input:
    file_bathy (str): Path to the bathymetry file.
    refined_bathymetry (bool, optional): Whether to use the refined ETOPOv2 dataset. Defaults to False.
    extent (list, optional): The geographical extent of the bathymetry data in the format [lon_min, lon_max, lat_min, lat_max]. Defaults to [-180, 180, -90, 90].

    Output:
    dpt1_mask (xarray.DataArray): Masked bathymetry data.
    zlon (xarray.DataArray): Longitude coordinates.
    zlat (xarray.DataArray): Latitude coordinates.
    """
    if not refined_bathymetry:
        [zlat, zlon, dpt1, var] = read_dpt(file_bathy) # load bathymetry file
        [lon_min, lon_max, lat_min, lat_max] = extent

        #convert to xarray
        dpt1 = xr.DataArray(dpt1, coords={'latitude':zlat, 'longitude': zlon}, dims=['latitude', 'longitude'])
        zlat = xr.DataArray(zlat, coords={'latitude':zlat}, dims=['latitude'])
        zlon = xr.DataArray(zlon, coords={'longitude':zlon}, dims=['longitude'])
        dpt1 = dpt1.sel(latitude = slice(lat_min, lat_max), longitude = slice(lon_min, lon_max))
        zlat = zlat.sel(latitude = slice(lat_min, lat_max))
        zlon = zlon.sel(longitude = slice(lon_min, lon_max))
        
    if refined_bathymetry:
        # load refined bathymetry ETOPOv2
        file_bathy = '/Users/tomasetl/Documents/code/ocean_source/pr_oceanswell/ocean_ampli/ETOPO_2022_v1_60s_N90W180_bed.nc'
        ds = xr.open_mfdataset(file_bathy, combine='by_coords')
        ds  = ds.rename({'lon':'longitude', 'lat': 'latitude'})
        #zlat = ds['latitude']
        #zlon = ds['longitude']
        z = ds['z']
        z *= -1 # ETOPOv2 to Depth
        z = z.where(z>0, other=np.nan)
        dpt1 = z.sel(latitude = slice(lat_min, lat_max), longitude = slice(lon_min, lon_max))
    
    ## Mask nan values
    dpt1_mask = dpt1.where(np.isfinite(dpt1))
    print(" dpt1 shape ", dpt1.shape)
    zlon = dpt1_mask.longitude
    zlat = dpt1_mask.latitude
    return dpt1_mask, zlon, zlat
def loop_SDF(paths, dpt1, zlon, zlat, date_vec=[2020, [], [], []], extent=[-180, 180, -90, 90],parameters= [2.8, 2830, 1/12, 0.2], **kwargs):
    """
    Input:
	paths = [file_bathy, ww3_local_path, longuet_higgins_file]:
					paths of additional files bathymetry, ww3 p2l file and Longuet-Higgins coefficients
	dpt1 : bathymetry grid in m (depth) with dimensions lon x lat
	zlon : longitude of bathymetry file (°)
	zlat : latitude of bathymetry file (°)
	date_vect : date vector [year, month, day, hour], with hour in [0, 3, 6, 9, 12, 15, 18, 21]
	extent : spatial extent format [lon_min, lon_max, lat_min, lat_max], default: [-180, 180, -90, 90]
	parameters : parameters vS of the crust, density of the crust, minimum frequency, maximum frequency,
								default: [2.8, 2830, 1/12, 0.2]
	**kwargs : catalog of other potential arguments
						plot : Bool default: True
									plot_hourly: Bool,  plot maps every 3-hours default: False
									plot_daily: Bool, plot maps every day, default: False
									plot_monthly: Bool, plot maps every month, default : True
									plot_yearly: Bool, plot map for the year average, default: False
						save: Bool, save 3-hourly matrix, default: False
    Output:
			None
			Saves in netcdf format power spectrum of the vertical displacement for Rayleigh waves,
			if save argument True.
			Plots in PNG source maps of Rayleigh waves at given intervals depending on plot variables. 
    """
    file_bathy = paths[0]
    ww3_local_path = paths[1]
    path_longuet_higgins = paths[2]
    
    # Constants
    radius = 6.371*1e6 # radius of the earth in meters
    lg10 = log(10) # log of 10
    #
    vs_crust = parameters[0]
    rho_s = parameters[1]
    f1 = parameters[2]
    f2 = parameters[3]
    ## Initialize variables
    if 'plot_type' in kwargs:
        plot_type = kwargs['plot_type']
        plot = True
        if plot_type == 'hourly':
            plot_hourly = True
            F_hourly =  np.zeros(dpt1.shape)
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
    else:
        plot = True
    
    if 'save' in kwargs:
        save = kwargs['save']
    else:
        save = False

    
    ## Adapt latitude and longitude to values in parameters file
    lon_min, lon_max, lat_min, lat_max = extent[0], extent[1], extent[2], extent[3]
    dpt1 = dpt1.sel(latitude = slice(lat_min, lat_max), longitude = slice(lon_min, lon_max))
    zlat = zlat.sel(latitude = slice(lat_min, lat_max))
    zlon = zlon.sel(longitude = slice(lon_min, lon_max))
    
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
            TOTAL_month = np.zeros(dpt1.shape)  # Initiate monthly source of Rayleigh wave matrix
            daymax = monthrange(iyear,imonth)[1]
            filename_p2l = '%s/LOPS_WW3-GLOB-30M_%d%02d_p2l.nc'%(ww3_local_path, iyear, imonth)
            # if imonth < 10:
            #     filepath = '%s/LOPSWW3-GLOB-30M_%d0%d'%(ww3_local_path, iyear, imonth)  # path to wave model files
            # else:
            #     filepath = '%s/WW3-GLOB-30M_%d%d'%(ww3_local_path, iyear, imonth)  # path to wave model files
            # filename_p2l = '%s_p2l.nc'%(filepath)
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
                SDF_day = np.empty(dpt1.shape)
                if isinstance(HOUR, int):
                    HOUR = np.array([HOUR])
                elif not len(HOUR):
                    HOUR = np.arange(0,24,3)
                else:
                    HOUR = np.array(HOUR)
                for ih in HOUR:
                    
                    ## Open F_p3D 
                    (lati, longi, freq_ocean, p2l, unit1) = read_p2l(filename_p2l, [iyear, imonth, iday, ih], [lon_min, lon_max], [lat_min, lat_max])
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
                        
                        #print('integral over frequency range f1 = %.3f, f2 = %.3f'%(f1, f2))
                        Fp = p2l.sel(frequency = freq_ocean[index_freq], latitude = slice(lat_min, lat_max), longitude = slice(lon_min, lon_max))
                        C = site_effect(dpt1, freq_seismic, vs_crust, path_longuet_higgins, zlat, zlon)  # computes Longuet-Higgins site effect given the bathymetryint(Fp.shape)
                        if C.shape == Fp.shape:
                            SDF_f = 2*np.pi/(rho_s**2*(vs_crust*1e3)**5)*C.data*Fp.data
                        else:
                            Fp = Fp.interp(latitude = zlat, longitude = zlon)
                            SDF_f = 2*np.pi/(rho_s**2*(vs_crust*1e3)**5)*C.data*Fp.data
                        if SDF_f.shape != C.shape:
                            print('SDF shape', SDF_f.shape)
                            return

                        ## Loop over frequencies of interest
                        M = np.zeros((n_freq, len(dpt1), len(dpt1[0])))
                        for i, fq in enumerate(freq_seismic):
                            SDF_f[i, :, :] = fq.data*SDF_f[i, :, :]
                            M[i, :, :] = df.data[i]*SDF_f[i, :, :]
                        M = xr.DataArray(M, coords={'frequency': freq_seismic, 'latitude': zlat,'longitude': zlon}, dims=["frequency", "latitude", "longitude"])
                        SDF = M.sum(dim = 'frequency', skipna = True)
                         
                    ## Single frequency
                    elif f1 == f2:
                        index_freq = np.squeeze(np.argmin(abs(freq-f1)))
                        print('unique frequency ', f1)
                        Fp = p2l[:, :, index_freq]
                        C = site_effect(dpt1, f1, vs_crust, path_longuet_higgins)
                        SDF_f = 2*np.pi*f1/(rho_s**2*(vs_crust*1e3)**5)*Fp.data*C.data
                        SDF = SDF_f
                        
                    ## Exception in parametrization of frequencies
                    else:
                        print('two frequencies with f2 < f1 were given')
                        return

                    ## Save SDF to file
                    if save == True:
                        ## Save SDF to file
                        print(os.getcwd())
                        path_out = './SDF/'
                        try:
                            os.path.exists(path_out)
                        except:
                            os.makedirs(path_out)  
                        # Create netCDF 3-hourly file
                        ncfile = Dataset(path_out+"rayleigh_SDF_%d%02d%02d%02d.nc"%(iyear, imonth, iday, ih), mode='w',format='NETCDF4_CLASSIC')
                        lat_dim = ncfile.createDimension('latitude', len(zlat))  # latitude axis
                        lon_dim = ncfile.createDimension('longitude', len(zlon))  # longitude axis
                        time_dim = ncfile.createDimension('time', daymax*8)  # unlimited axis (can be appended to).
                        freq_dim = ncfile.createDimension('frequency', n_freq)
                        for dim in ncfile.dimensions.items():
                            print(dim)
                        ncfile.title='Rayleigh waves power spectrum of vertical displacement on %d-%02d-%02d-%02d'%(iyear, imonth, iday, ih)
                        ncfile.subtitle='Equivalent Force maps every 3 hours for the secondary microseismic peak'
                        lat = ncfile.createVariable('latitude', np.float32, ('latitude',))
                        lat.units = 'degrees_north'
                        lat.long_name = 'latitude'
                        lon = ncfile.createVariable('longitude', np.float32, ('longitude',))
                        lon.units = 'degrees_east'
                        lon.long_name = 'longitude'
                        time = ncfile.createVariable('time', np.float32, ('time',))
                        time.units = 'hours since 1990-01-01'
                        time.long_name = 'time'
                        freq_nc = ncfile.createVariable('frequency', np.float32, ('frequency',))
                        freq_nc.units = 'Hz'
                        freq_nc.long_name = 'frequency'
                        sdf_f = ncfile.createVariable('sdf_f', np.float32, ('frequency','latitude', 'longitude'))
                        sdf_f.units = 'm.s'
                        sdf_f.long_name = 'Power spectrum of vertical displacement'
                        lat[:] = zlat
                        lon[:] = zlon
                        freq_nc[:] = freq_seismic
                        sdf_f[:] = SDF_f
                        time[:] = date2num(datetime(iyear, imonth, iday, ih), units='hours since 1990-01-01', calendar='standard')
                        ncfile.close()
                        
                    ## Plot SDF
                    if plot_hourly == True:
                        plt.close('all') 
                        SDF_plot = xr.DataArray(SDF, 
                                                coords={'latitude': zlat,'longitude': zlon}, 
                                                dims=["latitude", "longitude"],
                                                name = 'Source of the power spectrum for the vertical displacement.\nRayleigh waves.\nFrequency %.3f-%.3f Hz.%d-%02d-%02dT%02d'%(f1, f2, iyear, imonth, iday, ih))
                        fig = plt.figure(figsize=(9,6))
                        fig.suptitle('Source of the power spectrum for the vertical displacement. Rayleigh waves.\nFrequency %.3f-%.3f Hz.\n%d-%02d-%02dT%02d'%(f1, f2, iyear, imonth, iday, ih))
                        ax = plt.axes(projection=ccrs.Robinson())
                        ax.coastlines()
                        gl = ax.gridlines()
                        gl.xformatter = LONGITUDE_FORMATTER
                        gl.yformatter = LATITUDE_FORMATTER
                        ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k', facecolor='linen')
                        SDF_plot.plot(ax=ax, transform=ccrs.PlateCarree(),  cbar_kwargs={'label':'SDF (m)', 'orientation': 'horizontal'}) #vmin=0, vmax=3e-17,
                        plt.savefig('rayleigh_SDF_%d%02d%02dT%02d.png'%(iyear, imonth, iday, ih), dpi = 300, bbox_inches='tight')

                    ## Sum SDF
                    if plot_daily == True:
                        SDF_daily += SDF
                    if plot_monthly == True:
                        SDF_monthly += SDF
                    if plot_yearly == True:
                        SDF_yearly += SDF
                   
                if plot_daily == True:
                    plt.close('all')
                    SDF_plot = xr.DataArray(SDF_daily,
                                            coords={'latitude': zlat,'longitude': zlon}, 
                                            dims=["latitude", "longitude"])
                    fig = plt.figure(figsize=(9,6))
                    fig.suptitle('Source of the power spectrum for the vertical displacement. Rayleigh waves.\nFrequency %.3f-%.3f Hz %d-%02d-%02d'%(f1, f2, iyear, imonth, iday))
                    ax = plt.axes(projection=ccrs.Robinson())
                    ax.coastlines()
                    gl = ax.gridlines()
                    gl.xformatter = LONGITUDE_FORMATTER
                    gl.yformatter = LATITUDE_FORMATTER
                    ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k', facecolor='linen')
                    SDF_plot.plot(ax=ax, transform=ccrs.PlateCarree(), vmin = 0, vmax = 1e-16, cbar_kwargs={'label':'SDF (m)', 'orientation': 'horizontal'})
                    plt.savefig('rayleigh_SDF_daily_%d%02d%02d.png'%(iyear, imonth, iday), dpi = 300, bbox_inches='tight')
                    plt.close('all')
                    SDF_daily = np.zeros((dpt1.shape))
                    
            if plot_monthly == True:
                    plt.close('all')
                    SDF_plot = xr.DataArray(SDF_monthly,
                                            coords={'latitude': zlat,'longitude': zlon}, 
                                            dims=["latitude", "longitude"])
                    fig = plt.figure(figsize=(9,6))
                    fig.suptitle('Source of the power spectrum for the vertical displacement. Rayleigh waves.\nFrequency %.3f-%.3f Hz %d-%02d'%(f1, f2, iyear, imonth))
                    ax = plt.axes(projection=ccrs.Robinson())
                    ax.coastlines()
                    gl = ax.gridlines()
                    gl.xformatter = LONGITUDE_FORMATTER
                    gl.yformatter = LATITUDE_FORMATTER
                    ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k', facecolor='linen')
                    SDF_plot.plot(ax=ax, transform=ccrs.PlateCarree(), vmin = 0, vmax = 1e-15, cbar_kwargs={'label':'SDF (m)','orientation': 'horizontal'})
                    plt.savefig('rayleigh_SDF_monthly_%d%02d.png'%(iyear, imonth), dpi = 300, bbox_inches='tight')
                    #plt.show()
                    plt.close('all')
                    SDF_monthly = np.zeros((dpt1.shape))

        if plot_yearly == True:
                    plt.close('all')
                    SDF_plot = xr.DataArray(SDF_yearly,
                                            coords={'latitude': zlat,'longitude': zlon}, 
                                            dims=["latitude", "longitude"])
                    fig = plt.figure(figsize=(9,6))
                    ax = plt.axes(projection=ccrs.Robinson())
                    ax.coastlines()
                    gl = ax.gridlines()
                    gl.xformatter = LONGITUDE_FORMATTER
                    gl.yformatter = LATITUDE_FORMATTER
                    ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k', facecolor='linen')
                    fig.suptitle('Source of the power spectrum for the vertical displacement.Rayleigh waves.\nFrequency %.3f-%.3f Hz %d'%(f1, f2, iyear))
                    SDF_plot.plot(ax=ax, transform=ccrs.PlateCarree(), vmin = 0, vmax = 1e-14, cbar_kwargs={'label':'SDF (m)','orientation': 'horizontal'})
                    plt.savefig('rayleigh_SDF_yearly_%d.png'%(iyear), dpi = 300, bbox_inches='tight')
                    #plt.show()
                    plt.close('all')
                    SDF_yearly = np.zeros((dpt1.shape))
        plt.close('all')
    print('Rayleigh source maps done!')
    
    
    
def spectrogram(path_netcdf, dates, lon_sta=-21.3268, lat_sta=64.7474, Q=200, U=1800):
    """
    Calculate spectrogram for a given path to NetCDF file, dates, and optional station coordinates and constants.

    :param path_netcdf: str, path to NetCDF file containing the power spectrum of vertical displacement (m.s)
    :param dates: list or array-like, dates for which to calculate the spectrogram
    :param lon_sta: float, optional, longitude of the station (default=-21.3268) II.BORG 
    :param lat_sta: float, optional, latitude of the station (default=64.7474) II.BORG
    :param Q: int, optional, attenuation factor constant (default=200)
    :param U: int, optional, group velocity of Rayleigh waves constant (default=3500 m/s)
    
    Plot spectrogram for a given path to NetCDF file, dates, and optional station coordinates and constants.
    :return: None
    """
    """
    Attenuation Quality Factor
    attenuation factor used in Ardhuin et al. (2011) for 
    KIP (Hawaii) few reflections old crust, Q=580
    BORG (Iceland) few reflection young crust, Q=200
    BKS (California USA), Q=88
    SSB (France) respectively, Q=260
    """
    ## Constants
    radius_earth = 6371e3  # m
    res_mod = 0.5  # resolution ww3 model
    
    # open NetCDF for dims
    dates = pd.DatetimeIndex(dates)
    year = dates[0].year
    month = dates[0].month
    day = dates[0].day
    hour = dates[0].hour
    
    file_netcdf = path_netcdf + 'rayleigh_SDF_%d%02d%02d%02d.nc'%(year, month, day, hour)
    print(file_netcdf)
    dset = xr.open_mfdataset(file_netcdf, combine='by_coords')
    zlon = dset.longitude
    zlat = dset.latitude
    freq = dset.frequency
    
    # calculate spherical surface elementary elements
    msin = np.array([np.sin(np.pi/2 - np.radians(zlat))]).T
    ones = np.ones((1, len(zlon)))
    dA = radius_earth**2*res_mod**2*np.dot(msin,ones)
    
    # Compute distance of each gridpoint to station
    geoid = Geod(ellps='WGS84')
    lat_grid, lon_grid = np.meshgrid(zlat, zlon)
    lon_STA = np.ones((lon_grid.shape))*lon_sta
    lat_STA = np.ones((lat_grid.shape))*lat_sta
    _, _, distance_in_m = geoid.inv(lon_STA, lat_STA, lon_grid, lat_grid)
    distance = distance_in_m*180/(np.pi*radius_earth)
    # Initiate spectrogram
    n_dates = len(dates)
    spectro = np.zeros((n_dates, len(freq)))
    
    ## Loop Over Dates
    for idate, date in tqdm(enumerate(dates)):

        ## Open the netCDF File
        year = date.year
        month = date.month
        day = date.day
        hour = date.hour

        file_netcdf = path_netcdf + 'rayleigh_SDF_%d%02d%02d%02d.nc'%(year, month, day, hour)
        ds = xr.open_mfdataset(file_netcdf, combine='by_coords')
        ## Compute RMS Displacement
        sdf_f = ds.sdf_f
        F_delta = np.zeros((sdf_f.shape))
        for ifreq, f in enumerate(freq):
            SDF_freq = sdf_f.sel(frequency = freq[ifreq]).data
            EXP = np.exp(-2*np.pi*f.data*radius_earth/(U*Q)*np.radians(distance))
            denominateur = 1/(radius_earth*np.sin(np.radians(distance)))
            facteur = EXP*denominateur
            F_delta[ifreq, :, :] = facteur.T*SDF_freq*dA
        disp_RMS = 10*np.log(np.sqrt(np.nansum(F_delta, axis = (1,2))))
        spectro[idate, :] = disp_RMS
    
    ## Plot
    plt.figure(figsize=(16,9))
    plt.title("spectrogram")
    plt.pcolormesh(dates, freq, spectro.T, cmap='Spectral_r', vmin=-140, vmax=-80)
    plt.xlabel('Date ')
    plt.ylabel('Frequency [Hz]')
    plt.gcf().autofmt_xdate()
    plt.ylim(0.1, 0.5)
    plt.colorbar(label='$m^2/Hz$ [dB]')
    plt.show()
    return


def loop_ww3_sources(paths, dpt1, zlon, zlat, date_vec=[2020, [], [], []], extent=[-180, 180, -90, 90],parameters= [1/12, 1/2], **kwargs):
    """
    Input:
	paths = [file_bathy, ww3_local_path]: paths of additional files bathymetry, ww3 p2l file
	dpt1 : bathymetry grid in m (depth) with dimensions lon x lat
	zlon : longitude of bathymetry file (°)
	zlat : latitude of bathymetry file (°)
	date_vect : date vector [year, month, day, hour], with hour in [0, 3, 6, 9, 12, 15, 18, 21]
	extent : spatial extent format [lon_min, lon_max, lat_min, lat_max], default: [-180, 180, -90, 90]
	parameters : parameters minimum frequency, maximum frequency,
								default: [1/12, 1/2]
	**kwargs : catalog of other potential arguments
						plot : Bool default: True
									plot_hourly: Bool,  plot maps every 3-hours default: False
									plot_daily: Bool, plot maps every day, default: False
									plot_monthly: Bool, plot maps every month, default : True
									plot_yearly: Bool, plot map for the year average, default: False
						save: Bool, save 3-hourly matrix, default: False
    Output:
			None
			Saves in netcdf format the equivalent vertical force for each frequency if save argument True.
			Plots in PNG source maps of P/S waves at given intervals depending on plot variables.
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
    if 'plot_type' in kwargs:
        plot_type = kwargs['plot_type']
        plot = True
        if plot_type == 'hourly':
            plot_hourly = True
            F_hourly =  np.zeros(dpt1.shape)
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
    else:
        plot = True
    
    if 'save' in kwargs:
        save = kwargs['save']
    else:
        save = False

    ## Adapt latitude and longitude to values in parameters file
    lon_min = extent[0]
    lon_max = extent[1]
    lat_min = extent[2]
    lat_max = extent[3]
    
    dpt1 = dpt1.sel(latitude = slice(lat_min, lat_max), longitude = slice(lon_min, lon_max))
    zlat = zlat.sel(latitude = slice(lat_min, lat_max))
    zlon = zlon.sel(longitude = slice(lon_min, lon_max))
    
    ## Surface Element
    msin = np.array([np.sin(np.pi/2 - np.radians(zlat))]).T
    ones = np.ones((1, len(zlon)))
    dA = radius**2*res_mod**2*np.dot(msin,ones)
    
    ## Open Amplification Coefficient
    amplification_coeff = xr.open_dataarray('./C.nc')
    
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
            TOTAL_month = np.zeros(dpt1.shape)  # Initiate monthly source of Rayleigh wave matrix
            daymax = monthrange(iyear,imonth)[1]
            filename_p2l = '%s/LOPS_WW3-GLOB-30M_%d%02d_p2l.nc'%(ww3_local_path, iyear, imonth)
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
                    (lati, longi, freq_ocean, p2l, unit1) = read_p2l(filename_p2l, [iyear, imonth, iday, ih], [lon_min, lon_max], [lat_min, lat_max])
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
                        Fp = Fp.assign_coords({"freq": freq_seismic})
                        Fp = Fp.swap_dims({"frequency": "freq"})
                        Fp = Fp.drop('frequency')
                        Fp = Fp.rename({"freq": "frequency"})
                        
                        ## Amplification coefficient
                        amplification_coeff = amplification_coeff.sel(frequency = freq_seismic, method='nearest')
                        amplification_coeff = amplification_coeff.sel(latitude = slice(lat_min, lat_max), longitude = slice(lon_min, lon_max))
                        amplification_coeff_square = amplification_coeff**2
                        Fp_array = Fp.values
                        amplification_coeff_square_array = amplification_coeff_square.values
                        F_f = Fp_array * amplification_coeff_square_array
                        F = F_f
                        for ifq, fq in enumerate(freq_seismic.values):
                            df_fq = df.values[ifq]
                            F_fi = F_f[ifq, :, :]
                            F_fi *= dA
                            F[ifq, :, :] = F_fi*df_fq
                        F = 2*np.pi*np.sqrt(F.sum(axis = 0))
                    ## Single frequency
                    elif f1 == f2:
                        index_freq = np.squeeze(np.argmin(abs(freq-f1)))
                        print('unique frequency ', f1)
                        return
                        
                    ## Exception in parametrization of frequencies
                    else:
                        print('two frequencies with f2 < f1 were given')
                        return

                    ## Save F to file
                    if save == True:
                        print(os.getcwd())
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
                        ncfile.title='Equivalent Vertical Force on %d-%02d-%02d-%02d'%(iyear, imonth, iday, ih)
                        ncfile.subtitle='Equivalent Force every 3 hours for the secondary microseismic peak'
                        lat = ncfile.createVariable('latitude', np.float32, ('latitude',))
                        lat.units = 'degrees_north'
                        lat.long_name = 'latitude'
                        lon = ncfile.createVariable('longitude', np.float32, ('longitude',))
                        lon.units = 'degrees_east'
                        lon.long_name = 'longitude'
                        time = ncfile.createVariable('time', np.float32, ('time',))
                        time.units = 'hours since 1990-01-01'
                        time.long_name = 'time'
                        freq_nc = ncfile.createVariable('frequency', np.float32, ('frequency',))
                        freq_nc.units = 'Hz'
                        freq_nc.long_name = 'frequency'
                        F_freq = ncfile.createVariable('F_f', np.float32, ('frequency','latitude', 'longitude'))
                        F_freq.units = 'N.s^{1/2}'
                        F_freq.long_name = 'Equivalent Vertical Force spectrum'
                        lat[:] = zlat
                        lon[:] = zlon
                        freq_nc[:] = freq_seismic
                        F_freq[:] = F_f
                        time[:] = date2num(datetime(iyear, imonth, iday, ih), units='hours since 1990-01-01', calendar='standard')
                        ncfile.close()
                        
                    ## Plot F
                    if plot_hourly == True:
                        plt.close('all') 
                        F_plot = xr.DataArray(F, 
                                                coords={'latitude': zlat,'longitude': zlon}, 
                                                dims=["latitude", "longitude"],
                                                name = 'Equivalent Force. Rayleigh waves.\nFrequency %.3f-%.3f Hz.%d-%02d-%02dT%02d'%(f1, f2, iyear, imonth, iday, ih))
                        fig = plt.figure(figsize=(9,6))
                        fig.suptitle('Equivalent Force. Rayleigh waves.\nFrequency %.3f-%.3f Hz.%d-%02d-%02dT%02d'%(f1, f2, iyear, imonth, iday, ih))
                        ax = plt.axes(projection=ccrs.Robinson())
                        ax.coastlines()
                        gl = ax.gridlines()
                        gl.xformatter = LONGITUDE_FORMATTER
                        gl.yformatter = LATITUDE_FORMATTER
                        ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k', facecolor='linen')
                        F_plot.plot(ax=ax, transform=ccrs.PlateCarree(),  cbar_kwargs={'label':'F (N)', 'orientation': 'horizontal'}, vmin=0, vmax=3e10)
                        plt.savefig('F_R_%d%02d%02dT%02d.png'%(iyear, imonth, iday, ih), dpi = 300, bbox_inches='tight')

                    ## Sum F
                    if plot_daily == True:
                        F_daily += F
                    if plot_monthly == True:
                        F_monthly += F
                    if plot_yearly == True:
                        F_yearly += F
                   
                if plot_daily == True:
                    plt.close('all') 
                    F_plot = xr.DataArray(F_daily, 
                        coords={'latitude': zlat,'longitude': zlon}, 
                        dims=["latitude", "longitude"],
                        name = 'Equivalent Force. Rayleigh waves.Frequency %.3f-%.3f Hz.%d-%02d-%02d'%(f1, f2, iyear, imonth, iday))
                    fig = plt.figure(figsize=(9,6))
                    fig.suptitle('Equivalent Force. Rayleigh waves.\nFrequency %.3f-%.3f Hz.%d-%02d-%02d'%(f1, f2, iyear, imonth, iday))
                    ax = plt.axes(projection=ccrs.Robinson())
                    ax.coastlines()
                    gl = ax.gridlines()
                    gl.xformatter = LONGITUDE_FORMATTER
                    gl.yformatter = LATITUDE_FORMATTER
                    ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k', facecolor='linen')
                    F_plot.plot(ax=ax, transform=ccrs.PlateCarree(),  cbar_kwargs={'label':'F (N)', 'orientation': 'horizontal'}, vmin=0, vmax=3e10),
                    plt.savefig('F_R_%d%02d%02d.png'%(iyear, imonth, iday), dpi = 300, bbox_inches='tight')
                    F_daily = np.zeros((dpt1.shape))
                    
            if plot_monthly == True:
                plt.close('all')
                F_plot = xr.DataArray(F_monthly, 
                    coords={'latitude': zlat,'longitude': zlon}, 
                    dims=["latitude", "longitude"],
                    name = 'Equivalent Force. %s waves. Frequency %.3f-%.3f Hz.%d-%02d'%(wave_type,f1, f2, iyear, imonth))
                fig = plt.figure(figsize=(9,6))
                fig.suptitle('Equivalent Force. %s waves.\nFrequency %.3f-%.3f Hz.%d-%02d'%(wave_type, f1, f2, iyear, imonth))
                ax = plt.axes(projection=ccrs.Robinson())
                ax.coastlines()
                gl = ax.gridlines()
                gl.xformatter = LONGITUDE_FORMATTER
                gl.yformatter = LATITUDE_FORMATTER
                ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k', facecolor='linen')
                F_plot.plot(ax=ax, transform=ccrs.PlateCarree(),  cbar_kwargs={'label':'F (N)', 'orientation': 'horizontal'}, vmin=0, vmax=1e11)
                plt.savefig('F_R_%d%02d.png'%(iyear, imonth), dpi = 300, bbox_inches='tight')
                F_monthly = np.zeros((dpt1.shape))
                    

        if plot_yearly == True:
            plt.close('all')
            F_plot = xr.DataArray(F_yearly,
                                coords={'latitude': zlat,'longitude': zlon}, 
                                dims=["latitude", "longitude"],
                                name = 'Equivalent Force. %s waves.\nFrequency %.3f-%.3f Hz.%d'%(wave_type,f1, f2, iyear))
            fig = plt.figure(figsize=(9,6))
            fig.suptitle('Equivalent Force. %s waves.\nFrequency %.3f-%.3f Hz.%d'%(wave_type, f1, f2, iyear))
            ax = plt.axes(projection=ccrs.Robinson())
            ax.coastlines()
            gl = ax.gridlines()
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k', facecolor='linen')
            F_plot.plot(ax=ax, transform=ccrs.PlateCarree(),  cbar_kwargs={'label':'F (N)', 'orientation': 'horizontal'}, vmin=0, vmax=1e12)
            plt.savefig('F_%s_%d.png'%(wave_type, iyear), dpi = 300, bbox_inches='tight')
            F_daily = np.zeros((dpt1.shape))
            F_yearly = np.zeros((dpt1.shape))
        plt.close('all')
    print('Equivalent Vertical Force source maps done!')