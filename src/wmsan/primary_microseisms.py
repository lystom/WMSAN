#!/usr/bin/env python3

# Preamble
__author__ = "Fabrice Ardhuin" # mod. by Lisa Tomasetto 05/2026
__copyright__ = "Copyright 2026, CNES"
__credits__ = ["Fabrice Ardhuin","Lisa Tomasetto"] 
__version__ = "2026.1.0"
__maintainer__ = "Lisa Tomasetto"
__email__ = "lisa.tomasetto@partenaire-exterieur.ifremer.fr"

""" Functions """

import numpy as np
import xarray as xr
from datetime import datetime
from netCDF4 import Dataset
import pandas as pd
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
from pyproj import Geod
import fnmatch

from constants import R_E, g

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
plt.rcParams['font.family'] = "sans-serif"

#########################################################################
geoid = Geod(ellps='WGS84')


###################################################################################################################
############################################### SUBFUNCTIONS ######################################################
###################################################################################################################
def seismic_define_C2(readc, CgR):
    """
    Defines seismic response function (C2, nC2, c).
    
    Parameters
    ----------
    readc : int
        0 = analytic fit (Kedar et al. 2008 / Longuet-Higgins 1950)
        1 = read from Rayleigh_source.txt
    CgR : float
        Group speed. If < 0, read from Rayleigh_Cg.txt
    
    Returns
    -------
    C2 : ndarray
        Seismic response function (c^2)
    nC2 : int
        Number of discretized points
    c : ndarray
        
    """
    # --- Default version: analytic fit ---
    nC2 = 601
    x0 = np.linspace(0, 6, nC2)  # sig*h/beta, discretized
    
    a, b, c = 1.0, 10.0, 0.4
    c1 = 0.6 * (a + c * x0) / (a + b * (x0 - 0.85)**2) \
         + 0.28 / (2.0 + 1.0 * x0) \
         - 0.008 * (x0 - 2.5)

    a, b, c, d = 6.0, 30.0, 0.5, 10.0
    c2 = (0.30 * (a + c * x0) / (a + b * (x0 - 2.77)**2)
           + 0.40 / (1.0 + x0)
           + 0.006 * (x0 - 5)) * np.maximum(np.tanh((x0 - 1) * d), 0) \
         + 0.02 * np.exp(-(1.5 * (x0 - 3.6))**2) \
         - np.maximum(0.02 * (x0 - 5), 0)

    C2 = c1**2 + c2**2
    c0 = np.sqrt(c1**2)

    # --- Version 2: read from file ---
    if readc == 1:
        nmodes = 1
        all_c = np.loadtxt("Rayleigh_source.txt")
        nC, nmodein = all_c.shape[0], all_c.shape[1] - 1

        if nmodes > nmodein:
            nmodes = nmodein

        sighoverbetain = np.pi * all_c[:, 0]
        nC2 = 1001
        x0 = np.linspace(0, 10, nC2)
        C2 = np.zeros_like(x0)

        for i in range(nmodes):
            interp_func = interp1d(sighoverbetain, all_c[:, i+1],
                                   bounds_error=False, fill_value=np.nan)
            c = interp_func(x0)
            mask = np.isfinite(c)
            C2[mask] += c[mask]**2

        C2[0] = C2[1]
        c0 = np.sqrt(C2)

    # --- Group speed definition ---
    if CgR < 0:
        all_Cg = np.loadtxt("Rayleigh_Cg.txt")
        nCg, nmodein = all_Cg.shape[0], all_Cg.shape[1] - 1

        sighoverbetain = np.pi * all_Cg[:, 0]
        nC2 = 1001
        x0 = np.linspace(0, 10, nC2)

        # First mode (reference)
        interp_Cg = interp1d(sighoverbetain, all_Cg[:, 1],
                             bounds_error=False, fill_value=np.nan)
        CgRcmax = interp_Cg(x0)

        interp_c0 = interp1d(sighoverbetain, all_c[:, 1],
                             bounds_error=False, fill_value=np.nan)
        c0 = interp_c0(x0)

        # Loop over additional modes
        for i in range(1, nmodein):
            interp_c = interp1d(sighoverbetain, all_c[:, i+1],
                                bounds_error=False, fill_value=np.nan)
            interp_Cgi = interp1d(sighoverbetain, all_Cg[:, i+1],
                                  bounds_error=False, fill_value=np.nan)

            c = interp_c(x0)
            Cgi = interp_Cgi(x0)

            mask = np.isfinite(c) & np.isfinite(c0) & (c > c0)
            CgRcmax[mask] = Cgi[mask]
            c0[mask] = c[mask]

        # Fill NaNs with previous value
        isnan_mask = ~np.isfinite(CgRcmax)
        if np.any(isnan_mask):
            first_valid = np.argmax(np.isfinite(CgRcmax))
            CgRcmax[isnan_mask] = CgRcmax[first_valid]

    else:
        CgRcmax = np.zeros(nC2) + CgR

    return C2, nC2, c

def site_effect(z, f, zlat, zlon, vs_crust=2800, path='../../data/longuet_higgins.txt'):
    """ Bathymetry secondary microseismic excitation coefficients (Rayleigh waves).
    
    Args:
        z (np.ndarray): thickness of water layer, in m.
        f (np.ndarray): seismic frequency, in Hertz.
        vs_crust (float, optional): shear waves velocity in the crust (sea bed), in m/s.
        path (str, optional): the path to the Longuet Higgins file containing tabulated values of site effect coefficient for the 4th first modes.

    Returns:
        C (np.ndarray): Numerical value of the site effect in the shape (f.shape, z.shape).
    """

    df = pd.read_csv('%s'%path, sep='\t', header =0, usecols=[0, 1, 2, 3, 4, 5, 6, 7], names = ['fh1', 'c1', 'fh2', 'c2', 'fh3', 'c3', 'fh4', 'c4'])
    fc1 = interp1d(df.fh1, df.c1, kind='slinear', bounds_error=False, fill_value=0)
    fc2 = interp1d(df.fh2, df.c2, kind='slinear', bounds_error=False, fill_value=0)
    fc3 = interp1d(df.fh3, df.c3, kind='slinear', bounds_error=False, fill_value=0)
    fc4 = interp1d(df.fh4, df.c4, kind='slinear', bounds_error=False, fill_value=0)
    try :
        n = len(f)
        x = z.shape[1]
        y = z.shape[0]
        C = np.empty((n, y, x))
        for i, fq in enumerate(f):
            fh_v = 2*np.pi*fq*z/(vs_crust)
            C[i, :, :] = fc1(fh_v)**2 + fc2(fh_v)**2 + fc3(fh_v)**2 + fc4(fh_v)**2
        C = np.squeeze(C)
        ## C to xarray
        C = xr.DataArray(C, dims=('frequency','latitude', 'longitude'), coords={'frequency': f,'latitude': zlat, 'longitude': zlon})
    except:
        ## single frequency, f is a float
        fh_v = 2*np.pi*f*z/(vs_crust)
        C = fc1(fh_v)**2 + fc2(fh_v)**2 + fc3(fh_v)**2 + fc4(fh_v)**2
        C = np.squeeze(C)
        #C = xr.DataArray(C, dims=('latitude', 'longitude'), coords={'latitude': zlat, 'longitude': zlon})
    return C


#def dispNewtonTH(f, dep, eps=1e-6, max_iter=50):
#    """
#    Inverts the linear dispersion relation (2*pi*f)^2 = g*k*tanh(k*dep)
#    to get k from f and dep. Fully vectorized: `f` and `dep` can be
#    scalars or arrays of any (broadcastable) shape, e.g. f of shape (nf,)
#    and dep of shape (ny, nx) can be combined via
#    `dispNewtonTH(f[:, None, None], dep[None, :, :])` to get k of shape
#    (nf, ny, nx) without any Python loop.
#
#    If deep water (kh >= 6), use the deep-water approximation directly.
#    If shallow/intermediate water (kh < 6), use Newton iteration.
#
#    Parameters
#    ----------
#    f : array_like
#        Frequency in Hz (any shape).
#    dep : array_like
#        Water depth (m) (any shape broadcastable with f).
#    eps : float
#        Convergence tolerance on the iteration.
#    max_iter : int
#        Maximum number of Newton iterations.
#
#    Returns
#    -------
#    k : ndarray
#        Wavenumber, broadcast shape of f and dep.
#    """
#    f = np.asarray(f, dtype=float)
#    dep = np.asarray(dep, dtype=float)
#
#    sig = 2 * np.pi * f
#    # broadcast sig and dep together
#    sig, dep_b = np.broadcast_arrays(sig, dep)
#    dep_b = dep_b.copy()
#
#    Y = dep_b * sig**2 / g   # squared dimensionless frequency
#    X = np.sqrt(Y)           # initial guess valid for deep water
#
#    mask = X < 6
#
#    if np.any(mask):
#        Xm = X[mask]
#        Ym = Y[mask]
#        for _ in range(max_iter):
#            th = np.tanh(Xm)
#            f_val = Xm * th - Ym
#            df_val = th + Xm * (1 - th**2)
#            dX = f_val / df_val
#            Xm = Xm - dX
#            if np.max(np.abs(dX)) < eps:
#                break
#        X[mask] = Xm
#
#    # avoid division by zero for dep == 0
#    with np.errstate(divide='ignore', invalid='ignore'):
#        k = np.where(dep_b > 0, X / dep_b, 0.0)
#
#    return k
def dispNewtonTH(f, dep):
    """
    Inverts the linear dispersion relation (2*pi*f)^2 = g*k*tanh(k*dep)
    to get k from f and dep.
    
    Parameters
    ----------
    f : array_like
        Frequency in Hz
    dep : float
        Water depth (m)
        
    Returns
    -------
    k : ndarray
        Wavenumber corresponding to frequency and depth
    """
    g = 9.81
    eps = 1e-6

    f = np.atleast_1d(f)
    sig = 2 * np.pi * f
    Y = dep * sig**2 / g   # squared dimensionless frequency

    X = np.sqrt(Y)         # initial guess that fits deep water
    
    # mask: only apply Newton iteration where needed
    mask = X < 6

    if np.any(mask):
        F = np.ones_like(X[mask])
        while np.max(np.abs(F)) > eps:
            H = np.tanh(X[mask])
            F = Y[mask] - X[mask] * H
            FD = -H - X[mask] / np.cosh(X[mask])**2
            X[mask] = X[mask] - F / FD

    return X / dep

def open_bottom_topography_spectrum(path):
    with open(path, "r") as f:
        A = np.fromfile(f, sep=" ")
    nkbx = int(A[0])
    nkby = int(A[1])
    dkbx = float(A[2])
    dkby = float(A[3])
    nvals = nkbx * nkby
    flat = A[4:4 + nvals]
    if flat.size != nvals:
        raise ValueError(f"Expected {nvals} values for botspec, found {flat.size}")

    botspec = flat.reshape((nkby, nkbx), order="F").T  # -> shape (nkbx, nkby)
    kbxmax = ((nkbx - nkbx%2)/2.0) * dkbx
    kbymax = ((nkby - nkby%2)/2.0) * dkby

    return botspec, kbxmax, kbymax

def open_wave_spectrum(wave_spectrum):
    ds = xr.open_dataset(wave_spectrum)
    freq = np.array(ds.variables['frequency'][:])
    theta = np.array(ds.variables['direction'][:])
    #print('theta:',theta)
    times = ds.variables['time'][:]
    lonp  = ds.variables['longitude'][0,0].values # assumes fixed position
    latp  = ds.variables['latitude'][0,0].values
    #dptpall = ds.variables['dpt'][:]
    efthall = ds.variables['efth']
    nf = len(freq)
    nd = len(theta)
    nt = len(times)
    
    xfr = np.exp(np.log(freq[nf - 1] / freq[0]) / (nf - 1))  # geometric progression factor
    df = freq * 0.5 * (xfr - 1.0 / xfr)                       # frequency intervals in wave model
    dth = (2*np.pi/nd)   # direction increment in rad
    
    omega = 2*np.pi*freq

    return times, lonp, latp, omega, efthall

def open_bathy(file_bathy = '../../data/LOPS_WW3-GLOB-30M_dataref_dpt.nc', refined_bathymetry=False, extent=[-180, 180, -90, 90]):
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
    ## check file bathymetry name

    ## DEFAULT
    if fnmatch.fnmatch(file_bathy, '*/LOPS_WW3-GLOB-30M_dataref_dpt.nc'):
        ds = xr.open_mfdataset(file_bathy, combine='by_coords')
        print("Use default bathymetry or download refined.")

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

    ## ETOPO
    elif fnmatch.fnmatch(file_bathy, '*/ETOPO_20??_*_bed.nc'):
        ds = xr.open_mfdataset(file_bathy, combine='by_coords')
        print("Use refined bathymetry.ETOPO.")
        ds  = ds.rename({'lon':'longitude', 'lat': 'latitude'})
        if lon_min > lon_max:
            ## work on the pacific ocean
            ds = ds.assign_coords(longitude=((360 + (ds.longitude % 360)) % 360))
            ds = ds.roll(longitude=int(len(ds['longitude']) / 2),roll_coords=True)
            lon_min = ((360 + (lon_min % 360)) % 360)
            lon_max = ((360 + (lon_max % 360)) % 360)
        ds = ds.sel(latitude = slice(lat_min, lat_max), longitude = slice(lon_min, lon_max))
        dpt1 = ds['z']
        dpt1 *= -1 # ETOPOv2 to Depth
        dpt1 = dpt1.where(dpt1>0, other=np.nan)
        ## Mask nan values    
        dpt1_mask = dpt1.where(np.isfinite(dpt1))
        zlon = dpt1_mask.longitude
        zlat = dpt1_mask.latitude
        return dpt1_mask, zlon, zlat

    ## GEBCO
    elif fnmatch.fnmatch(file_bathy, '*/GEBCO_20??_*.nc'):
        ds = xr.open_mfdataset(file_bathy, combine='by_coords')
        print("Use refined bathymetry.GEBCO.")
        ds  = ds.rename({'lon':'longitude', 'lat': 'latitude'})
        if extent[0] > extent[1]:
            ## work on the pacific ocean
            ds = ds.assign_coords(longitude=((360 + (ds.longitude % 360)) % 360))
            ds = ds.roll(longitude=int(len(ds['longitude']) / 2),roll_coords=True)
            lon_min = ((360 + (lon_min % 360)) % 360)
            lon_max = ((360 + (lon_max % 360)) % 360)
        ds = ds.sel(latitude = slice(lat_min, lat_max), longitude = slice(lon_min, lon_max))
        dpt1 = ds['elevation']
        dpt1 *= -1 # GEBCO to Depth
        dpt1 = dpt1.where(dpt1>0, other=np.nan)
        ## Mask nan values    
        dpt1_mask = dpt1.where(np.isfinite(dpt1))
        zlon = dpt1_mask.longitude
        zlat = dpt1_mask.latitude
        return dpt1_mask, zlon, zlat

    ## refined bathymetry, no file given
    elif refined_bathymetry:
        print("Use refined bathymetry.")
        try: ## ETOPO
            file_bathy = '../../data/ETOPO_2022_v1_60s_N90W180_bed.nc'
            ds = xr.open_mfdataset(file_bathy, combine='by_coords')
            print("ETOPO")
            ds  = ds.rename({'lon':'longitude', 'lat': 'latitude'})
            if extent[0] > extent[1]:
                ## work on the pacific ocean
                ds = ds.assign_coords(longitude=((360 + (ds.longitude % 360)) % 360))
                ds = ds.roll(longitude=int(len(ds['longitude']) / 2),roll_coords=True)
                lon_min = ((360 + (lon_min % 360)) % 360)
                lon_max = ((360 + (lon_max % 360)) % 360)
            ds = ds.sel(latitude = slice(lat_min, lat_max), longitude = slice(lon_min, lon_max))
            dpt1 = ds['z']
            dpt1 *= -1 # ETOPOv2 to Depth
            dpt1 = dpt1.where(dpt1>0, other=np.nan)
            ## Mask nan values    
            dpt1_mask = dpt1.where(np.isfinite(dpt1))
            zlon = dpt1_mask.longitude
            zlat = dpt1_mask.latitude
            return dpt1_mask, zlon, zlat
        
        except:
            try: # GEBCO
                file_bathy = '../../data/GEBCO_2026_sub_ice.nc'
                ds = xr.open_mfdataset(file_bathy, combine='by_coords')
                print("GEBCO")
                ds  = ds.rename({'lon':'longitude', 'lat': 'latitude'})
                if extent[0] > extent[1]:
                    ## work on the pacific ocean
                    ds = ds.assign_coords(longitude=((360 + (ds.longitude % 360)) % 360))
                    ds = ds.roll(longitude=int(len(ds['longitude']) / 2),roll_coords=True)
                    lon_min = ((360 + (lon_min % 360)) % 360)
                    lon_max = ((360 + (lon_max % 360)) % 360)
                ds = ds.sel(latitude = slice(lat_min, lat_max), longitude = slice(lon_min, lon_max))
                dpt1 = ds['elevation']
                dpt1 *= -1 # GEBCO to Depth
                dpt1 = dpt1.where(dpt1>0, other=np.nan)
                ## Mask nan values    
                dpt1_mask = dpt1.where(np.isfinite(dpt1))
                zlon = dpt1_mask.longitude
                zlat = dpt1_mask.latitude
                return dpt1_mask, zlon, zlat        

            except:
                print("Refined bathymetry GEBCO not found. \nYou can download it from:\n https://www.gebco.net/\nSave in ../data/")
            print("Refined bathymetry ETOPOv2 not found. \nYou can download it from:\n https://www.ngdc.noaa.gov/thredds/catalog/global/ETOPO2022/60s/60s_bed_elev_netcdf/catalog.html?dataset=globalDatasetScan/ETOPO2022/60s/60s_bed_elev_netcdf/ETOPO_2022_v1_60s_N90W180_bed.nc\nSave in ../data/")
            return None, None, None
    return dpt1_mask, zlon, zlat

def read_ef_from_url(start, end, prefix = 'CCI_WW3-GLOB-30M_', url='https://data-ww3.ifremer.fr/PROJECT/CCI/RUNS/GLOB-30M/', lon = (-180, 180), lat = (-90, 90)):
    ## if default url then use the following url
    if url == 'https://data-ww3.ifremer.fr/PROJECT/CCI/RUNS/GLOB-30M/':
        url = url+str(start[0])+ '/FIELD_NC/' +prefix+str(start[0])+str(start[1]).zfill(2)+'_ef.nc'
    else:
        url = url+prefix+str(start[0])+str(start[1]).zfill(2)+'_ef.nc'
    try:
        nc_ds.close()
        extract_ds.close()
    except:
        pass
    (lat_min, lat_max) = lat
    (lon_min, lon_max) = lon
    start = datetime(start[0], start[1], start[2], 0, 0, 0)
    end = datetime(end[0], end[1], end[2], 0, 0, 0)
    #load netcdf from url as netCDF4 dataset
    ncfile = Dataset(url+'#mode=bytes')
    #load netCDF4 dataset as Xarray dataset
    nc_ds = xr.open_dataset(xr.backends.NetCDF4DataStore(ncfile))
    # extract from Jan 1st to Jan 5th included
    extract_ds=nc_ds.sel(f=slice(0, 0.1))
    print(extract_ds.keys())
    extract_ds = extract_ds.sel(time=slice(start, end), latitude=slice(lat_min, lat_max), longitude=slice(lon_min, lon_max))

    return extract_ds.latitude.values, extract_ds.longitude.values, extract_ds.f[:].values, extract_ds.time[:].values, extract_ds.ef
    # end of read_ef_from_url function
########################################################################################################
######################################## MAIN FUNCTION ################################################
########################################################################################################
def compute_Fp1_ef_map(
        start, end,
        bottom_topography_spectrum,
        dpt, lon_min, lon_max, lat_min, lat_max,
        rhow=1026, ifmax=None, nd=24):
    """
    Computes the wave-bottom coupling spectrum Fp1(f, lat, lon) over a grid,
    vectorized over the (lat, lon) grid points for each frequency, following
    the matrix-based style of `loop_SDF` (subfunctions_rayleigh_waves.py).

    Returns
    -------
    Fp1 : ndarray, shape (nt, nf, ny, nx)
    botspec_map : ndarray, shape (nf, ny, nx)
    k_map : ndarray, shape (nf, ny, nx)
    freq : ndarray
    df : ndarray
    times : xr.DataArray
    """
    botspec, kbxmax, kbymax = open_bottom_topography_spectrum(bottom_topography_spectrum)
    nkbx, nkby = botspec.shape
    dkbx = 2 * kbxmax / (nkbx - nkbx % 2) if nkbx > 1 else 1.0
    dkby = 2 * kbymax / (nkby - nkby % 2) if nkby > 1 else 1.0

    dtor = np.pi / 180.0

    ## Open ef from online source
    lat, lon, freq, times, efall = read_ef_from_url(start=start, end=end, lat=(lat_min, lat_max), lon=(lon_min, lon_max))

    nt, nf, ny, nx = np.shape(efall)
    if ifmax is None:
        ifmax = nf

    xfr = np.exp(np.log(freq[nf - 1] / freq[0]) / (nf - 1))
    df = freq * 0.5 * (xfr - 1.0 / xfr)

    theta = np.linspace(0, 360 - 360 / nd, nd)

    depth_sub = xr.open_dataarray(dpt)  # Assuming dpt_map is a path to a NetCDF file containing the depth data
    depth_sub = depth_sub.sel(latitude=slice(lat_min, lat_max), longitude=slice(lon_min, lon_max))
    ## size is (1, nys, nxs) because there is a singleton time dimension in the depth data
    depth_sub = depth_sub.isel(time=0)  # Remove the singleton time dimension
    valid = depth_sub.values > 1
    
    # valid = valid.values  # Convert to numpy array for indexing

    # --- wavenumber map (vectorized over sub-grid, looped over frequency) ---
    k_map = np.ones((nf, ny, nx))
    k_sub = np.ones((nf,) + depth_sub.shape)
    for i in range(nf):
        #k_sub[i][valid] = dispNewtonTH(np.full(np.sum(valid), freqp[i]),
        #                                depth_sub.values[valid])[:, 0] if False else \
        #                   np.array([dispNewtonTH(freqp[i], d)[0] for d in depth_sub.values[valid]])
        ## compute dispNewton for valid grid, dispNewton works with matrix 
        k_sub[i][valid] = np.array([dispNewtonTH(freq[i], d)[0] for d in depth_sub.values[valid]])
    # place back into full map
    for i in range(nf):
        k_map[i][np.ix_(np.where((lat >= lat_min) & (lat <= lat_max))[0], np.where((lon >= lon_min) & (lon <= lon_max))[0])] = k_sub[i]

    # --- direction-averaged bottom spectrum, vectorized over sub-grid per (freq, dir) ---
    botspec_map = np.zeros((nf, ny, nx))
    botspec_sub = np.zeros((nf,) + depth_sub.shape)
    for i in range(ifmax):
        print('Computing bottom spectra for frequency:', i, freq[i])
        k0 = k_sub[i]  # (nys, nxs)
        for j in range(nd):
            kbx = k0 * np.sin(theta[j] * dtor)
            kby = k0 * np.cos(theta[j] * dtor)

            kbotxi = (kbxmax + kbx) / dkbx
            kbotyi = (kbymax + kby) / dkby

            ibk = np.clip(np.floor(kbotxi).astype(int), 0, nkbx - 2)
            jbk = np.clip(np.floor(kbotyi).astype(int), 0, nkby - 2)
            xbk = kbotxi - np.floor(kbotxi)
            ybk = kbotyi - np.floor(kbotyi)

            botspeci = (
                (botspec[ibk, jbk] * (1 - ybk) + botspec[ibk, jbk + 1] * ybk) * (1 - xbk)
                + (botspec[ibk + 1, jbk] * (1 - ybk) + botspec[ibk + 1, jbk + 1] * ybk) * xbk
            )
            botspec_sub[i] += np.where(valid, botspeci / nd, 0.0)
        botspec_map[i][np.ix_(np.where((lat >= lat_min) & (lat <= lat_max))[0], np.where((lon >= lon_min) & (lon <= lon_max))[0])] = botspec_sub[i]

    # --- alpha (slope-coupling) coefficient ---
    alphas_sub = np.full((nf,) + depth_sub.shape, -1.0)
    for i in range(ifmax):
        
        om = 2 * np.pi * freq[i]
        k0 = k_sub[i]
        kh = k0 * depth_sub
        shallow = valid & (kh < 7)

        sinh2kh = np.sinh(2 * kh)
        dkdD = -2 * k0**2 / (2 * kh + sinh2kh)
        dkDdD = k0 * np.sinh(2 * kh) / (2 * kh + sinh2kh)

        Cg = om / k0 * 0.5 * (1 + (2 * kh) / sinh2kh)
        dCgdD = (om / k0 * (dkDdD / sinh2kh
                 - 2 * kh * dkDdD * np.cosh(2 * kh) / sinh2kh**2)
                 - dkdD * om / k0**2 * 0.5 * (1 + (2 * kh) / sinh2kh))
        alpha1 = -(k0 + depth_sub * dkdD) * np.tanh(kh) / k0
        alpha2 = -0.5 * dCgdD / (Cg * k0)
        alpha3 = -dkdD / (k0**2)
        alphas_sub[i] = np.where(shallow, alpha1 + alpha2 + alpha3, -1.0)

    # --- wave-bottom pressure coupling map ---
    Eftop1f_sub = np.zeros((nf,) + depth_sub.shape)
    for i in range(ifmax):
        k0 = k_sub[i]
        kh = k0 * depth_sub
        Eftop1f_sub[i] = np.where(
            valid,
            botspec_sub[i] * (rhow * g * k0 * alphas_sub[i] / np.cosh(kh))**2 * 3,
            0.0
        )

    # --- couple with 1D wave energy spectrum map, all time steps (vectorized) ---
    Fp1 = np.zeros((nt, nf, ny, nx))
    for it in range(nt):
        if np.mod(it, 6) == 0:
            print('time:', it, times[it])
        ef = 10**(efall.sel(time=times[it], f=freq[0:ifmax], latitude=slice(lat_min, lat_max), longitude=slice(lon_min, lon_max))) - 1e-12
        Fp1[it, 0:ifmax, :, :] = Eftop1f_sub[0:ifmax] * ef

    return Fp1, botspec_map, k_map, freq, df, times


def compute_F_delta_ef_map(
        start, end,
        bottom_topography_spectrum,
        dpt, lon_min, lon_max, lat_min, lat_max, CgR, Q,
        lono, lato, lon, lat, statname,
        rhow=1026, rhos=2600, betas=2800,
        site_effect_path='../../data/longuet_higgins.txt',
        ifmax=None):
    """
    Computes the primary microseism source spectrum F_delta(t, f) over a
    grid, fully vectorized over the (lat, lon) sub-grid at each frequency,
    using `site_effect` (matrix version, as in `loop_SDF`) in place of the
    analytic C2 fit.
    """
    Fp1, botspec_map, k_map, freq, df, times = compute_Fp1_ef_map(
        start=start, end=end, bottom_topography_spectrum=bottom_topography_spectrum, dpt=dpt,
        lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max, rhow=rhow, ifmax=ifmax)

    nt, nf, ny, nx = Fp1.shape
    if ifmax is None:
        ifmax = nf

    dtor = np.pi / 180.0
    omega = 2 * np.pi * freq

    lon_sub = lon[(lon >= lon_min) & (lon <= lon_max)]
    lat_sub = lat[(lat >= lat_min) & (lat <= lat_max)]
    lon_grid, lat_grid = np.meshgrid(lon_sub, lat_sub)  # (nys, nxs)

    depth_sub = xr.open_dataset(dpt)
    depth_sub = depth_sub.sel(latitude=slice(lat_min, lat_max), longitude=slice(lon_min, lon_max))['dpt']
    depth_sub = depth_sub.isel(time=0)  # Remove the singleton time dimension
    valid = depth_sub.values > 1

    # --- great-circle distance (alpha, rad), vectorized like `spectrogram` ---
    lon_STA = np.full(lon_grid.shape, lono)
    lat_STA = np.full(lat_grid.shape, lato)
    _, _, distance_in_m = geoid.inv(lon_STA, lat_STA, lon_grid, lat_grid)
    alpha_map = distance_in_m / R_E
    alpha_map = np.where(alpha_map <= 0, 1e-6, alpha_map)

    # --- surface element matrix dA, as in loop_SDF/spectrogram ---
    res_lon = abs(lon[1] - lon[0]) * dtor
    res_lat = abs(lat[1] - lat[0]) * dtor
    dA = R_E**2 * res_lon * res_lat * np.cos(lat_grid * dtor)

    factor1 = 2 * np.pi * (1 / rhos)**2 / (betas**5 * R_E)

    # --- site effect coefficient over the whole sub-grid & frequency band at once ---
    C_map = site_effect(depth_sub, freq[0:ifmax], lat_sub, lon_sub,
                         vs_crust=betas, path=site_effect_path)
    C_map = np.asarray(C_map)  # shape (ifmax, nys, nxs)

    coeff_all = np.zeros((nf,) + depth_sub.shape)
    for i in range(ifmax):
        coeff = (factor1 * omega[i] * C_map[i]) / np.sin(alpha_map)

        b = np.exp(-omega[i] * (2.0 * np.pi) * (R_E / (abs(CgR) * Q)))
        attenuation = (
            np.exp(-omega[i] * alpha_map * (R_E / (abs(CgR) * Q))) / (1.0 - b)
            + np.exp(-omega[i] * (2.0 * np.pi - alpha_map) * (R_E / (abs(CgR) * Q))) / (1.0 - b)
        )

        coeff_all[i] = np.where(valid, coeff * attenuation * dA, 0.0)

    F_delta = np.zeros((nt, nf))
    map_source = None
    for it in range(nt):
        source = coeff_all * Fp1[it]
        F_delta[it, 0:ifmax] = np.nansum(source, axis=(1, 2))
        if it == 9 * 24:
            map_source = source
    F_delta = xr.DataArray(F_delta, coords=[times, freq], dims=['time', 'frequency'])

    return F_delta, map_source, botspec_map, k_map, freq, df, times

if __name__ == '__main__':


    # For the following file, you can get it here: https://data-ww3.ifremer.fr/PROJECT/CCI/RUNS/GLOB-30M/2023/FIELD_NC/CCI_WW3-GLOB-30M_202308_ef.nc
    wave_spectrum_1D='/home/ltomaset/Documents/WMSAN_gitlab/microseisms_LOPS/data/waves/CCI_WW3-GLOB-30M_202308_ef.nc'
    depth_file='/home/ltomaset/Documents/WMSAN_gitlab/microseisms_LOPS/data/depth/LOPS_WW3-GLOB-30M_dataref_dpt.nc'
    bottom_topography_spectrum='/home/ltomaset/Documents/WMSAN_gitlab/microseisms_LOPS/data/bottom/spectrum_Ireland_shallow_rocks.bsp'

    ## dpt_map is the matrix from the depth map file used for the computation
    dpt_map, zlon, zlat = open_bathy(depth_file)
    longitude = dpt_map['longitude'].values
    latitude = dpt_map['latitude'].values
    ## focus on North Atlantic
    lon_min = -100
    lon_max = -10
    lat_min = 0
    lat_max = 60

    start=[2025,1,20]
    end=[2025,1,30]

    ## Open depth and plot it
    plt.figure(figsize=(10, 6))
    plt.pcolormesh(dpt_map, shading='auto')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Depth Map')
    plt.colorbar(label='Depth')
    plt.show()

    ## Open bottom topography spectrum and plot it
    botspec_map, kbxmax, kbymax = open_bottom_topography_spectrum(bottom_topography_spectrum)
    plt.figure(figsize=(10, 6))
    plt.pcolormesh(botspec_map, shading='auto')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Bottom Topography Spectrum')
    plt.colorbar(label='Spectrum')
    plt.show()

    exit()

    F_delta, map_source, botspec_map, k_map, freq, df, times = compute_F_delta_ef_map(start=start,
                                                                                        end=end,
                                                                                        bottom_topography_spectrum=bottom_topography_spectrum,
                                                                                        dpt=depth_file,
                                                                                        lon_min=lon_min,
                                                                                        lon_max=lon_max,
                                                                                        lat_min=lat_min,
                                                                                        lat_max=lat_max,
                                                                                        CgR=1800,Q=88, lono=4.542, lato=45.279,
                                                                                        lon=np.arange(-180, 180, 0.5),
                                                                                        lat=np.arange(-78, 80.5, 0.5),
                                                                                        statname='G.SSB')

    ## Plot F_delta as a function of time (x) and frequency (y)

    plt.figure(figsize=(10, 6))
    plt.pcolormesh(times, freq, 10*np.log10(F_delta.T), shading='auto')
    plt.xlabel('Time')
    plt.ylabel('Frequency')
    plt.title('F_delta as a function of time and frequency')
    plt.colorbar(label='F_delta (dB)')
    plt.show()  

    ### Plot map_source as a function of longitude (x) and latitude (y)
    #plt.figure(figsize=(10, 6))
    #plt.pcolormesh(map_source, shading='auto')
    #plt.xlabel('Longitude')
    #plt.ylabel('Latitude')
    #plt.title('Map Source')
    #plt.colorbar(label='Map Source')
    #plt.show()

    ## save F_delta as netcdf 
    ## xarray dataset

    F_delta_da = xr.DataArray(F_delta, coords=[times, freq], dims=['time', 'frequency'])
    F_delta_da.to_netcdf(f'F_delta_G.SSB_{start[0]}{start[1]:02d}{start[2]:02d}_{end[0]}{end[1]:02d}{end[2]:02d}.nc')