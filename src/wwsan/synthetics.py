#!/usr/bin/env python3

__author__ = "Lisa Tomasetto"
__copyright__ = "Copyright 2024, UGA"
__credits__ = ["Lisa Tomasetto"]
__version__ = "0.1"
__maintainer__ = "Lisa Tomasetto"
__email__ = "lisa.tomasetto@univ-grenoble-alpes.fr"
__status__ = "Developpement"

"""Functions to generate synthetic cross-correlations from secondary microseisms sources."""

##########################################################################

## Librairies
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("ggplot")
import obspy

import h5py 
import xarray as xr
import os
import multiprocessing

from scipy.fftpack import fftfreq,ifft, fft
import scipy.signal as signal
from datetime import datetime
from pyproj import Geod
from tqdm import tqdm
from obspy.core.utcdatetime import UTCDateTime
from obspy.taup import TauPyModel
from obspy.taup.taup_geo import calc_dist as calc_dist

multiprocessing.set_start_method('fork')
## Constants
radius_earth = 6371e3  # Earth's Radius in meters

## Functions
def apply_delay_f_domain(s,d=0.):
    """
    Apply a delay d to the input signal s using Fourier domain.
    
    :param s: The input signal to apply delay.
    :type s: numpy.ndarray
    
    :param d: The delay to apply to the input signal (default 0).
    :type d: float
    
    :return: The delayed input signal.
    :rtype: numpy.ndarray
    """
    ## apply a delay d to the input signal s
    d_phase = 2*np.pi*d*np.arange(len(s))/(len(s))  # phase associated to d
    d_phase[int(np.round(len(d_phase)/2))::]-=2*np.pi*d  # for the hermitian symmetry of the fft
    return np.real(np.fft.ifft(np.fft.fft(s)*np.exp(-1j*d_phase)))

def get_synthetic_info(path_file_axisem='../../data/NOISE_vertforce_dirac_0-ak135f_1.s_3600s.h5', comp='Z'):
    """
        Function to get synthetic info from a specified h5 file and return fe, time, and N.
        Parameters:
        - path_file_axisem: str, default='../../data/NOISE_vertforce_dirac_0-ak135f_1.s_3600s.h5'
        - comp: str, default='Z'
        Returns:
        - fe: float
        - time: numpy.ndarray
        - N: int
    """
    h5_file = h5py.File(path_file_axisem)
    fe = h5_file['_metadata']['fe'][()]
    dist = 0.1
    trace = h5_file['L']['SYNTH%03d.00'%(dist*10)][comp][:].astype(np.single)
    time = np.arange(0, len(trace)*1/fe, 1/fe).astype(np.single)
    N = len(trace)
    h5_file.close()
    return fe, time, N

def open_axisem(dist, path_file_axisem='../../data/NOISE_vertforce_dirac_0-ak135f_1.s_3600s.h5', comp='Z'):
    """ Reads an AxiSEM .h5 archive in a 1D model given its path and a distance,
    and returns sampling frequency, time vector and trace at the given distance.
    Input:
    path_file_axisem: path of axisem .h5 archive
    dist : distance of interest
     
    Output:
    fe_iasp: sampling frequency 
    time_iasp: time vector of the archive
    trace_synth: synthetic trace at the given distance."""
    dist = np.round(dist, 1)
    h5_file = h5py.File(path_file_axisem)
    trace_synth = h5_file['L']['SYNTH%03d.00'%(dist*10)][comp][:].astype(np.single)
    h5_file.close()
    return trace_synth

def taper_axisem_archive(time, distance, archive_name='../../data/NOISE_vertforce_dirac_0-ak135f_1.s_3600s.h5', umin = 2.5, umax = 3.5):
    R = radius_earth*1e-3  # Radius of the Earth in km
    dt = time[1] - time[0]
    fe = 1/dt
    tapered_archive = np.zeros((len(distance), len(time)))
    for i, dist in enumerate(distance):
        dist_in_km = dist*np.pi*R/180
        tmin = dist_in_km/umax
        tmax = dist_in_km/umin
        if tmax >= np.max(time):
            tmax = np.max(time)
        if tmin >= np.max(time):
            continue
        dt_width = tmax - tmin + 50
        ## Taper
        tukey = signal.windows.tukey(int(round(dt_width*fe)), alpha=0.1)
        dirac = np.zeros(len(time))
        index = np.argmin(np.abs(time - (tmax + tmin)//2))
        dirac[index] = 1
        taper = signal.fftconvolve(dirac, tukey, mode='same')
        W = open_axisem(dist, archive_name)
        tapered_archive[i, :] = W*taper
    ## Save tapered archive as h5 file
    h5_name_tapered = archive_name.replace('s.h5', '_tapered.s.h5')
    h5_file = h5py.File(h5_name_tapered, 'w')
    h5_file.create_dataset('WAVEFORMS', data=tapered_archive)
    h5_file.create_dataset('time', data=time)
    h5_file.create_dataset('distance', data=distance)
    h5_file.close()
    return tapered_archive

def create_spectrum_archive(time, distance, tapered_archive):
    dt = time[1] - time[0]
    fe = 1/dt
    N = len(time)
    freq = fftfreq(N, 1/fe)
    spectrum = np.zeros((len(distance), len(freq))).astype(complex)
    for dist in tqdm(distance):
        index = np.argmin(np.abs(distance - dist))
        W = tapered_archive[index, :]
        spectrum[index, :] = fft(W)
    
    # Save as h5
    h5_name_spectrum = '../../data/spectrum_archive_tapered.h5'
    h5_file = h5py.File(h5_name_spectrum, 'w')
    h5_file.create_dataset('SPECTRUM', data=spectrum)
    h5_file.create_dataset('frequency', data=freq)
    h5_file.create_dataset('distance', data=distance)
    h5_file.close()
    
def open_archive(h5_name = 'spectrum_archive_tapered.h5'):
    h5_file = h5py.File(h5_name, 'r')
    S = h5_file['SPECTRUM'][:]
    distance = h5_file['distance'][:]
    freq = h5_file['frequency'][:]
    h5_file.close()
    fe = 2*np.max(freq)
    return fe, freq, distance, S

def open_spectrum_axisem(path_file_axisem='./spectrum_vertforce_iasp91_1.s_256c_3600.s.h5', comp='Z'):
    """ Reads an AxiSEM .h5 archive in a 1D model given its path and a distance,
    and returns sampling frequency, time vector and trace at the given distance.
    Input:
    path_file_axisem: path of axisem .h5 archive
    dist : distance of interest
     
    Output:
    fe_iasp: sampling frequency 
    time_iasp: time vector of the archive
    dist_iasp: distance vector of the archive
    spectrum_synth: synthetic spectra """
    h5_file = h5py.File(path_file_axisem, mode = 'r')
    fe_iasp = h5_file['_metadata']['fe'][()].astype(np.single)
    time_iasp = h5_file['_metadata']['time'][:].astype(np.single)
    N = len(time_iasp)
    freq_iasp = h5_file['_metadata']['freq'][:].astype(np.single)
    index_freq = np.squeeze(np.argwhere(np.logical_and(freq_iasp>=0, freq_iasp<=2.)))
    freq_iasp = freq_iasp[index_freq]
    dist_iasp = h5_file['_metadata']['dist'][()].astype(np.single)
    spectrum_synth = np.zeros((len(dist_iasp), len(freq_iasp))).astype(np.csingle)
    for i, distance in enumerate(dist_iasp):
        spectrum_synth[i,:] = h5_file['L']['SYNTH%03d.00'%(distance*10)][comp][index_freq].astype(np.csingle)
    h5_file.close() 
    spectrum_synth = xr.DataArray(spectrum_synth, coords={'distance':dist_iasp, 'frequency':freq_iasp})
    return fe_iasp, freq_iasp, time_iasp, dist_iasp, spectrum_synth

## Create Date Vect
def create_date_vect(dates):
    """Create date vector"""
    date_vect = np.zeros((len(dates), 4))
    for i in range(len(dates)):
        date = dates[i]
        year, month, day, hour = date.year, date.month, date.day, date.hour
        date_vect[i] = [year, month, day, hour]
    return date_vect

def open_model(path_file_WW3, date_vect, N, fe, lon_slice=slice(-180, 180), lat_slice=slice(-78, 80)):
    """ Reads the WW3 model for ambient noise sources at a specific time and date. 
    Returns the PSD of the model on the whole grid.
    Input:
    path_file_WW3: path of WW3 model archive in N.s^{1/2}
    date_vect: date and time vector [YEAR, MONTH, DAY, HOUR]

    Output:
    Returns Hermitian PSD of the given model in N^2.s with dimensions (dim_lon, dim_lat, 2*dim_freq+1).
    """
    year = date_vect[0]
    month = date_vect[1]
    day = date_vect[2]
    hour = date_vect[3]

    ## Open WW3 model 
    ds = xr.open_dataset(path_file_WW3)
    ww3_data = ds.F_f.sel(longitude=lon_slice, latitude=lat_slice)

    ww3_data = ww3_data.dropna(dim='longitude', how='all').dropna(dim='latitude', how='all')
    del ds
    ww3_data = ww3_data.where(np.isfinite(ww3_data))

    ## Spectrum Output
    if N%2 == 0:
        lenght_spectrum = int(N)  # must be Hermitian
    else:
        lenght_spectrum = int(N-1)  # must be Hermitian
    frequency = np.squeeze(fftfreq(lenght_spectrum, 1/fe))
    
    ## Interpolate over frequency range
    force_0 = xr.DataArray(np.zeros((len(ww3_data.longitude), len(ww3_data.latitude), 1)), coords={'longitude':ww3_data.longitude, 'latitude': ww3_data.latitude, 'frequency': [0]})
    force_2 = xr.DataArray(np.zeros((len(ww3_data.longitude), len(ww3_data.latitude), 1)), coords={'longitude':ww3_data.longitude, 'latitude': ww3_data.latitude, 'frequency': [0.7]})
    ww3_data = xr.concat([force_0, ww3_data, force_2], dim='frequency')
    force_spectrum = ww3_data.interp(frequency=frequency[np.logical_and(frequency>=0., frequency <= 2.)], method="linear", kwargs={"fill_value": 0.})
    del ww3_data

    return force_spectrum**2  # in N^2.s

def distance_to_station(lon, lat, lon_s=0, lat_s=90, radius_earth=6371e3):
    """ Computes the distance of every point of the model to station of coordinates (lonS, latS)
    Input:
    lon: longitudes of the grid
    lat: latitude of the grid
    lon_s: station longitude 
    lat_s: station latitude
    radius_earth : the radius of the Earth
    
    Output:
    matrix of dimensions dim(lon) x dim(lat) 
    """
    geoid = Geod(ellps="WGS84")  # set reference ellipsoid
    (lat_grid, lon_grid) = np.meshgrid(lat, lon)  # grid used

    lonSTA = np.ones(lon_grid.shape)*lon_s
    latSTA = np.ones(lat_grid.shape)*lat_s
    (_, _, distance_in_metersS) = geoid.inv(lon_grid, lat_grid, lonSTA, latSTA, radians=False)
    distanceS = 180.0*distance_in_metersS/(radius_earth*np.pi)
    distanceS = xr.DataArray(distanceS, coords={'longitude': lon, 'latitude': lat}, name = 'distance', attrs={'units': '°'})
    return distanceS

def matrix_GF(spectrum_axi, lon, lat, N, distance_s, comp = 'Z', conjugate = False):
    """
    Generates a synthetic seismogram (Green's Functions) matrix in frequency domain
    using the given lon, lat, N, path_file_axisem, distance_s, and optional conjugate flag.
    Returns the synthetic seismogram matrix.

    Input:
    lon: longitude of the grid
    lat: latitude of the grid
    N: number of samples
    path_file_axisem: path to the axisem tapered file
    distance_s: distance matrix
    conjugate: conjugate flag

    Output:
    synth: synthetic seismogram
    """
    S_synth = np.zeros((len(lon), len(lat), N//2)).astype(complex)
    distance_s = np.round(distance_s, decimals=1)
    min_dist, max_dist = np.min(distance_s.data), np.max(distance_s.data)
    distance = np.arange(min_dist, max_dist+0.1, 0.1)
    distance_synth = 0.1*np.arange(0, len(spectrum_axi), 1)
    for dist in distance:
        dist = np.round(dist, decimals=1)
        index = np.squeeze(np.argwhere(abs(distance_s.data - dist) < 0.09))
        try:
            trace = spectrum_axi[np.squeeze(np.argmin(abs(distance_synth - dist))),0:N//2].astype(complex)
        except:
            raise
        if conjugate:
            S_synth[index, :] = np.squeeze(np.conj(trace))
        else:
            S_synth[index, :] = np.squeeze((trace))
    return S_synth

def compute_model_chunk(lon_inf, lon_sup, lat_inf, lat_sup, N, fe, date_vect, spectrum_axi, file_model, lon_staA, lat_staA, lon_staB, lat_staB, comp):
            """
            Compute correlation function between stations A and B for a chunk of the www3 model source.   

            Parameters:
            lon_inf (float): Lower longitude bound.
            lon_sup (float): Upper longitude bound.
            lat_inf (float): Lower latitude bound.
            lat_sup (float): Upper latitude bound.
            N (int): Number of points.
            fe (float): Sampling frequency.
            date_vect (array): Vector of dates.
            file_model (str): Model file.
            lon_staA (float): Longitude of station A.
            lat_staA (float): Latitude of station A.
            lon_staB (float): Longitude of station B.
            lat_staB (float): Latitude of station B.
            path_file_axi (str): Path to Axisem file.

            Returns:
            np.ndarray: Computed correlation function as a single precision complex array.
            """
            ## Open model
            psd_model = open_model(path_file_WW3=file_model, date_vect=date_vect, N=N, fe=fe, lon_slice=slice(lon_inf, lon_sup), lat_slice=slice(lat_inf, lat_sup))
            psd_model.data = psd_model.data.astype(complex)
            lon = psd_model.longitude
            lat = psd_model.latitude
            ## Distances
            distance_staA = distance_to_station(lon, lat, lon_s = lon_staA, lat_s = lat_staA)
            distance_staB = distance_to_station(lon, lat, lon_s = lon_staB, lat_s = lat_staB)

            ## Green's Functions spectrum
            psd_model.values *= matrix_GF(spectrum_axi, lon=lon, lat=lat, N=N, distance_s=distance_staA, conjugate=True, comp=comp).astype(complex)
            
            ## Green's Functions spectrum
            psd_model.values *= matrix_GF(spectrum_axi, lon=lon, lat=lat, N=N, distance_s=distance_staB, conjugate=False, comp=comp).astype(complex)
            
            del distance_staA, distance_staB, lon, lat 
            ## Sum along coordinates latitude and longitude
            return psd_model.sum(dim=['longitude', 'latitude']).data.astype(complex)


def ccf_computation(coords_staA, coords_staB, path_model, date_vect, spectrum_axi, fe=4., N=14400, extent = [-180, 181, -80, 81], comp='Z'):
    """
    A function to compute the cross-correlation function between two seismic stations. 
    It takes the coordinates of the stations, the path to the AxiSEM archive, the path to the model, 
    a vector of dates, a normalization factor, and a component parameter. 
    It returns the computed cross-correlation function and the corresponding time array.
    """

    ## Coordinates of stations
    lon_staA = coords_staA[0]
    lat_staA = coords_staA[1]

    lon_staB = coords_staB[0]
    lat_staB = coords_staB[1]

    ## Open WW3 PSD
    YEAR = date_vect[0]
    MONTH = date_vect[1]
    DAY = date_vect[2]
    HOUR = date_vect[3]

    file_model = path_model + 'F_%d%02d%02d%02d.nc'%(YEAR, MONTH, DAY, HOUR)
    paramlist = []

    ## Define longitude and latitude slices
    number_of_cpu = os.cpu_count()
    lon_inf, lon_sup = extent[0], extent[1]
    lat_inf, lat_sup = extent[2], extent[3]
    
    corr_f = np.zeros((N)).astype(complex)
    res = compute_model_chunk(lon_inf, lon_sup, lat_inf, lat_sup, N, fe, date_vect, spectrum_axi, file_model, lon_staA, lat_staA, lon_staB, lat_staB, comp)
    
    corr_f[0:N//2] = res
    if N%2 == 0:
        corr_f[N//2::] = np.flip(np.conjugate(corr_f[0:N//2])).astype(complex)
    else:
        corr_f[N//2+1::] = np.flip(np.conjugate(corr_f[0:N//2])).astype(complex)
    ## Correlation in Time Domain
    corr = ifft(corr_f.data).real
    corr = np.fft.fftshift(corr)
    corr *= 1e-40
    time_corr = np.arange(-N//2, N//2)*1/fe

    corr = xr.DataArray(corr, dims=['time'], coords={'time': time_corr}, name='synthetic correlation')
    del corr_f    
    return corr, time_corr