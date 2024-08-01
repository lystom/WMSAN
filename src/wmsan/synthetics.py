#!/usr/bin/env python3

#__author__ = "Lisa Tomasetto"
#__copyright__ = "Copyright 2024, UGA"
#__credits__ = ["Lisa Tomasetto"]
#__version__ = "0.1"
#__maintainer__ = "Lisa Tomasetto"
#__email__ = "lisa.tomasetto@univ-grenoble-alpes.fr"

"""Functions to generate synthetic cross-correlations from secondary microseisms sources.

It contains fourteen functions:

- `apply_delay_f_domain(s,d)`: Apply a delay d to the input signal s using Fourier domain.

- `get_synthetic_info(path_file_axisem, comp)`: Get synthetic info from a specified h5 file and return fe, time, and N.

- `open_axisem(dist, path_file_axisem, comp)`: Reads an AxiSEM .h5 archive in a 1D model given its path and a distance.

- `taper_axisem_archive(time, distance, archive_name, umin, umax)`: Taper an AxiSEM .h5 archive in a 1D model given its path and a distance and minimum and maximum Rayleigh wave velocity.

- `taper_axisem_archive_body_waves(time, distance, archive_name, phase1, phase2, model, **kwargs)`: Taper an AxiSEM .h5 archive around two body wave phases in a 1D model given its path and a distance.

- `create_spectrum_archive(time, distance, tapered_archive, h5_name_spectrum)`: Create a Green's function archive in frequency domain from an AxiSEM archive.

- `open_archive(h5_name)`: Open a Green's function archive in frequency domain.

- `open_spectrum_axisem(path_file_axisem, comp)`: Reads an AxiSEM .h5 archive in a 1D model given its path and a distance, and returns sampling frequency, time vector and trace.

- `create_date_vect(dates)`: Create a date vector from a list of dates.

- `open_model(path_file_WW3, date_vect, N, fe, lon_slice, lat_slice)`: Open WW3 model for ambient noise sources at a specific time and date.

- `distance_to_station(lon, lat, lon_s, lat_s, radius_earth)`: Computes the distance of every point of the model to station of coordinates (lonS, latS).

- `matrix_GF(spectrum_axi, lon, lat, N, distance_s, comp, conjugate)`: Generates a synthetic seismogram (Green's Functions) matrix to a specific station in frequency domain using the given lon, lat, N, path_file_axisem, distance_s, and optional conjugate flag.

- `compute_model_chunk(lon_inf, lon_sup, lat_inf, lat_sup, N, fe, date_vect, spectrum_axi, file_model, lon_staA, lat_staA, lon_staB, lat_staB, comp)`: Compute correlation function in frequency domain between stations A and B for a chunk of the www3 model source.

- `ccf_computation(coords_staA, coords_staB, path_model, date_vect, spectrum_axi, fe, N, extent, comp)`: Compute the cross-correlation function between two seismic stations.

"""

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
    """Apply a delay d to the input signal s using Fourier domain.
    
    Args:
        s (numpy.ndarray): The input signal to apply delay.
        d (float, optional): The delay to apply to the input signal.

    Returns:    
        s_delayed (numpy.ndarray): The delayed input signal.
    """
    
    ## apply a delay d to the input signal s
    d_phase = 2*np.pi*d*np.arange(len(s))/(len(s))  # phase associated to d
    d_phase[int(np.round(len(d_phase)/2))::]-=2*np.pi*d  # for the hermitian symmetry of the fft
    return np.real(np.fft.ifft(np.fft.fft(s)*np.exp(-1j*d_phase)))

def get_synthetic_info(path_file_axisem='../../data/NOISE_vertforce_dirac_0-ak135f_1.s_3600s.h5', comp='Z'):
    """Function to get synthetic info from a specified h5 file and return fe, time, and N.
    
    Args:
        path_file_axisem (str, optional): path of axisem .h5 archive.
        comp (str, optional): component.
    
    Returns:
        fe (float): sampling frequency
        time (numpy.ndarray): time vector
        N (int): number of points
    """
    try:
        h5_file = h5py.File(path_file_axisem)
        fe = h5_file['_metadata']['fe'][()]
        dist = 0.1
        trace = h5_file['L']['SYNTH%03d.00'%(dist*10)][comp][:].astype(np.single)
    except:
        print('File not found', path_file_axisem, "Download the file from: https://zenodo.org/records/11126562 \n", "Save in ../data/")
        return None, None, None  
    time = np.arange(0, len(trace)*1/fe, 1/fe).astype(np.single)
    N = len(trace)
    h5_file.close()
    return fe, time, N

def open_axisem(dist, path_file_axisem='../../data/NOISE_vertforce_dirac_0-ak135f_1.s_3600s.h5', comp='Z'):
    """Reads an AxiSEM .h5 archive in a 1D model given its path and a distance,
    and returns sampling frequency, time vector and trace at the given distance.
    
    Args:
        path_file_axisem (str, optional): path of axisem .h5 archive
        dist (float): distance of interest
     
    Returns:
        fe_iasp (float): sampling frequency 
        time_iasp (numpy.ndarray): time vector of the archive
        trace_synth (numpy.ndarray): synthetic trace at the given distance.
    """
    dist = np.round(dist, 1)
    try:
        h5_file = h5py.File(path_file_axisem)
        trace_synth = h5_file['L']['SYNTH%03d.00'%(dist*10)][comp][:].astype(np.single)
        h5_file.close()
    except:
        raise
        print('File not found', path_file_axisem, "Download the file from: https://zenodo.org/records/11126562 \n", "Save in ../data/")
        return None
    return trace_synth

def taper_axisem_archive(time, distance, archive_name='../../data/NOISE_vertforce_dirac_0-ak135f_1.s_3600s.h5', umin = 2.5, umax = 3.5):
    """Taper an AxiSEM .h5 archive in a 1D model given its path and distances as well as minimum and maximum Rayleigh wave velocity,and saves and returns the tapered archive.
    
    Args:
        time (numpy.ndarray): time vector
        distance (numpy.ndarray): distance vector
        archive_name (str, optional): path of axisem .h5 archive, default='../../data/NOISE_vertforce_dirac_0-ak135f_1.s_3600s.h5'
        umin (float, optional): minimum Rayleigh wave velocity
        umax (float, optional): maximum Rayleigh wave velocity
     
    Returns:
        tapered_archive (numpy.ndarray): tapered archive
    """
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

def taper_axisem_archive_body_waves(time, distance, archive_name='../../data/NOISE_vertforce_dirac_0-ak135f_1.s_3600s.h5', phase1 = 'P', phase2 = 'PP', model='ak135', **kwargs):
    """Taper an AxiSEM .h5 archive around two body wave phases in a 1D model given its path and a distance,
    and returns the tapered archive.
    
    Args:
        time (numpy.ndarray): time vector.
        distance (numpy.ndarray): distance vector.
        archive_name (str, optional): path of axisem .h5 archive
        phase1 (str, optional): first body wave phase.
        phase2 (str, optional): second body wave phase.
        model (str, optional): model to compute arrival times with TauP.
     
    Returns:
        fe_iasp (float): sampling frequency.
        time_iasp (numpy.ndarray): time vector of the archive.
        trace_synth (numpy.ndarray): synthetic trace at the given distance.
    """
    
    R = radius_earth*1e-3  # Radius of the Earth in km
    dt = time[1] - time[0]
    fe = 1/dt
    tapered_archive = np.zeros((len(distance), len(time)))
    
    ## IASP travel times TauP
    model_1D = TauPyModel(model='%s'%model)
    tP_IASP = []
    tPP_IASP = []
    dist = np.arange(0, 180.1, 0.1)
    for i, d in enumerate(dist):
        try:
            arr_phase1 = model_1D.get_travel_times(source_depth_in_km=0, distance_in_degree=d, phase_list=['%s'%phase1])
            t_phase1 = arr_phase1[0].time
        except:
            t_phase1 = np.nan
        
        try:    
            arr_phase2 = model_1D.get_travel_times(source_depth_in_km=0, distance_in_degree=d, phase_list=['%s'%phase2])
            t_phase2 = arr_phase2[0].time
        except:
            t_phase2 = np.nan
        dist_in_km = dist*np.pi*R/180
        # Window around TauP
        tukey_window = signal.windows.tukey(int(50*fe), alpha=0.3, sym=True)
        dirac = np.zeros(len(time))
        if t_phase1 >= 0:
            index = [int(t_phase1*fe), int(t_phase2*fe)]
        else:
            index = int(t_phase2*fe)
        dirac[index] = 1
        window = signal.fftconvolve(dirac, tukey_window, mode='same')
        W = open_axisem(d, archive_name)
        tapered_archive[i, :] = W*window
        
    ## Save tapered archive as h5 file
    h5_name_tapered = archive_name.replace('s.h5', '_tapered_%s_%s.h5'%(phase1, phase2))
    h5_file = h5py.File(h5_name_tapered, 'w')
    h5_file.create_dataset('WAVEFORMS', data=tapered_archive)
    h5_file.create_dataset('time', data=time)
    h5_file.create_dataset('distance', data=distance)
    h5_file.close()
    return tapered_archive

def create_spectrum_archive(time, distance, tapered_archive, h5_name_spectrum = 'spectrum_archive_tapered.h5'):
    """Create a Green's function archive in frequency domain from an AxiSEM archive.
    This function creates a spectrum archive by computing the Fourier transform of the input archive for each distance value.
    The resulting spectrum is saved as an HDF5 file with the specified name. The HDF5 file contains three datasets:
    
    'SPECTRUM': The spectrum data.
    
    'frequency': The frequency values.
    
    'distance': The distance values.
    
    Args:
        time (numpy.ndarray): An array of time values.
        distance (numpy.ndarray): An array of distance values.
        tapered_archive (numpy.ndarray): The Green's Function archive in time domain file path.
        h5_name_spectrum (str, optional): The name of the output HDF5 file.

    Returns:
        None (NoneType): The function saves the spectrum archive as an HDF5 file.
    """
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
    h5_file = h5py.File(h5_name_spectrum, 'w')
    h5_file.create_dataset('SPECTRUM', data=spectrum)
    h5_file.create_dataset('frequency', data=freq)
    h5_file.create_dataset('distance', data=distance)
    h5_file.close()
    
def open_archive(h5_name = 'spectrum_archive_tapered.h5'):
    """Open an AxiSEM .h5 archive in a 1D model given its path,
    and returns sampling frequency, time vector and trace.
    
    Args:
        h5_name (str, optional): path of axisem .h5 archive
        
    Returns:
        fe (float): sampling frequency
        freq (numpy.ndarray): frequency vector of the archive
        distance (numpy.ndarray): distance vector of the archive
        S (numpy.ndarray): synthetic spectra   
    """
    
    h5_file = h5py.File(h5_name, 'r')
    S = h5_file['SPECTRUM'][:]
    distance = h5_file['distance'][:]
    freq = h5_file['frequency'][:]
    h5_file.close()
    fe = 2*np.max(freq)
    return fe, freq, distance, S

def open_spectrum_axisem(path_file_axisem='./spectrum_vertforce_iasp91_1.s_256c_3600.s.h5', comp='Z'):
    """Reads an AxiSEM .h5 archive in a 1D model given its path and component,
    and returns sampling frequency, time vector and trace.
    
    Args:
        path_file_axisem (str, optional): path of axisem .h5 archive
        comp (float): component
     
    Returns:
        fe_iasp (float): sampling frequency 
        time_iasp (numpy.ndarray): time vector of the archive
        dist_iasp (numpy.ndarray): distance vector of the archive
        spectrum_synth (numpy.ndarray): synthetic spectra
    """
    
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
    """Creates a date vector from a list of dates.

    Args:
        dates (list): A list of datetime objects representing the dates.

    Returns:
        date_vect (numpy.ndarray): A 2D numpy array where each row represents a date and contains the year, month, day, and hour.
    """
    date_vect = np.zeros((len(dates), 4))
    for i in range(len(dates)):
        date = dates[i]
        year, month, day, hour = date.year, date.month, date.day, date.hour
        date_vect[i] = [year, month, day, hour]
    return date_vect

def open_model(path_file_WW3, date_vect, N, fe, lon_slice=slice(-180, 180), lat_slice=slice(-78, 80)):
    """Open WW3 model for ambient noise sources at a specific time and date. 
    Returns the PSD of the given model in N^2.s with dimensions (dim_lon, dim_lat, 2*dim_freq+1).

    Args:
        path_file_WW3 (str): Path of WW3 model archive in N.s^{1/2}.
        date_vect (numpy.ndarray): Date and time vector [YEAR, MONTH, DAY, HOUR].
        N (int): Length of the spectrum.
        fe (float): Sampling frequency.
        lon_slice (slice, optional): Slice of longitude to select.
        lat_slice (slice, optional): Slice of latitude to select.

    Returns:
        force_spectrum (xarray.DataArray): Hermitian PSD of the given model in N^2.s with dimensions (dim_lon, dim_lat, 2*dim_freq+1).
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
    force_2 = xr.DataArray(np.zeros((len(ww3_data.longitude), len(ww3_data.latitude), 1)), coords={'longitude':ww3_data.longitude, 'latitude': ww3_data.latitude, 'frequency': [2.]})
    ww3_data = xr.concat([force_0, ww3_data, force_2], dim='frequency')
    force_spectrum = ww3_data.interp(frequency=frequency[np.logical_and(frequency>=0., frequency <= 2.)], method="linear", kwargs={"fill_value": 0.})
    del ww3_data
    return force_spectrum**2  # in N^2.s

def distance_to_station(lon, lat, lon_s=0, lat_s=90, radius_earth=6371e3):
    """Computes the distance of every point of the model to station of coordinates (lonS, latS)
    
    Args:
        lon (numpy.ndarray): longitudes of the grid
        lat (numpy.ndarray): latitude of the grid
        lon_s (float): station longitude 
        lat_s (float): station latitude
        radius_earth (float): the radius of the Earth
    
    Returns:
        distanceS (np.ndarray): distance to station S matrix of dimensions dim(lon) x dim(lat) 
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
    """Generates a synthetic seismogram (Green's Functions) matrix in frequency domain
    using the given lon, lat, N, path_file_axisem, distance_s, and optional conjugate flag.
    Returns the synthetic seismogram matrix.

    Args:
        spectrum_axi (numpy.ndarray): synthetic seismogram in frequency domain
        lon (numpy.ndarray): longitude of the grid
        lat (numpy.ndarray): latitude of the grid
        N (int): number of samples
        distance_s (numpy.ndarray): distance matrix
        comp (str, optional): component
        conjugate (bool, optional): conjugate flag

    Returns:
        synth (numpy.ndarray): synthetic seismogram
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
    """Computes correlation function between stations A and B for a chunk of the www3 model source.   
    
    Args:
        lon_inf (float): Lower longitude bound.
        lon_sup (float): Upper longitude bound.
        lat_inf (float): Lower latitude bound.
        lat_sup (float): Upper latitude bound.
        N (int): Number of points.
        fe (float): Sampling frequency.
        date_vect (array): Vector of dates.
        spectrum_axi (np.ndarray): Axisem archive spectrum.
        file_model (str): WW3 PSD file.
        lon_staA (float): Longitude of station A.
        lat_staA (float): Latitude of station A.
        lon_staB (float): Longitude of station B.
        lat_staB (float): Latitude of station B.
        comp (str): Component.

    Returns:
        corr_f (np.ndarray): Computed correlation function as a single precision complex array.
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
    """Computes the cross-correlation function between two seismic stations.
    It takes the coordinates of the stations, the path to the AxiSEM archive, the path to the model, 
    a vector of dates, a normalization factor, and a component parameter. 
    It returns the computed cross-correlation function and the corresponding time array.

    Args:
        coords_staA (tuple): Coordinates of station A (longitude, latitude).
        coords_staB (tuple): Coordinates of station B (longitude, latitude).
        path_model (str): Path to the WW3 PSD model.
        date_vect (list): Vector of dates (year, month, day, hour).
        spectrum_axi (array): AxiSEM spectrum.
        fe (float, optional): Sampling frequency.
        N (int, optional): Number of points.
        extent (list, optional): Longitude and latitude slices.
        comp (str, optional): Component parameter.

    Returns:
        corr (xarray.DataArray): Synthetic correlation in time domain.
        time_corr (numpy.ndarray): Time array.
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