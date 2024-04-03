#!/usr/bin/env python3

__author__ = "Lisa Tomasetto"
__copyright__ = "Copyright 2024, UGA"
__credits__ = ["Lisa Tomasetto"]
__version__ = "1.0"
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
    # apply a delay d to the input signal s
    d_phase = 2*np.pi*d*np.arange(len(s))/(len(s)) #phase associated to d
    d_phase[int(np.round(len(d_phase)/2))::]-=2*np.pi*d #for the hermitian symmetry of the fft
    return np.real(np.fft.ifft(np.fft.fft(s)*np.exp(-1j*d_phase)))


## Create Date Vect
def create_date_vect(dates):
    """Create date vector"""
    date_vect = np.zeros((len(dates), 4))
    for i in range(len(dates)):
        date = dates[i]
        year, month, day, hour = date.year, date.month, date.day, date.hour
        date_vect[i] = [year, month, day, hour]
    return date_vect

### Create GF Archive via syngine
def create_station_file():
    """Create station file"""
    distance = np.arange(0, 181, 1)
    n = len(distance)
    datacenter = "CUSTOM"
    net = "AA"
    sta = ["STA%03d"%i for i in range(n)]
    lat = [round(90-distance[i], 1) for i in range(n)]
    lon = np.zeros(n)
    loc = "--"
    elev = np.zeros(n)
    depth = np.zeros(n)
    
    # Write station file .txt
    f = open('station_file.txt', 'w')
    for i in range(n):
        f.write("CUSTOM   AA   %s   --   %s   %s   %s   %s\n"%(sta[i], lat[i], lon[i], elev[i], depth[i]))
    f.close()
    
def create_syngine_archive(station_file_path ='station_file.txt', sourcelongitude=0, sourcelatitude=90, dt=0.25, comp='Z'):
    """Create syngine archive"""
    
    # check if file exists
    distance = np.arange(0, 181, 1)
    # Open Station Info
    # Receivers
    network = []
    station = []
    lat = []
    lon = []
    elev = []
    depth = []

    f = open(station_file_path)
    lines = f.readlines()
    for l in lines:
        l = l.split('   ')
        network.append(l[1])
        station.append(l[2])
        lat.append(float(l[4]))
        lon.append(float(l[5]))
        elev.append(float(l[6]))
        depth.append(float(l[7][:-2]))
    f.close()
    model = "iasp91_2s"

    N = len(station)
    M = np.zeros((N, 14762)) 
    
    # Get synthetics
    for i in tqdm(range(len(station))):
        sta, rcv_lat, rcv_lon, net = station[i], lat[i], lon[i], network[i]
        url = "http://service.iris.edu/irisws/syngine/1/query?"
        url += "model=%s&"%(model)
        url += "sourcelatitude=%d&sourcelongitude=%d&"%(sourcelatitude, sourcelongitude)
        url += "sourcedepthinmeters=0&sourceforce=-1e10,0,0&"
        url += "receiverlatitude=%d&receiverlongitude=%d&components=%s&"%(rcv_lat, rcv_lon, comp)
        url += "units=displacement&dt=%f&format=miniseed"%(dt)

        try:
            st_synth = obspy.read(url)
            st_synth.detrend("demean")
            st_filt = st_synth.copy()
            st_filt.filter("lowpass", freq=0.5, corners=3)
            tr = st_filt[0]
            dt = tr.stats.delta
            start = tr.stats.starttime
            end = tr.stats.endtime
            time = np.arange(0, tr.stats.npts)*dt
            M[i, :] = tr
            
        except:
            print('no data for %s'%station[i])
            raise  
    ## Save M as netcdf
    ds = xr.DataArray(M, coords={'distance': distance, 'time': time}, dims=['distance', 'time'])
    ds.to_netcdf('synthetics_iasp91_2s_Z.nc')
    
def create_spectrum_archive(archive_file = 'synthetics_iasp91_2s_Z.nc'):
    ## Open synthetics
    SYNTH = xr.open_dataarray(archive_file)
    distance = SYNTH.distance
    time = SYNTH.time.values
    dt = time[1]-time[0]
    fe = 1/dt
    ## thiner distance step
    new_distance = np.arange(0, 180.1, 0.1)
    model_obspy = TauPyModel("iasp91")
    
    WAVEFORMS = np.zeros((len(new_distance), len(time))).astype(float)
    SPECTRUM = np.zeros((len(new_distance), len(time))).astype(complex)

    for i, dist in tqdm(enumerate(new_distance)):
        st = SYNTH.sel(distance=dist, method='nearest')
        dist_original = st.distance 
        arrival_original = model_obspy.get_travel_times(source_depth_in_km=0, distance_in_degree=dist_original, phase_list=["PP"])[0].time    
        arrival_new = model_obspy.get_travel_times(source_depth_in_km=0, distance_in_degree=dist, phase_list=["PP"])[0].time
        delay = arrival_new - arrival_original
        new_st = st.values.copy()
        new_st = apply_delay_f_domain(new_st, delay/dt)
        
        ## Plot synthetic
        # plt.figure()
        # plt.title('%.1f degree'%dist)
        # plt.plot(time, st, 'b', label='synthetic', alpha=0.5)
        # plt.plot(time, new_st, 'r', label='delayed', alpha=0.5)
        # plt.axvline(arrival_original, color='b')
        # plt.axvline(arrival_new, color='r')
        # plt.xlim(min(arrival_new, arrival_original)-500, max(arrival_new, arrival_original)+500)
        # plt.savefig("waveform_%fdegree.png"%dist, dpi=300)
        WAVEFORMS[i, :] = new_st
        
        ## FFT
        fft_synth = fft(new_st)
        freq_vect = np.fft.fftfreq(len(time), dt)
        # # # Plot
        # plt.subplot(211)
        # plt.plot(freq_vect, np.abs(fft_synth))
        # plt.subplot(212)
        # plt.plot(freq_vect, np.angle(fft_synth))
        # plt.savefig("spectrum_%fdegree.png"%dist, dpi=300)
        
        plt.close('all')
        SPECTRUM[i, :] = fft_synth
    ## Save as h5 file
    h5_name = 'synthetics_iasp91_2s_Z_waveforms.h5'
    h5_file = h5py.File(h5_name, 'w')
    h5_file.create_dataset('WAVEFORMS', data=WAVEFORMS)
    h5_file.create_dataset('distance', data=new_distance)
    h5_file.create_dataset('time', data=time)
    h5_file.close()
    h5_name = 'synthetics_iasp91_2s_Z_spectrum.h5'
    h5_file = h5py.File(h5_name, 'w')
    h5_file.create_dataset('SPECTRUM', data=SPECTRUM)
    h5_file.create_dataset('distance', data=new_distance)
    h5_file.create_dataset('frequency', data=freq_vect)
    h5_file.close()

def open_archive(h5_name_spectrum = 'synthetics_iasp91_2s_Z_spectrum.h5', h5_name_waveforms = 'synthetics_iasp91_2s_Z_waveforms.h5'):
    h5_file = h5py.File(h5_name_waveforms, 'r')
    W = h5_file['WAVEFORMS'][:]
    time = h5_file['time'][:]
    distance = h5_file['distance'][:]
    h5_file.close()
    h5_file = h5py.File(h5_name_spectrum, 'r')
    S = h5_file['SPECTRUM'][:]
    distance = h5_file['distance'][:]
    freq = h5_file['frequency'][:]
    h5_file.close()
    fe = 2*np.max(freq)
    return fe, freq, time, distance, S

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
    ds = xr.open_dataset(path_file_WW3) #, chunks={'time': 1, 'lon': 30, 'lat': 40}, parallel=True,combine='by_coords')  # mfdataset  
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
            S_synth[index, :] = np.squeeze(trace)
        else:
            S_synth[index, :] = np.squeeze(np.conj(trace))
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
            print("PSD shape",psd_model.data.shape)
            lon = psd_model.longitude
            lat = psd_model.latitude
            ## Distances
            distance_staA = distance_to_station(lon, lat, lon_s = lon_staA, lat_s = lat_staA)
            distance_staB = distance_to_station(lon, lat, lon_s = lon_staB, lat_s = lat_staB)

            ## Green's Functions spectrum
            psd_model.values *= matrix_GF(spectrum_axi, lon=lon, lat=lat, N=N, distance_s=distance_staA, conjugate=False, comp=comp).astype(complex)
            
            ## Green's Functions spectrum
            psd_model.values *= matrix_GF(spectrum_axi, lon=lon, lat=lat, N=N, distance_s=distance_staB, conjugate=True, comp=comp).astype(complex)
            
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
    corr *= 1e-20
    time_corr = np.arange(-N//2, N//2)*1/fe

    corr = xr.DataArray(corr, dims=['time'], coords={'time': time_corr}, name='synthetic correlation')
    del corr_f    
    return corr, time_corr