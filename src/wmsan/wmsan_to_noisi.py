#!/usr/bin/env python3

# Preamble
__author__ = "Laura Ermert"
__credits__ = ["Laura Ermert"]
__version__ = "0.1"
__maintainer__ = "Laura Ermert, Lisa Tomasetto"
__email__ = "lisa.tomasetto@univ-grenoble-alpes.fr"
__status__ = "Development"
# Script by Laura Ermert
# mod. by Lisa Tomasetto 07/2024

##############################################################################
## Imports ##
import h5py
import numpy as np
import re
import os
import xarray as xr
from scipy.fftpack import next_fast_len
import matplotlib.pyplot as plt
from sklearn.neighbors import BallTree
from math import sin

import warnings ## to ignore warnings due to scikit learn see below
# https://stackoverflow.com/questions/40845304/runtimewarning-numpy-dtype-size-changed-may-indicate-binary-incompatibility
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

############################################################################
## Constants ##
R = 6371008.7714

############################################################################
## Functions ##

def create_sourcegrid_WW3(lat_min, lat_max, lon_min, lon_max, output_path):
    lon = np.arange(lon_min, lon_max, 0.5)
    lat = np.arange(lat_min, lat_max, 0.5)
    xx, yy = np.meshgrid(lon, lat)
    xx = xx.flatten()
    yy = yy.flatten()
    np.save(output_path + 'sourcegrid.npy', np.array([xx, yy]))
    
def get_nearest_neighbour_index(tree, query_point, n):
    distances, indices = tree.query([pt], k=n)
    return(indices)


def geo_nearest_neighbour(latlon_target, latlon_input, values, leaf_size=30):

    lat_in = latlon_input[0]
    lon_in = latlon_input[1]
    
    # transform into convenient array 
    # Create a meshgrid for x and y coordinates
    latlat, lonlon = np.meshgrid(lat_in, lon_in, indexing='ij')

    # Flatten the arrays
    lat_flat = latlat.flatten()
    lon_flat = lonlon.flatten()
    data_flat = values.flatten()
    
    lat = np.radians(lat_flat)
    lon = np.radians(lon_flat)
    grid_in = np.array([lat, lon]).T

    # construct a ball tree with haversine distance as the metric
    tree = BallTree(grid_in, leaf_size=leaf_size, metric='haversine')
    values_out = np.zeros(latlon_target[0].shape)

    for ix, query_point in enumerate(np.array(latlon_target).T):
        pt = np.radians(query_point)
        distance, index = tree.query([pt], k=1)
        value_nearest = data_flat[index[0][0]]
        values_out[ix] = value_nearest
    return(values_out)


def get_approx_surface_elements(lon):

    surfel = np.ones(len(lon)) #* 35_000 ** 2  # veeery rough
    # set to 1 to suppress multiplication by surface element ==> (LISA) done by muting  35_000 ** 2

    return(surfel)

def run(noise_source_file, grid_file, greens_function_file):
    fmin = 0.0    
    microseism_model = xr.load_dataset(noise_source_file)
    microseism_model = microseism_model.fillna(0.0)
    source_grid = np.load(grid_file)

    print(microseism_model)
    
    n_basis = len(microseism_model["frequency"])
    # here we use n_basis gaussians to make life simpler
    # frequency axis of the simulation
    with h5py.File(greens_function_file, "r") as wf:
        sampling_rate = wf["stats"].attrs["Fs"]
        n = next_fast_len(2 * wf["stats"].attrs['nt'] - 1)
        freqs_new = np.fft.rfftfreq(n, d=1. / sampling_rate)
    freqs_source = microseism_model.frequency.values

    # define arrays
    spect_basis =np.zeros((n_basis, len(freqs_new)), dtype=float)
    # # coefficient matrix
    n_loc = source_grid.shape[1]
    mod = np.zeros((n_loc, n_basis), dtype=float)

    # loop over (sparsely sampled) source model frequencies
    for ixf, ff in enumerate(freqs_source):
        print(f"Frequency {ixf+1} of {len(freqs_source)}")

        # define the basis function
        if ixf < n_basis - 1:
            gauss_sigma = (freqs_source[ixf + 1] - ff) / 1.5
        gaussian = 1./(gauss_sigma * np.sqrt(2. * np.pi)) * np.exp(0.5 * -np.power((freqs_new - ff) / gauss_sigma, 2.0))

        spect_basis[ixf, :] = gaussian
        
        # interpolate the map to the source grid
        temp = geo_nearest_neighbour([source_grid[1], source_grid[0]],
            [microseism_model.latitude.values, microseism_model.longitude.values],
            microseism_model.F_f.isel(frequency=ixf).values/1e6) ## Values to the square and multiply by 1e12 to normalise by the square of the GF force (not sure)  

        mod[:, ixf] = np.square(temp.copy())

    # Save to an hdf5 file
    noise_source_output_file = re.sub("nc", "h5", noise_source_file)
    with h5py.File(noise_source_output_file, 'w') as fh:
        fh.create_dataset('coordinates', data=np.load(grid_file))
        fh.create_dataset('frequencies', data=freqs_new)
        fh.create_dataset('surface_areas',
                          data=get_approx_surface_elements(source_grid[0]))
        fh.create_dataset("spectral_basis", data=spect_basis)
        fh.create_dataset("model", data=mod)

if __name__ == "__main__":
    # INPUT ######################################################################
    noise_source_file = "./rayleigh_waves/F/F_2008111506.nc"  # the source file from microseism_source.ipynb in N.s^{1/2}
    # To compute the source term we need the square of the values in N^2.s (already multiplied by the surface area in m^2)

    grid_file = "/Users/tomasetl/Documents/code/noisi/noisi/examples/wmsan/sourcegrid.npy"
    # the green's function file is needed to get the right frequency sampling
    # and the right number of time steps.
    greens_function_file = "/Users/tomasetl/Documents/code/noisi/noisi/examples/wmsan/greens/G.SSB..MXZ.h5"
    # END ##################################################################
    output_file = re.sub("nc", "h5", noise_source_file)
    if os.path.exists(output_file):
        raise ValueError("Output file must not exist yet: " + output_file)
    run(noise_source_file, grid_file, greens_function_file)