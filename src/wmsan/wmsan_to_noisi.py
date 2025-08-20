#!/usr/bin/env python3

# Preamble
#__author__ = "Laura Ermert"
#__credits__ = ["Laura Ermert"]
#__version__ = "2025.0.0"
#__maintainer__ = "Laura Ermert, Lisa Tomasetto"
#__email__ = "lisa.tomasetto@univ-grenoble-alpes.fr"
# Script by Laura Ermert
# mod. by Lisa Tomasetto 07/2024

""" Functions to convert WMSAN data to Noisi source distribution starting file format. 

It contains five functions:

- `create_sourcegrid_WW3(lat_min, lat_max, lon_min, lon_max, output_path)`: Create a grid of latitude and longitude values within the specified range and save it as a NumPy array.

- `get_nearest_neighbour_index(tree, query_point, n)`: Get the indices of the nearest neighbors to a given query point using a ball tree.

- `get_nearest_neighbour(latlon_target, latlon_input, values, leaf_size=30)`: Get the nearest neighbors to a given query point using a ball tree.

- `get_approx_surface_elements(lon)`: Returns an array of ones with the same length as the input array `lon`, representing an approximation of surface elements.

- `run(noise_source_file, grid_file, greens_function_file)`: Runs the computation to generate the secondary microseismic source distribution for noisi, given noise source file, grid file, and greens function file.

"""

##############################################################################
## Imports ##
import numpy as np
import xarray as xr
import h5py
import warnings
import re
import os

from scipy.fftpack import next_fast_len
from sklearn.neighbors import BallTree
 ## to ignore warnings due to scikit learn see below
# https://stackoverflow.com/questions/40845304/runtimewarning-numpy-dtype-size-changed-may-indicate-binary-incompatibility
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

############################################################################
## Constants ##
R = 6371008.7714

############################################################################
## Functions ##

def create_sourcegrid_WW3(lat_min, lat_max, lon_min, lon_max, output_path):
    """Create a grid of latitude and longitude values within the specified range and save it as a NumPy array.
    
    Args:
        lat_min (float): The minimum value of the latitude range.
        lat_max (float): The maximum value of the latitude range.
        lon_min (float): The minimum value of the longitude range.
        lon_max (float): The maximum value of the longitude range.
        output_path (str): The path to the directory where the grid will be saved.

    Returns:
        None (NoneType): The function saves the grid as a NumPy array in the specified path.
    """
    lon = np.arange(lon_min, lon_max, 0.5)
    lat = np.arange(lat_min, lat_max, 0.5)
    xx, yy = np.meshgrid(lon, lat)
    xx = xx.flatten()
    yy = yy.flatten()
    np.save(output_path + 'sourcegrid.npy', np.array([xx, yy]))
    
def get_nearest_neighbour_index(tree, query_point, n):
    """Get the indices of the nearest neighbors to a given query point using a ball tree.
    This function uses a ball tree to find the indices of the nearest neighbors to a given query point. It takes a ball tree object, a query point, and the number of nearest neighbors to find as input. It then queries the ball tree to find the distances and indices of the nearest neighbors to the query point. The function returns an array of indices of the nearest neighbors.
        
    Args:
        tree (BallTree): The ball tree object used for nearest neighbor search.
        query_point (array-like): The query point for which to find the nearest neighbors.
        n (int): The number of nearest neighbors to find.

    Returns:
        values_out (array-like): An array of indices of the nearest neighbors.
    """
    pt = np.radians(query_point)
    distances, indices = tree.query([pt], k=n)
    return(indices)


def geo_nearest_neighbour(latlon_target, latlon_input, values, leaf_size=30):
    """Find the nearest neighbors in a grid of latitude and longitude coordinates.
    This function takes in a target set of latitude and longitude coordinates, an input set of latitude and longitude coordinates, and a set of values associated with the input coordinates. It uses a BallTree to efficiently find the nearest neighbor in the input set for each point in the target set. The function returns an array of values corresponding to the nearest neighbors.
        
    Args:
        latlon_target (tuple of arrays): The target set of latitude and longitude coordinates.
        latlon_input (tuple of arrays): The input set of latitude and longitude coordinates.
        values (array-like): The set of values associated with the input coordinates.
        leaf_size (int, optional): The number of points to group together in the BallTree. Defaults to 30.

    Returns:
        values_out (array-like): An array of values corresponding to the nearest neighbors.
    """

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
    """Returns an array of ones with the same length as the input array `lon`, representing an approximation of surface elements.
    
    Args:
        lon (array-like): An array-like object of longitude coordinates.
        
    Returns:
        surfel (numpy.ndarray): A 1D array of ones with the same length as `lon`.
    """

    surfel = np.ones(len(lon))
    # set to 1 to suppress multiplication by surface element ==> (LISA) done by muting  35_000 ** 2

    return(surfel)

def run(noise_source_file, grid_file, greens_function_file):
    """Runs the computation to generate the WMSAN to NOISI secondary microseismic source distribution for the given noise source file, grid file, and greens function file.

    Args:
        noise_source_file (str): The path to the noise source file in N.s^{1/2}.
        grid_file (str): The path to the grid file.
        greens_function_file (str): The path to the greens function file.

    Returns:
        None (NoneType): This function does not return anything.
    """   
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