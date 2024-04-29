# WWSAN Python Package

## Description
This package is built to help computation of seismic ambient noise source maps and other products based on the WAVEWATCHIII hindcast output.

## Contents
- src/ : contains all Python scripts and subfunctions.
- data/ contains additional files used in computation
- notebooks/ : contains Jupyter Notebooks with detailed examples on how to use this package. Rayleigh waves and body waves are separated.

## Packages Needed
- [numpy](https://numpy.org/doc/stable/)
- [matplotlib](https://matplotlib.org/stable/)
- [cartopy](https://scitools.org.uk/cartopy/docs/latest/index.html)
- [xarray](https://docs.xarray.dev/en/stable/)
- [netCDF4](https://unidata.github.io/netcdf4-python/)
- [obspy](https://docs.obspy.org/)
- [datetime](https://docs.python.org/3/library/datetime.html)
- [scipy](https://scipy.org/)
- [pandas](https://pandas.pydata.org/pandas-docs/version/2.1.4/index.html)

# Architecture of WWSAN Python Package

![Scheme showing the different codes and Notebooks present in this repository and how they connect.](./package_archi.png)
