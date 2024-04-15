# Contents

This repository contains scripts to plot source maps of ambient noise in the secondary microseismic frequency band.
It uses the WAVEWATCHIII hindcast data to generate these maps for P and S waves.
Here is the detail of its content:
- microseismic_sources.ipynb : a Jupyter Notebook to generate and/or save source maps.
- subfunctions_source_DF.py: python background functions.
- read_hs_p2l.py : reads WW3 files.
- ww3.07121700.dpt: bathymetry file.
- readWW31.py : script reading the bathymetry file previously mentioned.
- cP.nc and cS.nc : amplification coefficients for P and S waves, from Gualtieri et al. (2014).

# Recommended Install

## for MacOSX
```wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ~/miniconda.sh```
## for Linux
```wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh```
```bash ~/miniconda.sh```
```source .bashrc```

## Install Environment 

```conda create --name notebook --file pkgs.txt```

or

```conda env create --file env_notebook.yml```

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
