# Contents
This repository provides scripts to compute Rayleigh source maps of ambient noise in the secondary microseismic band.
It uses the WAVEWATCHIII hindcast data to model source areas and generates spectrograms of the vertical displacement for a given station location.

Here is the detail of its content:
- this README.md 
- env_notebook.yml : an environnment YAML file to install python dependencies.
- pkgs.txt : an ASCII file to install python dependencies.
- rayleigh_source.ipynb : Jupyter Notebook to compute and/or save Rayleigh source maps.
- spectrogram.ipynb :Jupyter Notebook to plot spectrogram of the vertical displacement in the secondary microseismic band.
- amplification_coeffs.ipynb : Jupyter Notebook computing Longuet-Higgins amplification coefficients maps.
- source_microseism.py : background functions for source maps.
- longuet_higgins.txt : values of the amplification coefficients for Longuet-Higgins (1950).
- ww3.07121700.dpt : bathymetry file 
- readWW31.py : pythin script to read bathymetry file above
- read_hs_p2l : reads WW3 hindcast data

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

```conda env creat --file env_notebook.yml```

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
