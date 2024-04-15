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

# Architecture of WWSAN Python Package

![Scheme showing the different codes and Notebooks present in this repository and how they connect.](./package_archi.png)
