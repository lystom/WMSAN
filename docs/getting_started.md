# Getting Started

## Package Overview
The following scheme describes the different products one can compute using this package and how they interact.

![Scheme showing the different codes and Notebooks present in this repository and how they connect.](img/package_archi.png)

First, let's install the package!

## Installation

### PyPI
The package is available on [PyPI](https://pypi.org/project/wmsan).

#### Create an environment and install

- if you use [Conda](https://docs.anaconda.com/free/miniconda/#quick-command-line-install) environments:
    ```
    conda create --name wmsan 
    conda activate wmsan
    conda install pip
    python3 -m pip install wmsan
    ```
    to deactivate your environment:

    ```
    conda deactivate
    ```
- otherwise
    ```
    python3 -m venv venv
    source venv/bin/activate
    python3 -m pip install wmsan
    ```
    to deactivate your environment:
    ```
    deactivate
    ```

### From Source

1. Clone the repository 

    ``` 
    cd path_to_your_wmsan_directory/
    git clone https://gricad-gitlab.univ-grenoble-alpes.fr/tomasetl/ww3-source-maps.git 
    cd ww3-source-maps/
    ```

2. Create an environment and install 

- if you use [Conda](https://docs.anaconda.com/free/miniconda/#quick-command-line-install) environments:

    ```
    conda create --name wmsan 
    conda activate wmsan
    conda install pip
    pip install .
    ```
    to deactivate your environment:
    ```
    conda deactivate
    ```

- otherwise

    ```
    python3 -m venv venv
    source venv/bin/activate
    python3 -m pip install .
    ```
    to deactivate your environment:
    ```
    deactivate
    ```

### Dependencies
- [numpy](https://numpy.org/doc/stable/)
- [matplotlib](https://matplotlib.org/stable/)
- [cartopy](https://scitools.org.uk/cartopy/docs/latest/index.html)
- [xarray](https://docs.xarray.dev/en/stable/)
- [netCDF4](https://unidata.github.io/netcdf4-python/)
- [obspy](https://docs.obspy.org/)
- [datetime](https://docs.python.org/3/library/datetime.html)
- [scipy](https://scipy.org/)
- [pandas](https://pandas.pydata.org/pandas-docs/version/2.1.4/index.html)
- [dask](https://www.dask.org/)
- [ipykernel](https://pypi.org/project/ipykernel/)
- [pyproj](https://pyproj4.github.io/pyproj/stable/)
- [h5py](https://docs.h5py.org/en/stable/)
- [tqdm](https://tqdm.github.io/)
- [notebook](https://jupyter-notebook.readthedocs.io/en/stable/)

## Where Should I Start ?
Congratulations ! You successfully installed WMSAN.
Now feel free to browse the [Tutorial](user_guide.md) for examples. The Python functions' full descriptions can be found in the [API Overview](api_overview/api_overview.md). 