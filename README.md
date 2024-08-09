# WMSAN Python Package

## Description
This package is built to help computation of seismic ambient noise source maps and other products based on the WAVEWATCHIII hindcast output.

## Documentation
A detailed documentation is available [on this page](https://tomasetl.gricad-pages.univ-grenoble-alpes.fr/ww3-source-maps/). 

## Contents
```
ww3-source-maps/
|-- LICENSE
|-- pyproject.toml
|-- README.md
|-- mkdocs.yml
|-- docs/
|-- site/
|-- src/
│   └── wmsan/
│       ├── readWW31.py
│       ├── read_hs_p2l.py
│       ├── subfunctions_body_waves.py
│       ├── subfunctions_rayleigh_waves.py
│       └── synthetics.py
│       └── wmsan_to_noisi.py
│       └── temporal_variation.py
│       └── synthetic_CCF.ipynb
│
|-- notebooks/
|   └── body_waves/
│       ├── amplification_coeff.ipynb
│       └── microseismic_sources.ipynb 
│   └── rayleigh_waves/
│       ├── amplification_coeff.ipynb
│       ├── microseismic_sources.ipynb
│       ├── spectrograms.ipynb
│       ├── rayleigh_source.ipynb
│       └── synthetic_CCF.ipynb
│       └── wmsan_to_noisi.ipynb
│       └── temporal_variations.ipynb
|-- data/
│   ├── C.nc
│   ├── cP.nc
│   ├── cS.nc
│   ├── longuet_higgins.txt
│   ├── stations_pair.txt
│   └── ww3.07121700.dpt
```
- src/ : contains all Python scripts and subfunctions.
- notebooks/ : contains Jupyter Notebooks with detailed examples on how to use this package. Rayleigh waves and body waves are separated.
- data/: contains additional files used in computation.

## Installation

### PyPI

The package is available on [PyPI](https://pypi.org/).

- Create an environment and install

    - if you use [Conda](https://docs.anaconda.com/free/miniconda/#quick-command-line-install) environments:

            conda create --name wmsan 
            conda activate wmsan
            conda install pip
            python3 -m pip install wmsan

        to deactivate your environment:

            conda deactivate

    - otherwise

            python3 -m venv venv
            source venv/bin/activate
            python3 -m pip install wmsan
    
        to deactivate your environment:
    
            deactivate

### From Source
1. Clone the repository 
``` 
cd path_to_your_wmsan_directory/
git clone https://gricad-gitlab.univ-grenoble-alpes.fr/tomasetl/ww3-source-maps.git 
cd ww3-source-maps/
````
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

## Where should I start ?

![Table representing the differrent paths to Jupyter Notebooks examples and where to find what you wish to compute.](https://gricad-gitlab.univ-grenoble-alpes.fr/tomasetl/ww3-source-maps/-/raw/main/sumup.png) 

## Architecture of WMSAN Python Package

![Scheme showing the different codes and Notebooks present in this repository and how they connect.](https://gricad-gitlab.univ-grenoble-alpes.fr/tomasetl/ww3-source-maps/-/raw/main/package_archi.png)
