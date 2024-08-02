# Welcome to WSAN

## Description

WSAN for Wave Model Sources of Ambient Noise is a user-friendly Python package to help seismologists model their observations through maps of ambient noise sources from WAVEWATCHIII hindcast outputs, but also to compute spectrograms (e.g., Ardhuin et al., 2011; Stutzmann et al., 2012) and seismic noise correlations (Ermert et al., 2020). In particular, we provide Python scripts and Jupyter Notebooks to compute maps of secondary microseismic noise sources distribution, synthetic seismic cross-correlations, temporal variations of sources and synthetic seismic spectrograms.

## Project layout

    ww3-source-maps/
    |-- LICENSE
    |-- pyproject.toml
    |-- README.md
    |-- src/
    │   └── wmsan/
    │       ├── readWW31.py
    │       ├── read_hs_p2l.py
    │       ├── subfunctions_body_waves.py
    │       ├── subfunctions_rayleigh_waves.py
    │       └── synthetics.py
    │       └── wmsan_to_noisi.py
    │       └── temporal_variation.py
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

- src/ : contains all Python scripts and subfunctions.
- notebooks/ : contains Jupyter Notebooks with detailed examples on how to use this package. Rayleigh waves  and body waves are separated.
- data/ contains additional files used in computation.

[Getting Started](getting_started.md){: .btn}
[API Overview](api_overview/api_overview.md){: .btn}
[User Guide](user_guide.md){: .btn}
[Developer Guide](developer_guide.md){: .btn}
[About](about.md){: .btn}