# API Overview

In this page we provide a detailed description of each function in the package.
It is ordered by scripts.

## Read hs and p2l Files
Fonctions to read WAVEWATCHIII NetCDF4 files. Using xarray or basic Python packages.

[Read WW3 Files](read_hs_p2l.md){: .btn}

## Body Waves Modeling (P and SV)
Functions to compute body waves equivalent vertical force from secondary microseismic ocean sources and synthetic spectrograms.

[Body Waves](body_waves.md){: .btn}

## Rayleigh Waves Modeling
Functions to compute Rayleigh waves equivalent vertical force from secondary microseismic ocean sources, synthetic cross-correlation functions, and synthetic spectrograms.

[Rayleigh Waves](read_hs_p2l.md){: .btn}

## Temporal Variations of Microseismic Sources
Function to compute a time serie which correspond to the mean equivalent vertical force value over a given source area.

[Temporal Variations](temporal_variations.md){: .btn}

## Synthetics Cross-Correlation Functions
Functions to compute synthetic cross-correlations with a WW3 source distribution and a Green's functions archive. 

[Synthetics](synthetics.md){: .btn}

## WMSAN to noisi
Functions to create a source distribution file to use as starting model in [noisi](https://github.com/lermert/noisi).

[noisi](wmsan_to_noisi.md){: .btn}