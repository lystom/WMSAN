# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Calendar Versioning](https://calver.org/).

## [2025.0.0] - 2025-08-18

### Added 
- Temporal variations notebook for body waves
- Funding reference in About section + ANR logo

### Changed
- modified the code for Rayleigh waves amplification coefficients,
so that all parameters are in SI units (m/s instead of km/s).
- Revise Documentation to fix bugs and units in SI in all notebooks and descriptions.

### Fixed Bugs
- radians of the model resolution in the surface element for synthetic spectrograms computation, fixed.


## [2024.1.4] - 2025-03-03
### Added
- Condition on input latitudes : returns an error in absolute value higher than 90.

### Changed
- Rayleigh waves amplification coefficient in F_prox (microseismic_sources.ipynb) modified. We use C instead of C^2 previously. This matches with the cP coefficients use.
- fix bugs when reading refined bathymetry. 


## [2024.1.3] - 2024-12-18
### Added
- read from url function (to be implemented by default later)


### Changed
- dependency to dask in install, README, pyproject.toml
- modified try except in subfunction_rayleigh_waves.spectrogram

## [2014.1.2] - 2024-10-07

### Added

### Changed
- fix expression formula for spectrogram in documentation. 
- add exponential decay to temporal_variations.


## [2014.1.1] - 2024-09-09

### Added

### Changed

- Increase patch in version 2024.1.1
- modify documentation to rename equivalent vertical force to "proxy for the source force amplitude"
- update dependencies (remove dask)
- Spectrogram example to PPTF

### Deprecated


## [2014.1.0] - 2024-08-26

### Added

- New version with autocorr
- Release to host on Zenodo

### Changed

- Increase patch in version 2024.1.0


## [2014.0.6] - 2024-08-26

### Added

- Add functions to compute single station auto-correlations
- Add possibility to work on Pacific ocean.


### Changed

- Increase patch in version 2024.0.6 and comment

### Deprecated

- 

### Removed

- Unnecessary imports removed

### Fixed

- Modify the paths in the temporal variations notebook so that it can be run directly when cloning the repo.


### Security

- 

## [2014.0.5] - 2024-08-09

### Added

- Build documentation
- Add dependencies in README

### Changed

- Update docs with typo removed
- Update minimum python version and OS in pyproject.toml
- Update for Numpy 2.0.1
- Update minimum python version and OS

## [2014.0.4] - 2024-08-09

### Added

- Deploy documentation

### Changed

- Add documentation link at the beginning of README.md


### Fixed

- Correct typos for documentation

## [2014.0.3] - 2024-08-09

### Added

- Build docstrings for documentation.
- Update Documentation with mkdocs

## [2014.0.2] - 2024-08-09

### Added

- Initial commit
- README.md
- Add LICENSE.md
