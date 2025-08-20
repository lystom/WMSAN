# About 

## How to Cite WMSAN ?

- [Tomasetto, L., Boué, P., Ardhuin, F., Stutzman, E., Xu, Z., De Plaen, R., & Stehly, L. (2024). WMSAN Python Package: From Oceanic Forcing to Synthetic Cross-correlations of Microseismic Noise. EathArXiv.](https://doi.org/10.26443/seismica.v4i1.1483)

## Funding
This research was funded by the French National Research Agency (ANR) under the project TERRACORR (ANR-20-CE49-0003)
![ANR logo](https://gricad-gitlab.univ-grenoble-alpes.fr/tomasetl/ww3-source-maps/-/raw/main/ANR-logo-2021-complet.jpg)

## References

- [Ardhuin, F., Stutzmann, E., Schimmel, M., & Mangeney, A. (2011). Ocean wave sources of seismic noise. Journal of Geophysical Research: Oceans, 116(C9).](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2011JC006952)

- [Stutzmann, E., Ardhuin, F., Schimmel, M., Mangeney, A., & Patau, G. (2012). Modelling long-term seismic noise in various environments. Geophysical Journal International, 191(2), 707-722.](https://academic.oup.com/gji/article/191/2/707/644255)

- [Ermert, L., Igel, J., Sager, K., Stutzmann, E., Nissen-Meyer, T., & Fichtner, A. (2020). Introducing noisi: a Python tool for ambient noise cross-correlation modeling and noise source inversion. Solid Earth, 11(4), 1597-1615.](https://se.copernicus.org/articles/11/1597/2020/)

- [Igel, J. K., Bowden, D. C., & Fichtner, A. (2023). SANS: Publicly available daily multi‐scale seismic ambient noise source maps. Journal of Geophysical Research: Solid Earth, 128(1), e2022JB025114.](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2022JB025114)

- [Longuet-Higgins, M. S. (1950). A theory of the origin of microseisms. Philosophical Transactions of the Royal Society of London. Series A, Mathematical and Physical Sciences, 243(857), 1-35.](https://royalsocietypublishing.org/doi/10.1098/rsta.1950.0012)

- [Gimbert, F., & Tsai, V. C. (2015). Predicting short-period, wind-wave-generated seismic noise in coastal regions. Earth and Planetary Science Letters, 426, 280-292.](https://www.sciencedirect.com/science/article/abs/pii/S0012821X15003738)

- [Gualtieri, L., Stutzmann, É., Farra, V., Capdeville, Y., Schimmel, M., Ardhuin, F., & Morelli, A. (2014). Modelling the ocean site effect on seismic noise body waves. Geophysical Journal International, 197(2), 1096-1106.](https://doi.org/10.1093/gji/ggu042)

- [Boué, P., & Tomasetto, L. (2024). Opportune detections of global P-wave propagation from microseisms interferometry. Comptes Rendus. Géoscience, 356(S4), 1-16.](https://comptes-rendus.academie-sciences.fr/geoscience/articles/10.5802/crgeos.222/)

- [The WAVEWATCH III® Development Group (WW3DG), 2019: User manual and system documentation of WAVEWATCH III® version 6.07. Tech. Note 333, NOAA/NWS/NCEP/MMAB, College Park, MD, USA, 326 pp. + Appendices.](https://www.weather.gov/sti/coastalact_ww3)

- [Zhang, R., Boué, P., Campillo, M., & Ma, J. (2023). Quantifying P-wave secondary microseisms events: a comparison of observed and modelled backprojection. Geophysical Journal International, 234(2), 933-947.](https://doi.org/10.1093/gji/ggad103)

- [Nissen-Meyer, T., van Driel, M., Stähler, S. C., Hosseini, K., Hempel, S., Auer, L., ... & Fournier, A. (2014). AxiSEM: broadband 3-D seismic wavefields in axisymmetric media. Solid Earth, 5(1), 425-445.](https://se.copernicus.org/articles/5/425/2014/)

- [Tomasetto, L., Boué, P., Stehly, L., Ardhuin, F., & Nataf, H. C. (2024). On the stability of mantle‐sensitive P‐wave interference during a secondary microseismic event. Geophysical Research Letters, 51(8), e2023GL108018.](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2023GL108018)

## Change Logs

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Calendar Versioning](https://calver.org/).

### [2025.0.0] - 2025-08-18

#### Added 
- Temporal variations notebook for body waves
- Funding reference in About section + ANR logo

#### Changed
- modified the code for Rayleigh waves amplification coefficients,
so that all parameters are in SI units (m/s instead of km/s).
- Revise Documentation to fix bugs and units in SI in all notebooks and descriptions.

#### Fixed Bugs
- radians of the model resolution in the surface element for synthetic spectrograms computation, fixed.


### [2024.1.4] - 2025-03-03
#### Added
- Condition on input latitudes : returns an error in absolute value higher than 90.

#### Changed
- Rayleigh waves amplification coefficient in F_prox (microseismic_sources.ipynb) modified. We use C instead of C^2 previously. This matches with the cP coefficients use.
- fix bugs when reading refined bathymetry. 


### [2024.1.3] - 2024-12-18
#### Added
- read from url function (to be implemented by default later)


#### Changed
- dependency to dask in install, README, pyproject.toml
- modified try except in subfunction_rayleigh_waves.spectrogram

### [2014.1.2] - 2024-10-07

#### Added

#### Changed
- fix expression formula for spectrogram in documentation. 
- add exponential decay to temporal_variations.


### [2014.1.1] - 2024-09-09

#### Added

#### Changed

- Increase patch in version 2024.1.1
- modify documentation to rename equivalent vertical force to "proxy for the source force amplitude"
- update dependencies (remove dask)
- Spectrogram example to PPTF

#### Deprecated


### [2014.1.0] - 2024-08-26

#### Added

- New version with autocorr
- Release to host on Zenodo

#### Changed

- Increase patch in version 2024.1.0


### [2014.0.6] - 2024-08-26

#### Added

- Add functions to compute single station auto-correlations
- Add possibility to work on Pacific ocean.


#### Changed

- Increase patch in version 2024.0.6 and comment


#### Removed

- Unnecessary imports removed

#### Fixed

- Modify the paths in the temporal variations notebook so that it can be run directly when cloning the repo.
 

### [2014.0.5] - 2024-08-09

#### Added

- Build documentation
- Add dependencies in README

#### Changed

- Update docs with typo removed
- Update minimum python version and OS in pyproject.toml
- Update for Numpy 2.0.1
- Update minimum python version and OS

### [2014.0.4] - 2024-08-09

#### Added

- Deploy documentation

#### Changed

- Add documentation link at the beginning of README.md

#### Fixed

- Correct typos for documentation

### [2014.0.3] - 2024-08-09

#### Added

- Build docstrings for documentation.
- Update Documentation with mkdocs

### [2014.0.2] - 2024-08-09

#### Added

- Initial commit
- README.md
- Add LICENSE.md

## License
MIT License

Copyright (c) [2024] [Lisa Tomasetto]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Contributors
[Lisa Tomasetto](https://github.com/lystom)

[Raphaël De Plaen](https://orcid.org/0000-0003-3477-2001)

[Reza D. D. Esfahani](https://github.com/resfahani)

[Laura Ermert](https://lermert.github.io/)

[Lei Li](https://orcid.org/0000-0002-9012-4853)