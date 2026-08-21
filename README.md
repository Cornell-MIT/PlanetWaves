# *PlanetWaves*
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18906357.svg)](https://doi.org/10.5281/zenodo.18906357)
### Four-dimensional spectral wave model adapted from the University of Miami wave model ([UMWM](https://github.com/umwm/umwm)) to produce surface waves for a given bathymetry and a given wind climate on different planets. The model has been validated using buoy data at the Great Lakes and previous wave tank experiments for different surface pressure, temperature, and salinity.

## Getting started
For an example of the model output, run Titan_Ontario_Lacus.m in the applications/WaveModelling folder. This will produce plots of waves and wave spectrum in Ontario Lacus (Titan's largest south polar lake) for winds between 0 and 4 m/s.

## File Tree
```
umwm_titan
├── applications
│   └── WaveModelling : modelling waves in lakes across different planets/exoplanets
|   └── DeltaMorphology : modelling shoreline diffusivity of lakes on Titan
|   └── SedimentEntrainment : modelling nearshore sediment entrainment by waves on Titan
├── data  
│   ├── Earth
|   |   ├── find_fetch.py : find the fetch in all directions for a buoy given the lake bathymetry (produces WindFetchLS_45004.csv)
|   |   ├── WindFetchLS_45004.csv : direction wind originating from and resulting fetch for buoy 
│   │   ├── GreatLakes
│   │   │   └── LakeSuperior
│   │   │       ├── 45004_Buoy : Wave height and weather data at deepwater buoy 45004
│   │   │       └── BathyData : Bathymetric sounding data of Lake Superior
│   ├── Mars 
|   |     ├── Banfield2015_table1 : data from Banfield et al. 2015, Icarus (for model comparison)
|   |     └── M20_JezeroCrater_CTXDEM_20m : DTM mosaic of Jezero crater (for model comparison)
|   |     └── StevensRubins2022 : Data from Rubins et al. 2022, Icarus (for model (comparison)
│   └── Titan
│       ├── TAMnoTopo : Titan GCM (TAM) with no topography influence (old model) from Lora et al. 2015, Icarus
│       ├── TAMwTopo : Titan GCM (TAM) with topography influence (latest model) from Lora et al. 2022, Icarus
|       |   ├── TAM_LM_winds.mat : 10 years of GCM winds at Ligeia Mare
|       |   └── TAM_OL_winds.mat : 10 years of GCM winds at Ontario Lacus 
│       └── TitanLakes 
│           └── Bathymetries
│           |    ├── SAR_bathy_cleaned : bathymetry estimated from SAR darkness
│           |    ├── bathtub_bathy : bathymetry assumed to be due to constant slope and associated with largest fetch
│           └── shoreline : mapped shoreline coordinates from SAR mosaics
├── docs : documentation, data resources, readmes, and python venv requirements
├── figures : figures for understanding model
├── planetwaves : main model scripts
|   ├── intermediary_analysis : optional scripts run during the wave model
|   ├── post_analysis : optional scripts run after the wave model
|   ├── pre_analysis : optional scripts run before wave model 
|   └── unit_testing : scripts for unit testing model (IN PROGRESS)
└── validation : scripts comparing wave heights in model to buoy data on Earth and past empirical models
```


```
CHANGE LOG

7/29/2026: (U. Schneck) A limiter min_height was added to restrict the minimum height of retrieval for
the wind profile. For values very close to the boundary layer, the velocity profile steepens due
 to the log relationship, which can artificially boost/limit the forcing from the wind.
This created non-monotonic behavior for the high wavenumber section of the spectrum, which affected
the growth of MSS with increasing wind speed. 

```

## Citations

2\. [Detelich et al. (2026) Modeling the Seasonality of Wind-Driven Hydrocarbon Waves in Titan's Polar Lakes, Journal of Geophysical Research: Planets](10.1029/2026JE009693)

1\. [Schneck et al. (2026) Modeling Wind‐Driven Waves on Other Planets: Applications to Mars, Titan, and Exoplanets, Journal of Geophysical Research: Planets](10.1029/2025JE009490)





## Funding
This material is based upon work supported by the National Science Foundation Graduate Research Fellowship under grant 214106, the Future Investigators in NASA Earth and Space Sciences and Technology (FINESST) SMD's Graduate Student Research Fellowship under grant 171064, the Office of Naval Research grant N000142412598, and the NASA Cassini Data Analysis Program under grants 80NSSC18K1057 and 80NSSC20K0484.
