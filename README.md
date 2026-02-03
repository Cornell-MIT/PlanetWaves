# PLANETWAVES is a four-dimensional spectral wave model adapted from the University of Miami wave model ([UMWM](https://github.com/umwm/umwm)) to produce surface waves for a given bathymetry and a given wind climate on different planets. The model has been validated using buoy data at the Great Lakes and previous wave tank experiments for different surface pressure, temperature, and salinity.

```
umwm_titan
├── applications : model run scripts
│   └── WaveModelling : Scripts for modelling waves in lakes across different planets/exoplanets
|   └── DeltaMorphology : Scripts for modelling shoreline diffusivity of lakes on Titan
|   └── SedimentEntrainment : Scripts for modelling nearshore sediment entrainment by waves on Titan
├── data : 
│   ├── Earth
|   |   ├── analyze_buoy_data.m : script to collate buoy data with restrictions
|   |   ├── find_quiet_GREATLAKEs.m: function finds period of relatively uniform wind fields in buoy data
|   |   ├── compare_model.m : script to compare buoy data, PlanetWaves, UMWM, Pierson-Moskowitz, and JONSWAP estimates for wave heights
|   |   ├── find_fetch.py : script to find the fetch in all directions for a buoy given the lake bathymetry (produces WindFetchLS_45004.csv)
|   |   ├── WindFetchLS_45004.csv : direction wind originating from and resulting fetch for buoy 
│   │   ├── GreatLakes
│   │   │   └── LakeSuperior
│   │   │       ├── 45004_Buoy : Wave height and weather data at deepwater buoy 45004
│   │   │       └── BathyData : Bathymetric sounding data of Lake Superior
│   ├── Mars :
|   |     ├── Banfield2015_table1 : data from Banfield+2015 
|   |     └── M20_JezeroCrater_CTXDEM_20m : DTM mosaic of Jezero crater
|   |     └── StevensRubins2022 : Data from Stevens+2022 
│   └── Titan
│       ├── TAMnoTopo : Titan GCM with no topography influence (old model)
│       ├── TAMwTopo : Titan GCM with topography influence (latest model)
|       |   ├── TAM_LM_winds.mat : 10 years of GCM winds at Ligeia Mare
|       |   └── TAM_OL_winds.mat : 10 years of GCM winds at Ontario Lacus 
│       └── TitanLakes 
│           └── Bathymetries
│           |    ├── SAR_bathy_cleaned : bathymetry estimated from SAR darkness
│           |    ├── bathtub_bathy : bathymetry assumed to to be due to constant slope and associated with largest fetch
│           └── shoreline : mapped shoreline coordinates from SAR mosaics
├── docs : documentation, data resources, readmes, and python venv requirements
├── figures : figures from papers
|    ├── Schneck_2026 : figures from Schneck et al. 2026
|    └──  Detelich_2026 : figures from Detelich et al. 2026
├── planetwaves : main model scripts
|   ├── intermediary_analysis : optional scripts run during the wave model
|   ├── post_analysis : optional scripts run after the wave model
|   ├── pre_analysis : optional scripts run before wave model 
|   └── unit_testing : scripts for unit testing model
└── validation : scripts comparing wave heights in model to buoy data on Earth and past empirical models
```
