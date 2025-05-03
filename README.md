# PLANETWAVES is a four dimensional spectral wave model adapted from the University of Miami wave model ([UMWM](https://github.com/umwm/umwm)) to produce surface waves for a given bathymetry and a given wind climate on different planets. The model has been validated using buoy data at the Great Lakes and previous wave tank experiments for different surface pressure, temperature, salinity.

```
umwm_titan
├── applications : model run scripts
│   └── WaveModelling : Scripts for modelling waves in lakes across different planets/exoplanets
|   └── ShorelineSmoothing : Scripts for modelling shoreline diffusivity of lakes
|   └── SedimentEntrainment : Scripts for modelling nearshore sediment entrainment by waves
|   └── past_runs : useful wave model runs that took a long time to finish and were saved for future analysis 
├── data : 
│   ├── Earth
│   │   ├── GreatLakes
│   │   │   └── LakeSuperior
│   │   │       ├── 45004_Buoy : Wave height and weather data at deepwater buoy 45004
│   │   │       └── BathyData : Bathymetric sounding data of Lake Superior
│   ├── Mars :
|   |     ├── Banfield+2015 wave tank experiments
|   |     └── M20_JezeroCrater_CTXDEM_20m : DTM mosaic of Jezero crater
|   |     └── StevensRubins2022 : Data from Stevens+2022 
│   └── Titan
│       ├── TAMnoTopo : Titan GCM with no topography influence
│       ├── TAMwTopo : Titan GCM with topography influence
|       |   ├── TAM_LM_winds.mat : 10 years of GCM winds at Ligeia Mare
|       |   └── TAM_OL_winds.mat : 10 years of GCM winds at Ontario Lacus 
│       └── TitanLakes 
│           └── Bathymetries
│               ├── SAR_bathy_cleaned : bathymetry estimated from SAR darkness
│               ├── bathtub_bathy : bathymetry assumed to to be due to constant slope and associated with largest fetch
│               └── shoreline : mapped shoreline coordinates from SAR mosaics
├── docs

├── figures : figures for papers
├── planetwaves : main model scripts
|   ├── intermediary_analysis : optional scripts run during the wave model
|   ├── post_analysis : optional scripts run after the wave model
|   ├── pre_analysis : optional scripts run before wave model 
|   └── unit_testing : scripts for unit testing model
└── validation : scripts comparing wave heights in model to buoy data on Earth and past empirical models
```
## To do
```
(1) Unit testing for UMWM-PlanetWaves comparison with Lake Superior 
(2) Spatially varying wind climate
(2) Temporally varying wind climate
(3) Radar scattering from MSS
(4) Parallelize code for higher grid resolution

(a) Check generic runs work (ignores unneccessary packages or adjusts functionality for matlab versions) for Taylor's comments
(b) update docs
(c) migrate to final repo location
```
