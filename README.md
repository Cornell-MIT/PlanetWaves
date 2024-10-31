# PLANETWAVES is a four dimensional spectral wave model adapted from the University of Miami wave model ([UMWM](https://github.com/umwm/umwm)) to produce surface waves for a given bathymetry and a given wind climate on different planets. The model has been validated using buoy data at the Great Lakes and previous wave tank experiments for different surface pressure, temperature, salinity.

```
umwm_titan
├── data : 
│   ├── Earth
│   │   ├── GreatLakes
│   │   │   └── LakeSuperior
│   │   │       ├── 45004_Buoy : Wave height and weather data at deepwater buoy 45004
│   │   │       └── BathyData : Bathymetric sounding data of Lake Superior
│   ├── Mars : Banfield+2015 wave tank experiments
│   └── Titan
│       ├── TAMnoTopo : Titan GCM with no topography influence
│       ├── TAMwTopo : Titan GCM with topography influence
│       └── TitanLakes 
│           └── Bathymetries
│               ├── SAR_bathy_cleaned : bathymetry estimated from SAR darkness
│               ├── bathtub_bathy : bathymetry assumed to to be due to constant slope and associated with largest fetch
│               └── shoreline : mapped shoreline coordinates from SAR mosaics
├── docs
├── examples : model run scripts
│   └── model_logs : log of past runs
|   └── results : saved results from run
├── figures : figures for paper
└── planetwaves : main model scripts
    ├── intermediary_analysis : optional scripts run during the wave model
    ├── post_analysis : optional scripts run after the wave model
    ├── pre_analysis : optional scripts run before wave model 
    └── unit_testing : scripts for unit testing model
```
## To do
```
(1) Dynamic drag coefficient from wave maturity?
(2) Spatially varying wind climate
(3) Temporally varying wind climate

```
