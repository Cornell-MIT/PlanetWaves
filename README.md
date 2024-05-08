# PLANETWAVES is a four dimensional spectral wave model adapted from the University of Miami wave model ([UMWM](https://github.com/umwm/umwm)) to produce surface waves for a given bathymetry and a given wind climate on different planets. The model has been validated using buoy data at the Great Lakes and previous wave tank experiments for different surface pressure, temperature, salinity.

```
umwm_titan
├── data : 
│   ├── Earth
│   │   ├── GreatLakes
│   │   │   ├── LakeErie
│   │   │   │   ├── 45132_Buoy : Wave height and weather data at deepwater buoy 45132
│   │   │   │   └── BathyData : Bathymetric sounding data of Lake Erie  
│   │   │   └── LakeSuperior
│   │   │       ├── 45004_Buoy : Wave height and weather data at deepwater buoy 45004
│   │   │       └── BathyData : Bathymetric sounding data of Lake Superior
│   │   └── OpenOcean
│   │       └── 41139_Buoy : Wave height and weather data of deepwater (open ocean) buoy 41139
│   ├── Mars
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
├── figures : figures for paper
└── planetwaves : main model scripts
    ├── plot_runs : scripts for plotting output of model
    └── unit_testing : scripts for unit testing model
```
## To do
```
(1) Add shallow buoy for Lake Superior
(2) Add shallow and deep buoy for Lake Erie(?)
(3) Finish comparison with Banfield+2015 for changing environmental conditions
(4) Run for Titan lakes
```
