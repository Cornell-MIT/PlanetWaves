## PLANETWAVES
### Four-dimensional spectral wave model adapted from the University of Miami wave model ([UMWM](https://github.com/umwm/umwm)) to produce surface waves for a given bathymetry and a given wind climate on different planets. The model has been validated using buoy data at the Great Lakes and previous wave tank experiments for different surface pressure, temperature, and salinity.

## Papers
[2] Detelich, C. E., Schneck, U. G., Palermo, R. V., Steckloff, J., Ashton, A. A., Perron, J. T., “Seasonality of Wind-Driven Hydrocarbon Waves on Titan’s Polar Lakes”. Manuscript submitted.

[1] Schneck, U. G., Detelich, C. E., Curcic, M., Ashton, A. A., Hayes, A. H., Perron, J. T. “Modeling Wind-Driven Waves on Other Planets: Applications to Mars, Titan, and Exoplanets”. Manuscript under revision.


## File Tree
```
umwm_titan
├── applications : model run scripts
│   └── WaveModelling : Scripts for modelling waves in lakes across different planets/exoplanets
|   └── DeltaMorphology : Scripts for modelling shoreline diffusivity of lakes on Titan
|   └── SedimentEntrainment : Scripts for modelling nearshore sediment entrainment by waves on Titan
├── data : 
│   ├── Earth
|   |   ├── find_fetch.py : script to find the fetch in all directions for a buoy given the lake bathymetry (produces WindFetchLS_45004.csv)
|   |   ├── WindFetchLS_45004.csv : direction wind originating from and resulting fetch for buoy 
│   │   ├── GreatLakes
│   │   │   └── LakeSuperior
│   │   │       ├── 45004_Buoy : Wave height and weather data at deepwater buoy 45004
│   │   │       └── BathyData : Bathymetric sounding data of Lake Superior
│   ├── Mars :
|   |     ├── Banfield2015_table1 : data from Banfield et al. 2015, Icarus
|   |     └── M20_JezeroCrater_CTXDEM_20m : DTM mosaic of Jezero crater
|   |     └── StevensRubins2022 : Data from Rubins et al. 2022, JGR: Planets
│   └── Titan
│       ├── TAMnoTopo : Titan GCM with no topography influence (old model)
│       ├── TAMwTopo : Titan GCM with topography influence (latest model), includes 10 years from Lora et al. 2022, Icarus
|       |   ├── TAM_LM_winds.mat : 10 years of GCM winds at Ligeia Mare
|       |   └── TAM_OL_winds.mat : 10 years of GCM winds at Ontario Lacus 
│       └── TitanLakes 
│           └── Bathymetries
│           |    ├── SAR_bathy_cleaned : bathymetry estimated from SAR darkness
│           |    ├── bathtub_bathy : bathymetry assumed to to be due to constant slope and associated with largest fetch
│           └── shoreline : mapped shoreline coordinates from SAR mosaics
├── docs : documentation, data resources, readmes, and python venv requirements
├── figures : figures for understanding model
├── planetwaves : main model scripts
|   ├── intermediary_analysis : optional scripts run during the wave model
|   ├── post_analysis : optional scripts run after the wave model
|   ├── pre_analysis : optional scripts run before wave model 
|   └── unit_testing : scripts for unit testing model (WORK IN PROGRESS)
└── validation : scripts comparing wave heights in model to buoy data on Earth and past empirical models
```

## Funding
This material is based upon work supported by the National Science Foundation Graduate Fellowship under Grant No. 2141064, Future Investigators in NASA Earth and Space Sciences and Technology (FINESST) SMD's Graduate Student Research Fellowship under grant 171064, Naval Research grant N000142412598, and NASA Cassini Data Analysis Program under grants 80NSSC18K1057 and 80NSSC20K0484 
