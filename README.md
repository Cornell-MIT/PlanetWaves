# umwm_titan

To erase old logs (Ubuntu 22.04):
> python -c 'import makeWaves; makeWaves.remove_old_logs()'

## MATLAB
Model in makeWaves.m

Model specifics in test_runs.m

Makes plots of sig wave height from saved prev runs saved in '\Titan' in plot_sigH

EX.) H = plot_sigH('cutoff_5',{'New_1','New_2','New_3'}) where cutoff_5 is folder in current directory

## Python (in progress)
To work in virtual enviroment (miniconda)
> conda env create -f environment.yml

> conda activate wave_env

> conda env update --prefix ./env --file environment.yml  --prune

## Model Run Benchmarking

1. **(_US: ongoing_) Compare effect of changing cut off frequency** 
2. Compare to wind tank experiments (Banfield et al. 2015)
3. Comparing to closed lakes
   - ID periods of low variation in wind speed, direction, and gusts for (<ins>find_quiet_GREATLAKES.m</ins>)
     - Lake Erie
       1. Deep water buoy
       2. Shallow water buoy
     - **(_US: ongoing_)  Lake Superior** 
       1. **(_US: ongoing_) Deep water buoy**
       2. Shallow water buoy
   - Compare to model for
     - Deep uniform bathymetry
     - Realistic lake bathymetry
4. Comparing to open ocean
   - ID periods of low variation in wind speed, direction, and gusts near Hawai’i (Huppert et al. 2020) or Gulf of Mexico? (Lin et al. 2011)
     - Compare to model for deep uniform bathymetry



![Picture1](https://github.com/Cornell-MIT/umwm_titan/assets/24469269/d3ab52df-0260-4a08-a9b3-866d85b00e2b)

  ![ls_fetch](https://github.com/Cornell-MIT/umwm_titan/assets/24469269/f64b8d6c-16b9-465f-981e-3c8428b7cbea)








## Changes to model

1. Add Alves-Banner-Young curve to model
2. Add maturity threshold with minimum number of time steps
3. Set batch jobs to make model run in background with updates
