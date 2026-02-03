TO RUN:

1. Run UMWM with input values in UMWM_inputs to get wind and wave heights for Earth model
	>> source umwm_venv/bin/activate
	change wind speed in the input doc (umwm/namelist/main.nml)
	>> ./umwm
	>> cd output
	>> python3 plot_waves.py
	 {wind speed: sig wave height}
	>> add to excel sheet umwm_wind_waveheights.xlsx
2. Run Earth_LakeSuperior.m to create lakesuperior_run.mat for wind and wave heights for PlanetWaves
3. Fill in values for UMWM and PlanetWaves in umwm_wind_waveheights
4. Run compare_models.m for the final plot comparing empirical, data, and model results with error estimation


Note to self: Streamline this process. Super clunky right now. Maybe have Earth_LakeSuperior.m output into the excel doc? Put everything in a txt doc? The fortran model can be output to text.

========
To get fortran UMWM code running:
>> sudo apt install libnetcdff-dev
MakeFile
 NETCDFINC = -I/usr/include
NETCDFLIB = -L/usr/lib/x86_64-linux-gnu -lnetcdf -lnetcdff
Src/MakeFile
FC=gfortran 
CPPFLAGS=
>> make clean
>> make umwm
