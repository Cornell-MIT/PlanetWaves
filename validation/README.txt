TO RUN:

1. Run UMWM with input values in UMWM_inputs to get wind and wave heights for Earth model
2. Run Earth_LakeSupeerior.m to create lakesuperior_run.mat for wind and wave heights for PlanetWaves
3. Fill in values for UMWM and PlanetWaves in umwm_wind_waveheights
4. Run compare_models.m for the final plot comparing empirical, data, and model results with error estimation


Note to self: Streamline this process. Super clunky right now. Maybe have Earth_LakeSuperior.m output into the excel doc? Put everything in a txt doc? The fortran model can be output to text.