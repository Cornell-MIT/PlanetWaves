Notes from meeting with Taylor (1/8/2023) : 

Was able to run the model using the smoothed LM bathymetry but it was really slow so he coarsened it by a factor of 3 in each direction (9x fewer gridpoints)
and used the MATLAB profiler found that line 488 and line 544 take up over half the runtime. Could be a good target for optimizing.

To run UMWM_Titan_2c.m and/or UMWM_Titan_FetchLaws.m requires a directory named 'Titan' for results to be saved and loaded from.

--------------------------------------------------------------------------------------------------------------------------------------------------------------------
Initial Notes from Alex (12/16/2022):


MAIN SCRIPTS : 
|
|_UMWM_Titan_Ligeia_Noplot_3m_v1.m : wave model using LM bathymetry map for uniform 2.5 and 3 m/s winds <-- most commented/recently updated
|
|_UMWM_Titan_FetchLaws.m : wave model for small patch of a deep sea for winds 0.4 - 10 m/s. This script was used to determined the threshold for wave generation and 
			       generate plots for significant wave height vs. wind speed


SUBROUTINES : (must be in the same directory as the model script
|
|_WAVEKgT.m :
|
|_smooth_nu.m :


MAIN OUTPUTS : 
|
|_1. Map of signifigant wave heights
|
|_2. Wave field energy spectrum 
|
|_3. Magnitudes of dissipation terms in the model