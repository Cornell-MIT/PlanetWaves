Scripts for sediment entrainment under tides and waves

MAIN SCRIPTS
(1) RUN_tidal_entrainment.m : predictions for normalized shear from tidal flows in straits and nearshore on Shields diagram
(2) RUN_wave_entrainment.m  : predictions for depth of entrainment under waves on Titan (bedload)
			- relies on wave model for deepwater wave properties from calc_Titan_OL_waves.m
			  where past runs of this script are stored in past_runs/Titan_DeepwaterWaves.mat
			- for some plots with Ontario bathymetry, relies on map_entrainment_depth_in_Ontario_lacus
ADDITIONAL SCRIPTS
(1) calc_max_settling_time.m :                for a given depth, this script will make predictions for how long it will take grains of different
				              densities in Titan's lakes of varying composition to settle out of suspension (based on Dietrich's
				              still-water estimations)
(2) calc_hydraulic_smoothness_sensitivity.m : calculates the percent difference on the velocity profile with depth for different 
				              assumptions of hydraulic roughness