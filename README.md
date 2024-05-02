PLANETWAVES is a four dimensional spectral wave model adapted from the University of Miami wave model (UMWM) to produce surface waves for a given bathymetry and a given wave climate on different planets. The model has been validated using buoy data at the Great Lakes.
=============================================================================================================================================================================================================================================================================

To run:
Create a copy of EXAMPLE_RUN.m and edit as needed. The result includes the significant wave height along the grid, the wave energy spectrum in space and frequency/direction space as well as the mean slope of the liquid surface along the grid.

=============================================================================================================================================================================================================================================================================
EXAMPLE_RUN.m creates model parameters as groups:
(1) Planet Condition : description of the surface, atmospheric, and liquid properties for the planet
(2) Model Geometry   : description of the model geometry including bathymetry, time evolution, tuning parameters, and buoy information
(3) Near surface wind conditions : description of the wind direction and magnitude. A wind direction of 0 degrees corresponds to wind traveling from left to right on the grid and increases clockwise such that a wind direction of 90 degrees 
corresponds to wind traveling from top of grid to bottom of grid.
(4) Unidirectional currents : description of eastward and northward currents within the liquid that are not being driven by the wind (e.g. tides)
(5) Housekeeping : option to show intermediary plots of model working (e.g. frequency dependence of different components of model at different depths) and to save data from each time step

=============================================================================================================================================================================================================================================================================
Full description of model parameters:

MAKEWAVES calculates E(x,y,k,theta) for wave field using an energy balance between wind input and multiple dissipation terms including turbulent dissipation (Sdt), bottom friction (Sbf), wave breaking (Sds), and spilling breakers (Ssb) as well as a non-linear
interaction term that shifts energy conservatively within the wave spectrum
The equation to solve is:
   E_{n+1} = E_{n} + del*(-Cg*cos(theta)*dE/dx - Cg*sin(theta)*dE/dy] + Sin + Snl - Sds


 Dimensions of (x,y,k,theta) are (m,n,o,p)


   Arguments:
       planet
           rho_liquid      : liquid density [kg/m3]
           nu_liquid       : liquid kinematic viscocity [m2/s]
           nua             : atmospheric gas viscocity [m2/s]
           gravity         : gravitational acceleration [m/s2]
           surface_temp    : surface temperature [K]
           surface_press   : surface pressure [Pa]
           surface_tension : surface tension of liquid [N/m]
           name            : planet name ['string'] e.g. 'Titan'
       model
           m               : number of grid cells in x-dimension
           n               : number of grid cells in y-dimension
           o               : number of frequency bins
           p               : number of angular (theta) bins. Must be factorable by 8 for octants to satisfy the Courant condition of numerical stability
           long            : longitude location of grid to plot 
           lat             : latitude location of grid to plot
           gridX           : size of grid cell in x-dimension [m]
           gridY           : size of grid cell in y-dimension [m]
           mindelt         : minimum time step to minimize oscillations of delt to very small values [s]
           maxdelt         : maximum time step to prevent large values as wind ramps up at startup [s]
           tolH            : minimum change in SigH to stop model at (otherwise runs to full Tsteps) 
           cutoff_freq     : cutoff frequency bin seperating the diagnostic from advecting wavenumbers
           min_freq        : minimum frequency to model [Hz]
           max_freq        : maximum frequency to model [Hz]
           bathy_map       : m x n array of depth [m] (+ values = subsurface, - values = subaerial)
           time_step       : maximum size of time step [s]
           num_time_steps  : length of model run in terms of # of time steps
           tune_A1         : tuning parameter for input term (Donelan et al. 2012, eqn. 12)
           tune_mss_fac    : A3 tuning parameter in spilling breaker term (Donelan et al. 2012, eqn. 16)
           tune_Sdt_fac    : fraction of Sdt term going into spectrum (A4 in Donelan et al. 2012, eqn. 20)
           tune_Sbf_fac    : fraction of Sbf term going into spectrum (Gf in Donelan et al. 2012, eqn. 22)
           tune_cotharg    : tuneable constant in argument in coth() term for dissipation term (is 1 but later improved fit using 0.2 in Donelan et al. 2012, eqn. 17)
           tune_n          : exponent n term in dissipation term (Donelan et al. 2012, eqn. 15)
           z_data          : elevation where wind measurements are taken (e.g. u10 would correspond to z_data = 10 m) [m]
       wind
           dir             : direction of incoming near-surface wind (CCW from East) [radians]
           speed           : magnitude of incoming near-surface wind [m/s]
       uniflow
           East            : Eastward unidirectional current [m/s]
           North           : Northward unidirectional current [m/s]
       Etc
           showplots       : 0 = no plots made, 1 = plots intermediary steps of model
           savedata        : 1  = save data for time steps (will slow down model run), 0 = skip saving


   Returns:
       sigH                : significant wave height at specified (Model.lat, Model.lon) coordinates [m]
       htgrid              : significant wave height for each grid cell [m]
       E_each              : wave energy spectrum (x,y) in space and in (frequency,direction) space
       ms                  : mean slope of liquid surface

=============================================================================================================================================================================================================================================================================
Explanation of warnings:

(1) 'Numerical ringing. Suggested fix: (1) Re-run with a larger cutoff frequency (2) Re-run with a larger maximum time step.'
The model parameters and run information is saved to a *_RunLog.txt file named as: DATE_wind_speed_RunLog.txt. Previous mat files will be deleted if the model is re-run but the log files will be retained for reference. The log file includes the 
fraction of the wave model completed for a particular wind speed "u = _WINDSPEED_ : fraction time completed: _FRACTION_" as well as how the amplitude of the significant wave height at specified location as changed with new step in time 
" t_n-1/t_n : _FRACTION_" where the fraction represents the amplitude in previous time step divided by the amplitude in the current time step. For stable growth, this value should be less than or equal to 1. The model will throw a warning 
if there is numerical ringing resulting in a fraction larger than 1, which is usually a problem for weaker wind speeds at the beginning of the time evolution. This can generally be resolved by increasing the cutoff frequency (Model.cutoff_freq) 
or increasing the size of the maximum time step (Model.time_step).
