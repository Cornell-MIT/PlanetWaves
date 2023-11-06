    #!/usr/bin/env python3

import time
import os
import shutil
import re
import cmath


import datetime
import numpy as np

class ModelPlanet:
    def __init__(self):
        self.rho_liquid = None
        self.nu_liquid = None
        self.nua = None
        self.gravity = None
        self.surface_temp = None
        self.surface_press = None
        self.surface_tension = None
        self.name = None

    def default(self, planet_name):
        if planet_name == "Titan":
            self.rho_liquid = 540
            self.nu_liquid = 3e-7
            self.nua = 0.0126/1e4
            self.gravity = 1.352
            self.surface_temp = 92 
            self.surface_press = 1.5*101300
            self.surface_tension = 0.018
            self.name = planet_name
        elif planet_name == "Earth":
            pass
        else:
            raise ValueError(f"Invalid planet_name: {planet_name}, must be 'Earth' or 'Titan'")
            
    @classmethod
    def set_rho_liquid(self, _rho_liquid):
        pass

class ModelGeometry:
    def __init__(self):
        self.m = 10
        self.n = 10
        self.o = 25
        self.p = 288
        self.long = 5
        self.lat = 5
        self.gridX = 1000
        self.gridY = 1000
        self.mindelt = 0.0001
        self.maxdelt = 2000.0
        self.time_step = 10
        self.num_time_steps = 10
        self.tolH = None
        self.bathy_map = 100 * np.ones((self.m, self.n))

class ModelWind:
    def __init__(self):
        self.dir = 0
        self.speed = np.arange(5, 8, 1)

class ModelUniFlow:
    def __init__(self):
        self.East = 0 
        self.North = 0

class ModelHousekeeping:
    def __init__(self):
        self.showplots = 0
        self.savedata = 0
        self.showlog = 1


        
def make_waves(Planet, Model, Wind, Uniflow, Etc):
    start_time = time.time()

    sigH = None  
    htgrid = None  
    E_each = None  
    

    if Model.p % 8 != 0:
        raise AssertionError('Model input parameter p must be factorable by 8')
        
    current_time = datetime.datetime.now()
    date_string = current_time.strftime("%m%d%Y_%H%M%S")
    dfile = f"{date_string}_FetchLaws.txt"

    with open(dfile, 'w') as file:
        file.write(f"Run started at {current_time}\n")

        ModelResults = os.path.join(os.getcwd(), 'ModelResults')

        if not os.path.exists(ModelResults) and Etc.savedata:
            print(f"Results being saved to {ModelResults}")
            file.write(f"Results being saved to {ModelResults}\n")
            os.makedirs(ModelResults)
        else:
            oldmatfiles = os.path.join(ModelResults, '*.mat')
            for file in os.listdir(ModelResults):
                if file.endswith('.mat'):
                    fullmat = os.path.join(ModelResults, file)
                    print(f"Deleting previous .mat files {fullmat}")
                    file.write(f"Deleting previous .mat files {fullmat}\n")
                    os.remove(fullmat)

        print('================================================================')
        makeWaves_date = os.path.getmtime('makeWaves.py')
        makeWaves_date_str = str(makeWaves_date)
        print(f"Directional Wave Spectrum -- last updated: {makeWaves_date_str}")
        file.write(f"Directional Wave Spectrum -- last updated: {makeWaves_date_str}\n")
        speed_str_with_commas = ', '.join([f"{x:.1f}" for x in Wind.speed])
        print(f"Wind Speed(s) to Run: {speed_str_with_commas} m/s")
        file.write(f"Wind Speed(s) to Run: {speed_str_with_commas} m/s\n")
        print('================================================================')


        # Global constants
        kappa = 0.4  # Von-Karman constant
        i = 1j  # imaginary number i
        kgmolwt = 0.028  # gram molecular weight [Kgm/mol]
        RRR = 8.314  # Universal gas constant [J/K/mol]

        # Frequency and direction bins
        dr = cmath.pi / 180  # conversion from degrees to radians
        dthd = 360 / Model.p  # calculation of dthd in Python

        dely = Model.gridY * np.ones((Model.m, Model.n, Model.o, Model.p))
        delx = Model.gridX * np.ones((Model.m, Model.n, Model.o, Model.p))

        # Wind:
        gust = 0  # Gust factor applied to wind at each time step
        zref = 20  # Height of reference wind speed, normally at 20m [m]
        z = 10  # Height of measured wind speed [m]
        wfac = 0.035 # winddrift fraction of Uz ffor U10m
        
        
        # Ideal gas law: PV = nRT
        # Densities:
        rhoa = Planet.surface_press * kgmolwt / (RRR * Planet.surface_temp)  # air density [kg/m3]
        rhorat = rhoa / Planet.rho_liquid  # air-water density ratio

        # Wavenumber limits:
        kutoff = 1000  # wavenumber cut-off due to Kelvin-Hemholtz instabilities (Donelan 2012, pg. 3)
        kcg = np.sqrt(Planet.gravity * Planet.rho_liquid / Planet.surface_tension)  # wavenumber of slowest waves defined by the capillary-gravity waves, from Airy dispersion: omega^2 =gktanh{k)

        # modified wavenumbers
        kcgn = 1.0 * kcg  # n power min is not shifted with respect to kcg
        kcga = 1.15 * kcg  # a min is shifted above kcg.

        # frequency limits:
        f1 = 0.05  # minimum frequency
        f2 = 35  # maximum frequency

        # create frequency limits for spectrum
        dlnf = (np.log(f2) - np.log(f1)) / (Model.o - 1)  # frequency step size for log normal distribution
        f = np.exp(np.log(f1) + np.arange(Model.o) * dlnf)  # frequencies for spectrum
        dom = 2 * np.pi * dlnf * f  # angular frequency (w = 2pi*f)
        freqs = f  # save a copy of frequencies

        # Frequency bins:
        oa = 7  # top bin advected (make wavelengths [vector wn] at oa 1/20th of the grid spacing)
        ol = np.arange(0, oa+1)  # bins for long frequencies
        osh = np.arange(oa + 1, Model.o)  # bins for short frequencies

        # Compute diffusion values in 2 freqs and 2 directions:
        #   bf1 + bf2 = 1, and bt1 + bt2 = 1.
        bf1 = np.exp(-15.73 * dlnf * dlnf)  # part of non-linear source term for downshifting and spilling breakers, eqn. 21 Donelan 2012
        bf2 = np.exp(-15.73 * 4 * dlnf * dlnf)  # part of non-linear source term for downshifting and spilling breakers, eqn. 21 Donelan 2012
        bf1a = bf1 / (bf1 + bf2)  # normalization (eqn. 21, Donelan 2012)
        bf2 = bf2 / (bf1 + bf2)  # normalization (eqn. 21, Donelan 2012)
        bf1 = bf1a


        # Compute Snl_fac
        A = np.exp(dlnf)
        B = A * A
        fac = bf1 * (1 - 1 / B) + bf2 * (1 - 1 / B ** 2)  # eqn. 21, Donelan 2012
        Snl_fac = ((1.0942 / A) ** 1.9757) / fac  # eqn. 21, Donelan 2012

        # waveangle [radians]
        waveang = ((np.arange(Model.p) - Model.p / 2 + 0.5) * dthd * dr)  # wave angle [radians]
        th = ((np.arange(Model.p) - Model.p / 2 + 0.5) * dthd * dr)  # angle of wave propagation (phi in Donelan 2012)  (0.5 added to avoid division by zero in rotation term)
        cth = np.cos(th)  # cosine of angle in wave propagation direction 
        sth = np.sin(th)  # sine of angle in wave propagation direction
        dth = dthd * dr  # 2pi/p (small angle for integration [radians])
        
        # compute cos^2 for calculation of mss vs direction.
        tth = np.arange(Model.p) * dthd * dr  # angular difference between short waves and longer waves
        cth2 = np.cos(tth) ** 2  # cosine square of angular difference between short waves and longer waves
        # indices for refraction rotation:
        cw = np.concatenate(([Model.p], np.arange(1, Model.p)))  # clockwise indices
        ccw = np.concatenate((np.arange(2, Model.p + 1), [1]))  # counterclockwise indices
        # upwave indices for advective term:
        xp = np.concatenate(([1], np.arange(1, Model.m - 1)))
        yp = np.concatenate(([1], np.arange(1, Model.n - 1)))
        xm = np.concatenate((np.arange(2, Model.m), [Model.m]))
        ym = np.concatenate((np.arange(2, Model.n), [Model.n]))

        # Index being used for advection:
        cp = np.where(cth > 0)[0]
        cm = np.where(cth < 0)[0]
        sp = np.where(sth > 0)[0]
        sm = np.where(sth < 0)[0]

        # reshape matrices
        cth = np.tile(cth, (Model.o, 1, 1, 1))
        cth = np.moveaxis(cth, 0, 2)
        sth = np.tile(sth, (Model.o, 1, 1, 1))
        sth = np.moveaxis(sth, 0, 2)
        cth2 = np.tile(cth2, (Model.o, 1, 1, 1))
        cth2 = np.moveaxis(cth2, 0, 2)

    def close_diary():
        with open(dfile, 'a') as file:
            file.write(f"Run ended at {datetime.datetime.now()}\n")

    close_diary()


# Measure the elapsed time
    elapsed_time = time.time() - start_time
    print(f'Time elapsed: {elapsed_time} seconds')

    return sigH, htgrid, E_each


def main():
    
    planet_to_run = ModelPlanet()
    planet_to_run.default("Titan")
    
    model_to_run = ModelGeometry()
    
    wind_to_run = ModelWind()
    
    uniflow_to_run = ModelUniFlow()
    
    etc_to_run = ModelHousekeeping()
    
    
    DsigH, Dhtgrid, DE_each = make_waves(planet_to_run, model_to_run, wind_to_run, uniflow_to_run, etc_to_run)

    

    return


if __name__ == '__main__':
    main()
