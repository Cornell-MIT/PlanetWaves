import numpy as np
import math
import matplotlib
matplotlib.use('Agg') # can be commented out if not working from command line
from matplotlib import pyplot as plt
import sys

'''
    This script takes in a txt file of grain sizes, flow
    depths, and flow velocities and generates the Shields diagram for the 
    specified planet
    
    data should be ordered as: velocity | D50 | depth
    
    If you want to specify your own array of data within the script,
    the data is fed into the main() function and can be hard-coded
    from there.
    
    example of call:
    >> python shields_curve.py Titan mydata.txt
    
    author   :  Una Schneck (ugschneck@gmail.com)
    Created  :  7/30/2026
    
'''

# UNIVERSAL CONSTANTS
kappa  = 0.4;                                              # Von Karman coefficient

def MAKE_SHIELDS_DIAGRAM(planet,D50,u, depth):
    '''
        FUNCTION WILL DRAW POINTS ONTO A SHIELDS
        DIAGRAM WITH THE UPPER CURVE BEING THE
        THRESHOLD FOR SUSPENSION AND THE LOWER
        CURVE BEING THE THRESHOLD FOR ENTRAINMENT
        AND POINTS BEING THE SHIELDS VALUE FOR 
        THE GIVEN DATA
            INPUTS:
                planets : planetary conditions, including grain denisty, liquid density and viscosity, and gravity
                D50     : median grain size (m)
                u       : flow velocity (m/s)
                depth   : flow depth (m)
            OUTPUTS:
                Saved PDF of Shields diagram with points colored by transport regime

    '''
    # (1) Calculate Shields number for each grain
    my_Rep             = calculate_particle_Re(planet.rho_s, planet.rho, planet.nu, planet.g, D50);
    my_theta           = calc_shields_number(planet.rho_s, planet.rho, planet.g, depth, u,D50);
    # (2) Calculate the Shields curve for critical entrainment
    full_range_d50     = np.arange(1e-6, 0.3, 1e-6)
    Rep_all            = calculate_particle_Re(planet.rho_s, planet.rho, planet.nu, planet.g, full_range_d50);
    crit_theta         = calc_critical_shields(Rep_all)
    # (3) Calculate the curve for critical suspension
    crit_sus           = calc_suspension_threshold(planet.rho_s, planet.rho, planet.nu, planet.g, full_range_d50);

    # interpolate data to classify its transport regime
    crit_theta_data   = np.interp(D50, full_range_d50, crit_theta)
    crit_sus_data     = np.interp(D50, full_range_d50, crit_sus)

    # classify transport regime
    no_entrainment    = my_theta <= crit_theta_data
    bedload           = (my_theta > crit_theta_data) & (my_theta <= crit_sus_data)
    suspended         = my_theta > crit_sus_data

    # don't want the suspension threshold to be below the entrainment threshold
    # so cut it off when they intersect
    mask_sus         = crit_sus > crit_theta;
    Rep_sus          = Rep_all[mask_sus];
    crit_sus         = crit_sus[mask_sus];
    
    # make a plot
    f = plt.figure(figsize=(8, 8))
    plt.plot(Rep_all**2, crit_theta, linewidth = 2, color='black')
    plt.plot(Rep_sus**2, crit_sus, linewidth = 2, color='black')
    # plot data and color by transport regime
    plt.scatter(my_Rep[no_entrainment]**2,my_theta[no_entrainment]**2,color='grey', marker='o', s=100, label='no motion')
    plt.scatter(my_Rep[bedload]**2,my_theta[bedload],color='red', marker='o', s=100, label='bedload')
    plt.scatter(my_Rep[suspended]**2,my_theta[suspended],color='blue',marker='o',s=100, label='suspended load')
        
    # make plot look nice
    plt.title(f'Shields Diagram for {planet.name}')
    plt.xlabel(r'Dimensionless grain size $D_*$')
    plt.ylabel(r'Dimensionless shear stress $\tau$')
    plt.grid(True)
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    #plt.show()
    f.savefig(f"ShieldsDiagram_{planet.name}.pdf", bbox_inches='tight')
    print(f'Figure saved (ShieldsDiagram_{planet.name}.pdf)')
    

def calc_suspension_threshold(rho_s, rho, nu, g, D50):
    '''
    FUNCTION TO CALCULATE THE THRESHOLD FOR SUSPENSION
    AKA WHEN GRAINS ARE TRANSPORTED WITHOUT OFTEN
    INTERACTING WITH THE BED
        INPUTS: 
            rho_s                : density of grains (kg/m3) 
            rho                  : density of liquid (kg/m3) 
            nu                   : liquid kinematic viscotiy (m2/s)
            g                    : surface gravity (m/s2)
            D50                  : median grain diameter (m)
        OUTPUTS
            suspension_threshold : theta along the 
                                   curve defining the
                                   suspension threshold
                                   in the Shield's diagram
    '''
    
    P                    = 0.8;                                           # Critical Rouse number for 100% suspension (P = settling/turbulence)
    R                    = (rho_s / rho) - 1;                             # relative density of grains
    D_star               = (R * g * D50**3) / (nu**2);                    # dimensionless grain size
    log_D_star           = np.log10(D_star)                               # convert to log-space because that's how it was originally fit by Dietrich
    log_W_star           = (                                              # fit for dimensionless settling velocity (Dietrich, 1982)
        -3.76715
        + 1.92944 * log_D_star
        - 0.09815 * log_D_star**2
        - 0.00575 * log_D_star**3
        + 0.00056 * log_D_star**4
    )
    W_star               = 10**log_W_star                                 # convert back to linear space
    
    settling_velocity    = (W_star * R * g * nu) ** (1 / 3)               # settling velocity of grain (m/s)
    
    u_star               = settling_velocity/(P*kappa);                   # when the settling velocity is equal to the grain's velocity grain remains in suspension
    
    suspension_threshold = (rho * u_star**2 / ((rho_s - rho) * g * D50))
    
    return suspension_threshold
    

def calc_critical_shields(Re_particle): 
    '''
    FUNCTION TO CALCULATE THE CRTICIAL SHIELDS
        INPUTS:
            Re_p     : particle Reynolds number (dimensionless)
        OUTPUT:
            theta_c  : critical Shields value for entrainment (dimensionless)
    '''

    critical_shields = 0.5 * (0.22 * Re_particle**(-0.6) + 0.06 * 10**(-7.7 * Re_particle**(-0.6)));
     
    return critical_shields

    
def calculate_particle_Re(rho_s, rho, nu, g, D50):
    '''
    FUNCTION TO CALCULATE THE PARTICLE REYNOLDS NUMBER
        INPUTS: 
            rho_s : density of grains (kg/m3) 
            rho   : density of liquid (kg/m3) 
            nu    : liquid kinematic viscotiy (m2/s)
            g     : surface gravity (m/s2)
            D50   : median grain diameter (m)
        OUTPUT:
            Re_p  : particle Reynolds number (dimensionless)
                
    '''
    
    spec_weight = (rho_s / rho) - 1;
    Re_p        = np.sqrt(spec_weight * g * D50**3) / nu;
    
    return Re_p

def calc_shields_number(rho_s, rho, g, depth, u,D50):
    '''
    FUNCTION TO CALCULATE THE SHIELDS NUMBER FOR 
    GRAINS IN A FLOW
        INPUTS:
            rho_s : density of grains (kg/m3)
            rho   : density of liquid (kg/m3)
            g     : surface gravity (m/s2)
            depth : depth of flow (m)
            u     : flow velocity (m/s)
        OUTPUT
            theta : shields number (dimensionless)   
    '''

    z0     = D50 / 30                                          # boundary layer thickness for hydraulically rough flow (good assumption for rivers)
    u_star = (kappa * abs(u)) / np.log((0.37 * depth) / z0)    # Law of the Wall
    theta  = (rho * u_star**2) / ((rho_s - rho) * g * D50)     # dimensionless shear stress aka Shields number

    return theta
    
class Planet:
    '''
    Class contains information on planet of interest
        (1) Grain density 
        (2) Liquid density
        (3) Liquid kinematic viscosity
        (4) Gravity
        (5) Planet name (string)
    '''
    def __init__(self,rho_s,rho,nu,g,pname):
        self.rho_s = rho_s;
        self.rho   = rho;
        self.nu    = nu;
        self.g     = g;
        self.name  = pname;


def main(planet_name,data_filename) -> None:
   
   # populate planetary information into class (change values as desired 
   # or add more alternatives as needed)
    if planet_name == 'Titan':
        ice_grains                   = 940;       # kg/m3
        Ontario_Lacus_liquid_density = 588.15;    # kg/m3
        Ontario_Lacus_liquid_visc    = 8.08e-7;   # m2/s
        Titan_gravity                = 1.352;     # m/s2
        planet                       = Planet(ice_grains, Ontario_Lacus_liquid_density, Ontario_Lacus_liquid_visc, Titan_gravity,'Titan')
    elif planet_name == 'Earth':
        quartz_grains                = 2620;      # kg/m3
        water_18C_density            = 998.57;    # kg/m3
        water_18C_visc               = 1.0533e-6; # m2/s
        Earth_gravity                = 9.81;      # m/s2
        planet                       = Planet(quartz_grains, water_18C_density, water_18C_visc, Earth_gravity,'Earth')
    else:
        raise ValueError("Defined planets: 'Titan' or 'Earth'")
    
    # load in data
    data = np.loadtxt(data_filename,skiprows=1)   # skips first row of column names
    flow_velocity                   = data[:, 0]  # First column  : velocity
    grain_sizes                     = data[:, 1]  # Second column : grain size
    flow_depth                      = data[:, 2]  # Third column  : flow depth
    
    # make plot!
    MAKE_SHIELDS_DIAGRAM(planet,grain_sizes,flow_velocity,flow_depth);
    
    
if __name__ == '__main__':
    
    if len(sys.argv) != 3:
        print('NOT ENOUGH INPUTS\nexample use: python shields_curve.py Titan data.txt')
        sys.exit(1)
        
    planet_name       = sys.argv[1]
    data_filename     = sys.argv[2]
    
    main(planet_name, data_filename)
    
    
