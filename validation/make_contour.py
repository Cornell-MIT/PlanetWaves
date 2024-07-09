import numpy as np
import matplotlib.pyplot as plt

from scipy.io import loadmat
import matplotlib.font_manager as fm
import sys
'''
    make contour plots in polar coordinates for energy spectrograms
'''

def find_max_2d(array):
    max_val = float(0)
    for row in array:
        for val in row:
            if val > max_val:
                max_val = val
    return max_val

def find_min_2d(array):
    min_val = float('inf')
    for row in array:
        for val in row:
            if val < min_val:
                min_val = val
    return min_val


    
def main(matfile):

    '''
    matfile is saved such that spectrogram is a structure containing
    '''
    
    mat = loadmat(matfile)
    
    spectrogram = mat['spectrogram']

    
    spectrogram_dir = spectrogram['dir']
    spectrogram_freq = spectrogram['freq']
    spectrogram_energy9 = spectrogram['energy_9']
    spectrogram_energy5 = spectrogram['energy_5']
    
    spectrogram_dir = spectrogram_dir[0, 0]
    spectrogram_freq = spectrogram_freq[0, 0]
    spectrogram_energy9 = spectrogram_energy9[0, 0]
    spectrogram_energy9[spectrogram_energy9 == 0.1] = np.nan
    spectrogram_energy5 = spectrogram_energy5[0, 0]
    spectrogram_energy5[spectrogram_energy5 == 0.1] = np.nan
    
    [THETA, RHO] = np.meshgrid(spectrogram_dir, spectrogram_freq);
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    lvl = np.linspace(find_min_2d(spectrogram_energy9), find_max_2d(spectrogram_energy9), 11)
    cplot = ax.contour(THETA,RHO,(spectrogram_energy9),levels=lvl,colors=['red'])
    lvl2 = np.linspace(find_min_2d(spectrogram_energy5), find_max_2d(spectrogram_energy5), 11)
    cplot2 = ax.contour(THETA,RHO,(spectrogram_energy5),levels=lvl2,colors=['green'])
    #cbar = plt.colorbar(cplot)
    ax.set_yscale('log')
    #cbar.set_label('Energy Density',fontname='cmr10')
    ax.set_ylabel('Frequency [Hz]',fontname='cmr10')
    #plt.title(mytitle,fontname='cmr10')
    ax.set_rticks([0.01,0.1,1,10,40])
    ax.set_xlabel(r'$\theta$', fontname='cmr10')
    ax.grid(True)
    plt.rcParams['font.family'] = 'cmr10'
    plt.show()
   
    
if __name__ == "__main__":
    main(sys.argv[1])
