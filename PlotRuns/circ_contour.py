import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import matplotlib.font_manager as fm
import sys

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

def main(which_mat,wind_speed,plot_or_save='plot'):
    '''
    make circular polar plot of energy spectrogram
    e.g.
    > python circ_contour.py New_2_10_results.m -> plots
    > python circ_contour.py New_2_10_results.m save -> saves the plot as svg
    '''

    data = loadmat(which_mat)
    [THETA, RHO] = np.meshgrid(data['dir_bins'], data['freqs']);
    # Plot
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    lvl = np.linspace(find_min_2d(data['pp']), find_max_2d(data['pp']), 11)
    cplot = ax.contourf(THETA,RHO,np.transpose((data['pp'])),levels=lvl,cmap='YlOrRd')
    #cline = ax.contour(THETA, RHO, np.transpose(data['pp']), colors='grey')
    cbar = plt.colorbar(cplot)
    ax.set_yscale('log')
    cbar.set_label('Energy Density',fontname='cmr10')
    ax.set_ylabel('Frequency [Hz]',fontname='cmr10')
    mytitle = 'u = ' + wind_speed + ' m/s from 0 deg'
    plt.title(mytitle,fontname='cmr10')
    ax.set_rticks([0.01,0.1,1,10,40])
    ax.set_xlabel(r'$\theta$', fontname='cmr10')
    ax.grid(True)
    plt.rcParams['font.family'] = 'cmr10'

    if plot_or_save == 'save':
        name_fig = which_mat.removesuffix('_results.mat')
        plt.savefig(name_fig + ".svg", format="svg")
    else:
        plt.show()

if __name__ == "__main__":

    if len(sys.argv) == 4:
        which_mat = str(sys.argv[1])
        wind_speed = str(sys.argb[2])
        plot_or_save = str(sys.argv[3])
        main(which_mat,wind_speed,plot_or_save)
    elif len(sys.argv) == 3:
        which_mat = str(sys.argv[1])
        wind_speed = str(sys.argv[2])
        main(which_mat,wind_speed)
    else:
        raise Exception('Requires 2-3 arguments')

