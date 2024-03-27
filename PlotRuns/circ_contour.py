import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import matplotlib.font_manager as fm

#-- Generate Data -----------------------------------------
data = loadmat('New_2_10_results.mat')

[THETA, RHO] = np.meshgrid(data['dir_bins'], data['freqs']);
#-- Plot... ------------------------------------------------
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
lvl = np.linspace(0, 10, 11)
cplot = ax.contourf(THETA,RHO,np.transpose((data['pp'])),levels=lvl,cmap='YlOrRd')
#cline = ax.contour(THETA, RHO, np.transpose(data['pp']), colors='grey')
cbar = plt.colorbar(cplot)
ax.set_yscale('log')
cbar.set_label('Energy Density',fontname='cmr10')

ax.set_ylabel('Frequency [Hz]',fontname='cmr10')

plt.title('u = 5 m/s from 0 deg',fontname='cmr10')
ax.set_rticks([0.01,0.1,1,10,50])
ax.set_xlabel(r'$\theta$', fontname='cmr10')
ax.grid(True)
plt.rcParams['font.family'] = 'cmr10'
#plt.show()

plt.savefig("New_2_10.svg", format="svg")
