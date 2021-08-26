"""
Code to evaluate the effect of the mesh size in the VAFR
@author: marc
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import uvrct as uv

n_points = 20
VAFR = np.ones(n_points)
cells=np.arange(1, n_points+1)**2

for i in range(n_points):
    # Set the number of cells in r and z
    ncr = ncz = int(np.sqrt(cells[i]))
    
    # Define list of variables required for the case
    case_variables = [17, 0.28, 0.33, 0.0425, 0.0115, ncr, ncz]

    # Compute the volume averaged fluence rate (VAFR) for the given case
    VAFR[i] = uv.compute_vafr(case_variables, model='LSI', nls=1000, optimize=False, ppm=0)

#plt.style.use('ggplot')
plt.style.use('seaborn-pastel')
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['figure.autolayout'] = True

plt.figure()
plt.plot(cells, VAFR,'o-')
plt.xlabel('Number of cells')
plt.ylabel('VARI [W/m$^2$]')
plt.show()