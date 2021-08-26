"""
Code to evaluate the effect of the mesh size in the VAFR
@author: marc
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import uvrct as uv

n_points = 100
cells, VAFR  = np.ones(n_points), np.ones(n_points)

for i in range(n_points):
    # Set the number of cells in r
    ncr = i+1
    
    # Define list of variables required for the case
    case_variables = P, eff, L, r_outer, r_inner, ncr, ncz = [17, 0.28, 0.33, 0.0425, 0.0115, ncr, 0]

    # Set the number of cells in z
    #case_variables[-1] = round(L/(r_outer/ncr))  
    case_variables[-1]  = ncr
    cells[i] = case_variables[-1]*ncr
    print(f'i is {i}, cells is {cells[i]}')
    # Compute the volume averaged fluence rate (VAFR) for the given case
    VAFR[i] = uv.compute_vafr(case_variables, model='LSI', nls=1000, optimize=False, ppm=300)

#plt.style.use('ggplot')
plt.style.use('seaborn-pastel')
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['figure.autolayout'] = True

# plt.figure()
# plt.plot(cells, VAFR,'o-')
# plt.xlabel('Number of cells')
# plt.ylabel('VAFR [W/m$^2$]')
# plt.show()

VAFR_err = np.ones(len(VAFR)-1)
 
for i in range(len(VAFR)-1):
    
    VAFR_err[i] = np.abs((VAFR[i+1]-VAFR[i])/VAFR[i])

plt.figure()
plt.semilogy(cells[1:], VAFR_err*100,'o-')
plt.xlabel('Number of cells')
plt.ylabel('VARI rate of change [%]')
plt.show()    

    