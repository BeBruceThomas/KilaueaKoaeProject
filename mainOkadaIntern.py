#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Kilauea_Project
@author: bruce.eo.thomas
"""

"""
Okada Program
Run in Terminal
"""


import os
cwd = os.getcwd()
files = os.listdir(cwd)
print("Files in '%s': %s" % (cwd, files))


# Moduls imported
import numpy as np
import matplotlib.pyplot as plt
from scripts_okada.okada_wrapper.okada_wrapper import dc3d0wrapper, dc3dwrapper
"""
# nead scipy version 0.11.0 at least, not possible to use here because version 0.9.0 and upgrade doesn't work
import scipy 
print(scipy.__version__)
from scipy import optimize
from scipy.optimize import minimize
"""

# Choose the path to access data 
path = "/gps/Bruce/KilaueaKoaeProject"

os.chdir("/gps/Bruce/KilaueaKoaeProject")

# Load dataForOkada
from data import dataForOkada


from scripts.geographic import Geographic
geo = Geographic("geographic")

#--------------------------------------------------------------------------
# Fit an Okada : surface deformation due to a finite rectangular source.
#--------------------------------------------------------------------------



# Data
#--------------------------------------------------------------------------

m2017u = np.zeros((2, 66))
"""
with open(path+'/data/coords/m2017u_xmodif.dat') as f:
    fdata = [line.rstrip() for line in f]
with open(path+'/data/coords/m2017u_ymodif.dat') as f:
    fdata = [line.rstrip() for line in f]
""""""
m2017u_x = open(path+"/data/coords/m2017u_xmodif.dat", "r")
m2017u_y = open(path+"/data/coords/m2017u_ymodif.dat", "r")
for i in range(66):
    m2017u[0][i] = m2017u_x.readline()
    m2017u[1][i] = m2017u_y.readline()
m2017u_x.close()
m2017u_y.close()
"""
sx = [19.3431837944,
19.3422506702,
19.3411507623,
19.3402777396,
19.3441368616,
19.3450544805,
19.3457967819,
19.346701934,
19.3476266088,
19.3486451246,
19.3495574716,
19.3503529068,
19.3511659453,
19.3513182194,
19.3513537119,
19.3520332381,
19.3522362342,
19.3530807396,
19.354069174,
19.3548843478,
19.3557074447,
19.3525978368,
19.3565835199,
19.3575105878,
19.358452022,
19.3589565966,
19.3602352638,
19.361223611,
19.3621760582,
19.3627013622,
19.3632852774,
19.3638188487,
19.364510132,
19.3655268379,
19.366478685,
19.3673549415,
19.3683213761,
19.369053562,
19.3698776561,
19.3704110906,
19.3728909219,
19.3738895233,
19.3752569547,
19.3396913395,
19.339250081,
19.3385274397,
19.3375389476,
19.3366788491,
19.3359081524,
19.3350272435,
19.3346049714,
19.3336801721,
19.3326583054,
19.3317480469,
19.3305714852,
19.3296479408,
19.3287879814,
19.3278723596,
19.3269731617,
19.3261243646]

sy=[-155.2748394013,
-155.2750032434,
-155.2750414596,
-155.2749086251,
-155.2748356402,
-155.2746789446,
-155.2750771334,
-155.2753894104,
-155.2758572172,
-155.2759574908,
-155.2761424814,
-155.2766807601,
-155.2770896893,
-155.2772099524,
-155.2776624273,
-155.2774815836,
-155.2773982494,
-155.2777013839,
-155.2779278233,
-155.2785125307,
-155.2787198439,
-155.2774483646,
-155.2791579015,
-155.2795261023,
-155.2796584979,
-155.2797728311,
-155.2803391777,
-155.280484322,
-155.2807421608,
-155.2811062492,
-155.280864248,
-155.2812150923,
-155.2812197932,
-155.2817508325,
-155.2819278833,
-155.2823120674,
-155.2824961558,
-155.2830700192,
-155.2829301108,
-155.2836030452,
-155.2845846893,
-155.2828341525,
-155.2831964825,
-155.2757858956,
-155.2767236917,
-155.2776667973,
-155.2785306116,
-155.2788573605,
-155.2794388173,
-155.2795961193,
-155.2804422771,
-155.2806599467,
-155.2805933755,
-155.2802468515,
-155.2800880759,
-155.2799819712,
-155.2803848798,
-155.2801726072,
-155.280320725,
-155.2803817922]


sitex = [0]*60
sitey = [0]*60
for isite in range(60):
    result = geo.from_latlon(sx[isite], sy[isite])
    sitex[isite] = result[0]
    sitey[isite] = result[1]     
        
print(min(sitex), max(sitex), min(sitey), max(sitey))        

# Functions
#--------------------------------------------------------------------------


def okada_SWZR_fit():
    """
    Evaluates the misfit of an okada solution defined by the passed parameters to the slip (and errors) globally defined.
    """
    
    nsite = len(site_neu_err[0])
    
    slip_weights = np.zeros((3, nsite))
    for i in range(3):
        for j in range(nsite):
            slip_weights[i][j] = 1 / (site_neu_err[i][j]**2)
    
    # only for z first
    calc_slip = np.zeros((3, nsite))
    slip_misfit = np.zeros((3, nsite))
    
    for i in range(3):
        for j in range(nsite):
            site_slip = calc_SWZR_okada(dataForOkada.okada_start, site_neu_posn)
            calc_slip[i][j] = site_slip
            slip_misfit[i][j] = site_neu_slip[i][j] - calc_slip[i][j]
    
    misfit = np.zeros((1, nsite))
    
    for isite in range(15):            
        misfit[0][isite] = slip_misfit[0][isite] * slip_weights[2][isite] / (sum(slip_weights[2]))
    return misfit 


def get_params():
    """
    okada_pr are the okada parameters that matlab program gives us as the okada fit doesn't work on python
    """
    poisson_ratio = okada_pr[10]
    mu = 30
    lmda = (2 * mu * poisson_ratio) / (1 - 2 * poisson_ratio)
    alpha = (lmda + mu) / (lmda + 2 * mu)
    
    x0 = [ okada_pr[0], okada_pr[1], - okada_pr[2] ]
    depth = okada_pr[2]
    dip = okada_pr[4]
    strike_width = [ -okada_pr[5]/2, okada_pr[5]/2 ]
    dip_width = [ -okada_pr[6]/2, okada_pr[6]/2 ]
    dislocation = [ okada_pr[7], okada_pr[8], okada_pr[9] ]
    
    return alpha, x0, depth, dip, strike_width, dip_width, dislocation


def calc_SWZR_okada():
    """
    """
    alpha, x0, depth, dip, strike_width, dip_width, dislocation = get_params()
    n0 = len(sitex)
    x = [0] * n0
    for i in range(n0):
        x[i] = sitex[i]
    n1 = len(sitey)
    y = [0] * n1
    for i in range(n0):
        y[i] = sitey[i]
    ux = np.zeros((n0, n1))
    for i in range(n0):
        for j in range(n1):
            success, u, grad_u = dc3dwrapper(alpha, 
                                             [x[i], y[j], -1.0],
                                             depth, 
                                             dip,
                                             strike_width, 
                                             dip_width,
                                             dislocation
                                            )
            assert(success == 0)
            #ux[i, j] = u[0]
            #ux[i, i] = u[0]
            ux[i] = u[0]

    levels = np.linspace(-0.5, 0.5, 21)
    cntrf = plt.contourf(x, y, ux.T, levels = levels)
    plt.contour(x, y, ux.T, colors = 'k', levels = levels, linestyles = 'solid')
    plt.xlabel('x')
    plt.ylabel('y')
    cbar = plt.colorbar(cntrf)
    #tick_locator = plt.ticker.MaxNLocator(nbins=5)
    #cbar.locator = tick_locator
    #cbar.update_ticks()
    #cbar.set_label('$u_{\\textrm{x}}$')
    #plt.savefig("strike_slip.png")
    plt.show()

    



# Main Okada 
#--------------------------------------------------------------------------

# Okada parameters

okada_prA = dataForOkada.okada_intern[0]

len_bounds = len(okada_prA)

okada_pr = np.zeros(len_bounds)

for i in range(len_bounds):
    okada_pr[i] = okada_prA[i]



okada_params = okada_pr
calc_SWZR_okada()  



# Fit Okada 

"""
fun = okada_SWZR_fit()

pmin = lower_bounds # mimimum bounds
pmax = upper_bounds # maximum bounds
p_guess = (pmin + pmax)/2  
bounds = np.c_[pmin, pmax]  # [[pmin[0],pmax[0]], [pmin[1],pmax[1]]]

sol = scipy.optimize.minimize(fun, p_guess, args=(), bounds=bounds)   
    
if not sol.success:
    raise RuntimeError("Failed to solve")
okada_params = sol.x   
"""


"""
nsite = 15
calc_slip = np.zeros((1, nsite))
#for isite in range(nsite):
    #site_slip = calc_SWZR_okada(okada_params, site_neu_posn)
#test_dc3d()
    #calc_slip[0][isite] = site_slip
"""


"""
h_okada_vert=quiver(sitex-x0,sitey-y0,0*stepE',calc_slip(3,:)',1)
h_okada_horiz=quiver(sitex-x0,sitey-y0,calc_slip(2,:)',calc_slip(1,:)',1)
set(h_okada_horiz,'Color','k')
set(h_okada_vert,'Color',[.7 .7 .7])
"""
 





