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




#--------------------------------------------------------------------------
# Fit an Okada : surface deformation due to a finite rectangular source.
#--------------------------------------------------------------------------



# Data
#--------------------------------------------------------------------------


site_neu_posn = np.zeros((3, 15))

site_neu_posnX = open(path+"/data/site_neu/site_neu_posnX.dat", "r")
for i in range(15):
    site_neu_posn[0][i] = site_neu_posnX.readline()
site_neu_posnX.close()
site_neu_posnY = open(path+"/data/site_neu/site_neu_posnY.dat", "r")
for i in range(15):
    site_neu_posn[1][i] = site_neu_posnY.readline()
site_neu_posnY.close()
site_neu_posnZ = open(path+"/data/site_neu/site_neu_posnZ.dat", "r")
for i in range(15):
    site_neu_posn[2][i] = site_neu_posnZ.readline()
site_neu_posnZ.close()

site_neu_slip = np.zeros((3, 15))

site_neu_slipX = open(path+"/data/site_neu/site_neu_slipX.dat", "r")
for i in range(15):
    site_neu_slip[0][i] = site_neu_slipX.readline()
site_neu_slipX.close()
site_neu_slipY = open(path+"/data/site_neu/site_neu_slipY.dat", "r")
for i in range(15):
    site_neu_slip[1][i] = site_neu_slipY.readline()
site_neu_slipY.close()
site_neu_slipZ = open(path+"/data/site_neu/site_neu_slipZ.dat", "r")
for i in range(15):
    site_neu_slip[2][i] = site_neu_slipZ.readline()
site_neu_slipZ.close()

site_neu_err = np.zeros((3, 15))

site_neu_errX = open(path+"/data/site_neu/site_neu_errX.dat", "r")
for i in range(15):
    site_neu_err[0][i] = site_neu_errX.readline()
site_neu_errX.close()
site_neu_errY = open(path+"/data/site_neu/site_neu_errY.dat", "r")
for i in range(15):
    site_neu_err[1][i] = site_neu_errY.readline()
site_neu_errY.close()
site_neu_errZ = open(path+"/data/site_neu/site_neu_errZ.dat", "r")
for i in range(15):
    site_neu_err[2][i] = site_neu_errZ.readline()
site_neu_errZ.close()

print(site_neu_posn)



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


def test_dc3d():
    """
    """
    alpha, x0, depth, dip, strike_width, dip_width, dislocation = get_params()
    n = [10, 10]
    x = np.linspace(lower_bounds[0], upper_bounds[0], n[0])
    print(x)
    y = np.linspace(lower_bounds[1], upper_bounds[1], n[1])
    print(y)
    ux = np.zeros((n[0], n[1]))
    for i in range(n[0]):
        for j in range(n[1]):
            success, u, grad_u = dc3dwrapper(alpha, 
                                             [x[i], y[j], - 1.0],
                                             depth, 
                                             dip,
                                             strike_width, 
                                             dip_width,
                                             dislocation
                                            )
            assert(success == 0)
            ux[i, j] = u[0]
            #print(u[0])

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
    

def calc_SWZR_okada():
    """
    """
    alpha, x0, depth, dip, strike_width, dip_width, dislocation = get_params()
    x = site_neu_posn[0]
    #x = mux
    n0 = len(x)
    print(x)
    y = site_neu_posn[1]
    #y = muy
    n1 = len(y)
    print(y)
    ux = np.zeros((n0, n1))
    for i in range(n0):
        #for j in range(10):
        success, u, grad_u = dc3dwrapper(alpha, 
                                         [x[i], y[i], -1.0],
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
        print(u[0])

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


    
def calc_slip():
    """
    """
    alpha, x0, depth, dip, strike_width, dip_width, dislocation = get_params()
    n = [100, 100]
    x = np.linspace(260038, 261038, n[0])
    print(x)
    y = np.linspace(2138491, 2143935, n[1])
    print(y)
    ux = np.zeros((n[0], n[1]))
    for i in range(n[0]):
        for j in range(n[1]):
            success, u, grad_u = dc3dwrapper(alpha, 
                                             [x[i], y[j], - 1.0],
                                             depth, 
                                             dip,
                                             strike_width, 
                                             dip_width,
                                             dislocation
                                            )
            assert(success == 0)
            ux[i, j] = u[0]
            #print(u[0])

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

upper_boundsA = dataForOkada.okada_initial_params[0]
okada_startA  = dataForOkada.okada_initial_params[1]
lower_boundsA = dataForOkada.okada_initial_params[3]
okada_prA = dataForOkada.okada_initial_params[2]

len_bounds = len(okada_prA)

upper_bounds = np.zeros(len_bounds)
okada_start = np.zeros(len_bounds)
lower_bounds = np.zeros(len_bounds)
okada_pr = np.zeros(len_bounds)

for i in range(len_bounds):
    upper_bounds[i] = upper_boundsA[i]
    okada_start[i] = okada_startA[i]
    lower_bounds[i] = lower_boundsA[i]
    okada_pr[i] = okada_prA[i]


okada_params = okada_pr

#test_dc3d()
#calc_SWZR_okada() 
calc_slip()



"""
# Create a sample fault and print out some information about it: use of CLAWPACK version for Okada.
fault, subfault = okaC.set_params(okada_start)
print ("This sample fault has %s meter of slip over a %s by %s km patch" % (subfault.slip,subfault.length/1e3,subfault.width/1e3))
print ("With shear modulus %4.1e Pa the seismic moment is %4.1e" % (subfault.mu, subfault.Mo()))
print ("   corresponding to an earthquake with moment magnitude %s" % fault.Mw())
print ("The depth at the top edge of the fault plane is %s km" % (subfault.depth/1e3)) 
"""

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
h_okada_vert=quiver(sitex-x0,sitey-y0,0*stepE',calc_slip(3,:)',1)
h_okada_horiz=quiver(sitex-x0,sitey-y0,calc_slip(2,:)',calc_slip(1,:)',1)
set(h_okada_horiz,'Color','k')
set(h_okada_vert,'Color',[.7 .7 .7])
"""
 





