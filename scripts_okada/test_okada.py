#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Kilauea_Project
@author: bruce.eo.thomas
"""

"""

"""


from okada_wrapper.okada_wrapper import dc3d0wrapper, dc3dwrapper
from numpy import linspace, zeros, log
from matplotlib.pyplot import contourf, contour,\
    xlabel, ylabel, title, colorbar, show, savefig
import matplotlib
import time

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = ['Computer Modern Roman']
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
matplotlib.rcParams['lines.linewidth'] = 1

def get_params():
    source_depth = 1.0
    obs_depth = 1.0
    poisson_ratio = 0.25
    mu = 1.0
    dip = 0.0
    lmda = 2 * mu * poisson_ratio / (1 - 2 * poisson_ratio)
    alpha = (lmda + mu) / (lmda + 2 * mu)
    return source_depth, obs_depth, poisson_ratio, mu, dip, alpha

def test_dc3d():
    source_depth, obs_depth, poisson_ratio, mu, dip, alpha = get_params()
    n = (1000, 1000)
    x = linspace(-4.0, 4.0, n[0])
    y = linspace(-4.0, 4.0, n[1])
    ux = zeros((n[0], n[1]))
    for i in range(n[0]):
        for j in range(n[1]):
            success, u, grad_u = dc3dwrapper(alpha,
                                               [x[i], y[j], -obs_depth],
                                               source_depth, dip,
                                               [-2.0, 2.0], [-1.0, 1.0],
                                               [0.0, 0.0, 1.0])
            assert(success == 0)
            ux[i, j] = u[0]

    levels = linspace(-0.5, 0.5, 21)
    cntrf = contourf(x, y, ux.T, levels = levels)
    contour(x, y, ux.T, colors = 'k', levels = levels, linestyles = 'solid')
    xlabel('x')
    ylabel('y')
    cbar = colorbar(cntrf)
    tick_locator = matplotlib.ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.set_label('$u_{\\textrm{x}}$')
    savefig("strike_slip.png")
    show()
    
if __name__ == '__main__':
    test_dc3d()
