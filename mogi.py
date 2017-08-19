#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 12:40:08 2017

@author: gps
"""

"""
Code for computing mogi model as a function of horintal distance from 
source (radius) and plotting

Nick Voss 
USF geodesy lab

October 2015

Bruce Thomas modification june 2017
"""

import numpy as np
import matplotlib.pylab as plt

m2017u_y = [2140372.66466,
2140269.57525,
2140147.84576,
2140051.00059,
2140478.18352,
2140579.56606,
2140662.30476,
2140762.95566,
2140865.98344,
2140978.89289,
2141080.16454,
2141168.98073,
2141259.56711,
2141276.59354,
2141281.14994,
2141356.13711,
2141378.49758,
2141472.42172,
2141582.17563,
2141673.24251,
2141764.6638,
2141418.6039,
2141862.27082,
2141965.4271,
2142069.84716,
2142125.87269,
2142268.23351,
2142377.86567,
2142483.67936,
2142542.34682,
2142606.66268,
2142666.2272,
2142742.77337,
2142856.08138,
2142961.71682,
2143059.2705,
2143166.53099,
2143248.39674,
2143339.44709,
2143399.44476,
2143675.37971,
2143783.5128,
2143935.42013,
2139987.28731,
2139939.72837,
2139861.02239,
2139752.77188,
2139657.99369,
2139573.46698,
2139476.15004,
2139430.56812,
2139328.47529,
2139215.2412,
2139113.97658,
2138983.4869,
2138881.08446,
2138786.42731,
2138684.755,
2138585.40038,
2138491.50559
]

poisson_ratio = 0.25
mu = 30.0
lmda = 2 * mu * poisson_ratio / (1 - 2 * poisson_ratio)
alpha = (lmda + mu) / (lmda + 2 * mu)

#x = np.arange(0, 14, .01)#set up array of horizontal values to compute 
x = m2017u_y
n = len(x)
#alpha = 0.2 #km radius of body
d = 1.0 #km depth
delV = 0.2 #km^3 change in volume
v = .25 #poissons ratio

#calculate horizontal displacement
h = [0] * n
for i in range(n):
    h[i] = (3.0/4.0)*(delV/np.pi)*((x[i]*(10**(-3)))/(((x[i]*(10**(-3)))**2+d**2)**(3/2))) #units of km
h = [i*100000 for i in h] # convert km to cm

#calculate vertical displacement 
v = [0] * n
for i in range(n):
    v[i] = (3.0/4.0)*(delV/np.pi)*(d/(((x[i]*10**(-3))**2+d**2)**(3/2)))
v = [i*100000 for i in v] #convert km to cm
# plot results
ax = plt.subplot(211)
plt.plot(x,v)
plt.title('Displacement due to a Mogi Source')
plt.ylabel('Vertical Displacement in (cm)')

plt.subplot(212, sharex = ax)
plt.plot(x,h)
plt.xlabel('Horizontal distance from center of sphere (m)')
plt.ylabel('Horizontal displacement (cm)')
plt.savefig('mogiPython.png')
plt.show()