#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Kilauea_Project
@author: bruce.eo.thomas
"""

"""
Main Program for Internship Kilauea 2017
Run in Console 
"""


# Moduls imported
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# Name variables
from collections import defaultdict

# Get all the files in the current working directory: optionnal, just to check if the directory is the good one 
cwd = os.getcwd()
files = os.listdir(cwd)
print("Files in '%s': %s" % (cwd, files))


# Choose the path to access data: have to find a solution to change only in the main only or directly emter in the interface 
path = "C:/Users/bruce/Desktop/KilaueaKoaeProject/KilaueaKoaeProject"


import coords
# Load BI_linefile
from data import bi


# Main program : all the run is done here
if __name__ == "__main__":
    
    ux03, uy03 = coords.m2003u[1], coords.m2003u[2]
    uy06 = coords.m2006u[2]
    uy09 = coords.m2009u[2]
    uy11 = coords.m2011u[2]
    uy17 = coords.m2017gamitu[2]
    
    line = [0] * coords.nsites
    name = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66]

    
    
    #--------------------------------------------------------------------------
    # Compare 2017gamit / 2017lgo
    #-------------------------------------------------------------------------- 
    
    dde = plt.figure()
    axdde = dde.add_subplot(111)   
          
    plt.axvline(x=2140269.57525397, color='k') # fault in 2
    plt.axvline(x=2141378.49758047, color='k') # fault in 17
    
    plt.plot(uy17, coords.m2017ddu[1], marker='+', linestyle='', color='b', label="2017e")     
    for i in range(coords.nsites):
        plt.annotate(int(name[i]), (uy17[i], coords.m2017ddu[1][i]))   

    plt.xlabel('Northing UTM (m) : --> South to North -->')
    plt.ylabel('Difference in Easting UTM (m)')
    plt.title('Difference in Easting between GAMIT and LGO for 2017')
    plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
    
    
    ddn = plt.figure()
    axddn = ddn.add_subplot(111)   
          
    plt.axvline(x=2140269.57525397, color='k') # fault in 2
    plt.axvline(x=2141378.49758047, color='k') # fault in 17
    
    plt.plot(uy17, coords.m2017ddu[2], marker='+', linestyle='', color='b', label="2017n")     
    for i in range(coords.nsites):
        plt.annotate(int(name[i]), (uy17[i], coords.m2017ddu[2][i]))   

    plt.xlabel('Northing UTM (m) : --> South to North -->')
    plt.ylabel('Difference in Northing UTM (m)')
    plt.title('Difference in Northing between GAMIT and LGO for 2017')
    plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
    
    
    ddh = plt.figure()
    axddh = ddh.add_subplot(111)   
          
    plt.axvline(x=2140269.57525397, color='k') # fault in 2
    plt.axvline(x=2141378.49758047, color='k') # fault in 17
    
    plt.plot(uy17, coords.m2017ddl[3], marker='+', linestyle='', color='b', label="2017h")     
    for i in range(coords.nsites):
        plt.annotate(int(name[i]), (uy17[i], coords.m2017ddl[3][i]))   

    plt.xlabel('Northing UTM (m) : --> South to North -->')
    plt.ylabel('Difference in Altitude (m)')
    plt.title('Difference in Altitude between GAMIT and LGO for 2017')
    plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
    
    
    
    #--------------------------------------------------------------------------
    # View in 3D: color for each year: could try only the faults
    #--------------------------------------------------------------------------    
  
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
    """   
    # Loop on each year
    for year in range(coords.len_mYears):                
        x = coords.mYearsXY[year][1]
        y = coords.mYearsXY[year][2]
        z = coords.mYearsZ[year][3]
        cmap = mpl.cm.autumn
        ax1.scatter(x, y, z, c=cmap(year / float(coords.len_mYears)), marker='o')
        plt.xlim([260754, 260755])
        plt.ylim([2141355, 2141356])
        #plt.zlim([1023, 1024])
    """             
    # Just 2017                   
    x = coords.m2017gamitu[1]
    y = coords.m2017gamitu[2]
    z = coords.m2017gamitl[3]
    cmap = mpl.cm.autumn
    ax1.scatter(x, y, z, c='red', marker='o')    
    
    ax1.set_xlabel('Easting UTM (m)')
    ax1.set_ylabel('Northing UTM (m)')
    ax1.set_zlabel('Altitude (m)')
    
    
    """    
    #--------------------------------------------------------------------------
    # View in 3D: color for each site.
    #--------------------------------------------------------------------------  
    
    #can also do it only for each site, have to put orientation

    
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection='3d')
        
    # Loop on each site
    for site in range(1, 67):                
                       
        x = coords.se['se_%02d' % site][1]
        y = coords.se['se_%02d' % site][2]
        z = coords.se['se_%02d' % site][3]
        
        cmap = mpl.cm.autumn
        ax2.scatter(x, y, z, c=cmap(site / float(67)), marker='o')
        
        #ax.scatter(xt, yt, zt, c='b', marker='^')
        
    ax2.set_xlabel('X Label')
    ax2.set_ylabel('Y Label')
    ax2.set_zlabel('Z Label')
    """
    

    #--------------------------------------------------------------------------
    # Cumulative Distance from line.
    #--------------------------------------------------------------------------    
        
    # dx utm 
    d1_03_17 = coords.m2017_2003[1]
    d1_06_17 = coords.m2017_2006[1]
    d1_09_17 = coords.m2017_2009[1]
    d1_11_17 = coords.m2017_2011[1]
    d1_03_06 = coords.m2006_2003[1]
    d1_03_09 = coords.m2009_2003[1]
    
    # dy utm
    d2_03_17 = coords.m2017_2003[2]
    d2_06_17 = coords.m2017_2006[2]
    d2_09_17 = coords.m2017_2009[2]
    d2_11_17 = coords.m2017_2011[2]
    d2_03_06 = coords.m2006_2003[2]
    d2_03_09 = coords.m2009_2003[2]
    
    # dz utm
    d3_03_17 = coords.m2017_2003[3]
    d3_06_17 = coords.m2017_2006[3]
    d3_09_17 = coords.m2017_2009[3]
    d3_11_17 = coords.m2017_2011[3]
    d3_03_06 = coords.m2006_2003[3]
    d3_03_09 = coords.m2009_2003[3]
    
   
    
    # Xcart
    X = plt.figure()
    axX = X.add_subplot(111)
          
    plt.axhline(y=0.0, color='k')
    plt.axvline(x=2140269.57525397, color='k') # fault in 2
    plt.axvline(x=2141378.49758047, color='k') # fault in 17
    
    plt.plot(uy03, d1_03_17, marker='+', linestyle='', color='b', label="2017-2003")     
    for i in range(coords.nsites):
        plt.annotate(int(name[i]), (uy03[i], d1_03_17[i]))     
       
    plt.plot(uy06, d1_06_17, marker='+',  linestyle='', color='g', label="2017-2006")  
    for i in range(coords.nsites):
        plt.annotate(int(name[i]), (uy06[i], d1_06_17[i]))    
    """
    plt.plot(uy03, d1_03_06, marker='+', linestyle='', color='m', label="2006-2003") 
    for i in range(coords.nsites):
        plt.annotate(int(name[i]), (uy03[i], d1_03_06[i]))    
    """
    plt.plot(uy09, d1_09_17, marker='+', linestyle='', color='r', label="2017-2009")     
    for i in range(coords.nsites):
        plt.annotate(int(name[i]), (uy09[i], d1_09_17[i])) 
    """    
    plt.plot(uy11, d1_11_17, marker='+', linestyle='', color='c', label="2017-2011") 
    for i in range(coords.nsites):
        plt.annotate(int(name[i]), (uy11[i], d1_11_17[i]))    
    """
    plt.plot(uy03, d1_03_09, marker='+', linestyle='', color='y', label="2009-2003") 
    for i in range(coords.nsites):
        plt.annotate(int(name[i]), (uy03[i], d1_03_09[i]))   
          
    plt.xlabel('Northing UTM (m) : --> South to North -->')
    plt.ylabel('Cumulative Displacement in Easting UTM (m)')
    plt.title('Cumulative Displacement in Easting : years')
    plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
       
    
    # Ycart
    Y = plt.figure()
    axY = Y.add_subplot(111)   
          
    plt.axhline(y=0.0, color='k')
    plt.axvline(x=2140269.57525397, color='k') # fault in 2
    plt.axvline(x=2141378.49758047, color='k') # fault in 17
    
    plt.plot(uy03, d2_03_17, marker='+', linestyle='', color='b', label="2017-2003")     
    for i in range(coords.nsites):
        plt.annotate(int(name[i]), (uy03[i], d2_03_17[i]))   
      
    plt.plot(uy06, d2_06_17, marker='+',  linestyle='', color='g', label="2017-2006")  
    for i in range(coords.nsites):
        plt.annotate(int(name[i]), (uy06[i], d2_06_17[i]))    
    """
    plt.plot(uy03, d2_03_06, marker='+', linestyle='', color='m', label="2006-2003") 
    for i in range(coords.nsites):
        plt.annotate(int(name[i]), (uy03[i], d2_03_06[i]))    
    """
    plt.plot(uy09, d2_09_17, marker='+', linestyle='', color='r', label="2017-2009")     
    for i in range(coords.nsites):
        plt.annotate(int(name[i]), (uy09[i], d2_09_17[i])) 
    """    
    plt.plot(uy11, d2_11_17, marker='+', linestyle='', color='c', label="2017-2011") 
    for i in range(coords.nsites):
        plt.annotate(int(name[i]), (uy11[i], d2_11_17[i]))    
    """
    plt.plot(uy03, d2_03_09, marker='+', linestyle='', color='y', label="2009-2003") 
    for i in range(coords.nsites):
        plt.annotate(int(name[i]), (uy03[i], d2_03_09[i]))  
    
    plt.xlabel('Northing UTM (m) : --> South to North -->')
    plt.ylabel('Cumulative Displacement in Northing UTM (m)')
    plt.title('Cumulative Displacement in Northing : years')
    plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
    
    
    # Zcart
    Z = plt.figure()
    axZ = Z.add_subplot(111)
       
    plt.axhline(y=0.0, color='k')
    plt.axvline(x=2140269.57525397, color='k') # fault in 2
    plt.axvline(x=2141378.49758047, color='k') # fault in 17
    
    plt.plot(uy03, d3_03_17, marker='+', linestyle='', color='b', label="20017-2003")     
    for i in range(coords.nsites):
        plt.annotate(int(name[i]), (uy03[i], d3_03_17[i])) 
       
    plt.plot(uy06, d3_06_17, marker='+',  linestyle='', color='g', label="2017-2006")  
    for i in range(coords.nsites):
        plt.annotate(int(name[i]), (uy06[i], d3_06_17[i]))    
    """
    plt.plot(uy03, d3_03_06, marker='+', linestyle='', color='m', label="2006-2003") 
    for i in range(coords.nsites):
        plt.annotate(int(name[i]), (uy03[i], d3_03_06[i]))    
    """
    plt.plot(uy09, d3_09_17, marker='+', linestyle='', color='r', label="2017-2009")     
    for i in range(coords.nsites):
        plt.annotate(int(name[i]), (uy09[i], d3_09_17[i]))  
    """    
    plt.plot(uy11, d3_11_17, marker='+', linestyle='', color='c', label="2017-2011") 
    for i in range(coords.nsites):
        plt.annotate(int(name[i]), (uy11[i], d3_11_17[i]))    
    """
    plt.plot(uy03, d3_03_09, marker='+', linestyle='', color='y', label="2009-2003") 
    for i in range(coords.nsites):
        plt.annotate(int(name[i]), (uy03[i], d3_03_09[i]))     
    
    plt.xlabel('Northing UTM (m) : --> South to North -->')
    plt.ylabel('Cumulative Displacement in Altitude (m)')
    plt.title('Cumulative Displacement in Altitude : years')
    plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
    
    
    plt.grid()
    plt.show()
    

    
    
    #--------------------------------------------------------------------------
    # Plot Kilauea map and vectors. 
    #--------------------------------------------------------------------------
    

    plt.figure(figsize = (10, 8))
    # Plot faults?  
    plt.plot(bi.fx, bi.fy, 'k-')       
    # Plot Big Island
    plt.plot(bi.cx, bi.cy, 'k-')
    
    # Define limit axes
    axes = plt.gca()
    
    x0, y0, xe, ye = 254000, 2138000, 264000, 2148000
    # limits for the two big faults
    x0f1, y0f1, xef1, yef1 = 260000, 2140000, 262000, 2141000
    x0f2, y0f2, xef2, yef2 = 260600, 2141000, 260900, 2141500
    
    plt.xlim([x0, xe])
    plt.ylim([y0, ye])
    
    plt.plot(ux03, uy03, 'bo')
    
    for i in range(66):
        ux = ux03[i]
        uy = uy03[i]
        dux = d1_03_17[i]
        duy = d2_03_17[i]
        dh = d3_03_17[i]
        # Magenta arrows for elevation
        plt.arrow(ux, uy, 0, dh*1000, head_width=10, head_length=10, fc='m', ec='m', clip_on=False)
        # Red arrows for motions
        plt.arrow(ux, uy, dux*1000, duy*1000, head_width=10, head_length=10, fc='r', ec='r', clip_on=False)    
      
    plt.title("GPS Vectors")
