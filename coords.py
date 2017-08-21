#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Kilauea_Project
@author: bruce.eo.thomas
"""

"""
Script for the GPS data for each year.
You will have to change the PATHS to data (excel or txt) !  
Could be done automatic.
"""

import os
import matplotlib.pyplot as plt
# Get all the files in the current working directory: optionnal, just to check if the directory is the good one 
cwd = os.getcwd()
files = os.listdir(cwd)
print("Files in '%s': %s" % (cwd, files))

import numpy as np
# Excel modul 
import xlrd
from openpyxl import load_workbook
# Name variables
from collections import defaultdict

nsites = 66

#--------------------------------------------------------------------------
# Open useful data.
#-------------------------------------------------------------------------- 

# Open Excel
path = "/gps/Bruce/KilaueaKoaeProject/data/data_koae/excel"
k2003 = xlrd.open_workbook(path+"/KOAE_2003.xlsx")
k2004 = xlrd.open_workbook(path+"/KOAE_2004.xlsx")
k2006 = xlrd.open_workbook(path+"/KOAE_2006.xlsx")
k2007 = xlrd.open_workbook(path+"/KOAE_2007.xlsx")
k2008 = xlrd.open_workbook(path+"/KOAE_2008.xlsx")
k2009 = xlrd.open_workbook(path+"/KOAE_2009.xlsx")
k2011 = xlrd.open_workbook(path+"/KOAE_2011.xlsx")
k2017 = xlrd.open_workbook(path+"/KOAE_2017.xlsx")
#lendata = load_workbook(filename=path+"/KOAE_2003.xlsx", read_only=True)

# Open all sheets
sheets2003 = k2003.sheet_names()
sheets2004 = k2004.sheet_names()
sheets2006 = k2006.sheet_names()
sheets2007 = k2007.sheet_names()
sheets2008 = k2008.sheet_names()
sheets2009 = k2009.sheet_names()
sheets2011 = k2011.sheet_names()
sheets2017 = k2017.sheet_names()

# 2003
s2003 = k2003.sheet_by_name(sheets2003[2])
# 2004
s2004 = k2004.sheet_by_name(sheets2004[2])
# 2006
s2006 = k2006.sheet_by_name(sheets2006[2])
# 2007
s2007 = k2007.sheet_by_name(sheets2007[2])
# 2008
s2008 = k2008.sheet_by_name(sheets2008[3])
# 2009
s2009 = k2009.sheet_by_name(sheets2009[2])
# 2011
s2011 = k2011.sheet_by_name(sheets2011[2])
# 2017
s2017gamit = k2017.sheet_by_name(sheets2017[2])
s2017lgo = k2017.sheet_by_name(sheets2017[4])

#--------------------------------------------------------------------------
# Convert in Python transposed matrix.
#-------------------------------------------------------------------------- 

"""
Exemple of a matrix per year
For one year
array([[ id site ...],
       [ x ...],
       [ y ...],
       [ z ...]])
"""

"""
mm = defaultdict(int)
for year in [2003, 2017]:
    print('s%04d' % year)
    mm['mm%04d' % year] = np.zeros((12,nsites))
    for i in range(12):
        for j in range(nsites):
            mm['mm%04d' % year][i][j] = mm['s%04d' % year].cell_value(j,i)
    print(mm['mm%04d' % year])
"""



# Copy in Python matrix

m2003 = np.zeros((12,nsites))
for i in range(12):
    for j in range(nsites):
        m2003[i][j] = s2003.cell_value(j,i)

m2004 = np.zeros((12,nsites))
for i in range(12):
    for j in range(nsites):
        m2004[i][j] = s2004.cell_value(j,i)
        
m2006 = np.zeros((12,nsites))
for i in range(12):
    for j in range(nsites):
        m2006[i][j] = s2006.cell_value(j,i)
        
m2007 = np.zeros((12,nsites))
for i in range(12):
    for j in range(nsites):
        m2007[i][j] = s2007.cell_value(j,i)

m2008 = np.zeros((12,nsites))
for i in range(12):
    for j in range(nsites):
        m2008[i][j] = s2008.cell_value(j,i)

m2009 = np.zeros((12,nsites))
for i in range(12):
    for j in range(nsites):
        m2009[i][j] = s2009.cell_value(j,i)

m2011 = np.zeros((12,nsites))
for i in range(12):
    for j in range(nsites):
        m2011[i][j] = s2011.cell_value(j,i)
        
m2017gamit = np.zeros((12,nsites))
for i in range(12):
    for j in range(nsites):
        m2017gamit[i][j] = s2017gamit.cell_value(j,i)

m2017lgo = np.zeros((12,nsites))
for i in range(12):
    for j in range(nsites):
        m2017lgo[i][j] = s2017lgo.cell_value(j,i)



# Matrix cartesian X, Y, Z in meters

m2003c = np.zeros((4,nsites))
for j in range(nsites):
    m2003c[0][j] = m2003[0][j]
for i in range(1, 4):
    for j in range(nsites):
        m2003c[i][j] = m2003[i+3][j]

m2004c = np.zeros((4,nsites))
for j in range(nsites):
    m2004c[0][j] = m2004[0][j]
for i in range(1, 4):
    for j in range(nsites):
        m2004c[i][j] = m2004[i+3][j]

m2006c = np.zeros((4,nsites))
for j in range(nsites):
    m2006c[0][j] = m2006[0][j]
for i in range(1, 4):
    for j in range(nsites):
        m2006c[i][j] = m2006[i+3][j]

m2007c = np.zeros((4,nsites))
for j in range(nsites):
    m2007c[0][j] = m2007[0][j]
for i in range(1, 4):
    for j in range(nsites):
        m2007c[i][j] = m2007[i+3][j]

m2008c = np.zeros((4,nsites))
for j in range(nsites):
    m2008c[0][j] = m2008[0][j]
for i in range(1, 4):
    for j in range(nsites):
        m2008c[i][j] = m2008[i+3][j]

m2009c = np.zeros((4,nsites))
for j in range(nsites):
    m2009c[0][j] = m2009[0][j]
for i in range(1, 4):
    for j in range(nsites):
        m2009c[i][j] = m2009[i+3][j]

m2011c = np.zeros((4,nsites))
for j in range(nsites):
    m2011c[0][j] = m2011[0][j]
for i in range(1, 4):
    for j in range(nsites):
        m2011c[i][j] = m2011[i+3][j]
        
m2017gamitc = np.zeros((4,nsites))
for j in range(nsites):
    m2017gamitc[0][j] = m2017gamit[0][j]
for i in range(1, 4):
    for j in range(nsites):
        m2017gamitc[i][j] = m2017gamit[i+3][j]

m2017lgoc = np.zeros((4,nsites))
for j in range(nsites):
    m2017lgoc[0][j] = m2017lgo[0][j]
for i in range(1, 4):
    for j in range(nsites):
        m2017lgoc[i][j] = m2017lgo[i+3][j]



# Matrix latitude longitude in deg and h in meters 

m2003l = np.zeros((4,nsites))
for j in range(nsites):
    m2003l[0][j] = m2003[0][j]
for i in range(1, 4):
    for j in range(nsites):
        m2003l[i][j] = m2003[i+6][j]

m2004l = np.zeros((4,nsites))
for j in range(nsites):
    m2004l[0][j] = m2004[0][j]
for i in range(1, 4):
    for j in range(nsites):
        m2004l[i][j] = m2004[i+6][j]

m2006l = np.zeros((4,nsites))
for j in range(nsites):
    m2006l[0][j] = m2006[0][j]
for i in range(1, 4):
    for j in range(nsites):
        m2006l[i][j] = m2006[i+6][j]

m2007l = np.zeros((4,nsites))
for j in range(nsites):
    m2007l[0][j] = m2007[0][j]
for i in range(1, 4):
    for j in range(nsites):
        m2007l[i][j] = m2007[i+6][j]

m2008l = np.zeros((4,nsites))
for j in range(nsites):
    m2008l[0][j] = m2008[0][j]
for i in range(1, 4):
    for j in range(nsites):
        m2008l[i][j] = m2008[i+6][j]

m2009l = np.zeros((4,nsites))
for j in range(nsites):
    m2009l[0][j] = m2009[0][j]
for i in range(1, 4):
    for j in range(nsites):
        m2009l[i][j] = m2009[i+6][j]

m2011l = np.zeros((4,nsites))
for j in range(nsites):
    m2011l[0][j] = m2011[0][j]
for i in range(1, 4):
    for j in range(nsites):
        m2011l[i][j] = m2011[i+6][j]
        
m2017gamitl = np.zeros((4,nsites))
for j in range(nsites):
    m2017gamitl[0][j] = m2017gamit[0][j]
for i in range(1, 4):
    for j in range(nsites):
        m2017gamitl[i][j] = m2017gamit[i+6][j]

m2017lgol = np.zeros((4,nsites))
for j in range(nsites):
    m2017lgol[0][j] = m2017lgo[0][j]
for i in range(1, 4):
    for j in range(nsites):
        m2017lgol[i][j] = m2017lgo[i+6][j]

# Difference between gamit and lgo
m2017ddl = np.zeros((4,nsites))
for j in range(nsites):
    m2017ddl[0][j] = m2017gamit[0][j]
for i in range(1, 4):
    for j in range(nsites):
        m2017ddl[i][j] = m2017gamit[i+6][j] - m2017lgo[i+6][j] 



# Matrix utm x y in meters 

m2003u = np.zeros((3,nsites))
for j in range(nsites):
    m2003u[0][j] = m2003[0][j]
for i in range(1, 3):
    for j in range(nsites):
        m2003u[i][j] = m2003[i+9][j]

m2004u = np.zeros((3,nsites))
for j in range(nsites):
    m2004u[0][j] = m2004[0][j]
for i in range(1, 3):
    for j in range(nsites):
        m2004u[i][j] = m2004[i+9][j]

m2006u = np.zeros((3,nsites))
for j in range(nsites):
    m2006u[0][j] = m2006[0][j]
for i in range(1, 3):
    for j in range(nsites):
        m2006u[i][j] = m2006[i+9][j]

m2007u = np.zeros((3,nsites))
for j in range(nsites):
    m2007u[0][j] = m2007[0][j]
for i in range(1, 3):
    for j in range(nsites):
        m2007u[i][j] = m2007[i+9][j]

m2008u = np.zeros((3,nsites))
for j in range(nsites):
    m2008u[0][j] = m2008[0][j]
for i in range(1, 3):
    for j in range(nsites):
        m2008u[i][j] = m2008[i+9][j]

m2009u = np.zeros((3,nsites))
for j in range(nsites):
    m2009u[0][j] = m2009[0][j]
for i in range(1, 3):
    for j in range(nsites):
        m2009u[i][j] = m2009[i+9][j]

m2011u = np.zeros((3,nsites))
for j in range(nsites):
    m2011u[0][j] = m2011[0][j]
for i in range(1, 3):
    for j in range(nsites):
        m2011u[i][j] = m2011[i+9][j]
        
m2017gamitu = np.zeros((3,nsites))
for j in range(nsites):
    m2017gamitu[0][j] = m2017gamit[0][j]
for i in range(1, 3):
    for j in range(nsites):
        m2017gamitu[i][j] = m2017gamit[i+9][j]    

m2017lgou = np.zeros((3,nsites))
for j in range(nsites):
    m2017lgou[0][j] = m2017lgo[0][j]
for i in range(1, 3):
    for j in range(nsites):
        m2017lgou[i][j] = m2017lgo[i+9][j] 
        
# Difference between gamit and lgo
m2017ddu = np.zeros((3,nsites))
for j in range(nsites):
    m2017ddu[0][j] = m2017gamit[0][j]
for i in range(1, 3):
    for j in range(nsites):
        m2017ddu[i][j] = m2017gamit[i+9][j] - m2017lgo[i+9][j] 

        
        
mYearsXY = [m2003u, m2004u, m2006u, m2007u, m2008u, m2009u, m2011u, m2017gamitu]
mYearsZ = [m2003l, m2004l, m2006l, m2007l, m2008l, m2009l, m2011l, m2017gamitl]
len_mYears = len(mYearsXY)
# use gamit for 2017 !!!
mYears_name = [2003, 2004, 2006, 2007, 2008, 2009, 2011, 2017]



#--------------------------------------------------------------------------
# Table per site (1 to 66) with coordinates for each year.
#--------------------------------------------------------------------------

"""
Exemple of a matrix for site evolution
For one site
array([[ year ...],
       [ x ...],
       [ y ...],
       [ z ...]])
"""

# Site evolution
se = defaultdict(int)

# Loop on each site

for site in range(1, nsites+1):
    # Create matrix changing name 
    se['se_%02d' % site] = np.empty((4,len_mYears))
    se['se_%02d' % site][:] = np.NAN   
    # Loop on each year of data
    for year in range(len_mYears):
        m = mYearsXY[year]
        # Loop to find the site
        for i in range(len(m[0])):
            se['se_%02d' % site][0][year] = mYears_name[year]
            # Check if data for this site
            if site == m[0][i]:
                se['se_%02d' % site][1][year] = m[1][i]
                se['se_%02d' % site][2][year] = m[2][i]
    for year in range(len_mYears):
        m = mYearsZ[year]
        # Loop to find the site
        for i in range(len(m[0])):
            # Check if data for this site
            if site == m[0][i]:
                se['se_%02d' % site][3][year] = m[3][i]
    
                
              
site = int(input('Choose a site : '))
f, axarr = plt.subplots(3, sharex=True)
axarr[0].set_title('Site'+str(site))
axarr[0].scatter(mYears_name, se['se_%02d' % site][1])
axarr[1].scatter(mYears_name, se['se_%02d' % site][2])
axarr[2].scatter(mYears_name, se['se_%02d' % site][3])


"""
for key, value in sorted(se.items()):
    print (key, value)
"""


#--------------------------------------------------------------------------
# Table of displacements for each point. 
#--------------------------------------------------------------------------

"""
This is between first year of measure and last year of measure.
Not always the same years for each site, so there will be many matrix like that.

Exemple of a matrix between 2 years.
For 2 years
array([[ id site ...],
       [ dx ...],
       [ dy ...],
       [ dz ...]])
"""

# Define matrix of years 
m2017_2003 = np.empty((4,nsites))
m2017_2003[:] = np.NAN
m2017_2006 = np.empty((4,nsites))
m2017_2006[:] = np.NAN
m2017_2009 = np.empty((4,nsites))
m2017_2009[:] = np.NAN
m2017_2011 = np.empty((4,nsites))
m2017_2011[:] = np.NAN
m2006_2003 = np.empty((4,nsites))
m2006_2003[:] = np.NAN
m2009_2003 = np.empty((4,nsites))
m2009_2003[:] = np.NAN

# Loop on each site
for site in range(1, nsites+1):
    # First lien for id sites
    m2017_2003[0][site-1] = site
    m2017_2006[0][site-1] = site
    m2017_2009[0][site-1] = site
    m2017_2011[0][site-1] = site
    m2006_2003[0][site-1] = site
    m2009_2003[0][site-1] = site
    # Test to know which year of start 
    if se['se_%02d' % site][1][0] != 'nan' and se['se_%02d' % site][1][-1] != 'nan':
        # m2017_2003
        m2017_2003[1][site-1] = - se['se_%02d' % site][1][0] + se['se_%02d' % site][1][-1]
        m2017_2003[2][site-1] = - se['se_%02d' % site][2][0] + se['se_%02d' % site][2][-1]
        m2017_2003[3][site-1] = - se['se_%02d' % site][3][0] + se['se_%02d' % site][3][-1]   
    if se['se_%02d' % site][1][2] != 'nan' and se['se_%02d' % site][1][-1] != 'nan':
        # m2017_2006
        m2017_2006[1][site-1] = - se['se_%02d' % site][1][2] + se['se_%02d' % site][1][-1]
        m2017_2006[2][site-1] = - se['se_%02d' % site][2][2] + se['se_%02d' % site][2][-1]
        m2017_2006[3][site-1] = - se['se_%02d' % site][3][2] + se['se_%02d' % site][3][-1] 
    if se['se_%02d' % site][1][5] != 'nan' and se['se_%02d' % site][1][-1] != 'nan':
        # m2017_2009
        m2017_2009[1][site-1] = - se['se_%02d' % site][1][5] + se['se_%02d' % site][1][-1]
        m2017_2009[2][site-1] = - se['se_%02d' % site][2][5] + se['se_%02d' % site][2][-1]
        m2017_2009[3][site-1] = - se['se_%02d' % site][3][5] + se['se_%02d' % site][3][-1] 
    if se['se_%02d' % site][1][6] != 'nan' and se['se_%02d' % site][1][-1] != 'nan':
        # m2017_2011
        m2017_2011[1][site-1] = - se['se_%02d' % site][1][6] + se['se_%02d' % site][1][-1]
        m2017_2011[2][site-1] = - se['se_%02d' % site][2][6] + se['se_%02d' % site][2][-1]
        m2017_2011[3][site-1] = - se['se_%02d' % site][3][6] + se['se_%02d' % site][3][-1] 
    if se['se_%02d' % site][1][0] != 'nan' and se['se_%02d' % site][1][2] != 'nan':
        # m2006_2003
        m2006_2003[1][site-1] = - se['se_%02d' % site][1][0] + se['se_%02d' % site][1][2]
        m2006_2003[2][site-1] = - se['se_%02d' % site][2][0] + se['se_%02d' % site][2][2]
        m2006_2003[3][site-1] = - se['se_%02d' % site][3][0] + se['se_%02d' % site][3][2] 
    if se['se_%02d' % site][1][0] != 'nan' and se['se_%02d' % site][1][5] != 'nan':
        # m2009_2003
        m2009_2003[1][site-1] = - se['se_%02d' % site][1][0] + se['se_%02d' % site][1][5]
        m2009_2003[2][site-1] = - se['se_%02d' % site][2][0] + se['se_%02d' % site][2][5]
        m2009_2003[3][site-1] = - se['se_%02d' % site][3][0] + se['se_%02d' % site][3][5] 


