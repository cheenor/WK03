# -*- coding: utf-8 -*-
"""
Created on Tue Jul 03 14:34:21 2018

@author: jhchen
"""

from netCDF4 import Dataset
fpath='L:/DATA/ERA_Interim/X2.5/201505.nc'
f=Dataset(fpath,'a')
lon=f.variables['longitude'][:]
lat=f.variables['latitude'][:]
print f.variables