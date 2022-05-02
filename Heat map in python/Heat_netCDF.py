#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Heat_netCDF.py   Apr 2022

This file
- reads the .nc file of Tmax downloaded for a 0.01º resolution from the CHIRTS project
- reads the .tif file of Relative Humidity  downloaded from the SB project
- combines both in order to construct a heat preassure map following the National Ocean
and Atmospheric Administration’s (NOAA) guidelines by:
    1. Using Steadman's equation
    2. Using Rothfusz's equation
    3. Applying the first and the second adjustments
- identifies the days with extreme heat events per location following the
  ISO occupational heat stress criteria and the US National Weather Service’s 
  definition.

@author: naiacasina
"""
import netCDF4 as nc
import pickle
import matplotlib.pyplot as plt
import datetime as dt  
import pandas as pd
import numpy as np
# To unmask netcdf data
import numpy.ma as ma
# To read .tif RH files
import rasterio
from dateutil import rrule        
from datetime import datetime
import math
import pickle
import os

# Set working directory
os.chdir('/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Heat/Results/')       
# Read pickle
df = pd.read_pickle('my_df.pickle') 
heat_data = df.to_numpy()

fn = '/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Heat/Tmax_2016.nc'
nc_fid = nc.Dataset(fn)

# Extract data from NetCDF file
lats = nc_fid.variables['latitude'][:]  # extract the data
lons = nc_fid.variables['longitude'][:]
time = nc_fid.variables['time'][:]
Tmax = nc_fid.variables['Tmax'][:] 

# Unmask
Tmax.mask = ma.nomask
# Create list with RH name files
a = '20160101'    # Start date
b = '20161231'    # End date
date_list = []

for dt in rrule.rrule(rrule.DAILY,
                      dtstart=datetime.strptime(a, '%Y%m%d'),
                      until=datetime.strptime(b, '%Y%m%d')):
    date_list.append(dt.strftime('RH.%Y.%m.%d.tif'))

T_sh = Tmax.shape

heat_data = np.zeros(shape=(T_sh[1],T_sh[2]))

HImax_cels1 =  np.zeros(shape=(T_sh[1],T_sh[2]))
HImax_cels1[:] = np.NAN

HImax_cels2 =  np.zeros(shape=(T_sh[1],T_sh[2]))
HImax_cels2[:] = np.NAN


WBGT =  np.zeros(shape=(T_sh[1],T_sh[2]))
WBGT[:] = np.NAN

# Compute first value for HImax and WBGT
i = 183
fp = r"/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Heat/RH/"+date_list[i]
raster = rasterio.open(fp)
RH = raster.read(1)
for j in range(0,T_sh[1]):
    for k in range(0,T_sh[2]):
        RH_v = RH[j,k]
        Tmax_v = Tmax[i,j,k]*9/5+32    # convert into Fah to use Rothfusz eq
        if ((Tmax[i,j,k]!=-9999)and(not math.isnan(RH_v))):
            # Steadman’s equation
            HImax = (0.5*(Tmax_v+61+((Tmax_v-68)*1.2)+(0.094*RH_v))+Tmax_v)/2
            if (HImax>80): # Rothfusz equation
                HImax = -42.379 + 2.04901523*Tmax_v + 10.14333127*RH_v -  0.22475541*Tmax_v*RH_v - 0.00683783*Tmax_v**2 - 0.05481717*RH_v**2 + 0.00122874*Tmax_v**2*RH_v + 0.00085282*Tmax_v*RH_v**2 - 0.00000199*Tmax_v*RH_v**2
            # First adjustment
            if ((80<Tmax_v<112)and(RH_v<13)): 
                HImax = HImax - (0.25*(13-RH_v)*np.sqrt(17-np.abs(Tmax_v-95)))/(17)
            # Second adjustment
            elif ((80<Tmax_v<87)and(RH_v>85)):
                HImax = HImax + ((RH_v-85)/10)*(87-Tmax_v)/5
            # Fah to celsius
            HImax_cels1[j,k] = (HImax-32)*5/9
            WBGT[j,k] = -0.0034*HImax**2 + 0.96*HImax - 34


for i in range(184, 213):
    fp = r"/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Heat/RH/"+date_list[i]
    raster = rasterio.open(fp)   # open raster data for relative humidity
    RH = raster.read(1)       # take first band
    for j in range(0,T_sh[1]):
        for k in range(0,T_sh[2]):
            RH_v = RH[j,k]
            Tmax_v = Tmax[i,j,k]*9/5+32   # Convert to Fah
            if ((Tmax[i,j,k]!=-9999)and(not math.isnan(RH_v))):
                # Steadman’s equation
                HImax = (0.5*(Tmax_v+61+((Tmax_v-68)*1.2)+(0.094*RH_v))+Tmax_v)/2
                if (HImax>80): # Rothfusz equation
                    HImax = -42.379 + 2.04901523*Tmax_v + 10.14333127*RH_v -  0.22475541*Tmax_v*RH_v - 0.00683783*Tmax_v**2 - 0.05481717*RH_v**2 + 0.00122874*Tmax_v**2*RH_v + 0.00085282*Tmax_v*RH_v**2 - 0.00000199*Tmax_v*RH_v**2
                # First adjustment
                if ((80<Tmax_v<112)and(RH_v<13)): 
                    HImax = HImax - (0.25*(13-RH_v)*np.sqrt(17-np.abs(Tmax_v-95)))/(17)
                # Second adjustment
                elif ((80<Tmax_v<87)and(RH_v>85)):
                    HImax = HImax + ((RH_v-85)/10)*(87-Tmax_v)/5
                # Fah to celsius
                HImax_cels2[j,k] = (HImax-32)*5/9
                if ((WBGT[j,k]>30)or((HImax_cels1[j,k]>40.6)and(HImax_cels2[j,k]>40.6))):
                    # Identify number of extreme heat events per location
                    heat_data[j,k]+= 1
                WBGT[j,k] = -0.0034*HImax**2 + 0.96*HImax - 34
                HImax_cels1[j,k] = HImax_cels2[j,k]
    print(i)


os.chdir('/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Heat/Results/')       
df = pd.DataFrame(heat_data)
df.to_pickle('my_df.pickle')

df = pd.read_pickle('my_df.pickle')    


# Load data
heat_data_colab = np.load('heat_colab', allow_pickle=True)

# Merge colab with Spyder processed data
heat_map = np.add(heat_data_colab,heat_data)

# Quick view
from rasterio.plot import show
show(heat_map)


# Build raster from array
from rasterio.crs import CRS
new_transform = raster.transform

with rasterio.open("/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Rasters/heat.tif", "w",
                   driver = "GTiff",
                   height = heat_map.shape[0],
                   width = heat_map.shape[1],
                   count = 1,
                   dtype = heat_map.dtype,
                   crs = CRS.from_epsg(4326),
                   transform = new_transform) as dst:
    dst.write(heat_map,1)

          
