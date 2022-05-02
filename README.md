# Spatial Economics: Project

## Code in R
Here lie the main files to reproduce all the results in the document: reading of the data, creation of the rasters, matching of the rasters, main computations, creation of the map.

## Heat map in Python
Since I could not get access to proper heat preassure maps (that combined properly both maximum temperature and relative humidity, along with a fine resolution), I construct one in Heat_netCDF.py. The file:

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
  
  ## Results
  This folder contains two subfolders: Outcomes and Rasters. The _outcomes_ are the main figures/pdf files germane to the results commented in the term paper. The _rasters_ folder contains the final layers superposed to get the metric in .tif format. 
