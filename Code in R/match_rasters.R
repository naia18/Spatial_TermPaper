# ========================================================================
# Match_rasters.R    -   Naia Ormaza Zulueta   -  Apr 2022
# In this file:
# - We gather all the variables of interest for the analysis 
# - Read them into R through the terra package to check extents and resolutions
# - Match the resolution.
# - Match the extent.
# - Check CRS.
# - Create a raster file for each and save.
# ========================================================================

# Clear the environment
rm(list=ls()) 
# ------ Load libraries ------
packages <- c("terra", "raster", "tidyverse", "rasterVis", "ncdf4", 
              "lattice", "foreign", "rworldmap")
lapply(packages, require, character=TRUE)


# ------ Open (already processed) raster data ------
airQuality <- rast("/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Air Quality/Raster/V5GL02.HybridPM25c_0p10.Global.202001-202012.tif")
water <- rast("/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Water/WS_blue_monthly_rasters/WSbl_monthly_30m/ws_avg/hdr.adf")
biodiv <- rast("/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Biodiversity/lbii.asc")
pop <- rast("/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/World-Pop/gpw-v4-population-count-rev11_2020_2pt5_min_tif/gpw_v4_population_count_rev11_2020_2pt5_min.tif")
heat <- rast("/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Rasters/heat.tif")
flood <- rast("/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Floods/floodMapGL_rp10y/floodMapGL_rp10y.tif")

rast.list <- list(airQuality, water, biodiv, heat, flood, pop)
names(rast.list) <- c("airQuality", "water", "biodiv", "heat", "flood", "pop")

# ------ Aggregation/Disaggregation ------
# Methods of interpolation: bilinear for continuous data and NN for categorical
newres <- 0.1
airQuality <- aggregate(airQuality, fact=newres/res(airQuality), fun=mean)
biodiv <- aggregate(biodiv,newres/res(biodiv),fun=mean)
heat <- aggregate(heat, newres/res(heat), fun=mean)
flood <- aggregate(flood, newres/res(flood), fun=mean)
pop <- aggregate(pop,newres/res(pop),fun=sum)
water <- disagg(water,fact=5, method="bilinear")

# ------ Match extents ------
heat <- resample(heat,pop)
flood <- resample(flood,pop)
water <- resample(water, pop)
biodiv <- resample(biodiv, pop)
airQuality <- resample(airQuality,pop)

# ------ Write raster files ------
writeRaster(airQuality, filename="/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Rasters/airQuality.tif")
writeRaster(water, filename="/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Rasters/water.tif", overwrite=TRUE)
writeRaster(biodiv, filename="/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Rasters/biodiv.tif", overwrite=TRUE)
writeRaster(heat, filename="/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Rasters/heat.tif", overwrite=TRUE)
writeRaster(pop, filename="/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Rasters/pop.tif", overwrite=TRUE)
writeRaster(flood, filename= "/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Rasters/flood.tif", overwrite=TRUE)


