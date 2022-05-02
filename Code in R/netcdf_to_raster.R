# ========================================================================
# ncdf_to_raster.R    -   Naia Ormaza Zulueta   -  Apr 2022
# This file converts air quality pm2.5 data in netcdf format to a raster
# file GTiff format
# ========================================================================

library(raster)
library(rasterVis)
library(ncdf4)
library(lattice)

ncpath <- "/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Air Quality/Global 2/Annual/"
ncname <- "V5GL02.HybridPM25c_0p10.Global.202001-202012"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")

ncfile = ncdf4::nc_open(ncfname)
names(ncfile$var)


# set input path
input_nc <-  ncfname
varname <- 'GWRPM25'
nc2raster <- raster(input_nc,varname = varname,band = 1)

# quick view for the dataset
png("F:\\plot2019.png",
    height = 15,
    width = 20,
    units = 'cm',
    res = 1000)
print(levelplot(nc2raster))
dev.off()

nc2raster <- stack(input_nc,varname = varname)

outpath <- "/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Air Quality/Raster/"
outname <- "airQuality"
output <- paste(outpath, outname, ".tif", sep="")

# write raster to file
writeRaster(nc2raster,output,format = 'GTiff',overwrite = TRUE)

