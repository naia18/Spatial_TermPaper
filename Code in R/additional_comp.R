# ========================================================================
# additional_comp.R    -   Naia Ormaza Zulueta   -  Apr 2022
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
packages <- c("terra", "raster", "tidyverse", "lattice", "foreign", "ncdf4",
              "dplyr", "stringr")
lapply(packages, require, character=TRUE)

# Load main dataframe
load(file="/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/R/Data/stack_layers.Rdata")
load(file="/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/R/Data/stack_layers_two.Rdata")

hrv.0 <- df[df$HRV==0,]
hrv.1 <- df[df$HRV==1,]
hrv.2 <- df[df$HRV==2,]
hrv.3 <- df[df$HRV==3,]
hrv.4 <- df[df$HRV==4,]
hrv.5 <- df[df$HRV==5,]

pop.0 <- sum(hrv.0$gpw_v4_population_count_rev11_2020_2pt5_min)
pop.1 <- sum(hrv.1$gpw_v4_population_count_rev11_2020_2pt5_min)
pop.2 <- sum(hrv.2$gpw_v4_population_count_rev11_2020_2pt5_min)
pop.3 <- sum(hrv.3$gpw_v4_population_count_rev11_2020_2pt5_min)
pop.4 <- sum(hrv.4$gpw_v4_population_count_rev11_2020_2pt5_min)
pop.5 <- sum(hrv.5$gpw_v4_population_count_rev11_2020_2pt5_min)
pop_total <- pop.0 + pop.1 + pop.2 + pop.3 + pop.4 + pop.5
pop_2 <- pop.2 + pop.3 + pop.4 + pop.5
pop_3 <- pop.3 + pop.4 + pop.5
pop_4 <- pop.4 + pop.5

# SOUTH AND SOUTHEAST ASIA
df_sea <- df[(df$y<45)&(df$y>-10)&(df$x>45)&(df$x<140),]  # Delimit extent
hrv_sea.0 <- df_sea[df_sea$HRV==0,]
hrv_sea.1 <- df_sea[df_sea$HRV==1,]
hrv_sea.2 <- df_sea[df_sea$HRV==2,]
hrv_sea.3 <- df_sea[df_sea$HRV==3,]
hrv_sea.4 <- df_sea[df_sea$HRV==4,]
hrv_sea.5 <- df_sea[df_sea$HRV==5,]

pop_sea.0 <- sum(hrv_sea.0$gpw_v4_population_count_rev11_2020_2pt5_min)
pop_sea.1 <- sum(hrv_sea.1$gpw_v4_population_count_rev11_2020_2pt5_min)
pop_sea.2 <- sum(hrv_sea.2$gpw_v4_population_count_rev11_2020_2pt5_min)
pop_sea.3 <- sum(hrv_sea.3$gpw_v4_population_count_rev11_2020_2pt5_min)
pop_sea.4 <- sum(hrv_sea.4$gpw_v4_population_count_rev11_2020_2pt5_min)
pop_sea.5 <- sum(hrv_sea.5$gpw_v4_population_count_rev11_2020_2pt5_min)
pop_total_sea <- pop_sea.0 + pop_sea.1 + pop_sea.2 + pop_sea.3 + pop_sea.4 + pop_sea.5
pop_2_sea <- pop_sea.2 + pop_sea.3 + pop_sea.4 + pop_sea.5
pop_3_sea <- pop_sea.3 + pop_sea.4 + pop_sea.5
pop_4_sea <- pop_sea.4 + pop_sea.5

# SSA
df_ssa <- df[(df$y<25)&(df$y>-35)&(df$x>-30)&(df$x<60),]  # Delimit extent
hrv_ssa.0 <- df_ssa[df_ssa$HRV==0,]
hrv_ssa.1 <- df_ssa[df_ssa$HRV==1,]
hrv_ssa.2 <- df_ssa[df_ssa$HRV==2,]
hrv_ssa.3 <- df_ssa[df_ssa$HRV==3,]
hrv_ssa.4 <- df_ssa[df_ssa$HRV==4,]
hrv_ssa.5 <- df_ssa[df_ssa$HRV==5,]

pop_ssa.0 <- sum(hrv_ssa.0$gpw_v4_population_count_rev11_2020_2pt5_min)
pop_ssa.1 <- sum(hrv_ssa.1$gpw_v4_population_count_rev11_2020_2pt5_min)
pop_ssa.2 <- sum(hrv_ssa.2$gpw_v4_population_count_rev11_2020_2pt5_min)
pop_ssa.3 <- sum(hrv_ssa.3$gpw_v4_population_count_rev11_2020_2pt5_min)
pop_ssa.4 <- sum(hrv_ssa.4$gpw_v4_population_count_rev11_2020_2pt5_min)
pop_ssa.5 <- sum(hrv_ssa.5$gpw_v4_population_count_rev11_2020_2pt5_min)
pop_total_ssa <- pop_ssa.0 + pop_ssa.1 + pop_ssa.2 + pop_ssa.3 + pop_ssa.4 + pop_ssa.5
pop_2_ssa <- pop_ssa.2 + pop_ssa.3 + pop_ssa.4 + pop_ssa.5
pop_3_ssa <- pop_ssa.3 + pop_ssa.4 + pop_ssa.5
pop_4_ssa <- pop_ssa.4 + pop_ssa.5



# Summary statistics
prop2 <- (pop_2_sea+pop_2_ssa)/pop_2
prop3 <- (pop_3_sea+pop_3_ssa)/pop_3
prop4 <- (pop_4_sea+pop_4_ssa)/pop_4

# Second approach computations
hrv_two.0 <- df_two[df_two$HRV==0,]
hrv_two.1 <- df_two[df_two$HRV==1,]
hrv_two.2 <- df_two[df_two$HRV==2,]
hrv_two.3 <- df_two[df_two$HRV==3,]
hrv_two.4 <- df_two[df_two$HRV==4,]
hrv_two.5 <- df_two[df_two$HRV==5,]

pop_two.0 <- sum(hrv_two.0$gpw_v4_population_count_rev11_2020_2pt5_min)
pop_two.1 <- sum(hrv_two.1$gpw_v4_population_count_rev11_2020_2pt5_min)
pop_two.2 <- sum(hrv_two.2$gpw_v4_population_count_rev11_2020_2pt5_min)
pop_two.3 <- sum(hrv_two.3$gpw_v4_population_count_rev11_2020_2pt5_min)
pop_two.4 <- sum(hrv_two.4$gpw_v4_population_count_rev11_2020_2pt5_min)
pop_two.5 <- sum(hrv_two.5$gpw_v4_population_count_rev11_2020_2pt5_min)
pop_total_two <- pop_two.0 + pop_two.1 + pop_two.2 + pop_two.3 + pop_two.4 + pop_two.5
pop_2_two <- pop_two.2 + pop_two.3 + pop_two.4 + pop_two.5



# open a netCDF file
ncin <- nc_open("/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Human Development Index/doi_10.5061_dryad.dk1j0__v2/GDP_PPP_30arcsec_v3.nc")
print(ncin)

# get longitude and latitude
lon <- ncvar_get(ncin,"longitude")
nlon <- dim(lon)
head(lon)

lat <- ncvar_get(ncin,"latitude")
nlat <- dim(lat)
head(lat)

print(c(nlon,nlat))

# get the PM25
dname <- "GDP_PPP"
hdi_array <- ncvar_get(ncin,dname)
dlname <- ncatt_get(ncin,dname,"long_name")
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")
dim(hdi_array)

# Take last slice of HDI
hdi_slice <- hdi_array[,,26]

