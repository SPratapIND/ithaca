# Review era5-land data ----
source('source/partition_evap.R')
source("source/cdo_functions.R")

library(pRecipe)

## geospatial 
library(raster)
library(ncdf4)
library(sp)
library(sf)
library(stars)


## Load data ----
fnames_nc <- list.files(path = PATH_EVAP_SIM, pattern = ".nc$", 
                        full.names = T)

fnames_product <- fnames_nc[grep("era5-land", fnames_nc)]

fname <- fnames_product[grep("yearly", fnames_product)]
## 1. Check name: Data product, variables, unit, scale, beginning of time, end of time, spatial resolution, temporal resolution ----
fname

## 2. Check unit conversion/scaling (nc history)----
# multipliers "-mulc"
# muldpm
cmd_ncdump <- paste0("ncdump -h ", fname)
system(cmd_ncdump)

# variable name
cmd_cdo_variable_name <- paste0("cdo vardes ", fname)
system(cmd_cdo_variable_name)

# summary
cdo_info_fnc(fname)

## 3. plot map to see if it looks okay ----
data <- brick(fname) 
pRecipe::plot_map(data[[12]])

## 4. Check projection ----
# projection = lonlat
# lon: from east to west and -180 to 180 (1440)
# lat:from north to south and 90 to -90 (720)
crs(data)
extent(data)

# OG grid location for comparison
cmd_cdo_griddes <- paste0("cdo -griddes ", fname)
system(cmd_cdo_griddes)

grid_location <- "/mnt/shared/data_projects/ithaca/misc/grid_ithaca"
cmd_cat_grid <- paste0("cat ", grid_location)
system(cmd_cat_grid)

## 5. Check time units ----
# daily: hours should always be 00:00:00
# monthly: day should always be 01
# yearly: month and day should be 01-01
cdo_sinfo_fnc(fname)
data@z$Date

## 6. Check data type: Data type F32z ----
# scroll to top of output
cdo_sinfo_fnc(fname)

## 7. Compare values to at least one other publication ----

cdo_fldmean_fnc(inputfile_name = fname, outputfile_name = "~/Review/era5land_review_sel_fldmean.nc")
cdo_info_fnc(inputfile_name = "~/Review/era5land_review_sel_fldmean.nc")

cdo_timmean_fnc(inputfile_name = fname, outputfile_name = "~/Review/era5land_review_sel_time_mean.nc")
cdo_info_fnc(inputfile_name = "~/Review/era5land_review_sel_time_mean.nc")
