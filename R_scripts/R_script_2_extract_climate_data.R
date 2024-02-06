# ============================================================================ #
# R script 2: extracting climate data from tiff and ncdf4 files
# ============================================================================ #



# Article title: Landscape-Level Analysis of a Fig-Pollinator-Antagonist Community: Spatial and Temporal Variation in a Fig Wasp Community and its Response to Biotic and Abiotic Factors
# By Finn Piatscheck, Justin van Goor, Derek Houston and John Nason.



# The script allows to extract time- and site-specific weather/climate data from estimates provided by:
# CHIRPS (daily precipitation estimates): https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/tifs/p05/
# MODIS11A1 (daily land surface temperature estimates): https://lpdaac.usgs.gov/products/mod11a1v006/ (obtained with the tool: https://lpdaacsvc.cr.usgs.gov/appeears/task/point)
# TerraClimate (monthly climate estimates): http://thredds.northwestknowledge.net:8080/thredds/catalog/TERRACLIMATE_ALL/data/catalog.html
# Daymet (daily estimates): https://thredds.daac.ornl.gov/thredds/catalog/ornldaac/1328/catalog.html (obtained with package daymetr)

# The data for MODIS11A1 has been previously downloaded, as well as study site-nearby weather stations (obtained from NOAA: https://www.ncdc.noaa.gov/cdo-web/datatools/findstation)
# Many Mexican weather stations cover only periods until mid or end of 2012 (at the time this script was written), therefore, 2011 data were used to compare weather station measurements and weather estimates.



# ---------------------------------------------------------------------------- #
# Set directory and load packages
# ---------------------------------------------------------------------------- #



# R version --------------------

# R version used for these analyses: 4.3.1.


# Set working directory --------------------

# To replicate the following analyses, create a folder named "working_directory".
# It should have 4 sub-folders: "raw_data", "generated_data", "tables" and "figures".
# workdir <- "PathToWorkingDirectory/" # Replace here with the path to the directory in which you placed the raw data files.
setwd(workdir)


# Load packages --------------------

library(reshape2) # version 1.4.4
library(dplyr) # version 1.1.3
library(maptools) # version 1.1-8
library(daymetr) # version 1.7.1
library(climateR) # version 0.3.1.3
library(AOI) # version 0.3.0
library(sf) # version 10.0-14
library(ggplot2) # version 3.4.4
library(stringr) # version 1.5.1



# ---------------------------------------------------------------------------- #
# Get meteorological station weather data
# ---------------------------------------------------------------------------- #



# To compare global estimates to actual meteorological stations, we downloaded weather station data that were near our study site.
# https://www.ncdc.noaa.gov/cdo-web/datatools/findstation
# Stations were chosen if they were near study sites and have data to analyze.
# Most stations only have precipitation, minimum and maximum temperature data.
# Punta Prieta weather station was changed to another because it was containing too many missing data.
# Loreto DGE station will be moved slightly to the west because it is near the water.

# Extracting and formatting weather station near sites data
NOAA_weather_stations <- read.csv(paste(workdir,"raw_data/","NOAA_BC_weather_stations_near_9_sites_y2011_(PuntaPrieta_changed)_02.JUN.2020.csv", sep = ""))

# Formating 
NOAA_weather_stations$STATION <- NULL; NOAA_weather_stations$PRCP_ATTRIBUTES <- NULL; NOAA_weather_stations$TMAX_ATTRIBUTES <- NULL; NOAA_weather_stations$TMIN_ATTRIBUTES <- NULL
colnames(NOAA_weather_stations) <- c("weather_station","latitude","longitude","elevation","date","prcp","tmax","tmin")
NOAA_weather_stations$date <- as.Date(NOAA_weather_stations$date, format = c("%Y-%m-%d"))
tmpsplit <- colsplit(NOAA_weather_stations$weather_station, ", ", c('a', 'b'))
NOAA_weather_stations$weather_station <- tmpsplit$a; rm(tmpsplit)
NOAA_weather_stations$weather_station[NOAA_weather_stations$weather_station == "CHAPALA"]<-"Chapala"
NOAA_weather_stations$weather_station[NOAA_weather_stations$weather_station == "BAHIA DE LOS ANGELES"]<-"Bahia de los Angeles"
NOAA_weather_stations$weather_station[NOAA_weather_stations$weather_station == "SAN FRANCISCO DE LA SIERRA"]<-"San Francisco de la Sierra"
NOAA_weather_stations$weather_station[NOAA_weather_stations$weather_station == "SANTA AGUEDA"]<-"Santa Agueda"
NOAA_weather_stations$weather_station[NOAA_weather_stations$weather_station == "OJO DE AGUA"]<-"Ojo de Agua"
NOAA_weather_stations$weather_station[NOAA_weather_stations$weather_station == "LORETO DGE"]<-"Loreto DGE"
NOAA_weather_stations$weather_station[NOAA_weather_stations$weather_station == "SAN RAMON"]<-"San Ramon"
NOAA_weather_stations$weather_station[NOAA_weather_stations$weather_station == "LOS ROBLES"]<-"Los Robles"
NOAA_weather_stations$weather_station[NOAA_weather_stations$weather_station == "SAN BARTOLO"]<-"San Bartolo"
NOAA_weather_stations$weather_station <- factor(NOAA_weather_stations$weather_station, levels = c("Chapala","Bahia de los Angeles","San Francisco de la Sierra","Santa Agueda","Ojo de Agua","Loreto DGE","San Ramon","Los Robles","San Bartolo")) # Stations ordered from north (first) to south (last)

# Save the weather stations data table
write.csv(NOAA_weather_stations, paste(workdir,"generated_data/","NOAA_weather_stations_meteorological_data_2011-2012.csv",sep = ""), row.names=FALSE)

# Extract weather stations coordinates
weather_stations_coord <- NOAA_weather_stations %>% 
  group_by(weather_station)  %>% 
  summarise(latitude=unique(latitude),longitude=unique(longitude))
weather_stations_coord <- as.data.frame(weather_stations_coord)

# Fix LORETO DGE longitude
# The coordinates are from the NOAA weather station catalog but that station is too close to the water, often resulting in missing data when trying to extract estimates at that location
# The longitude needs to be moved slightly to the left
weather_stations_coord[6,3] <- -111.3500
weather_stations_coord

# Save formatted table
write.csv(weather_stations_coord, paste(workdir,"generated_data/","NOAA_weather_stations_coord.csv",sep = ""), row.names=FALSE)



# ---------------------------------------------------------------------------- #
# Import main study sites coordinates and nearby weather stations coordinates 
# ---------------------------------------------------------------------------- #



# Main study sites coordinates --------------------

Fpet_sites_coord <-  read.csv(paste(workdir,"raw_data/","Fpet_sites_coord.csv",sep=""))

# Fpet_sites_coord <- data.frame(latitude = c(29.26359723,28.29038826,27.56491059,27.0995915,26.35797318,25.913455,25.38127688,24.03565692,23.7377989),
#                                longitude = c(-114.0216665,-113.1110026,-113.0711841,-112.4968451,-111.8027891,-111.349716,-111.3151591,-110.1232369,-109.8303991),
#                                pop = c("158","172","112","113","95","179","201","96","70"))

# Formating
Fpet_sites_coord$pop <- factor(Fpet_sites_coord$pop, levels = c("158","172","112","113","95","179","201","96","70")) # Sites ordered from north (first) to south (last)


# Weather stations near main study site data --------------------

weather_stations_coord <-  read.csv(paste(workdir,"generated_data/","NOAA_weather_stations_coord.csv",sep=""))

# weather_stations_coord <- data.frame(weather_station = c("Chapala","Bahia de los Angeles","San Francisco de la Sierra","Santa Agueda","Ojo de Agua","Loreto DGE","San Ramon","Los Robles","San Bartolo"),
#                                      latitude = c(29.4833,28.6,27.5833,27.25,26.3167,26,25.2667,24.0333,23.7333),
#                                      longitude = c(-114.35,-113.55,-113.0167,-112.35,-111.9833,-111.35,-111.2833,-110.1167,-109.8333))

# LORETO DGE longitude was moved slightly to the east and changed to "-111.3500" to avoid being in the water and failing to be compared with other data sets.
weather_stations_coord$weather_station <- factor(weather_stations_coord$weather_station, levels = c("Chapala","Bahia de los Angeles","San Francisco de la Sierra","Santa Agueda","Ojo de Agua","Loreto DGE","San Ramon","Los Robles","San Bartolo")) # Stations ordered from north (first) to south (last)


# Quick spatial look at main sites and weather stations --------------------

data(wrld_simpl)
plot(wrld_simpl, xlim=c(-117,-107), ylim=c(22,31), axes=TRUE, col="light yellow") # Rough representation of Baja California
box()
points(Fpet_sites_coord$longitude, Fpet_sites_coord$latitude, col='red', cex=1)
points(weather_stations_coord$longitude, weather_stations_coord$latitude, col='orange', cex=0.5) 
# Main study sites are red points, weather stations are smaller orange points.
rm(wrld_simpl)



# ---------------------------------------------------------------------------- #
# Access daily weather and monthly climate summaries at main study sites
# ---------------------------------------------------------------------------- #


# ____________________________________________________________________________ #
# DAYMET V4 (with daymetr)


Daymet_main_sites_list <- download_daymet_batch(file_location = paste(workdir,"raw_data/","Fpet_sites_coord.csv",sep=""),
                                                  start = 2012,
                                                  end = 2014,
                                                  internal = TRUE)

# Format
Daymet_main_sites_df <- Daymet_main_sites_list
loop_df <- Daymet_main_sites_df[[1]]$data[FALSE,]
for (i in 1:9){
  Daymet_main_sites_df[[i]]$data$site <- Daymet_main_sites_df[[i]]$site[1]
  Daymet_main_sites_df[[i]]$data$latitude <- Daymet_main_sites_df[[i]]$latitude[1]
  Daymet_main_sites_df[[i]]$data$longitude <- Daymet_main_sites_df[[i]]$longitude[1]
  Daymet_main_sites_df[[i]]$data$altitude <- Daymet_main_sites_df[[i]]$altitude[1]
  loop_df <- rbind(loop_df,Daymet_main_sites_df[[i]]$data)
}
loop_df$date <- as.Date(paste(loop_df$year, loop_df$yday, sep = "-"), "%Y-%j")
Daymet_main_sites_df <- loop_df[,c(10:14,4:9)]
Daymet_main_sites_df <- Daymet_main_sites_df[,c(1:6,9:11)] # select variables of interest
colnames(Daymet_main_sites_df) <- c("site", "latitude","longitude","altitude","date","DM_prec","DM_tmax","DM_tmin","DM_vp")
head(Daymet_main_sites_df); rm(loop_df)
Daymet_Fpet <- Daymet_main_sites_df # Rename

write.csv(Daymet_Fpet, paste(workdir,"generated_data/","Daymet_Fpet_2012-2014.csv", sep=""), row.names = FALSE)
rm(Daymet_main_sites_df,Daymet_main_sites_list)


# ____________________________________________________________________________ #
# CHIRPS (with climateR)


# Create AOI objects for climateR functions
# Main study sites
Fpet_sites_coord_sf <- st_as_sf(x = Fpet_sites_coord, 
         coords = c("longitude", "latitude"),
         crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
Fpet_sites_coord_aoi <- aoi_get(Fpet_sites_coord_sf)


# Extract data
CHIRPS_main_sites_stack <- getCHIRPS(Fpet_sites_coord_aoi, startDate = "2012-01-01", endDate = "2014-12-31" )
CHIRPS_main_sites_list <- extract_sites(r = CHIRPS_main_sites_stack, pts = Fpet_sites_coord_sf, id = "pop")

# Format 
CHIRPS_main_sites_df <- as.data.frame(CHIRPS_main_sites_list)
colnames(CHIRPS_main_sites_df) <- c("date","158","172","112","113","95","179","201","96","70")
CHIRPS_Fpet <- melt(CHIRPS_main_sites_df, id = "date")
colnames(CHIRPS_Fpet) <- c("date","site","CHIRPS_prec")

# Save formatted table
write.csv(CHIRPS_Fpet, paste(workdir,"generated_data/","CHIRPS_Fpet_2012-2014.csv", sep=""), row.names = FALSE)
rm(CHIRPS_main_sites_stack,CHIRPS_main_sites_list,CHIRPS_main_sites_df)

  
# ____________________________________________________________________________ #
# TerraClimate (with climateR)



# # Create AOI object for climatR functions
# Fpet_sites_coord_sf <- st_as_sf(x = Fpet_sites_coord, 
#          coords = c("longitude", "latitude"),
#          crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# Fpet_sites_coord_aoi <- aoi_get(Fpet_sites_coord_sf)

# Extract data
TerraClimate_main_sites_stack <- getTerraClim(Fpet_sites_coord_aoi, startDate = "2012-01-01", endDate = "2014-12-31" )
TerraClimate_main_sites_list <- extract_sites(r = TerraClimate_main_sites_stack, pts = Fpet_sites_coord_sf, id = "pop")

# Format 
# We are only interested in:
# - ppt (Precipitation, monthly total)
# - tmax (Max Temperature)
# - tmin (Min Temperature)
# - vap (Vapor pressure)
# - ws (Wind speed)

TerraClimate_main_sites_list <- TerraClimate_main_sites_list[c("ppt","tmax","tmin","vap","ws")]
TerraClimate_main_sites_df <- TerraClimate_main_sites_list
# loop_df <- data.frame(date = "", variable ="", value ="")[FALSE,]
loop_df <- melt(as.data.frame(TerraClimate_main_sites_df[["ppt"]]), id = "date")
colnames(loop_df) <- c("date", "site", "ppt")
for (i in c("tmax","tmin","vap","ws")) {
  tmp_melt <- melt(as.data.frame(TerraClimate_main_sites_df[[i]]), id = "date")
  colnames(tmp_melt) <- c("date", "site", i)
  loop_df <- cbind(loop_df,tmp_melt)
}
# The resulting data frame has repeated dates and sites. We will remove these rows.
# They are identical, it was checked with the function setdiff().
TerraClimate_main_sites <- loop_df[,c(1:3,6,9,12,15)]; rm(loop_df); rm(tmp_melt)
# We are not interested in the time of recording.
TerraClimate_main_sites$date <- as.Date(do.call(rbind, strsplit(as.character(TerraClimate_main_sites$date)," "))[,1])
# Rename variables
colnames(TerraClimate_main_sites) <- c("date","site","TC_prec","TC_tmax","TC_tmin","TC_vp","TC_ws")
# TerraClimate data include January 2015 (even if we didn't ask for year 2015).
# We remove this last row, as it will complicate the merging of all the data.
TerraClimate_main_sites_no2015 <- subset(TerraClimate_main_sites, !date =="2015-01-01")
TC_Fpet <- TerraClimate_main_sites_no2015 # Rename

# Save formatted table
write.csv(TC_Fpet, paste(workdir,"generated_data/","TC_Fpet_2012-2014.csv", sep=""), row.names = FALSE)
rm(TerraClimate_main_sites_stack,TerraClimate_main_sites_list,TerraClimate_main_sites_df,TerraClimate_main_sites,TerraClimate_main_sites_no2015)


# ____________________________________________________________________________ #
# MODIS MODA11A1 v0-61 LST (obtained with AppEears)


# The data were downloaded previously with the tool AppEears
# AppEears also allows to download DAYMET v4 data
# Alternatively, check the R packages MODISTools or MODIStsp

# Main sites
MODIS_LST_Fpet <- read.csv(paste(workdir,"raw_data/","Fpet-sites-MOD11A1v0-61-MOD11A1-061-results.csv",sep = ""))
MODIS_LST_Fpet <- MODIS_LST_Fpet[,c(1:4,14)]
MODIS_LST_Fpet$Date <- as.Date(MODIS_LST_Fpet$Date)
MODIS_LST_Fpet$MOD11A1_061_LST_Day_1km[MODIS_LST_Fpet$MOD11A1_061_LST_Day_1km == 0] <- NA
MODIS_LST_Fpet$MOD11A1_061_LST_Day_1km <- MODIS_LST_Fpet$MOD11A1_061_LST_Day_1km  - 273.15
colnames(MODIS_LST_Fpet) <- c("site","latitude","longitude","date","MODIS_LST")

# Save formatted table
write.csv(MODIS_LST_Fpet, paste(workdir,"generated_data/","MODIS_LST_Fpet_2011-2014.csv", sep=""), row.names = FALSE)



# ---------------------------------------------------------------------------- #
# Access daily weather and monthly climate summaries at nearby meteorological stations
# ---------------------------------------------------------------------------- #



# NOAA Weather Stations found nearest main study sites in Baja California.
# Data source: https://www.ncdc.noaa.gov/cdo-web/datatools/findstation
# The stations were chosen by eye. Most stations have only data up to 2012.
# Therefor, we download data for 2011 and compare with remote-sensing estimates


# ____________________________________________________________________________ #
# DAYMET V4 (with daymetr)


Daymet_weather_stations_list <- download_daymet_batch(file_location = paste(workdir,"generated_data/","NOAA_weather_stations_coord.csv",sep=""),
                                                 start = 2011,
                                                 end = 2011,
                                                 internal = TRUE)
Daymet_weather_stations_df <- Daymet_weather_stations_list
loop_df <- Daymet_weather_stations_df[[1]]$data[FALSE,]
for (i in 1:9){
  Daymet_weather_stations_df[[i]]$data$weather_station <- Daymet_weather_stations_df[[i]]$site[1]
  Daymet_weather_stations_df[[i]]$data$latitude <- Daymet_weather_stations_df[[i]]$latitude[1]
  Daymet_weather_stations_df[[i]]$data$longitude <- Daymet_weather_stations_df[[i]]$longitude[1]
  Daymet_weather_stations_df[[i]]$data$altitude <- Daymet_weather_stations_df[[i]]$altitude[1]
  loop_df <- rbind(loop_df,Daymet_weather_stations_df[[i]]$data)
}
loop_df$date <- as.Date(paste(loop_df$year, loop_df$yday, sep = "-"), "%Y-%j")
Daymet_weather_stations_df <- loop_df[,c(10:14,4:9)]; rm(loop_df)
colnames(Daymet_weather_stations_df) <- c("weather_station", "latitude","longitude","altitude","date","DM_prec","DM_rad","DM_ws","DM_tmax","DM_tmin","DM_vp")
Daymet_weather_stations <- Daymet_weather_stations_df # Rename

# Save formatted table
write.csv(Daymet_weather_stations, paste(workdir,"generated_data/","Daymet_weather_stations_2011.csv", sep=""), row.names = FALSE)
rm(Daymet_weather_stations_list,Daymet_weather_stations_df)

# ____________________________________________________________________________ #
# CHIRPS (with climateR)


# Weather stations
weather_stations_coord_sf <- st_as_sf(x = weather_stations_coord, 
                                      coords = c("longitude", "latitude"),
                                      crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
weather_stations_coord_aoi <- aoi_get(weather_stations_coord_sf)

CHIRPS_weather_stations_stack <- getCHIRPS(weather_stations_coord_aoi, startDate = "2011-01-01", endDate = "2011-12-31" )
CHIRPS_weather_stations_list <- extract_sites(r = CHIRPS_weather_stations_stack, pts = weather_stations_coord_sf, id = "weather_station")

# Format 
CHIRPS_weather_stations_df <- as.data.frame(CHIRPS_weather_stations_list)
colnames(CHIRPS_weather_stations_df) <- c("date","Chapala","Bahia de los Angeles","San Francisco de la Sierra","Santa Agueda","Ojo de Agua","Loreto DGE","San Ramon","Los Robles","San Bartolo")
CHIRPS_weather_stations <- melt(CHIRPS_weather_stations_df, id = "date")
colnames(CHIRPS_weather_stations) <- c("date","weather_station","CHIRPS_prec")

# Save formatted table
write.csv(CHIRPS_weather_stations, paste(workdir,"generated_data/","CHIRPS_weather_stations_2011.csv", sep=""), row.names = FALSE)
rm(CHIRPS_weather_stations_stack,CHIRPS_weather_stations_list,CHIRPS_weather_stations_df)



# ____________________________________________________________________________ #
# TerraClimate (with climateR)


# # Create AOI object for climatR functions

# # Weather stations
# weather_stations_coord_sf <- st_as_sf(x = weather_stations_coord, 
#                                       coords = c("longitude", "latitude"),
#                                       crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# weather_stations_coord_aoi <- aoi_get(weather_stations_coord_sf)


TerraClimate_weather_stations_stack <- getTerraClim(weather_stations_coord_aoi, startDate = "2011-01-01", endDate = "2011-12-31" )
TerraClimate_weather_stations_list <- extract_sites(r = TerraClimate_weather_stations_stack, pts = weather_stations_coord_sf, id = "weather_station")

# Format 
# We are only interested in:
# - ppt (Precipitation, monthly total)
# - tmax (Max Temperature)
# - tmin (Min Temperature)
# - vap (Vapor pressure)
# - ws (Wind speed)

TerraClimate_weather_stations_list <- TerraClimate_weather_stations_list[c("ppt","tmax","tmin","vap","ws")]
TerraClimate_weather_stations_df <- TerraClimate_weather_stations_list
# loop_df <- data.frame(date = "", variable ="", value ="")[FALSE,]
loop_df <- melt(as.data.frame(TerraClimate_weather_stations_df[["ppt"]]), id = "date")
colnames(loop_df) <- c("date", "site", "ppt")
for (i in c("tmax","tmin","vap","ws")) {
  tmp_melt <- melt(as.data.frame(TerraClimate_weather_stations_df[[i]]), id = "date")
  colnames(tmp_melt) <- c("date", "site", i)
  loop_df <- cbind(loop_df,tmp_melt)
}
# The resulting data frame has repeated dates and sites. We will remove these rows.
# They are identical, it was checked with the function setdiff().
TerraClimate_weather_stations <- loop_df[,c(1:3,6,9,12,15)]; rm(loop_df); rm(tmp_melt)
# We are not interested in the time of recording.
TerraClimate_weather_stations$date <- as.Date(do.call(rbind, strsplit(as.character(TerraClimate_weather_stations$date)," "))[,1])
# Rename variables
colnames(TerraClimate_weather_stations) <- c("date","weather_station","TC_prec","TC_tmax","TC_tmin","TC_vp","TC_ws")
# Rename stations
levels(TerraClimate_weather_stations$weather_station) <- c("Chapala","Bahia de los Angeles","San Francisco de la Sierra","Santa Agueda","Ojo de Agua","Loreto DGE","San Ramon","Los Robles","San Bartolo")
# TerraClimate data include January 2012 even if we didn't ask for year 2012).
# We remove this last row, as it will complicate the merging of all the data.
TerraClimate_weather_stations_no2012 <- subset(TerraClimate_weather_stations, !date =="2012-01-01")
TC_weather_stations <- TerraClimate_weather_stations_no2012 # Rename

# Save formatted table
write.csv(TC_weather_stations, paste(workdir,"generated_data/","TC_weather_stations_2011.csv", sep=""), row.names = FALSE)
rm(TerraClimate_weather_stations,TerraClimate_weather_stations_list,TerraClimate_weather_stations_df,TerraClimate_weather_stations_stack,TerraClimate_weather_stations_no2012)


# ____________________________________________________________________________ #
# MODIS MODA11A1 v0-61 LST (obtained with AppEears)


# Weather stations
MODIS_LST_weather_stations <- read.csv(paste(workdir,"raw_data/","weather-stations-2011-MOD11A1v0-61-MOD11A1-061-results.csv",sep = ""))
MODIS_LST_weather_stations <- MODIS_LST_weather_stations[,c(1:4,14)]
MODIS_LST_weather_stations$Date <- as.Date(MODIS_LST_weather_stations$Date)
MODIS_LST_weather_stations$MOD11A1_061_LST_Day_1km[MODIS_LST_weather_stations$MOD11A1_061_LST_Day_1km == 0] <- NA
MODIS_LST_weather_stations$MOD11A1_061_LST_Day_1km <- MODIS_LST_weather_stations$MOD11A1_061_LST_Day_1km  - 273.15 # Change to Celcius
colnames(MODIS_LST_weather_stations) <- c("weather_station","latitude","longitude","date","MODIS_LST")

# Save formatted table
write.csv(MODIS_LST_weather_stations, paste(workdir,"generated_data/","MODIS_LST_weather_stations_2011.csv",sep=""), row.names = FALSE)



# ---------------------------------------------------------------------------- #
# Combine into unique data sets for study sites and weather stations
# ---------------------------------------------------------------------------- #



# Study sites, 2012-2014 --------------------

# Format and merge
Daymet_Fpet$date <- as.Date(Daymet_Fpet$date)
CHIRPS_Fpet$date <- as.Date(CHIRPS_Fpet$date)
MODIS_LST_Fpet$date <- as.Date(MODIS_LST_Fpet$date)
TC_Fpet$date <- as.Date(TC_Fpet$date)


# Merge data
Fpet_all_estimates <-  merge(Daymet_Fpet, 
                             merge(CHIRPS_Fpet,
                                   merge(MODIS_LST_Fpet[,c(1,4:5)],TC_Fpet, by = c("site", "date"), all = TRUE), by = c("site", "date"), all = TRUE), by = c("site", "date"), all = TRUE)
Fpet_all_estimates$site <- factor(Fpet_all_estimates$site, levels=c("158","172","112","113","95","179","201","96","70"))
head(Fpet_all_estimates); dim(Fpet_all_estimates) # number of rows should be 9864


# Save
write.csv(Fpet_all_estimates, paste(workdir,"generated_data/","Fpet_all_estimates_2012-2014.csv", sep=""), row.names = FALSE)


# Weather stations 2011 --------------------

# Many weather stations contain data downloaded for some or all days of 2012.
# We remove these rows.
NOAA_weather_stations <- subset(NOAA_weather_stations, format(as.Date(date),"%Y")==2011)

# Format and merge
NOAA_weather_stations$date <- as.Date(NOAA_weather_stations$date)
Daymet_weather_stations$date <- as.Date(Daymet_weather_stations$date)
CHIRPS_weather_stations$date <- as.Date(CHIRPS_weather_stations$date)
MODIS_LST_weather_stations$date <- as.Date(MODIS_LST_weather_stations$date)
TC_weather_stations$date <- as.Date(TC_weather_stations$date)

weather_station_all_estimates <-  merge(NOAA_weather_stations, 
                                        merge(Daymet_weather_stations[,c(1,5,6,9,10)],
                                              merge(CHIRPS_weather_stations,
                                                    merge(MODIS_LST_weather_stations[,c(1,4:5)], TC_weather_stations[,c(1:5)], 
                                                                                 by = c("weather_station", "date"), all = TRUE), by = c("weather_station", "date"), all = TRUE), by = c("weather_station", "date"), all = TRUE), by = c("weather_station", "date"), all = TRUE)
head(weather_station_all_estimates); dim(weather_station_all_estimates) # number of rows should be 3285
# Save
write.csv(weather_station_all_estimates, paste(workdir,"generated_data/","weather_station_all_estimates_2011.csv", sep=""), row.names = FALSE)



# ---------------------------------------------------------------------------- #
# Get annual temperature and precipitation for a transect in Baja California 
# ---------------------------------------------------------------------------- #



# ____________________________________________________________________________ #
# Create transect of coordinates through study sites


# We create a series of coordinates across Baja California which passes through our study sites.
# The following will generate coordinates values for a "transect" which will pass by each sites in Baja (detour around La Paz).
# Here I calculated the distance (with an online tool) between all sites, summed the total distance, divided each value by the total, multiplied by 350 and rounded which gave me the number of points between each sites.
# We create a loop around La Paz to avoid crossing the water and be able to get annual summaries.
latitutes_values1 <- seq(23.73780, 24.03566, length.out=18)
latitutes_values2 <- seq(24.03566, 24.1, length.out=18)
latitutes_values2l1 <- seq(24.1, 24.5, length.out=23)          # Loop around La Paz
latitutes_values2l2 <- seq(24.5, 25.38128, length.out=44)      # Loop around La Paz
latitutes_values3 <- seq(25.38128, 25.913455, length.out=22)
latitutes_values4 <- seq(25.913455, 26.35797, length.out=30)
latitutes_values5 <- seq(26.35797, 27.09959, length.out=74)
latitutes_values6 <- seq(27.09959, 27.56491, length.out=31)
latitutes_values7 <- seq(27.56491, 28.29039, length.out=33)
latitutes_values8 <- seq(28.29039, 29.26360, length.out=57)
latitutes_values <- c(latitutes_values1, latitutes_values2, latitutes_values2l1, latitutes_values2l2, latitutes_values3, latitutes_values4, latitutes_values5, latitutes_values6, latitutes_values7, latitutes_values8)

longitude_values1 <- seq(-109.8304, -110.1232, length.out=18)
longitude_values2 <- seq(-110.1232, -110.55, length.out=18)
longitude_values2l1 <- seq(-110.55, -110.9, length.out=23)
longitude_values2l2 <- seq(-110.9, -111.3151591, length.out=44)
longitude_values3 <- seq(-111.3151591, -111.349716, length.out=22)
longitude_values4 <- seq(-111.349716, -111.8027891, length.out=30)
longitude_values5 <- seq(-111.8027891, -112.4968451, length.out=74)
longitude_values6 <- seq(-112.4968451, -113.0712, length.out=31)
longitude_values7 <- seq(-113.0712, -113.1110, length.out=33)
longitude_values8 <- seq(-113.1110, -114.0217, length.out=57)
longitude_values <- c(longitude_values1, longitude_values2, longitude_values2l1, longitude_values2l2, longitude_values3, longitude_values4, longitude_values5, longitude_values6, longitude_values7, longitude_values8)

virtualoop <- seq(301, 650, length.out=350) # Fake pop names, starting at 300 to not overlap with actual site names.

# Import main study sites coordinates (if not already loaded in environment).
Fpet_sites_coord <- read.csv(paste(workdir,"raw_data/","Fpet_sites_coord.csv",sep=""))
# Formatting.
Fpet_sites_coord$pop <- factor(Fpet_sites_coord$pop, levels = c("158","172","112","113","95","179","201","96","70")) # Sites ordered from north (first) to south (last)

# Plot a simple map and visualize 
data(wrld_simpl)
plot(wrld_simpl, xlim=c(-117,-107), ylim=c(22,31), axes=TRUE, col="light yellow")
box()
points(longitude_values, latitutes_values, col='orange', cex=0.5, lwd=2)
points(Fpet_sites_coord$longitude, Fpet_sites_coord$latitude, col='red', cex=1,lwd=3)
# We can see the line passing trough the study sites.

# Create transect coordinates data frame
transect <- as.data.frame(cbind(latitutes_values, longitude_values, virtualoop))
colnames(transect) <- c("latitude", "longitude", "pop")
transect <- transect[,c(3,1,2)]
transect$pop <- as.factor(transect$pop)
transect <- rbind(Fpet_sites_coord, transect)                         


# ____________________________________________________________________________ #
# Obtain annual temperature and precipitation summaries for all these points 


# Create AOI object for climatR functions
transect_sf <- st_as_sf(x = transect,
         coords = c("longitude", "latitude"),
         crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
transect_aoi <- aoi_get(transect_sf)

# Extract monthly ClimateTerra data
transect_stack_2012 <- getTerraClim(transect_aoi, startDate = "2012-01-01", endDate = "2012-12-01" )
transect_stack_2012_prec_annual <- sum(transect_stack_2012$ppt)
transect_stack_2012_tmax_annual <- max(transect_stack_2012$tmax)
transect_2012_prec_annual <- extract_sites(r = transect_stack_2012_prec_annual, pts = transect_sf, id = "pop")
transect_2012_tmax_annual <- extract_sites(r = transect_stack_2012_tmax_annual, pts = transect_sf, id = "pop")

transect_stack_2013 <- getTerraClim(transect_aoi, startDate = "2013-01-01", endDate = "2013-12-01" )
transect_stack_2013_prec_annual <- sum(transect_stack_2013$ppt)
transect_stack_2013_tmax_annual <- max(transect_stack_2013$tmax)
transect_2013_prec_annual <- extract_sites(r = transect_stack_2013_prec_annual, pts = transect_sf, id = "pop")
transect_2013_tmax_annual <- extract_sites(r = transect_stack_2013_tmax_annual, pts = transect_sf, id = "pop")

transect_stack_2014 <- getTerraClim(transect_aoi, startDate = "2014-01-01", endDate = "2014-12-01" )
transect_stack_2014_prec_annual <- sum(transect_stack_2014$ppt)
transect_stack_2014_tmax_annual <- max(transect_stack_2014$tmax)
transect_2014_prec_annual <- extract_sites(r = transect_stack_2014_prec_annual, pts = transect_sf, id = "pop")
transect_2014_tmax_annual <- extract_sites(r = transect_stack_2014_tmax_annual, pts = transect_sf, id = "pop")

# Create a data frame with all coordinates and annual summaries
transect_2012_tmax_annual_t <- as.data.frame(t(transect_2012_tmax_annual[[1]]))
transect_2012_tmax_annual_t$pop <- rownames(transect_2012_tmax_annual_t)
transect_2012_tmax_annual_t <- transect_2012_tmax_annual_t[-1,]
colnames(transect_2012_tmax_annual_t) <- c("TC_tmax_2012","pop")
rownames(transect_2012_tmax_annual_t) <- NULL

transect_2013_tmax_annual_t <- as.data.frame(t(transect_2013_tmax_annual[[1]]))
transect_2013_tmax_annual_t$pop <- rownames(transect_2013_tmax_annual_t)
transect_2013_tmax_annual_t <- transect_2013_tmax_annual_t[-1,]
colnames(transect_2013_tmax_annual_t) <- c("TC_tmax_2013","pop")
rownames(transect_2013_tmax_annual_t) <- NULL

transect_2014_tmax_annual_t <- as.data.frame(t(transect_2014_tmax_annual[[1]]))
transect_2014_tmax_annual_t$pop <- rownames(transect_2014_tmax_annual_t)
transect_2014_tmax_annual_t <- transect_2014_tmax_annual_t[-1,]
colnames(transect_2014_tmax_annual_t) <- c("TC_tmax_2014","pop")
rownames(transect_2014_tmax_annual_t) <- NULL

transect_2012_prec_annual_t <- as.data.frame(t(transect_2012_prec_annual[[1]]))
transect_2012_prec_annual_t$pop <- rownames(transect_2012_prec_annual_t)
transect_2012_prec_annual_t <- transect_2012_prec_annual_t[-1,]
colnames(transect_2012_prec_annual_t) <- c("TC_prec_2012","pop")
rownames(transect_2012_prec_annual_t) <- NULL

transect_2013_prec_annual_t <- as.data.frame(t(transect_2013_prec_annual[[1]]))
transect_2013_prec_annual_t$pop <- rownames(transect_2013_prec_annual_t)
transect_2013_prec_annual_t <- transect_2013_prec_annual_t[-1,]
colnames(transect_2013_prec_annual_t) <- c("TC_prec_2013","pop")
rownames(transect_2013_prec_annual_t) <- NULL

transect_2014_prec_annual_t <- as.data.frame(t(transect_2014_prec_annual[[1]]))
transect_2014_prec_annual_t$pop <- rownames(transect_2014_prec_annual_t)
transect_2014_prec_annual_t <- transect_2014_prec_annual_t[-1,]
colnames(transect_2014_prec_annual_t) <- c("TC_prec_2014","pop")
rownames(transect_2014_prec_annual_t) <- NULL


# Combine all
transect$id_merge <- 1:nrow(transect)
transect_annual_data <- merge(transect, merge(transect_2012_tmax_annual_t, 
                                              merge(transect_2013_tmax_annual_t,transect_2014_tmax_annual_t, 
                                              by = "pop", all = TRUE), by = "pop", all = TRUE), by = "pop", all = TRUE)
transect_annual_data <- merge(transect_annual_data, merge(transect_2012_prec_annual_t, 
                                              merge(transect_2013_prec_annual_t,transect_2014_prec_annual_t, 
                                                    by = "pop", all = TRUE), by = "pop", all = TRUE), by = "pop", all = TRUE)
transect_annual_data <- transect_annual_data[order(transect_annual_data$id_merge),]                    
transect_annual_data$id_merge <- NULL
rownames(transect_annual_data) <- NULL


write.csv(transect_annual_data, paste(workdir,"generated_data/","transect_annual_data_2012-2014.csv", sep=""), row.names = FALSE)
rm(transect_stack_2012,transect_stack_2013,transect_stack_2014,
   transect_stack_2012_prec_annual,transect_stack_2013_prec_annual,transect_stack_2014_prec_annual,
   transect_stack_2012_tmax_annual,transect_stack_2013_tmax_annualtransect_stack_2014_tmax_annual,
   transect_2012_prec_annual,transect_2013_prec_annual,transect_2014_prec_annual,
   transect_2012_tmax_annual,transect_2013_tmax_annual,transect_2014_tmax_annual)



# ---------------------------------------------------------------------------- #
# End of R script 2
# ---------------------------------------------------------------------------- #



# Clear working environment.
rm(list = ls())
dev.off()
