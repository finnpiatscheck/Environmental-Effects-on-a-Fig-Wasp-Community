# ============================================================================ #
# R script 3: reproduction of climate data analysis
# ============================================================================ #



# Article title: Landscape-Level Analysis of a Fig-Pollinator-Antagonist Community: Spatial and Temporal Variation in a Fig Wasp Community and its Response to Biotic and Abiotic Factors
# By Finn Piatscheck, Justin van Goor, Derek Houston and John Nason.



# The following script allows to download, extract and analyse the daily weather estimates and monthly climate summaries gathered with script "extract_climate_data.R".
# The goal here is to visualize, compare, plot and compute correlation among these variables.
# The script allows also the generation of weather summaries prior to collecting period necessary for downstream analyses.



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

library(ggplot2) # version 3.4.4
library(scales) # version 1.3.0
library(reshape2) # version 1.4.4
library(corrplot) # version 0.92
library(dplyr) # version 1.3.3
library(Hmisc) # version 5.1-1
# library(data.table) # version 1.14.8 
# data.table will be loaded later.



# ---------------------------------------------------------------------------- #
# Import climate data (extracted and formatted in extract_climate_data.R)
# ---------------------------------------------------------------------------- #



# The data have been already extracted with "R_script_X1_extract_climate_data.R"


# Main study sites and weather stations coordinates --------------------

Fpet_sites_coord <- read.csv(paste(workdir,"raw_data/","Fpet_sites_coord.csv", sep = ""))
weather_stations_coord <- read.csv(paste(workdir,"generated_data/","NOAA_weather_stations_coord.csv",sep=""))


# Weather estimates/climate summaries merged data sets --------------------

# Study sites (2012-2014)
Fpet_all_estimates <-  read.csv(paste(workdir,"generated_data/","Fpet_all_estimates_2012-2014.csv",sep = ""))
Fpet_all_estimates$TC_vp <- Fpet_all_estimates$TC_vp * 1000 # TerraClimate vapor pressure data are in kPa, we change it to Pa.

# Weather station data (2011)
weather_station_all_estimates <- read.csv(paste(workdir,"generated_data/","weather_station_all_estimates_2011.csv",sep = ""))



# ---------------------------------------------------------------------------- #
# Weather stations and weather estimates comparison for the year 2011
# ---------------------------------------------------------------------------- #



# Format data --------------------
weather_station_tmax_melt <- melt(weather_station_all_estimates, id.vars = c("weather_station","date"), measure.vars = , c("tmax","MODIS_LST","DM_tmax","TC_tmax")) # Why is there a "," after measure.vars ???? double check later
weather_station_tmin_melt <- melt(weather_station_all_estimates, id.vars = c("weather_station","date"), measure.vars = , c("tmin","DM_tmin","TC_tmin"))
weather_station_prec_melt <- melt(weather_station_all_estimates, id.vars = c("weather_station","date"), measure.vars = , c("prcp","CHIRPS_prec","DM_prec","TC_prec"))
weather_station_tmax_melt$date <- as.Date(weather_station_tmax_melt$date)
weather_station_tmin_melt$date <- as.Date(weather_station_tmin_melt$date)
weather_station_prec_melt$date <- as.Date(weather_station_prec_melt$date)

# Plots --------------------

# Maximum temperature
ggplot(data = weather_station_tmax_melt) + 
  geom_line(aes(x = date, colour=variable, y = value)) +
  geom_line(data = na.omit(subset(weather_station_tmax_melt, variable=="TC_tmax")), aes(x = date, colour=variable, y = value)) +
  facet_wrap(~ weather_station, ncol=3)

# Minimum temperature
ggplot(data = weather_station_tmin_melt) + 
  geom_line(aes(x = date, colour=variable, y = value)) +
  geom_line(data = na.omit(subset(weather_station_tmin_melt, variable=="TC_tmin")), aes(x = date, colour=variable, y = value)) +
  facet_wrap(~ weather_station, ncol=3)

# Precipitation
ggplot(data = weather_station_prec_melt) + 
  geom_line(aes(x = date, colour=variable, y = value)) +
  geom_line(data = na.omit(subset(weather_station_prec_melt, variable=="TC_prec")), aes(x = date, colour=variable, y = value)) +
  facet_wrap(~ weather_station, ncol=3)


# ---------------------------------------------------------------------------- #
# Daily weather and monthly climate summaries plot for 2012 and 2014
# ---------------------------------------------------------------------------- #



# For the following we do not consider MODIS's LST variable.


# library(data.table) # version 1.14.6

# Colors for correlation plot
corrPalette <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", 
                                  "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                  "#4393C3", "#2166AC", "#053061"))


# ____________________________________________________________________________ #
# Compute correlation and plot correlation heatmap


# Daily estimates
daily_abiotic_variable_corr <- cor(na.omit(Fpet_all_estimates[,6:10], method="pearson"))
colnames(daily_abiotic_variable_corr) <- c("Daymet Prec", "Daymet Tmax", "Daymet Tmin", "Daymet VP", "CHIRPS Prec")
rownames(daily_abiotic_variable_corr) <- c("Daymet Prec", "Daymet Tmax", "Daymet Tmin", "Daymet VP", "CHIRPS Prec")
diag(daily_abiotic_variable_corr) = NA
daily_abiotic_variable_corr_sp <- rcorr(as.matrix(daily_abiotic_variable_corr), type = "spearman") # the function rcorr() resolve the issue above.

corrplot(daily_abiotic_variable_corr, method = "color", type = "lower", col=corrPalette(200), cl.pos = "r",
         tl.col="black", tl.offset = 0.2, tl.srt = 70 , order = NULL, diag = TRUE, na.label = "NA", 
         addCoef.col = "darkgrey", sig.level = 0.05, insig = "pch", p.mat = daily_abiotic_variable_corr_sp$P)

# Monthly estimates
monthly_abiotic_variable_corr <- cor(na.omit(Fpet_all_estimates[,11:15]), method="spearman") 
colnames(monthly_abiotic_variable_corr) <- c("TC Prec.", "TC Tmax", "TC Tmin", "TC VP", "TC wind")
rownames(monthly_abiotic_variable_corr) <- c("TC Prec.", "TC Tmax", "TC Tmin", "TC VP", "TC wind")
diag(monthly_abiotic_variable_corr) = NA
monthly_abiotic_variable_corr_sp <- rcorr(as.matrix(monthly_abiotic_variable_corr), type = "spearman") # the function rcorr() resolve the issue above.

corrplot(monthly_abiotic_variable_corr, method = "color", type = "lower", col=corrPalette(200), cl.pos = "r",
         tl.col="black", tl.offset = 0.2, tl.srt = 70 , order = NULL, diag = TRUE, na.label = "NA", 
         addCoef.col = "darkgrey", sig.level = 0.05, insig = "pch", p.mat = monthly_abiotic_variable_corr_sp$P)

# With the information above, we decide to ignore LST, minimum temperature and vapor pressure.


# ____________________________________________________________________________ #
# Plots


# Format
Fpet_tmax_melt <- melt(Fpet_all_estimates, id.vars = c("site","date"), measure.vars = , c("DM_tmax","MODIS_LST","TC_tmax"))
Fpet_prec_melt <- melt(Fpet_all_estimates, id.vars = c("site","date"), measure.vars = , c("DM_prec","CHIRPS_prec","TC_prec"))
Fpet_ws_melt <- melt(Fpet_all_estimates, id.vars = c("site","date"), measure.vars = , c("TC_ws"))
Fpet_vp_melt <- melt(Fpet_all_estimates, id.vars = c("site","date"), measure.vars = , c("DM_vp","TC_vp"))
Fpet_tmax_melt$date <- as.Date(Fpet_tmax_melt$date)
Fpet_prec_melt$date <- as.Date(Fpet_prec_melt$date)
Fpet_ws_melt$date <- as.Date(Fpet_ws_melt$date)
Fpet_vp_melt$date <- as.Date(Fpet_vp_melt$date)

# Plot
# Maximum temperature
ggplot(data = Fpet_tmax_melt) + 
  geom_line(aes(x = date, colour=variable, y = value)) +
  geom_line(data = na.omit(subset(Fpet_tmax_melt, variable=="TC_tmax")), aes(x = date, colour=variable, y = value)) +
  facet_wrap(~ site, ncol=3)

# Precipitation
ggplot(data = Fpet_prec_melt) + 
  geom_line(aes(x = date, colour=variable, y = value)) +
  geom_line(data = na.omit(subset(Fpet_prec_melt, variable=="TC_prec")), aes(x = date, colour=variable, y = value)) +
  facet_wrap(~ site, ncol=3)

# Wind speed
ggplot(data = na.omit(Fpet_ws_melt)) + 
  geom_line(aes(x = date, colour=variable, y = value)) +
  facet_wrap(~ site, ncol=3)

# Vapor pressure
ggplot(data = Fpet_vp_melt) + 
  geom_line(aes(x = date, colour=variable, y = value)) +
  geom_line(data = na.omit(subset(Fpet_vp_melt, variable=="TC_vp")), aes(x = date, colour=variable, y = value)) +
  facet_wrap(~ site, ncol=3)


# ---------------------------------------------------------------------------- #
# Create weather summaries prior to collecting period
# ---------------------------------------------------------------------------- #



# First, we replace TerraClimate NA values with the preceding observation.
# This will allow to average and weight when computing climate summaries during wasp reproductive activity.
# We use the package data.table for this. Code from: https://stackoverflow.com/questions/7735647/replacing-nas-with-latest-non-na-value
# Alternatively, the package zoo allows to do this.

library(data.table) # Careful, there will be conflicting functions with the package reshape2

TC_prec_cp <- data.table(TC_prec_cp = Fpet_all_estimates$TC_prec)
TC_prec_cp[, y_forward_fill := TC_prec_cp[1], .(cumsum(!is.na(TC_prec_cp)))]
TC_tmax_cp <- data.table(TC_tmax_cp = Fpet_all_estimates$TC_tmax)
TC_tmax_cp[, y_forward_fill := TC_tmax_cp[1], .(cumsum(!is.na(TC_tmax_cp)))]
TC_tmin_cp <- data.table(TC_tmin_cp = Fpet_all_estimates$TC_tmin)
TC_tmin_cp[, y_forward_fill := TC_tmin_cp[1], .(cumsum(!is.na(TC_tmin_cp)))]
TC_ws_cp <- data.table(TC_ws_cp = Fpet_all_estimates$TC_ws)
TC_ws_cp[, y_forward_fill := TC_ws_cp[1], .(cumsum(!is.na(TC_ws_cp)))]
TC_vp_cp <- data.table(TC_vp_cp = Fpet_all_estimates$TC_vp)
TC_vp_cp[, y_forward_fill := TC_vp_cp[1], .(cumsum(!is.na(TC_vp_cp)))]

Fpet_all_estimates$TC_prec_cp <- as.numeric(TC_prec_cp[[2]])
Fpet_all_estimates$TC_tmax_cp <- as.numeric(TC_tmax_cp[[2]])
Fpet_all_estimates$TC_tmin_cp <- as.numeric(TC_tmin_cp[[2]])
Fpet_all_estimates$TC_ws_cp <- as.numeric(TC_ws_cp[[2]])
Fpet_all_estimates$TC_vp_cp <- as.numeric(TC_vp_cp[[2]])

# We estimated that the reproductive period of the wasps is within 5-8 weeks before we collected their offspring.
# This estimate is based on observation on the length of fig maturation, as well as destructive observations on fig before ripening.

Fpet_collecting_dates <- data.frame(
  site=c("158","172","112","113","95","179","201","96","70"),
  F2012=c("03-11-2012","07-11-2012","10-11-2012","06-12-2012","29-11-2012","25-11-2012","13-11-2012","24-11-2012","19-11-2012"),
  S2013=c("18-05-2013","21-05-2013","22-05-2013","28-05-2013","15-06-2013",NA,"12-06-2013","28-06-2013","13-06-2013"),
  F2013=c("05-11-2013","07-11-2013","09-11-2013","12-11-2013","17-11-2013",NA,NA,"30-11-2013","27-11-2013"),
  S2014=c("04-07-2014","18-05-2014","20-05-2014","25-05-2014","14-06-2014","11-06-2014","19-06-2014","31-05-2014","09-06-2014"))
Fpet_collecting_dates$F2012 <- as.Date(Fpet_collecting_dates$F2012, format = "%d-%m-%Y")
Fpet_collecting_dates$S2013 <- as.Date(Fpet_collecting_dates$S2013, format = "%d-%m-%Y")
Fpet_collecting_dates$F2013 <- as.Date(Fpet_collecting_dates$F2013, format = "%d-%m-%Y")
Fpet_collecting_dates$S2014 <- as.Date(Fpet_collecting_dates$S2014, format = "%d-%m-%Y")
head(Fpet_collecting_dates)

# # We want the range of dates of one month (30 days), 6 weeks (42 days) prior to sampling.
# Fpet_6weeksprior <- cbind(Fpet_collecting_dates[,1], Fpet_collecting_dates[,2:5] - 42); colnames(Fpet_6weeksprior) <- colnames(Fpet_collecting_dates)
# Fpet_10weeksprior <- cbind(Fpet_6weeksprior[,1], Fpet_6weeksprior[,2:5] - 30); colnames(Fpet_10weeksprior) <- colnames(Fpet_collecting_dates)
# Fpet_priorsampling_dates.list <- list(Fpet_6weeksprior, Fpet_10weeksprior); rm(Fpet_6weeksprior,Fpet_10weeksprior)
# Fpet_priorsampling_dates.list


# We want the range of dates of 3 weeks (21 days), 5 weeks (35 days) prior to sampling.
Fpet_5weeksprior <- cbind(Fpet_collecting_dates[,1], Fpet_collecting_dates[,2:5] - 35); colnames(Fpet_5weeksprior) <- colnames(Fpet_collecting_dates)
Fpet_8weeksprior <- cbind(Fpet_5weeksprior[,1], Fpet_5weeksprior[,2:5] - 21); colnames(Fpet_8weeksprior) <- colnames(Fpet_collecting_dates)
Fpet_priorsampling_dates.list <- list(Fpet_5weeksprior, Fpet_8weeksprior); rm(Fpet_5weeksprior,Fpet_8weeksprior)
Fpet_priorsampling_dates.list

# Use loops to obtain the weather summary wanted by using the date range from the list.
weather_summaries_prior_sampling <- data.frame(Date=as.Date(character()),
                                               File=character(), 
                                               User=character(), 
                                               stringsAsFactors=FALSE) 
for(i in 2:5) {
  print(colnames(Fpet_priorsampling_dates.list[[1]])[[i]])
  for(j in 1:9) { 
    print(Fpet_priorsampling_dates.list[[1]][[j,1]])
    tmp <- Fpet_all_estimates[which(Fpet_all_estimates$site == Fpet_priorsampling_dates.list[[1]][[j,1]] & 
                                      Fpet_all_estimates$date > Fpet_priorsampling_dates.list[[2]][[j,i]] & 
                                      Fpet_all_estimates$date <= Fpet_priorsampling_dates.list[[1]][[j,i]]), ]
    tmp2 <- tmp %>%
      group_by(site) %>%
      summarise(latitude = mean(latitude), longitude = mean(longitude),  # Here, change sums or means for the parameters of choice (sums or mean for precipitation CHIRPS and Daymet)
                DM_tmax=mean(DM_tmax, na.rm = T), DM_tmin=mean(DM_tmin, na.rm = T), DM_prec=sum(DM_prec, na.rm = T),  
                DM_vp=mean(DM_vp, na.rm = T), CHIRPS_prec=sum(CHIRPS_prec, na.rm = T), MODIS_LST=mean(MODIS_LST, na.rm = T), # Remove vapor pressure
                TC_tmax=mean(TC_tmax_cp, na.rm = T), TC_tmin=mean(TC_tmin_cp, na.rm = T),# Do not sum for TerraClimate data.Inconsistent with other data?
                TC_prec=mean(TC_prec_cp, na.rm = T), TC_vp=mean(TC_vp_cp, na.rm = T), # Remove vapor pressure?
                TC_ws=mean(TC_ws_cp, na.rm = T))
    tmp2$season <- colnames(Fpet_priorsampling_dates.list[[1]])[[i]]
    tmp2 <- tmp2[,c(length(tmp2),1:length(tmp2)-1)]
    weather_summaries_prior_sampling <- rbind(weather_summaries_prior_sampling, tmp2)
    rm(tmp,tmp2)
  }
}

View(weather_summaries_prior_sampling)
write.csv(weather_summaries_prior_sampling, paste(workdir,"generated_data/","Fpet_weather_summaries_prior_sampling_2012-2014.csv", sep = ""), row.names = FALSE)


# ____________________________________________________________________________ #
# The following code is to test if the loop above does the intended computations:


# Test for site 158, fall 2012 collected 03-11-2012 --------------------

# 8 weeks prior sampling:
setdiff(as.Date("03-11-2012", format = "%d-%m-%Y") -56, Fpet_priorsampling_dates.list[[2]][[1,2]]) # If returns "numeric(0)" then it is the same date
# 5 weeks prior sampling:
setdiff(as.Date("03-11-2012", format = "%d-%m-%Y") -35, Fpet_priorsampling_dates.list[[1]][[1,2]]) # If returns "numeric(0)" then it is the same date
# Dates for site 158 collection in fall 2012 in the "Fpet_priorsampling_dates.list" are accurate.

# Now we subset the weather/climate data for these time period:
test_S158_F2012 <- Fpet_all_estimates[which(Fpet_all_estimates$site == "158" &
                                       Fpet_all_estimates$date > as.Date("03-11-2012", format = "%d-%m-%Y") -56 &
                                       Fpet_all_estimates$date <= as.Date("03-11-2012", format = "%d-%m-%Y") -35), ]

# We summarize:
test_S158_F2012_sum <- test_S158_F2012 %>%
  group_by(site) %>%
  summarise(latitude = mean(latitude), longitude = mean(longitude), 
            DM_tmax=mean(DM_tmax, na.rm = T), DM_tmin=mean(DM_tmin, na.rm = T), DM_prec=sum(DM_prec, na.rm = T),  
            DM_vp=mean(DM_vp, na.rm = T), CHIRPS_prec=sum(CHIRPS_prec, na.rm = T), MODIS_LST=mean(MODIS_LST, na.rm = T), 
            TC_tmax=mean(TC_tmax_cp, na.rm = T), TC_tmin=mean(TC_tmin_cp, na.rm = T),
            TC_prec=mean(TC_prec_cp, na.rm = T), TC_vp=mean(TC_vp_cp, na.rm = T),
            TC_ws=mean(TC_ws_cp, na.rm = T))

# We obtain one row for site 158 collection fall 2012, we compare it to the generated data set:
rbind(test_S158_F2012_sum,weather_summaries_prior_sampling[1,-1])
# The rows are identical, which was expected.


# Test for site 113, fall 2013, collected 2013-11-12 --------------------

# 8 weeks prior sampling:
setdiff(as.Date("12-11-2013", format = "%d-%m-%Y") -56, Fpet_priorsampling_dates.list[[2]][[4,4]]) # If returns "numeric(0)" then it is the same date
# 5 weeks prior sampling:
setdiff(as.Date("12-11-2013", format = "%d-%m-%Y") -35, Fpet_priorsampling_dates.list[[1]][[4,4]]) # If returns "numeric(0)" then it is the same date
# Dates for site 158 collection in fall 2012 in the "Fpet_priorsampling_dates.list" are accurate.

# Now we subset the weather/climate data for these time period:
test_S113_F2013 <- Fpet_all_estimates[which(Fpet_all_estimates$site == "113" &
                                              Fpet_all_estimates$date > as.Date("12-11-2013", format = "%d-%m-%Y") -56 &
                                              Fpet_all_estimates$date <= as.Date("12-11-2013", format = "%d-%m-%Y") -35), ]

# We summarize:
test_S113_F2013_sum <- test_S113_F2013 %>%
  group_by(site) %>%
  summarise(latitude = mean(latitude), longitude = mean(longitude), 
            DM_tmax=mean(DM_tmax, na.rm = T), DM_tmin=mean(DM_tmin, na.rm = T), DM_prec=sum(DM_prec, na.rm = T),  
            DM_vp=mean(DM_vp, na.rm = T), CHIRPS_prec=sum(CHIRPS_prec, na.rm = T), MODIS_LST=mean(MODIS_LST, na.rm = T), 
            TC_tmax=mean(TC_tmax_cp, na.rm = T), TC_tmin=mean(TC_tmin_cp, na.rm = T),
            TC_prec=mean(TC_prec_cp, na.rm = T), TC_vp=mean(TC_vp_cp, na.rm = T),
            TC_ws=mean(TC_ws_cp, na.rm = T))

# We obtain one row for site 158 collection fall 2012, we compare it to the generated data set:
rbind(test_S113_F2013_sum,weather_summaries_prior_sampling[21,-1])
# The rows are identical, which was expected.



# ---------------------------------------------------------------------------- #
# Plot annual data
# ---------------------------------------------------------------------------- #



# Import data
transect_annual_data <- read.csv(paste(workdir,"generated_data/","transect_annual_data_2012-2014.csv", sep = ""))


# Plots
ggplot(transect_annual_data, aes( x = latitude)) +
  geom_line(aes( y = TC_tmax_2012 ), linetype="solid", colour="black", size=1) +
  scale_x_reverse("",breaks=seq(23.5,30,0.5))

ggplot(transect_annual_data, aes( x = latitude)) +
  geom_line(aes( y = TC_tmax_2013 ), linetype="solid", colour="black", size=1) +
  scale_x_reverse("",breaks=seq(23.5,30,0.5))

ggplot(transect_annual_data, aes( x = latitude)) +
  geom_line(aes( y = TC_tmax_2014 ), linetype="solid", colour="black", size=1) +
  scale_x_reverse("Latitude",breaks=seq(23.5,30,0.5))

ggplot(transect_annual_data, aes( x = latitude)) +
  geom_line(aes( y = TC_prec_2012 ), linetype="solid", colour="black", size=1) +
  scale_x_reverse("",breaks=seq(23.5,30,0.5))

ggplot(transect_annual_data, aes( x = latitude)) +
  geom_line(aes( y = TC_prec_2013 ), linetype="solid", colour="black", size=1) +
  scale_x_reverse("",breaks=seq(23.5,30,0.5))

ggplot(transect_annual_data, aes( x = latitude)) +
  geom_line(aes( y = TC_prec_2014 ), linetype="solid", colour="black", size=1) +
  scale_x_reverse("Latitude",breaks=seq(23.5,30,0.5))



# ---------------------------------------------------------------------------- #
# End of R script 3
# ---------------------------------------------------------------------------- #



# Clear working environment 
rm(list = ls())
dev.off()
