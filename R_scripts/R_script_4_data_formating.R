# ============================================================================ #
# R script 4: formatting the raw data
# ============================================================================ #



# Article title: Landscape-Level Analysis of a Fig-Pollinator-Antagonist Community: Spatial and Temporal Variation in a Fig Wasp Community and its Response to Biotic and Abiotic Factors
# By Finn Piatscheck, Justin van Goor, Derek Houston and John Nason.



# The following script allows to reproduce the main data sets needed for the downstream analyses.
# The raw files contain all the raw observations made in the field at the tree level and in the field and lab at the fig (syconium) level.
# Some data sets were generated using other scripts: flowering_and_syconium_landscape.csv and Fpet_weather_summaries_prior_sampling_2012-2014.csv.
# This script computes some variables, and merges fig- and tree-level observations, as well as landscape indices and site- and time- specific climate estimates.
# Thus, it is recommended to just run the whole script at once to obtain the two data sets formatted in the wanted way.



# ---------------------------------------------------------------------------- #
# Load packages and set working directory
# ---------------------------------------------------------------------------- #



# R version --------------------

# R version used for these analyses: 4.3.1.


# Set working directory --------------------

# To replicate the following analyses, create a folder named "working_directory".
# It should have 4 sub-folders: "raw_data", "generated_data", "tables" and "figures".
# workdir <- "PathToWorkingDirectory/" # Replace here with the path to the directory in which you placed the raw data files.
setwd(workdir)



# Load packages --------------------

library(dplyr) # Version 1.1.3



# ---------------------------------------------------------------------------- #
# Import and format tree-level observations
# ---------------------------------------------------------------------------- #



# Import fig-level data --------------------

# The raw data contains empty cells and NAs. NAs were noted when measurements on fig and trees were missing, empty cells in the wasp count correspond to 0 wasp.
df_fig <- read.csv(paste(workdir,"raw_data/","syconium_level_variables.csv", sep = ""), h=T, na.strings=c("NA",""), stringsAsFactors = FALSE)
# Be careful when importing these data. Tree names are sometimes leading with a "0". They will be changed to single digit numbers if opened and saved with EXCEL.


# Format fig-level data --------------------

# Sums of female and male wasp counts are created for each species.
# Note that for Idarnes species counts, only females will be considered because males were too hard to identify at the species level.
# Thus, variables LO1_f, SO1_f and SO2_f will be used in the analyses.
# Males are however considered in the total parasite count.
df_fig$pollinators <- df_fig$pegoscapus_f+df_fig$pegoscapus_m
df_fig$heterandrium_1 <- df_fig$heterandrium_1_f+df_fig$heterandrium_1_mw+df_fig$heterandrium_1_mu
df_fig$heterandrium_2 <- df_fig$heterandrium_2_f+df_fig$heterandrium_2_mw+df_fig$heterandrium_2_mu
df_fig$ficicola <- df_fig$ficicola_f+df_fig$ficicola_m
df_fig$physothorax <- df_fig$physothorax_f+df_fig$physothorax_m
df_fig$parasites <- df_fig$LO1_f+df_fig$SO1_f+df_fig$SO2_f+df_fig$idarnes_m+df_fig$heterandrium_1+df_fig$heterandrium_2+df_fig$ficicola+df_fig$physothorax+df_fig$sycophila+df_fig$unidentified_idarnes

# Another column will sum all species counts.
df_fig$wasps_total <- df_fig$pollinators+df_fig$parasites+df_fig$unidentified_wasps+df_fig$unidentified_idarnes


# Another column will sum all species counts.
df_fig$wasps_total <- df_fig$pollinators+df_fig$parasites+df_fig$unidentified_wasps+df_fig$unidentified_idarnes



# Summarize data at the tree-level --------------------

# Summarize rows by season, site and tree, creates the data frame wanted to be merged with tree level variables.
# Summarize by sums or means. These summary will not be part of the analyses, but will be used for visual representation in figures.

df_to_merge <- df_fig %>% 
  group_by(season, site, tree, sampling_date) %>%
  summarise(pollinators = mean(na.omit(pollinators)), 
            LO1=mean(na.omit(LO1_f)), 
            SO1=mean(na.omit(SO1_f)), 
            SO2=mean(na.omit(SO2_f)), 
            heterandrium_1=mean(na.omit(heterandrium_1)), 
            heterandrium_2=mean(na.omit(heterandrium_2)), 
            ficicola=mean(na.omit(ficicola)), 
            physothorax=mean(na.omit(physothorax)),
            sycophila=mean(na.omit(sycophila)), 
            parasites=mean(na.omit(parasites)), 
            wasps_total=mean(na.omit(wasps_total)))
head(df_to_merge)



# ---------------------------------------------------------------------------- #
# Import and format tree-level observations
# ---------------------------------------------------------------------------- #



df_tree <- read.csv(paste(workdir,"raw_data/","tree_level_variables.csv", sep = ""), h=T, na.strings=c("NA",""), stringsAsFactors = FALSE)

# Compute the flowering asynchrony variable
df_tree$asynchrony <- 1/((df_tree$pre_female/df_tree$total_syconia)^2 + 
                          (df_tree$early_female/df_tree$total_syconia)^2 +
                          (df_tree$late_female/df_tree$total_syconia)^2 +
                          (df_tree$early_interphase/df_tree$total_syconia)^2 +
                          (df_tree$late_interphase/df_tree$total_syconia)^2 +
                          (df_tree$early_male/df_tree$total_syconia)^2 +
                          (df_tree$late_male/df_tree$total_syconia)^2)

# Compute crown volume
df_tree$tree_volume <- (pi/6)*df_tree$tree_diameter*df_tree$tree_diameter*df_tree$tree_height # Note that this formula is equivalent to the one used for fruit_volume: 4/3*pi*d/2*d/2*h/2

# Compute crop size 
df_tree$crop_size <- df_tree$tree_volume*df_tree$reproduction


# re-order
df_tree <- df_tree[, c(1:11,23,12:13,24,14:22)]

# Merge the data
df_tree$id_merge  <- 1:nrow(df_tree) # This allows to re-order the data frame in the original order
df_tree_fig <- merge(df_tree,df_to_merge, by = c("season", "site", "tree"), all = TRUE)
df_tree_fig <- df_tree_fig[order(df_tree_fig$id_merge), ]


# Add tree-level flowering and syconia landscape variables --------------------

# Note: these variable have been realized before with the function neighborVar detailed in the script "density neighborhood lanscape function" 2.

df_landscape <- read.csv(paste(workdir,"generated_data/","flowering_and_syconium_landscape.csv", sep = ""))
df_tree_fig$ID <- as.factor(paste(df_tree_fig$season, df_tree_fig$site, df_tree_fig$tree, sep = "-"))
setdiff(df_landscape[,1], df_tree_fig$ID) # Same IDs if character(0)
setdiff(df_tree_fig$ID, df_landscape[,1]) # Difference here
# The difference is due to the outliers that were eliminated in the script neighborhood_landscape_spatial_analysis.R.
# Thus we merge with the argument "all = TRUE".
# The outliers will be examined but will have NA values for statistical modeling. 
df_tree_fig_neig <- merge(df_tree_fig, df_landscape, by = "ID", sort = FALSE, all = TRUE)
df_tree_fig_neig <- df_tree_fig_neig[order(df_tree_fig_neig$id_merge), ]


# Add site- and time- specific prior sampling climate data
df_climdata <- read.csv(paste(workdir,"generated_data/","Fpet_weather_summaries_prior_sampling_2012-2014.csv", sep = ""))
df_tree_fig_neig_clim <- merge(df_tree_fig_neig, df_climdata[,c(-3,-4)], by = c("season", "site"), all = TRUE, sort = FALSE) # Make sure to remove the site coordinates, they are already in the previous dataset.
df_tree_fig_neig_clim <- df_tree_fig_neig_clim[order(df_tree_fig_neig_clim$id_merge), ]
df_tree_fig_neig_clim$id_merge <- NULL
df_tree_fig_neig_clim$ID <- NULL


# Rename the dataframe
df_tree_data <- df_tree_fig_neig_clim # This one will be scaled and prepared to be used for tree level analysis or to be merged with the fig data.


# Clear
rm(df_to_merge, df_tree, df_tree_fig, df_tree_fig_neig, df_tree_fig_neig_clim, df_landscape, df_climdata)

# Save
write.csv(df_tree_data, paste(workdir,"generated_data/","df_tree_data.csv", sep = ""), row.names = FALSE)
# WARNING: When opening the file in EXCEL, the "tree" variable will be missing the leading "0" of some tree names.



# ---------------------------------------------------------------------------- #
# Add tree-level observations to the fig data set
# ---------------------------------------------------------------------------- #


# We remove the wasp summary and sampling dates already present in the fig data set
df_tree_tmp <- df_tree_data[,!names(df_tree_data) %in% c("pollinators","LO1","SO1","SO2","heterandrium_1","heterandrium_2","ficicola","physothorax","sycophila","parasites","wasps_total","sampling_date")] 
df_fig$id_merge  <- 1:nrow(df_fig)
df_fig_data <- merge(df_fig, df_tree_tmp, by = c("season", "site", "tree"), all = FALSE)
df_fig_data <- df_fig_data[order(df_fig$id_merge), ]
df_fig_data$id_merge <- NULL 

# Difference of 2 in number of observation in both dataframes, let's find out which individuals didn't make it when merging the data set..
fullID1 <- as.factor(paste(df_fig$season, df_fig$site, df_fig$tree, sep = "-"))
fullID2 <- as.factor(paste(df_fig_data$season, df_fig_data$site, df_fig_data$tree, sep = "-"))
setdiff(fullID1, fullID2); setdiff(fullID2, fullID1);rm(fullID1, fullID2)
# Results should be "character(0)", two mislabeled trees were removed from the raw data as it was not known to whom they corresponded.

# Clear
rm(df_fig, df_tree_tmp)

# Save
write.csv(df_fig_data, paste(workdir,"generated_data/","df_fig_data.csv", sep = ""), row.names = FALSE)
# WARNING: When opening the file in EXCEL, the "tree" variable will be missing the leading "0" of some tree names.



# ---------------------------------------------------------------------------- #
# End of R script 4
# ---------------------------------------------------------------------------- #



# Clear working environment 
rm(list = ls())
