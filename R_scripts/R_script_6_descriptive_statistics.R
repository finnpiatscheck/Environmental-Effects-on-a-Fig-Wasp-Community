# ============================================================================ #
# R script 6: reproduction of the descriptive statistical analyses
# ============================================================================ #



# Article title: Landscape-Level Analysis of a Fig-Pollinator-Antagonist Community: Spatial and Temporal Variation in a Fig Wasp Community and its Response to Biotic and Abiotic Factors
# By Finn Piatscheck, Justin van Goor, Derek Houston and John Nason.


# This R script allows to reproduce all the descriptive statistics of wasp species counts and fig- and tree-level host-related ecological variables.
# Table 1 and supplementary tables S1,S2 and S3 can directly be obtained from this script in a text format (i.e,. csv).


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
library(reshape2) # Version 1.4.4
library(plotrix) # Version 3.8-2



# ---------------------------------------------------------------------------- #
# Import the data and formatting
# ---------------------------------------------------------------------------- #



# Open df_fig_data.csv and df_tree_data.csv which are formatted files combining fig- and tree-level observations.
# To read the raw data and format them, refer to the R script: R_script_Piatscheck.et.al.2020_Fpet_data_formating.R

df_fig_data <- read.csv(paste(workdir,"generated_data/","df_fig_data.csv", sep = ""))
# df_fig_data is the data set that contains both fig and tree- level observations. 
# However, we measured ecological variables on trees that were not containing E phase syconia.
# These observations are not in this data set.

df_tree_data <- read.csv(paste(workdir,"generated_data/","df_tree_data.csv", sep = ""))
# df_tree_data is a simplified version of the dataset, each row correspond to one tree observed (at a particular season).
# Fig-level variables are absent here.
# Wasp counts are averages (which can be changed to sums in R script data formating)
# The dataset contains, however, measurments made on tree not containing E phase syconia ("male" phase figs containing wasps), which are absent in the fig-level dataset.



# ____________________________________________________________________________ #
# Getting the dataframes formatted

df_fig_data$site <- factor(df_fig_data$site, levels = c("158", "172","112","113","95","179","201","96","70", "250")) # Keep the factors in this order, better for plotting, also site 250 here used, not in tree level analyses
df_fig_data$season <- factor(df_fig_data$season, levels = c("F2012", "S2013","F2013","S2014"))
df_tree_data$site <- factor(df_tree_data$site, levels = c("158", "172","112","113","95","179","201","96","70", "250")) # Keep the factors in this order, better for plotting, also site 250 here used, not in tree level analyses
df_tree_data$season <- factor(df_tree_data$season, levels = c("F2012", "S2013","F2013","S2014"))

# Replace outliers with NAs
df_fig_data[rownames(df_fig_data[which(df_fig_data$tree_volume > 6000),]),c(which(colnames(df_fig_data) == "tree_volume"), which(colnames(df_fig_data) == "crop_size"))] <- NA
df_tree_data[rownames(df_tree_data[which(df_tree_data$tree_volume > 6000),]),c(which(colnames(df_tree_data) == "tree_volume"), which(colnames(df_tree_data) == "crop_size"))] <- NA


# Add Idarnes variables
df_fig_data$idarnes_allsp <- df_fig_data$LO1_f + df_fig_data$SO1_f + df_fig_data$SO2_f + df_fig_data$idarnes_m + df_fig_data$unidentified_idarnes



# ---------------------------------------------------------------------------- #
# Descriptive statistics
# ---------------------------------------------------------------------------- #

# [Note: Here I will reorder these calculation to follow the order in the manuscript]

# ____________________________________________________________________________ #
# Field observations summary statistics


# Number of tree georeferenced
length(unique(as.factor(paste(df_tree_data$site, df_tree_data$tree, sep = "-"))))
# 947 trees in 9 sites

# Number of flowering tree recorded throughout the 4 field trips
length(which(df_tree_data$flowering == "Yes"))
# 868

# Number of tree sampled with D phases syconia
length(unique(df_fig_data$tree))
# 122


# Number of syconia sampled
length(df_fig_data$fig)
# 2304

# Number of adult fig wasp sampled
sum(df_fig_data$wasps_total)
# 216699


# ____________________________________________________________________________ #
# Host-related ecological variables


# Tree sizes --------------------

# Mean tree volume
mean(df_tree_data$tree_volume, na.rm = TRUE)
std.error(df_tree_data$tree_volume, na.rm = TRUE)
# 319m, SE = 13.03

# Mean height volume
mean(df_tree_data$tree_height, na.rm = TRUE)
# 6.94m

# Range height volume
range(df_tree_data$tree_height, na.rm = TRUE)
# 0.1m to 27m


# Tree reproduction --------------------

# Mean reproduction
mean(df_tree_data$reproduction, na.rm = TRUE)
# 0.11%: low reproduction is the most common reproductive effort


# Within-tree Asynchrony --------------------

# Asynchronous trees proportion
length(which(df_tree_data$asynchrony == 1))
length(which(df_tree_data$asynchrony > 1))
77/(658+77)*100
658/(658+77)*100
# 10.48% synchronous
# 89.52% asynchronus

# ____________________________________________________________________________ #
# Wasps count statistics


# All figs: 2367
# Early D phase with no wasps emerged removed: 2304 figs


# General --------------------

# Number of wasps collected
sum(df_fig_data$pollinators) + sum(df_fig_data$parasites)
# 215209 identified wasps
sum(df_fig_data$wasps_total) # Larger number due to unidentified wasps
# 216699 with identified and unidentified wasps
216699-215209
# 1490 unidentified wasps

# Average number of wasp per syconia
mean(df_fig_data$wasps_total)
sd(df_fig_data$wasps_total)
std.error(df_fig_data$wasps_total)
# 94.05, sd = 65.89, se = 1.37

# Average number of wasp per syconia
mean(df_fig_data$wasps_total)
# 94.05

# Wasp range
range(df_fig_data$wasps_total)
# 2 to 547


# Foundresses --------------------

# Average number of foundresses
mean(na.omit(df_fig_data$foundress_number))
# 1.32

# Syconia without foundresses
length(which(df_fig_data$foundress_number == 0))
# 171 
171/2304*100
# proportion = 7.42%

# Foundresses without offspring but parasites that developed
length(which(df_fig_data$foundress_number >0 & df_fig_data$pegoscapus_f == 0 & df_fig_data$pegoscapus_m == 0))
# 296
296/2304*100
# proportion = 12.85%


# Pegoscapus (pollinators) --------------------

# Proportion of pollinators 
sum(df_fig_data$pollinators)/(sum(df_fig_data$wasps_total)-sum(df_fig_data$unidentified_wasps))*100
# 40.98%

# Syconia with pollinators
length(which(df_fig_data$pollinators >0))
# 1810
1810/2304*100
# proportion = 78.56%

# Syconia with only pollinators 
length(which(df_fig_data$parasites == 0))
# 24
24/2304*100
# proportion = 1.04

# Syconia without pollinators
length(which(df_fig_data$pollinators == 0 & df_fig_data$parasites >0))
# 494
494/2304*100
# proportion = 21.44%


# Parasites --------------------

# No pollinators, no foundresses, only parasites
length(which(df_fig_data$pollinators == 0 & df_fig_data$parasites > 0 & df_fig_data$foundress_number == 0))
# 165
165/2304*100
# proportion = 7.16%

# Pollinators/Parasites
sum(df_fig_data$pollinators)/sum(df_fig_data$parasites)
# 0.7

# Syconia with Idarnes
length(which(df_fig_data$idarnes_allsp > 0))
# 2260
2260/2304*100
# proportion = 98.09%

# Syconia without parasites
length(which(df_fig_data$pollinators >0 & df_fig_data$parasites == 0))
# 24
length(which(df_fig_data$pollinators >0 & df_fig_data$parasites == 0))/2304*100
# 1.04%


# Idarnes flavicolis sp1 (LO1) --------------------

# Proportion of Idarnes flavicolis LO1
sum(df_fig_data$LO1_f)/(sum(df_fig_data$wasps_total)-sum(df_fig_data$unidentified_wasps))*100
# 17.63%

# Syconia with LO1s
length(which(df_fig_data$LO1_f >0 ))
# 1878
1878/2304*100
# proportion = 81.51%

# Syconia with only LO1s
length(which(df_fig_data$pegoscapus_f == 0 & df_fig_data$pegoscapus_m == 0 & df_fig_data$LO1_f > 0 & df_fig_data$SO1_f == 0 & df_fig_data$SO2_f == 0))
# Only LO1 in 165 syconia
165/2304*100
# 7.16%

# Syconia with no pollinators
length(which(df_fig_data$pegoscapus_f == 0 & df_fig_data$pegoscapus_m == 0 & df_fig_data$LO1_f > 0))
425/2304*100
# 18.45%

# Idarnes carme sp1 (SO1) --------------------

# Proportion of Idarnes species SO1
sum(df_fig_data$SO1_f)/(sum(df_fig_data$wasps_total)-sum(df_fig_data$unidentified_wasps))*100
# 6.13%

# Syconia with SO1s
length(which(df_fig_data$SO1_f >0 ))
# 1272
1272/2304*100
# proportion = 55.21%
length(which(df_fig_data$pegoscapus_f == 0 & df_fig_data$pegoscapus_m == 0 & df_fig_data$LO1_f == 0 & df_fig_data$SO1_f > 0 & df_fig_data$SO2_f == 0))
# only SO1 in 11 syconia


# Idarnes carme sp2 --------------------

# Proportion of Idarnes species SO2
sum(df_fig_data$SO2_f)/(sum(df_fig_data$wasps_total)-sum(df_fig_data$unidentified_wasps))*100
# 10.30%

# Syconia with SO2s
length(which(df_fig_data$SO2_f >0 ))
# 1453
1453/2304*100
# proportion = 63.06%
length(which(df_fig_data$pegoscapus_f == 0 & df_fig_data$pegoscapus_m == 0 & df_fig_data$LO1_f == 0 & df_fig_data$SO1_f == 0 & df_fig_data$SO2_f > 0))
# only SO2 in 17 syconia



# Heterandrium sp1 --------------------

# Proportion of Heterandrium 1
sum(df_fig_data$heterandrium_1)/(sum(df_fig_data$wasps_total)-sum(df_fig_data$unidentified_wasps))*100
# 3.87%

# Syconia with Heterandrium 1
length(which(df_fig_data$heterandrium_1 >0))
# 1165
1165/2304*100
# proportion = 50.56%


# Heterandrium sp2 --------------------

# Proportion of Heterandrium 2
sum(df_fig_data$heterandrium_2)/(sum(df_fig_data$wasps_total)-sum(df_fig_data$unidentified_wasps))*100
# 2.45%

# Syconia with Heterandrium 2
length(which(df_fig_data$heterandrium_2 >0))
# 1049
1049/2304*100
# proportion = 45.53%


# Ficicola sp --------------------

# Proportion of Ficicola
sum(df_fig_data$ficicola)/(sum(df_fig_data$wasps_total)-sum(df_fig_data$unidentified_wasps))*100
# 0.83%

# Syconia with Ficicola
length(which(df_fig_data$ficicola >0))
# 595
595/2304*100
# proportion = 25.82%


# Physothorax sp --------------------

# Proportion of Physothorax
sum(df_fig_data$physothorax)/(sum(df_fig_data$wasps_total)-sum(df_fig_data$unidentified_wasps))*100
# 1.45%

# Syconia with Physothorax
length(which(df_fig_data$physothorax >0))
# 742
742/2304*100
# proportion = 32.20%


# Sycophyla sp --------------------

# Proportion of Sycophyla
sum(df_fig_data$sycophila)/(sum(df_fig_data$wasps_total)-sum(df_fig_data$unidentified_wasps))*100
# 0.08%

# Syconia with Sycophyla
length(which(df_fig_data$sycophila >0))
# 85
85/2304*100
# proportion = 3.69%


# Other --------------------

# Syconia with both Ficicola and Physothorax 
length(which(df_fig_data$ficicola >0 | df_fig_data$physothorax >0))
# 987
987/2304*100
# proportion = 42.84%

# Proportion of parasites (should sum to 1 added to pollinators)
sum(df_fig_data$parasites)/(sum(df_fig_data$wasps_total)-sum(df_fig_data$unidentified_wasps))*100
# 58.68%

# Proportion of Idarnes species (should NOT sum to 1 added to pollinators)
sum(df_fig_data$idarnes_allsp)/(sum(df_fig_data$wasps_total)-sum(df_fig_data$unidentified_wasps))*100
# 50.00%

# Heterandrium 1 without Idarnes
length(which(df_fig_data$heterandrium_1 > 0 & df_fig_data$LO1_f == 0 & df_fig_data$SO1_f == 0 & df_fig_data$SO2_f == 0 & df_fig_data$idarnes_m == 0))
9/2304*100
# proportion = 0.39%


# Heterandrium 2 without Idarnes
length(which(df_fig_data$heterandrium_2 > 0 & df_fig_data$LO1_f == 0 & df_fig_data$SO1_f == 0 & df_fig_data$SO2_f == 0 & df_fig_data$idarnes_m == 0))
6/2304*100
# proportion = 0.26%


# Wasp proportion at a particular site and collecting trip --------------------

# Parasite proportion at site 158 fall 2012
sum(subset(df_fig_data, site == "158" & season == "F2012")$pollinators, na.rm = TRUE)
# 87
sum(subset(df_fig_data, site == "158" & season == "F2012")$parasites, na.rm = TRUE)
# 2996
2996/(2996+87)*100
# 97.18%

# Change the species, site and season to get different proportions.
# The example above shows one site (158, most northern site) during the first field trip that was particularly dominated by the parasites.
# A figure will later be made showing this proportions.


# ____________________________________________________________________________ #
# Summary statistic tables


# Warning: if plyr is loaded, the following won't do what is intended.

# Table 2 --------------------

# Number of flowering trees by site and season
reproducing_tree_summary <- df_tree_data %>% 
  group_by(season, site) %>% 
  summarize(trees=length(tree), flowering_trees = length(which(flowering == "Yes")))

# Re-arrange the table
reproducing_tree_summary <- cbind(dcast(reproducing_tree_summary, site~season, value.var="flowering_trees"), reproducing_tree_summary[1:9,3])
reproducing_tree_summary <- reproducing_tree_summary[,c(1,6,2:5)]

# Save the table  
write.csv(reproducing_tree_summary, paste(workdir,"tables/","reproducing_tree_summary.csv", sep = ""))     


# Supplementary table 1 --------------------

# Number of syconia collected by site and season
syconia_number <- df_fig_data %>% 
  group_by(season, site) %>%
  summarize(syconia_collected=length(fig))

# Wasp sum, mean and SD table
fig_summary_stats <- df_fig_data %>% 
  group_by(season, site) %>%
  summarize_each( 
    funs(sum = sum(., na.rm = TRUE), mean = mean(., na.rm = TRUE), sd = sd(., na.rm = TRUE)), 
    foundress_number, pollinators, LO1_f, SO1_f, SO2_f, heterandrium_1, heterandrium_2, ficicola, physothorax, sycophila, parasites, wasps_total) %>%
  mutate_if(is.numeric, round, 2)

# Parasite proportions
fig_summary_stats$parasites_prop <- paste(format(fig_summary_stats$parasites_sum/(fig_summary_stats$pollinators_sum+fig_summary_stats$parasites_sum)*100, digits = 2), "%")


# Re-arrange the table
fig_summary_stats <- fig_summary_stats[,c(1:3,15,27,4,16,28,5,17,29,6,18,30,7,19,31,8,20,32,9,21,33,10,22,34,11,25,35,12,24,36,14,39)]

# merge with number of syconia collected
fig_summary_stats <- cbind(syconia_number, fig_summary_stats[,3:length(fig_summary_stats)])

# Save the table
write.csv(fig_summary_stats, paste(workdir,"tables/","fig_summary_stats.csv", sep = ""))
# Note: some columns will be removed manually, dropped column will be added manually.


# Supplementary table 2 --------------------

tree_height_stats <- (df_tree_data %>% 
                        group_by(site) %>%
                        summarize_each(funs(mean = mean(., na.rm = TRUE), min = min(., na.rm = TRUE), max = max(., na.rm = TRUE)),
                                       tree_volume, tree_height)) %>% mutate_if(is.numeric, round, 2)

# Re-arrange the table
tree_height_stats$tree_height_range <- paste(tree_height_stats$tree_height_min, "-", tree_height_stats$tree_height_max)
tree_height_stats$tree_volume_range <- paste(tree_height_stats$tree_volume_min, "-", tree_height_stats$tree_volume_max)
tree_height_stats <- tree_height_stats[,c(1,3,8,2,9)]

# Save the table
write.csv(tree_height_stats, paste(workdir,"tables/","tree_height_stats.csv", sep = "")) 

# Supplementary table 3 --------------------

tree_summary_stats <- df_tree_data %>% 
                         group_by(season, site) %>%
                         summarize_each(funs(mean = mean(., na.rm = TRUE), sd = sd(., na.rm = TRUE)), 
                                        reproduction, crop_size, asynchrony, flowering_landscape_250_125, syconium_landscape_250_125) %>%
                         mutate_if(is.numeric, round, 2)
  
  
# Re-arrange the table
tree_summary_stats <- tree_summary_stats[,c(1:3,8,4,9,5,10,6,11,7,12)]

# Save the table
write.csv(tree_summary_stats, paste(workdir,"tables/","tree_summary_stats.csv", sep = "")) 



# ---------------------------------------------------------------------------- #
# End of script 4
# ---------------------------------------------------------------------------- #



# Clear working environment 
rm(list = ls())
