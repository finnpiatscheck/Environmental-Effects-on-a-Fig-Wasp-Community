# ============================================================================ #
# R script 5: data exploration and pre-analyse filtering
# ============================================================================ #



# Article title: Landscape-Level Analysis of a Fig-Pollinator-Antagonist Community: Spatial and Temporal Variation in a Fig Wasp Community and its Response to Biotic and Abiotic Factors
# By Finn Piatscheck, Justin van Goor, Derek Houston and John Nason.


# The following script is to explore and analyse the response and predictor variables that will be used later in the statistical analyses necessary to test our hypotheses.
# The plots realized in this script are NOT the ones presented in the article.
# For publication quality figures, refer to the R script: Appendix 7 figure reproduction



# ---------------------------------------------------------------------------- #
# Load packages and set working directory
# ---------------------------------------------------------------------------- #



# R version used for these analyses: 4.3.1.


# Set working directory --------------------

# To replicate the following analyses, create a folder named "working_directory".
# It should have 4 sub-folders: "raw_data", "generated_data", "tables" and "figures".
# workdir <- "PathToWorkingDirectory/" # Replace here with the path to the directory in which you placed the raw data files.
setwd(workdir)


# Load packages --------------------

library(reshape2) # Version 1.4.4
library(ggplot2) # Version 3.4.4
library(ggridges) # Version 0.5.4
library(vcdExtra) # Version 0.8-5
library(DataExplorer) # Version 0.8.2
library(RColorBrewer) # Version 1.1-3
library(Hmisc) # Version 5.1-1
library(corrplot) # Version 0.92
library(plyr) # Version 1.8.9



# ---------------------------------------------------------------------------- #
# Import the data and formatting
# ---------------------------------------------------------------------------- #



# Open df_fig_data which is a formatted file combining fig- and tree-level observations.
# To read the raw data and format them, refer to the R script: R_script_Piatscheck.et.al.2020_Fpet_data_formating.R

df_fig_data <- read.csv(paste(workdir,"generated_data/","df_fig_data.csv", sep = ""))
df_tree_data <- read.csv(paste(workdir,"generated_data/","df_tree_data.csv", sep = ""))

# While the statistical analyses necessary to test our hypotheses are done at the fig (syconium) level, we will look at tree-level variable with df_tree_data.csv.


# ____________________________________________________________________________ #
# Getting the data frame formatted


df_fig_data$site <- factor(df_fig_data$site, levels = c("158", "172","112","113","95","179","201","96","70", "250"))
df_fig_data$season <- factor(df_fig_data$season, levels = c("F2012", "S2013","F2013","S2014"))
df_tree_data$site <- factor(df_tree_data$site, levels = c("158", "172","112","113","95","179","201","96","70", "250"))
df_tree_data$season <- factor(df_tree_data$season, levels = c("F2012", "S2013","F2013","S2014"))


# ---------------------------------------------------------------------------- #
# Explore the data set
# ---------------------------------------------------------------------------- #



# Wasp response variables --------------------

# The variables of interests are:
# - each single wasp species counts
# - a two column matrix of both the counts of pollinators and parasites.



#Visualize
df_fig_data_reponse_plot <- subset(df_fig_data, select=c(season,site,pollinators,LO1_f,SO1_f,SO2_f,heterandrium_1,heterandrium_2,ficicola,physothorax,sycophila,wasps_total))
df_fig_data_reponse_plot_m <- melt(df_fig_data_reponse_plot, id.vars = c("season", "site"))

ggplot(df_fig_data_reponse_plot_m, aes(x = value, y = variable, fill = variable, colour = variable)) + 
  geom_density_ridges(alpha = 0.8) + 
  geom_boxplot(alpha = 0.1, outlier.alpha = 0.2) +
# facet_wrap(season ~ ., scales="free_x",  ncol=1) + # Add if wanting to visualize for each collection trip.
  xlab(NULL) +
  ylab(NULL) +
  xlim(c(0,100)) +
  scale_colour_manual(values = rainbow(10),  guide = "none") +
  scale_fill_manual(values = rainbow(10),  name="Wasp species", 
                    breaks=levels(df_fig_data_reponse_plot_m$variable), 
                    labels=c("Pegoscapus sp.","I. flavicolis sp.1","I. carme sp.1","I. carme sp.2","Heterandrium sp. 1","Heterandrium sp. 2","Ficicola sp.","Physothorax sp.","Sycophila sp.","Wasp Total"))

rm(df_fig_data_reponse_plot, df_fig_data_reponse_plot_m )


# Host-related response variables --------------------

# The variables of interests are:
# - foundress number
# - crop size (product of tree volume estimation and reproductive effort estimation)
# - asynchrony
# - flowering and syconium landscape (for which we have several variables)

# Format
df_tree_data_host_variable_plot <- subset(df_tree_data, select=c(season, site, tree_volume,crop_size,asynchrony,flowering_landscape_250_125,syconium_landscape_250_125)) # We chose flowering_landscape_250_50 an syconium_landscape_250_50 based on preliminary analyses and glmm best fits.
df_tree_data_host_variable_plot_m <- melt(df_tree_data_host_variable_plot, id.vars = c("season", "site"))

# Plot
ggplot(df_tree_data_host_variable_plot_m, aes(x = variable, y=value, fill=variable, colour=variable)) +
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5) +
  theme_bw() + 
  theme(legend.position = "bottom") +
  facet_wrap( ~ variable, scales="free", ncol = 3)
# The plot allows us to spot some outliers in the crop size variable.
# This are probably due to overestimated tree volumes (and consequently large crop sizes when these where flowering).
# These outliers will be problematic for downstream analyses, so we set these values to NA.
df_tree_data[which(df_tree_data$tree_volume > 6000),][,c(1:3)] # These rows are the trees that have outliers.
# Replace outliers with NA.
df_tree_data[rownames(df_tree_data[which(df_tree_data$tree_volume > 6000),]),c(which(colnames(df_tree_data) == "tree_volume"), which(colnames(df_tree_data) == "crop_size"))] <- NA

# We re-run the code above to plot the variables without the outliers.
df_tree_data_host_variable_plot <- subset(df_tree_data, select=c(season, site, tree_volume,crop_size,asynchrony,flowering_landscape_250_125,syconium_landscape_250_125)) # We chose flowering_landscape_250_50 an syconium_landscape_250_50 based on preliminary analyses and glmm best fits.
df_tree_data_host_variable_plot_m <- melt(df_tree_data_host_variable_plot, id.vars = c("season", "site"))
ggplot(df_tree_data_host_variable_plot_m, aes(x = variable, y=value, fill=variable, colour=variable)) +
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5) +
  theme_bw() + 
  theme(legend.position = "bottom") +
  facet_wrap( ~ variable, scales="free", ncol = 3)
# Removing these rows helps resolving an outliers problem identified in preliminary analyses. 
# While there seems still to be large values in syconium landscape (syconium_landscape_250_50), it is harder to judge what threshold should be used and we will let the variable as it is for now.


# Abiotic response variables --------------------

df_tree_data_abiotic_variable_plot <- subset(df_tree_data, select=c(season,site,DM_tmax,DM_tmin,DM_prec,DM_vp,CHIRPS_prec,MODIS_LST,TC_tmax,TC_tmin,TC_prec,TC_vp,TC_ws))
# These climate estimates are scaled at the site level, so we remove duplicated colums corresponding to each trees observed.
df_tree_data_abiotic_variable_plot <- df_tree_data_abiotic_variable_plot[!duplicated(df_tree_data_abiotic_variable_plot), ]
df_tree_data_abiotic_variable_plot_m <- melt(df_tree_data_abiotic_variable_plot, id.vars = c("season", "site"))

ggplot(df_tree_data_abiotic_variable_plot_m, aes(x = variable, y=value, fill=variable, colour=variable)) +
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5) +
  theme_bw() + 
  theme(legend.position = "bottom") +
  facet_wrap( ~ variable, scales="free", ncol = 4)

# Clear working environment.
rm(df_tree_data_host_variable_plot, df_tree_data_host_variable_plot_m, df_tree_data_abiotic_variable_plot, df_tree_data_abiotic_variable_plot_m)



# ---------------------------------------------------------------------------- #
# Exploring missing data
# ---------------------------------------------------------------------------- #



# The statistical modeling will be done at the fig level, at which the fig wasp community is counted.
# Thus, we explore missing data in this data set only.

# First we remove the ouliers identified above in the fig data set.
df_fig_data[rownames(df_fig_data[which(df_fig_data$tree_volume > 6000),]),c(which(colnames(df_fig_data) == "tree_volume"), which(colnames(df_fig_data) == "crop_size"))] <- NA

# With base R
p_missing <- unlist(lapply(df_fig_data, function(x) sum(is.na(x))))/nrow(df_fig_data)
sort(p_missing[p_missing > 0], decreasing = TRUE) # Code from https://data.library.virginia.edu/getting-started-with-multiple-imputation-in-r/

# With DataExplorer
plot_missing(df_fig_data[,-c(which(colnames(df_fig_data) == "season"),
                             which(colnames(df_fig_data) == "site"),
                             which(colnames(df_fig_data) == "tree"),
                             which(colnames(df_fig_data) == "sampling_date"),
                             which(colnames(df_fig_data) == "fig"),
                             which(colnames(df_fig_data) == "latitude"),
                             which(colnames(df_fig_data) == "longitude"),
                             which(colnames(df_fig_data) == "elevation"),
                             which(colnames(df_fig_data) == "site_latitude"),
                             which(colnames(df_fig_data) == "site_longitude"),
                             which(colnames(df_fig_data) == "site_elevation"))])

# Individual trees elevation has a lot of missing data (these are missing GPS data), but this variable is not used for the analyses.
# Crop size, (which is the product of reproduction and tree volume) and asynchrony can be problematic (more than 10% missing data).

# Problematic variables of interest:
# crop size (due to tree size and reproduction missing data)
# asynchrony

# Other variable where missing data could be imputed
# flowering landscape
# syconium landscape
# foundress number

# If we are using all variables and omitting missing data, we lose:
dim(df_fig_data) - dim(na.omit(subset(df_fig_data, select=c(foundress_number,crop_size,asynchrony,flowering_landscape_250_125,syconium_landscape_250_125,
                                                            DM_tmax,DM_tmin,DM_prec,DM_vp,CHIRPS_prec,MODIS_LST,TC_tmax,TC_tmin,TC_prec,TC_vp,TC_ws))))
# 969 observations over 2304 observation. This represent a lose of 42% of the observations.
rm(p_missing)



# ---------------------------------------------------------------------------- #
# Pre-model fitting tests
# ---------------------------------------------------------------------------- #



# ____________________________________________________________________________ #
# Zero inflation in the response variables


# Check for zero-inflation in the response variables with Poisson distribution.
df_fig_data_reponse_zero_t <- subset(df_fig_data, select=c(pollinators, LO1_f,SO1_f,SO2_f,heterandrium_1,heterandrium_2,ficicola,physothorax,sycophila,parasites,wasps_total))


# We use the command zero.test() of the vcdExtra package. 
# This performs a zero test that has a Chi-square distribution, thus p-values below 0.05 suggest zero-inflation.
# See documentation for details and references.
 
zero_res <- c()
for(i in 1:length(colnames(df_fig_data_reponse_zero_t))) {
  tmp <- paste(colnames(df_fig_data_reponse_zero_t)[i], zero.test(df_fig_data_reponse_zero_t[,i])$prob, sep="\t\t")
  zero_res <- paste(zero_res, tmp, sep="\n")
}
cat(paste0("variable\t\t", "p-value", zero_res)) 
# Only the total of the wasp count is not zero-inflated
# In our statistical models, we will first assume zero-inflation in the . 

# Clear working environment
rm(df_fig_data_reponse_zero_t, zero_res, tmp, i)


# ____________________________________________________________________________ #
# Deviation of environmental variables from theoretical distribution tests


df_fig_data_qq <- subset(df_fig_data, select=c(foundress_number,tree_volume,reproduction,crop_size,asynchrony,flowering_landscape_250_125,syconium_landscape_250_125,
                                               DM_tmax, DM_tmin,DM_prec,DM_vp,CHIRPS_prec,MODIS_LST,TC_tmax,TC_tmin,TC_prec,TC_vp,TC_ws))

# We first look at the environmental variables's distributions.
par(mfrow= c(5,4)); for (i in colnames(df_fig_data_qq)){
plot(density(na.omit(df_fig_data_qq[,i])))
} # Zoom to look at distributions.

# QQ plots.
plot_qq(df_fig_data_qq[,1:9])
plot_qq(df_fig_data_qq[,10:length(df_fig_data_qq)])

# We use Shapiro-Wilk's method to test deviation from a normal distribution.
# p-values below 0.05 indicate for a deviation from normality.
shapiro_res <- c()
for(i in 1:length(colnames(df_fig_data_qq))) {
  tmp <- paste(colnames(df_fig_data_qq)[i], shapiro.test(df_fig_data_qq[,i])$p.value, sep="\t\t")
  shapiro_res <- paste(shapiro_res, tmp, sep="\n")
}
cat(paste0("variable\t\t", "p-value", shapiro_res))

# All explanatory variables deviate from a normal distribution.
# This suggest that transformation is necessary for GLMMs, and/or Glaussian family is not recommended.
# Spearman's rank correlation test is also advised.

# Clear working environment 
rm(df_fig_data_qq, shapiro_res, tmp, i)
dev.off()



# ____________________________________________________________________________ #
# Spearman correlation between explanatory variables



# Host-related variables --------------------

# Select the variables
df_fig_data_hv_corr <- subset(df_fig_data, select=c(foundress_number,crop_size,asynchrony,flowering_landscape_250_125,syconium_landscape_250_125))
df_fig_data_hv_corr <- na.omit(df_fig_data_hv_corr) # Remove missing data, which reduces the data set size.


# Correlation calculation and significance test
variable_hv_corr <- cor(df_fig_data_hv_corr, method="spearman") 
colnames(variable_hv_corr) <- c("Foundresses", "Crop size", "Asynchrony", "Flowering landscape", "Syconium landscape")
rownames(variable_hv_corr) <- c("Foundresses", "Crop size", "Asynchrony", "Flowering landscape", "Syconium landscape")
diag(variable_hv_corr) = NA
# variable_hv_corr_test <- cor.mtest(df_fig_data_hv_corr, method="pearson") # I fail to tun this test with the spearman method. 
variable_hv_corr_test_sp <- rcorr(as.matrix(variable_hv_corr), type = "spearman") # the function rcorr() resolve the issue above.

# Visualize significant correlations
corrplot(variable_hv_corr, method = "color", type = "full", 
         tl.col="black", tl.srt=35, order = NULL, diag = TRUE, na.label = "NA",
         addCoef.col = "darkgrey", sig.level = 0.05, insig = "pch", p.mat = variable_hv_corr_test_sp$P) # change to p.mat = variable_hv_corr_test$P
# Significant correlation exists but are too low to consider removed one of the variable (i.e. correlation >0.7).


# Abiotic variables --------------------


df_fig_data_av_corr <- subset(df_fig_data, select=c(DM_tmax, DM_tmin,DM_prec,DM_vp,CHIRPS_prec,MODIS_LST,TC_tmax,TC_tmin,TC_prec,TC_vp,TC_ws))
df_fig_data_av_corr <- na.omit(df_fig_data_av_corr) # Remove missing data, which reduces the data set size.

# Correlation calculation and significance test
variable_av_corr <- cor(df_fig_data_av_corr, method="spearman") 
colnames(variable_av_corr) <- c("DM T.max.", "DM T.min.","DM Prec.","DM V.P.","CHIRPS Prec.","MODIS LST","TC T.max.","TC T.min.","TC Prec.","TC V.P.","TC W.S.")
rownames(variable_av_corr) <- c("DM T.max.", "DM T.min.","DM Prec.","DM V.P.","CHIRPS Prec.","MODIS LST","TC T.max.","TC T.min.","TC Prec.","TC V.P.","TC W.S.")
diag(variable_av_corr) = NA
variable_av_corr_test_sp <- rcorr(as.matrix(variable_av_corr), type = "spearman")

# Visualize significant correlations
corrplot(variable_av_corr, method = "color", type = "full", 
         tl.col="black", tl.srt=35, order = NULL, diag = TRUE, na.label = "NA",
         addCoef.col = "darkgrey", sig.level = 0.05, insig = "pch", p.mat = variable_av_corr_test_sp$P) # change to p.mat = variable_hv_corr_test$P
# Significant correlations exist and are strong (i.e., ~ and > 0.7)



# All explanatory variables --------------------

# Select the variables
df_fig_data_all_corr <- subset(df_fig_data, select=c(foundress_number,crop_size,asynchrony,flowering_landscape_250_125,syconium_landscape_250_125, 
                                                    DM_tmax, DM_tmin,DM_prec,DM_vp,CHIRPS_prec,MODIS_LST,TC_tmax,TC_tmin,TC_prec,TC_vp,TC_ws))
df_fig_data_all_corr <- na.omit(df_fig_data_all_corr) # Remove missing data, which reduces the data set size.


# Correlation calculation and significance test
variable_all_corr <- cor(df_fig_data_all_corr, method="spearman") 
colnames(variable_all_corr) <- c("Foundresses", "Crop size", "Asynchrony", "Flowering landscape", "Syconium landscape", "DM T.max.", "DM T.min.","DM Prec.","DM V.P.","CHIRPS Prec.","MODIS LST","TC T.max.","TC T.min.","TC Prec.","TC V.P.","TC W.S.")
rownames(variable_all_corr) <- c("Foundresses", "Crop size", "Asynchrony", "Flowering landscape", "Syconium landscape", "DM T.max.", "DM T.min.","DM Prec.","DM V.P.","CHIRPS Prec.","MODIS LST","TC T.max.","TC T.min.","TC Prec.","TC V.P.","TC W.S.")
diag(variable_all_corr) = NA
# variable_hv_corr_test <- cor.mtest(df_fig_data_hv_corr, method="pearson") # I fail to tun this test with the spearman method. 
variable_all_corr_test_sp <- rcorr(as.matrix(variable_all_corr), type = "spearman") # the function rcorr() resolve the issue above.

# Visualize significant correlations
corrplot(variable_all_corr, method = "color", type = "full", 
         tl.col="black", tl.srt=35, order = NULL, diag = TRUE, na.label = "NA",
         addCoef.col = "darkgrey", sig.level = 0.05, insig = "pch", p.mat = variable_all_corr_test_sp$P) # change to p.mat = variable_hv_corr_test$P
# Strong and significant correlations concern mainly abiotic variables.
# Now let's remember that some of these variables are expected to be correlated, and that daily and monthly climate estimates will be analyzed separately.

# ____________________________________________________________________________ #
# Variance inflation factors


# Before variable selection --------------------

# Variance inflation factors calculation.
# Values > 10 (or 3 depending on references) are problematic. 
vif <- diag(solve(cor(df_fig_data_all_corr, method="spearman")))
vif
# Some abiotic variables are an issue: vapor pressure, max. and min. temperatures, precipitations.
# We know from previous analysis that LST is probably not a good variable to use.
# Is is also highly correlated with precipitation variables.
# We are not much interested in T. minimum and it is correlated to maximum temperature.
# We know that vapor pressure is often correlated with T. minimum and precipitations.
# Also, we will analyse daily and monthly abiotic variables separatly for consistency.
# Thus, we remove the following variables:
df_fig_data_day_corr <- df_fig_data_all_corr[,-c(which(colnames(df_fig_data_all_corr) == "MODIS_LST"),
                                                 which(colnames(df_fig_data_all_corr) == "DM_tmin"),
                                                 which(colnames(df_fig_data_all_corr) == "TC_tmax"),
                                                 which(colnames(df_fig_data_all_corr) == "TC_tmin"),
                                                 which(colnames(df_fig_data_all_corr) == "TC_prec"),
                                                 which(colnames(df_fig_data_all_corr) == "TC_vp"),
                                                 which(colnames(df_fig_data_all_corr) == "TC_ws"))]
vif <- diag(solve(cor(df_fig_data_day_corr, method="spearman")))
vif
# Vapor pressure could be an issue. We know that precipitation variables will show multicolinearity.


df_fig_data_month_corr <- df_fig_data_all_corr[,-c(which(colnames(df_fig_data_all_corr) == "MODIS_LST"),
                                                 which(colnames(df_fig_data_all_corr) == "TC_tmin"),
                                                 which(colnames(df_fig_data_all_corr) == "DM_tmax"),
                                                 which(colnames(df_fig_data_all_corr) == "DM_tmin"),
                                                 which(colnames(df_fig_data_all_corr) == "DM_prec"),
                                                 which(colnames(df_fig_data_all_corr) == "DM_vp"),
                                                 which(colnames(df_fig_data_all_corr) == "CHIRPS_prec"))]

vif <- diag(solve(cor(df_fig_data_month_corr, method="spearman")))
vif
# Vapor pressure and T. maximum are problematic, as they are highly correlated in the monthly estimates data set.
# We therefore decide to ignore vapor pressure and remove CHIRPS to compute VIF (CHIRPS will be used as a variable instead of DM_prec):


# After variable selection --------------------

df_fig_data_day_corr <- df_fig_data_all_corr[,-c(which(colnames(df_fig_data_all_corr) == "MODIS_LST"),
                                                 which(colnames(df_fig_data_all_corr) == "DM_tmin"),
                                                 which(colnames(df_fig_data_all_corr) == "DM_vp"),
                                                 which(colnames(df_fig_data_all_corr) == "TC_tmax"),
                                                 which(colnames(df_fig_data_all_corr) == "TC_tmin"),
                                                 which(colnames(df_fig_data_all_corr) == "TC_prec"),
                                                 which(colnames(df_fig_data_all_corr) == "TC_vp"),
                                                 which(colnames(df_fig_data_all_corr) == "TC_ws"),
                                                 which(colnames(df_fig_data_all_corr) == "CHIRPS_prec"))]
vif <- diag(solve(cor(df_fig_data_day_corr, method="spearman")))
vif

# With CHIRPS precipitation
df_fig_data_day_corr <- df_fig_data_all_corr[,-c(which(colnames(df_fig_data_all_corr) == "MODIS_LST"),
																								 which(colnames(df_fig_data_all_corr) == "DM_tmin"),
																								 which(colnames(df_fig_data_all_corr) == "DM_vp"),
																								 which(colnames(df_fig_data_all_corr) == "TC_tmax"),
																								 which(colnames(df_fig_data_all_corr) == "TC_tmin"),
																								 which(colnames(df_fig_data_all_corr) == "TC_prec"),
																								 which(colnames(df_fig_data_all_corr) == "TC_vp"),
																								 which(colnames(df_fig_data_all_corr) == "TC_ws"),
																								 which(colnames(df_fig_data_all_corr) == "DM_prec"))]
vif <- diag(solve(cor(df_fig_data_day_corr, method="spearman")))
vif
# No worrying VIF values anymore.

df_fig_data_month_corr <- df_fig_data_all_corr[,-c(which(colnames(df_fig_data_all_corr) == "MODIS_LST"),
                                                   which(colnames(df_fig_data_all_corr) == "TC_tmin"),
                                                   which(colnames(df_fig_data_all_corr) == "TC_vp"),
                                                   which(colnames(df_fig_data_all_corr) == "DM_tmax"),
                                                   which(colnames(df_fig_data_all_corr) == "DM_tmin"),
                                                   which(colnames(df_fig_data_all_corr) == "DM_prec"),
                                                   which(colnames(df_fig_data_all_corr) == "DM_vp"),
                                                   which(colnames(df_fig_data_all_corr) == "CHIRPS_prec"))]

vif <- diag(solve(cor(df_fig_data_month_corr, method="spearman")))
vif
# No worrying VIF values anymore.
# The selected predictive variables here do not seem to present issues that would lead in multicolinearity problems.

# From the pre-model fitting analyses done in the script above, we conclude:
# - Host-related variables contain many missing data that should be imputed.
# - Non-Gaussian families and zero-inflation need to be concidered in model fitting.
# - Daymet precipitation and CHIRPS precipitation will be used alternatively in model fitting.
# - Some abiotic variables will be dropped:
#     - MODIS LST
#     - Vapor pressure



# Clear working environment
rm(df_fig_data_hv_corr, variable_hv_corr, variable_hv_corr_test_sp, 
   df_fig_data_av_corr, variable_av_corr, variable_av_corr_test_sp,
   df_fig_data_all_corr, variable_all_corr, variable_all_corr_test_sp, 
   df_fig_data_day_corr, df_fig_data_month_corr,
   vif)



# ---------------------------------------------------------------------------- #
# End of script 5
# ---------------------------------------------------------------------------- #



# Clear working environment 
rm(list = ls())
dev.off()
