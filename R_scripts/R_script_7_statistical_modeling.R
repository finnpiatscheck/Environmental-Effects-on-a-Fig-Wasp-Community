# ============================================================================ #
# R script for reproduction of the statistical analyses
# ============================================================================ #



# Article title: Landscape-Level Analysis of a Fig-Pollinator-Antagonist Community: Spatial and Temporal Variation in a Fig Wasp Community and its Response to Biotic and Abiotic Factors
# By Finn Piatscheck, Justin van Goor, Derek Houston and John Nason.


# This script allows to reproduce the statistical analyses performed in the study.
# The methods used here are multivariate modeling and generalized linear mixed models.
# A few outliers are removed and empty figs (i.e. D phase syconia collected too early with wasp on the brink of emerging but not emerged yet) were previously removed.
# The plots realized in this script are not the ones presented in the article.
# For publication quality figures, refer to the R script 8: figure reproduction.



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

library(mice) # Version 3.16.0
library(RRPP) # Version 1.4.0
library(lme4) # Version 1.1-34
library(glmmTMB) # Version 1.1.7
library(PLNmodels) # version 1.0.4
library(broom.mixed) # Version 0.2.9.4
library(tidyverse) # Version 2.0.0
library(performance) # Version 0.10.8
library(bbmle) # Version 1.0.25
library(DHARMa) # Version 0.4.6
library(corrplot) # Version 0.92



# ---------------------------------------------------------------------------- #
# Import data and formatting
# ---------------------------------------------------------------------------- #



# Open df_fig_data.csv and df_tree_data.csv which are formatted files combining fig- and tree-level observations.
# To read the raw data and format them, refer to the R script: R_script_Piatscheck.et.al.2020_Fpet_data_formating.R

df_fig_data <- read.csv(paste(workdir,"generated_data/","df_fig_data.csv", sep = ""))
# df_fig_data is the data set that contains both fig and tree- level observations. 
# The statistical analyses are realized at the fig level.
# We will not use the tree-level data set here.


df_fig_analysis <- df_fig_data
# This data set will be modified and variable transformed for downstream analyses



# ____________________________________________________________________________ #
# Data frame formatting


# Format variables
df_fig_analysis$site <- factor(df_fig_analysis$site, levels = c("158", "172","112","113","95","179","201","96","70", "250"))
df_fig_analysis$site <- droplevels(df_fig_analysis$site)
df_fig_analysis$season <- factor(df_fig_analysis$season, levels = c("F2012", "S2013","F2013","S2014"))
df_fig_analysis$tree <- factor(df_fig_analysis$tree)
df_fig_analysis$fig <- factor(df_fig_analysis$fig)
df_fig_analysis$figID <- as.factor(paste(df_fig_analysis$site, df_fig_analysis$tree, df_fig_analysis$fig, sep = "-")) 
df_fig_analysis$treeID <- as.factor(paste(df_fig_analysis$site, df_fig_analysis$tree, sep = "-")) 
df_fig_analysis$wasp_resp_mat <- cbind(df_fig_analysis$pollinators, df_fig_analysis$parasites)

# Replace outliers with NAs.
df_fig_analysis[rownames(df_fig_data[which(df_fig_data$tree_volume > 6000),]),c(which(colnames(df_fig_data) == "tree_volume"), which(colnames(df_fig_data) == "crop_size"))] <- NA
which(df_fig_data$tree_volume > 6000)
which(df_fig_analysis$tree_volume > 6000) # Outliers removed


# Subset abundances for fig wasp interaction analysis
wasp_abund <- df_fig_analysis[,c("pollinators","LO1_f","SO1_f","SO2_f","heterandrium_1","heterandrium_2","ficicola","physothorax","sycophila")] 
colnames(wasp_abund) <- c("Pegoscapus sp.","Idarnes flavicolis sp.","Idarnes carme sp. 1","Idarnes carme sp. 2","Heterandrium sp. 1","Heterandrium sp. 2","Ficicola sp.","Physothorax sp.","Sycophila sp.")



# ---------------------------------------------------------------------------- #
# Standardization
# ---------------------------------------------------------------------------- #



# ____________________________________________________________________________ #
# Standardization: zero-centered

df_fig_analysis_zero <- df_fig_analysis

# Standardizing variables around 0.
stdize_zero_centered = function(x, ...) {(x - mean(x, ...)) / sd(x, ...)}

colnames(df_fig_analysis_zero)
# Scaling explanatory variables around 0.
df_fig_analysis_zero[,c("foundress_number","crop_size","asynchrony","flowering_landscape_250_125","syconium_landscape_250_125","DM_tmax","DM_prec","CHIRPS_prec","TC_tmax","TC_prec","TC_ws")] <- lapply(df_fig_analysis_zero[,c("foundress_number","crop_size","asynchrony","flowering_landscape_250_125","syconium_landscape_250_125","DM_tmax","DM_prec","CHIRPS_prec","TC_tmax","TC_prec","TC_ws")] , stdize_zero_centered, na.rm = T)
# The wasp counts won't be standardized yet because response variables should not be transformed. 

# Check the variables' range.
numdata<-df_fig_analysis_zero[sapply(df_fig_analysis_zero, is.numeric)]
sapply(numdata[,c("foundress_number","crop_size","asynchrony","flowering_landscape_250_125","syconium_landscape_250_125","DM_tmax","DM_prec","CHIRPS_prec","TC_tmax","TC_prec","TC_ws")], range, na.rm = T)

# Rename data frame for analyses
df_fig_analysis <- df_fig_analysis_zero

# Clear working environment
rm(numdata, df_fig_analysis_zero)



# ____________________________________________________________________________ #
# Impute missing data with the method outlined in script 5

# PMM Imputatio,
df_fig_data_to_be_imputed <- subset(df_fig_analysis, select = c(season,site,treeID,wasps_total ,pollinators,parasites,foundress_number,crop_size,asynchrony,flowering_landscape_250_125,syconium_landscape_250_125,flowering_landscape_250_125,syconium_landscape_250_125,DM_tmax,DM_tmin,DM_prec,DM_vp,CHIRPS_prec,MODIS_LST,TC_tmax,TC_tmin,TC_prec,TC_vp,TC_ws))
# mice doesn't handle matrix type columns, so we will recreate it calling GlmmTMB


imputed_fig_data_pmm <- mice(df_fig_data_to_be_imputed, m=50, maxit = 50, method = 'pmm')

# Diagnostic
densityplot(imputed_fig_data_pmm) # Check if imputed data are not way off the distribution of the observed data
stripplot(imputed_fig_data_pmm) # Similar but with dots
# checking the summary.
imputed_fig_data_pmm



# ---------------------------------------------------------------------------- #
# Hypotheses 1 test
# ---------------------------------------------------------------------------- #



# Hypothesis 1 (concise): fig wasp community structure are significantly variable across space and time


# ____________________________________________________________________________ #
# Linear model with multivariate data using a randomized residual permutation procedure 


# Prepare data set for multivariate analyses.
df_fig_RRPP <- rrpp.data.frame(logY = log(as.matrix(df_fig_analysis[, c("pollinators", "LO1_f", "SO1_f", "SO2_f", "heterandrium_1", "heterandrium_2", "ficicola", "physothorax", "sycophila")]) + 1), 
                         season = df_fig_data$season,
                         site = df_fig_data$site)

# Resolve the issue due to lack of data in all the levels (i.e. lack of D phase syconia during some collection).
# Using interaction() to create groups and pairwise() will result in an error telling that not all levels have data.
# The following works.
gps.all <- as.factor(paste(df_fig_data$season,df_fig_data$site,sep = "-"))

# Run multivariate linear model
wasp_abund_fit <- lm.rrpp(logY ~ gps.all, data = df_fig_RRPP, iter = 9999, SS.type = "I", print.progress = TRUE)
anova(wasp_abund_fit)
summary(wasp_abund_fit)

# Supplementary Table 4 --------------------

write.csv(anova(wasp_abund_fit)$table, paste(workdir,"tables/","wasp_communities_multivar_anal_anova_table.csv", sep = ""), row.names = TRUE)


# Visualization
plot(wasp_abund_fit, type = "PC", pch=21, bg = gps.all, cex=2)
legend("topright", levels(gps.all), pch = 21, pt.bg = 1:24)

# Post-hoc pairwise comparison
wasp_abund_pw <- pairwise(wasp_abund_fit, groups = gps.all)
pw_RRPP <- summary(wasp_abund_pw)
pw_RRPP
# How many non-significant differences between pairs?
length(which(pw_RRPP$summary.table$`Pr > d`> 0.05))
# How many significant differences between pairs?
length(which(pw_RRPP$summary.table$`Pr > d`< 0.05))

# Optional: use residuals from a null model (same conclusion than using residuals from the fit model)
# wasp_abund_null <- lm.rrpp(logY ~ 1, data = df_fig_RRPP, iter = 999, SS.type = "I", print.progress = TRUE)
# summary(pairwise(wasp_abund_fit, wasp_abund_null, groups = gps.all, print.progress = TRUE))


# Supplementary table X --------------------

write.csv(pw_RRPP$summary.table, paste(workdir,"tables/","wasp_communities_multivar_anal_posthoc_pw.csv", sep = ""), row.names = TRUE)


# Conclusions --------------------

# Fig wasps communities are significantly different among collections (sites by seasons) 
# All but one pairwise wasp community means are significantly different.



# ---------------------------------------------------------------------------- #
# Hypotheses 2 test
# ---------------------------------------------------------------------------- #



# Hypotheses 2: fig wasp abundance will be larger in: 
# - a) syconia with larger number of foundresses within 
# - b) trees with larger crops  
# - c) a larger flowering tree landscape
# - d) a larger syconium landscape
# - e) lower maximum temperature
# - f) lower precipitation
# - g) lower wind speed


# Model error distribution fitting testing --------------------

# Investigation of a better fit with alternatives to Poisson family
hyp2_poisson <- glmmTMB(wasps_total ~ foundress_number + crop_size + flowering_landscape_250_125 + syconium_landscape_250_125 + DM_tmax + DM_prec + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula=~0, family = poisson)
hyp2_nbinom1 <- glmmTMB(wasps_total ~ foundress_number + crop_size + flowering_landscape_250_125 + syconium_landscape_250_125 + DM_tmax + DM_prec + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula=~0, family = nbinom1)
hyp2_nbinom2 <- glmmTMB(wasps_total ~ foundress_number + crop_size + flowering_landscape_250_125 + syconium_landscape_250_125 + DM_tmax + DM_prec + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula=~0, family = nbinom2)
compare_performance(hyp2_poisson, hyp2_nbinom1, hyp2_nbinom2, metrics= c("AIC","BIC","RMSE"), rank = FALSE)


# Residual simulation and evaluation with DHARMa
sims_poisson <- simulateResiduals(hyp2_poisson, n = 10000)
plot(sims_poisson)
testResiduals(sims_poisson, plot = T)
testQuantiles(sims_poisson, plot = T) 
testZeroInflation(sims_poisson, plot = T)

sims_nbinom1 <- simulateResiduals(hyp2_nbinom1, n = 10000)
plot(sims_nbinom1)
testResiduals(sims_nbinom1, plot = T)
testQuantiles(sims_nbinom1, plot = T) 
testZeroInflation(sims_nbinom1, plot = T)

sims_nbinom2 <- simulateResiduals(hyp2_nbinom2, n = 10000)
plot(sims_nbinom2)
testResiduals(sims_nbinom2, plot = T)
testQuantiles(sims_nbinom2, plot = T) 
testZeroInflation(sims_nbinom2, plot = T)

# Although GLMMs fitted with negative binomial with quadratic parameterization provides a better fit, model evaluation leans toward selecting negative binomial with linear parameterization family.

# Clear memory
rm(hyp2_mod1_poisson, hyp2_mod1_nbinom1, hyp2_mod1_nbinom2, sims_poisson, sims_nbinom1, sims_nbinom2)


# Generalized linear mixed models with glmmTMB to test hypotheses 2 --------------------

# No interactions (models 1)
hyp2_mod1_DM <- glmmTMB(wasps_total ~ foundress_number + crop_size + flowering_landscape_250_125 + syconium_landscape_250_125 + DM_tmax + DM_prec + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula=~0, family = nbinom1)
hyp2_mod1_DM_CH <- glmmTMB(wasps_total ~ foundress_number + crop_size + flowering_landscape_250_125 + syconium_landscape_250_125 + DM_tmax + CHIRPS_prec  + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula=~0, family = nbinom1)
hyp2_mod1_TC <- glmmTMB(wasps_total ~ foundress_number + crop_size + flowering_landscape_250_125 + syconium_landscape_250_125 + TC_tmax + TC_prec + TC_ws  + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula=~0, family = nbinom1)
compare_performance(hyp2_mod1_DM, hyp2_mod1_DM_CH, hyp2_mod1_TC, metrics= c("AIC","BIC","RMSE"), rank = FALSE)

# Residual simulation and evaluation with DHARMa
hyp2_mod1_DM_sims <- simulateResiduals(hyp2_mod1_DM, n = 10000)
plot(hyp2_mod1_DM_sims)
testResiduals(hyp2_mod1_DM_sims, plot = T)
testQuantiles(hyp2_mod1_DM_sims, plot = T) 
testZeroInflation(hyp2_mod1_DM_sims, plot = T)
# The model performs relatively well, though residual vs fitted values are not uniform around the regression line.
# check post-fit collinearity
check_collinearity(hyp2_mod1_DM) # No issues (all VIF values < 3)

hyp2_mod1_DM_CH_sims <- simulateResiduals(hyp2_mod1_DM_CH, n = 10000)
plot(hyp2_mod1_DM_CH_sims)
testResiduals(hyp2_mod1_DM_CH_sims, plot = T)
testQuantiles(hyp2_mod1_DM_CH_sims, plot = T) 
testZeroInflation(hyp2_mod1_DM_CH_sims, plot = T)
# The model performs relatively well, though residual vs fitted values are not uniform around the regression line, less well than with DM data.
# check post-fit collinearity
check_collinearity(hyp2_mod1_DM_CH) # No issues (all VIF values < 3)

hyp2_mod1_TC_sims <- simulateResiduals(hyp2_mod1_TC, n = 10000)
plot(hyp2_mod1_TC_sims)
testResiduals(hyp2_mod1_TC_sims, plot = T)
testQuantiles(hyp2_mod1_TC_sims, plot = T) 
testZeroInflation(hyp2_mod1_TC_sims, plot = T)
# The model performs poorly.
# check post-fit collinearity
check_collinearity(hyp2_mod1_TC) # We note a VIF of 2.3 of wind speed. 


# Remove wind speed?
hyp2_mod1_TC_noWs <- glmmTMB(wasps_total ~ foundress_number + crop_size + flowering_landscape_250_125 + syconium_landscape_250_125 + TC_tmax + TC_prec  + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula=~0, family = nbinom1)
compare_performance(hyp2_mod1_DM, hyp2_mod1_DM_CH, hyp2_mod1_TC, hyp2_mod1_TC_noWs, metrics= c("AIC","BIC","RMSE"), rank = FALSE)
compare_performance(hyp2_mod1_DM, hyp2_mod1_DM_CH, hyp2_mod1_TC, hyp2_mod1_TC_noWs, metrics= c("AIC","BIC","RMSE"), rank = TRUE)
# No improvements

hyp2_mod1_TC_noWs_sims <- simulateResiduals(hyp2_mod1_TC_noWs, n = 10000)
plot(hyp2_mod1_TC_noWs_sims)
testResiduals(hyp2_mod1_TC_noWs_sims, plot = T)
testQuantiles(hyp2_mod1_TC_noWs_sims, plot = T) 
testZeroInflation(hyp2_mod1_TC_noWs_sims, plot = T)
# The model performs relatively well

# Now we consider interactions (models 2)
hyp2_mod2_DM <- glmmTMB(wasps_total ~ foundress_number * DM_tmax * DM_prec + crop_size * flowering_landscape_250_125 * syconium_landscape_250_125  + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula=~0, family = nbinom1)
hyp2_mod2_DM_CH <- glmmTMB(wasps_total ~ foundress_number * DM_tmax * CHIRPS_prec + crop_size * flowering_landscape_250_125 * syconium_landscape_250_125  + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula=~0, family = nbinom1)
hyp2_mod2_TC <- glmmTMB(wasps_total ~ foundress_number * TC_tmax * TC_prec * TC_ws + crop_size * flowering_landscape_250_125 * syconium_landscape_250_125  + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula=~0, family = nbinom1)
hyp2_mod2_TC_noWs <- glmmTMB(wasps_total ~ foundress_number * TC_tmax * TC_prec + crop_size * flowering_landscape_250_125 * syconium_landscape_250_125  + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula=~0, family = nbinom1)

# Residual simulation and evaluation with DHARMa
hyp2_mod2_DM_sims <- simulateResiduals(hyp2_mod2_DM, n = 10000)
plot(hyp2_mod2_DM_sims)
testResiduals(hyp2_mod2_DM_sims, plot = T)
testQuantiles(hyp2_mod2_DM_sims, plot = T) 
testZeroInflation(hyp2_mod2_DM_sims, plot = T)
# The model performs relatively well.
# check post-fit collinearity
check_collinearity(hyp2_mod2_DM) # High collinearity but this is expected in interaction models. Probably not relevant here.

hyp2_mod2_DM_CH_sims <- simulateResiduals(hyp2_mod2_DM_CH, n = 10000)
plot(hyp2_mod2_DM_CH_sims)
testResiduals(hyp2_mod2_DM_CH_sims, plot = T)
testQuantiles(hyp2_mod2_DM_CH_sims, plot = T) 
testZeroInflation(hyp2_mod2_DM_CH_sims, plot = T)
# The model performs poorly, residuals very skewed.
check_collinearity(hyp2_mod2_DM_CH) # High collinearity, same as above.

hyp2_mod2_TC_sims <- simulateResiduals(hyp2_mod2_TC, n = 10000)
plot(hyp2_mod2_TC_sims)
testResiduals(hyp2_mod2_TC_sims, plot = T)
testQuantiles(hyp2_mod2_TC_sims, plot = T) 
testZeroInflation(hyp2_mod2_TC_sims, plot = T)
# The model performs relatively well. Residuals quite skewed though.
check_collinearity(hyp2_mod2_TC) # High collinearity, same as above.

hyp2_mod2_TC_noWs_sims <- simulateResiduals(hyp2_mod2_TC_noWs, n = 10000)
plot(hyp2_mod2_TC_noWs_sims)
testResiduals(hyp2_mod2_TC_noWs_sims, plot = T)
testQuantiles(hyp2_mod2_TC_noWs_sims, plot = T) 
testZeroInflation(hyp2_mod2_TC_noWs_sims, plot = T)
# The model performs relatively well. Residuals quite skewed though.
check_collinearity(hyp2_mod2_TC_noWs) # High collinearity, same as above.

# Compare all models
compare_performance(hyp2_mod1_DM, hyp2_mod1_DM_CH, hyp2_mod1_TC, hyp2_mod1_TC_noWs, hyp2_mod2_DM, hyp2_mod2_DM_CH, hyp2_mod2_TC, hyp2_mod2_TC_noWs, metrics= c("AIC","BIC","RMSE"), rank = FALSE)
compare_performance(hyp2_mod1_DM, hyp2_mod1_DM_CH, hyp2_mod1_TC, hyp2_mod1_TC_noWs, hyp2_mod2_DM, hyp2_mod2_DM_CH, hyp2_mod2_TC, hyp2_mod2_TC_noWs, metrics= c("AIC","BIC","RMSE"), rank = TRUE)

# Results of the models
summary(hyp2_mod1_DM)
summary(hyp2_mod1_DM_CH)
summary(hyp2_mod1_TC)
summary(hyp2_mod1_TC_noWs)
summary(hyp2_mod2_DM)
summary(hyp2_mod2_DM_CH)
summary(hyp2_mod2_TC) # Best fit model
summary(hyp2_mod2_TC_noWs)

# Rerun GLMMs with best fit with 50 completed data sets --------------------

hyp2_mod1_DM_imputdata <- with(data=imputed_fig_data_pmm, glmmTMB(wasps_total ~ foundress_number + crop_size + flowering_landscape_250_125 + syconium_landscape_250_125 + DM_tmax + DM_prec + (1|season) + (1|site/treeID), ziformula=~0, family = nbinom1))
hyp2_mod1_DM_imputdata_pool_fit <- pool(hyp2_mod1_DM_imputdata)
summary(hyp2_mod1_DM_imputdata_pool_fit)

hyp2_mod2_TC_imputdata <- with(data=imputed_fig_data_pmm, glmmTMB(wasps_total ~ foundress_number * TC_tmax * TC_prec * TC_ws + crop_size * flowering_landscape_250_125 * syconium_landscape_250_125  + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula=~0, family = nbinom1))
hyp2_mod2_TC_imputdata_pool_fit <- pool(hyp2_mod2_TC_imputdata)
summary(hyp2_mod2_TC_imputdata_pool_fit)

# Conclusions  --------------------

# Larger foundress number in figs is associated with larger production of wasp.
# Trees with larger crop size have a larger wasp production.
# Effect of the neighborhood seems to not impact wasp production.
# Maximum temperature could positively affect wasp count
# Precipitation and wind speed is inconsistent.


# ---------------------------------------------------------------------------- #
# Hypotheses 3 test
# ---------------------------------------------------------------------------- #



# Hypotheses 3: Pollinators will be more abundant than parasites in
# - a) syconia with larger number of foundresses within 
# - b) synchronous trees
# - c) a larger flowering tree landscape
# - d) a larger syconium landscape
# - e) higher maximum temperature
# - f) higher precipitation
# - g) higher wind speed


# Model error distribution fitting testing --------------------

# Investigation of a better fit with alternatives to Poisson family
hyp3_binomial <- glmmTMB(wasp_resp_mat ~ foundress_number + crop_size + flowering_landscape_250_125 + syconium_landscape_250_125 + DM_tmax + DM_prec + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula= ~0, family = binomial(link = "logit"))
hyp3_betabinom <- glmmTMB(wasp_resp_mat ~ foundress_number + crop_size + flowering_landscape_250_125 + syconium_landscape_250_125 + DM_tmax + DM_prec + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula= ~0, family = betabinomial(link = "logit"))
# hyp3_betabinom_cloglog <- glmmTMB(wasp_resp_mat ~ foundress_number + crop_size + flowering_landscape_250_125 + syconium_landscape_250_125 + DM_tmax + DM_prec + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula= ~0, family = betabinomial(link = "cloglog"))
# compare_performance(hyp3_binomial, hyp3_betabinom, hyp3_betabinom_cloglog, metrics= c("AIC","BIC","RMSE"), rank = FALSE)
# compare_performance(hyp3_binomial, hyp3_betabinom, hyp3_betabinom_cloglog, metrics= c("AIC","BIC","RMSE"), rank = TRUE)
compare_performance(hyp3_binomial, hyp3_betabinom, metrics= c("AIC","BIC","RMSE"), rank = FALSE)
compare_performance(hyp3_binomial, hyp3_betabinom, metrics= c("AIC","BIC","RMSE"), rank = TRUE)


# Residual simulation and evaluation with DHARMa
sims_binomial <- simulateResiduals(hyp3_binomial, n = 10000)
plot(sims_binomial)
testResiduals(sims_binomial, plot = T)
testZeroInflation(sims_binomial, plot = T)

sims_betabinom <- simulateResiduals(hyp3_betabinom, n = 10000)
plot(sims_betabinom)
testResiduals(sims_betabinom, plot = T)
testZeroInflation(sims_betabinom, plot = T)
# Both models fail the residual tests.
# We will continue with betabinomial error distribution, as it is appropriated for over dispersed data.


# Generalized linear mixed models with glmmTMB to test hypotheses 3 --------------------

# No interactions (models 1)
hyp3_mod1_DM <- glmmTMB(wasp_resp_mat ~ foundress_number + asynchrony + flowering_landscape_250_125 + syconium_landscape_250_125 + DM_tmax + DM_prec + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula= ~0, family = betabinomial(link = "logit"))
hyp3_mod1_DM_CH <- glmmTMB(wasp_resp_mat ~ foundress_number + asynchrony + flowering_landscape_250_125 + syconium_landscape_250_125 + DM_tmax + CHIRPS_prec + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula= ~0, family = betabinomial(link = "logit"))
hyp3_mod1_TC <- glmmTMB(wasp_resp_mat ~ foundress_number + asynchrony + flowering_landscape_250_125 + syconium_landscape_250_125 + TC_tmax + TC_prec + TC_ws + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula= ~0, family = betabinomial(link = "logit"))
compare_performance(hyp3_mod1_DM, hyp3_mod1_DM_CH, hyp3_mod1_TC, metrics= c("AIC","BIC","RMSE"), rank = FALSE)
compare_performance(hyp3_mod1_DM, hyp3_mod1_DM_CH, hyp3_mod1_TC, metrics= c("AIC","BIC","RMSE"), rank = TRUE)

# Based on wind speed causing potential multicollinearity, we try removing wind speed
hyp3_mod1_TC_noWs <- glmmTMB(wasp_resp_mat ~ foundress_number + asynchrony + flowering_landscape_250_125 + syconium_landscape_250_125 + TC_tmax + TC_prec  + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula= ~0, family = betabinomial(link = "logit"))
compare_performance(hyp3_mod1_DM, hyp3_mod1_DM_CH, hyp3_mod1_TC, hyp3_mod1_TC_noWs, metrics= c("AIC","BIC","RMSE"), rank = FALSE)
compare_performance(hyp3_mod1_DM, hyp3_mod1_DM_CH, hyp3_mod1_TC, hyp3_mod1_TC_noWs, metrics= c("AIC","BIC","RMSE"), rank = TRUE)


# Residual simulation and evaluation with DHARMa
hyp3_mod1_DM_sims <- simulateResiduals(hyp3_mod1_DM, n = 10000)
plot(hyp3_mod1_DM_sims)
testResiduals(hyp3_mod1_DM_sims, plot = T)
testZeroInflation(hyp3_mod1_DM_sims, plot = T)
# The model performs poorly.
check_collinearity(hyp3_mod1_DM) # 

hyp3_mod1_DM_CH_sims <- simulateResiduals(hyp3_mod1_DM_CH, n = 10000)
plot(hyp3_mod1_DM_CH_sims)
testResiduals(hyp3_mod1_DM_CH_sims, plot = T)
testZeroInflation(hyp3_mod1_DM_CH_sims, plot = T)

hyp3_mod1_TC_sims <- simulateResiduals(hyp3_mod1_TC, n = 10000)
plot(hyp3_mod1_TC_sims)
testResiduals(hyp3_mod1_TC_sims, plot = T)
testZeroInflation(hyp3_mod1_TC_sims, plot = T)

hyp3_mod1_TC_noWs_sims <- simulateResiduals(hyp3_mod1_TC_noWs, n = 10000)
plot(hyp3_mod1_TC_noWs_sims)
testResiduals(hyp3_mod1_TC_noWs_sims, plot = T)
testZeroInflation(hyp3_mod1_TC_noWs_sims, plot = T)
# All models are failing at least one test.

summary(hyp3_mod1_DM)
summary(hyp3_mod1_DM_CH)
summary(hyp3_mod1_TC)
summary(hyp3_mod1_TC_noWs)

check_collinearity(hyp3_mod1_TC)
check_collinearity(hyp3_mod1_TC_noWs)
# Wind speed can be an issue.


# Now we consider interactions (models 2)
hyp3_mod2_DM <- glmmTMB(wasp_resp_mat ~ foundress_number * DM_tmax * DM_prec + asynchrony + flowering_landscape_250_125 * syconium_landscape_250_125 + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula= ~0, family = betabinomial(link = "logit"))
hyp3_mod2_DM_CH <- glmmTMB(wasp_resp_mat ~ foundress_number * DM_tmax * CHIRPS_prec + asynchrony + flowering_landscape_250_125 * syconium_landscape_250_125 + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula= ~0, family = betabinomial(link = "logit"))
hyp3_mod2_TC <- glmmTMB(wasp_resp_mat ~ foundress_number * TC_tmax * TC_prec * TC_ws + asynchrony + flowering_landscape_250_125 * syconium_landscape_250_125  + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula= ~0, family = betabinomial(link = "logit"))
hyp3_mod2_TC_noWs <- glmmTMB(wasp_resp_mat ~ foundress_number * TC_tmax * TC_prec + asynchrony + flowering_landscape_250_125 * syconium_landscape_250_125   + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula= ~0, family = betabinomial(link = "logit"))
compare_performance(hyp3_mod2_DM, hyp3_mod2_DM_CH, hyp3_mod2_TC,hyp3_mod2_TC_noWs, metrics= c("AIC","BIC","RMSE"), rank = FALSE)
compare_performance(hyp3_mod2_DM, hyp3_mod2_DM_CH, hyp3_mod2_TC,hyp3_mod2_TC_noWs, metrics= c("AIC","BIC","RMSE"), rank = TRUE)


hyp3_mod2_DM_sims <- simulateResiduals(hyp3_mod2_DM, n = 10000)
plot(hyp3_mod2_DM_sims)
testResiduals(hyp3_mod2_DM_sims, plot = T)
testQuantiles(hyp3_mod2_DM_sims, plot = T) 
testZeroInflation(hyp3_mod2_DM_sims, plot = T)

hyp3_mod2_DM_CH_sims <- simulateResiduals(hyp3_mod2_DM_CH, n = 10000)
plot(hyp3_mod2_DM_CH_sims)
testResiduals(hyp3_mod2_DM_CH_sims, plot = T)
testQuantiles(hyp3_mod2_DM_CH_sims, plot = T) 
testZeroInflation(hyp3_mod2_DM_CH_sims, plot = T)

hyp3_mod2_TC_sims <- simulateResiduals(hyp3_mod2_TC, n = 10000)
plot(hyp3_mod2_TC_sims)
testResiduals(hyp3_mod2_TC_sims, plot = T)
testQuantiles(hyp3_mod2_TC_sims, plot = T) 
testZeroInflation(hyp3_mod2_TC_sims, plot = T)

hyp3_mod2_TC_noWs_sims <- simulateResiduals(hyp3_mod2_TC_noWs, n = 10000)
plot(hyp3_mod2_TC_noWs_sims)
testResiduals(hyp3_mod2_TC_noWs_sims, plot = T)
testQuantiles(hyp3_mod2_TC_noWs_sims, plot = T) 
testZeroInflation(hyp3_mod2_TC_noWs_sims, plot = T)


summary(hyp3_mod2_DM)
summary(hyp3_mod2_DM_CH)
summary(hyp3_mod2_TC)
summary(hyp3_mod2_TC_noWs)

compare_performance(hyp3_mod2_DM, hyp3_mod2_DM_CH, hyp3_mod2_TC,hyp3_mod2_TC_noWs, metrics= c("AIC","BIC","RMSE"), rank = FALSE)
compare_performance(hyp3_mod2_DM, hyp3_mod2_DM_CH, hyp3_mod2_TC,hyp3_mod2_TC_noWs, metrics= c("AIC","BIC","RMSE"), rank = TRUE)



# Using with Binomial error distribution

# No interactions (models 1)
hyp3_mod1_DM_b <- glmmTMB(wasp_resp_mat ~ foundress_number + asynchrony + flowering_landscape_250_125 + syconium_landscape_250_125 + DM_tmax + DM_prec + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula= ~0, family = binomial(link = "logit"))
hyp3_mod1_DM_CH_b <- glmmTMB(wasp_resp_mat ~ foundress_number + asynchrony + flowering_landscape_250_125 + syconium_landscape_250_125 + DM_tmax + CHIRPS_prec + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula= ~0, family = binomial(link = "logit"))
hyp3_mod1_TC_b <- glmmTMB(wasp_resp_mat ~ foundress_number + asynchrony + flowering_landscape_250_125 + syconium_landscape_250_125 + TC_tmax + TC_prec + TC_ws + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula= ~0, family = binomial(link = "logit"))
hyp3_mod1_TC_noWs_b <- glmmTMB(wasp_resp_mat ~ foundress_number + asynchrony + flowering_landscape_250_125 + syconium_landscape_250_125 + TC_tmax + TC_prec  + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula= ~0, family = binomial(link = "logit"))
compare_performance(hyp3_mod1_DM_b, hyp3_mod1_DM_CH_b, hyp3_mod1_TC_b, hyp3_mod1_TC_noWs_b, metrics= c("AIC","BIC","RMSE"), rank = FALSE)
compare_performance(hyp3_mod1_DM_b, hyp3_mod1_DM_CH_b, hyp3_mod1_TC_b, hyp3_mod1_TC_noWs_b, metrics= c("AIC","BIC","RMSE"), rank = TRUE)

hyp3_mod1_DM_b_sims <- simulateResiduals(hyp3_mod1_DM_b, n = 10000)
plot(hyp3_mod1_DM_b_sims)
testResiduals(hyp3_mod1_DM_b_sims, plot = T)
testQuantiles(hyp3_mod1_DM_b_sims, plot = T) 
testZeroInflation(hyp3_mod1_DM_b_sims, plot = T)
# Model failed

hyp3_mod1_DM_CH_b_sims <- simulateResiduals(hyp3_mod1_DM_CH_b, n = 10000)
plot(hyp3_mod1_DM_CH_b_sims)
testResiduals(hyp3_mod1_DM_CH_b_sims, plot = T)
testQuantiles(hyp3_mod1_DM_CH_b_sims, plot = T) 
testZeroInflation(hyp3_mod1_DM_CH_b_sims, plot = T)
# Model failed

hyp3_mod1_TC_b_sims <- simulateResiduals(hyp3_mod1_TC_b, n = 10000)
plot(hyp3_mod1_TC_b_sims)
testResiduals(hyp3_mod1_TC_b_sims, plot = T)
testQuantiles(hyp3_mod1_TC_b_sims, plot = T) 
testZeroInflation(hyp3_mod1_TC_b_sims, plot = T)
# Only model that doesn't break model assumptions

hyp3_mod1_TC_noWs_b_sims <- simulateResiduals(hyp3_mod1_TC_noWs_b, n = 10000)
plot(hyp3_mod1_TC_noWs_b_sims)
testResiduals(hyp3_mod1_TC_noWs_b_sims, plot = T)
testQuantiles(hyp3_mod1_TC_noWs_b_sims, plot = T) 
testZeroInflation(hyp3_mod1_TC_noWs_b_sims, plot = T)
# Model failed

# With interactions
hyp3_mod2_DM_b <- glmmTMB(wasp_resp_mat ~ foundress_number * DM_tmax * DM_prec + asynchrony + flowering_landscape_250_125 * syconium_landscape_250_125 + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula= ~0, family = binomial(link = "logit"))
hyp3_mod2_DM_CH_b <- glmmTMB(wasp_resp_mat ~ foundress_number * DM_tmax * CHIRPS_prec + asynchrony + flowering_landscape_250_125 * syconium_landscape_250_125 + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula= ~0, family = binomial(link = "logit"))
hyp3_mod2_TC_b <- glmmTMB(wasp_resp_mat ~ foundress_number * TC_tmax * TC_prec * TC_ws + asynchrony + flowering_landscape_250_125 * syconium_landscape_250_125  + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula= ~0, family = binomial(link = "logit"))
hyp3_mod2_TC_noWs_b <- glmmTMB(wasp_resp_mat ~ foundress_number * TC_tmax * TC_prec + asynchrony + flowering_landscape_250_125 * syconium_landscape_250_125  + (1|season) + (1|site/treeID), data=df_fig_analysis, ziformula= ~0, family = binomial(link = "logit"))
compare_performance(hyp3_mod2_DM_b, hyp3_mod2_DM_CH_b, hyp3_mod2_TC_b, hyp3_mod2_TC_noWs_b, metrics= c("AIC","BIC","RMSE"), rank = FALSE)
compare_performance(hyp3_mod2_DM_b, hyp3_mod2_DM_CH_b, hyp3_mod2_TC_b, hyp3_mod2_TC_noWs_b, metrics= c("AIC","BIC","RMSE"), rank = TRUE)

compare_performance(hyp3_mod2_DM_b, hyp3_mod2_DM_CH_b, hyp3_mod2_TC_b, hyp3_mod2_TC_noWs_b, metrics= c("AIC","BIC","RMSE"), rank = FALSE)
compare_performance(hyp3_mod2_DM, hyp3_mod2_DM_CH, hyp3_mod2_TC, hyp3_mod2_TC_noWs, hyp3_mod2_DM_b, hyp3_mod2_DM_CH_b, hyp3_mod2_TC_b, hyp3_mod2_TC_noWs_b, metrics= c("AIC","BIC","RMSE"), rank = FALSE)
compare_performance(hyp3_mod2_DM, hyp3_mod2_DM_CH, hyp3_mod2_TC, hyp3_mod2_TC_noWs, hyp3_mod2_DM_b, hyp3_mod2_DM_CH_b, hyp3_mod2_TC_b, hyp3_mod2_TC_noWs_b, metrics= c("AIC","BIC","RMSE"), rank = TRUE)


# Interactions (models 2)
hyp3_mod2_DM_b_sims <- simulateResiduals(hyp3_mod2_DM_b, n = 10000)
plot(hyp3_mod2_DM_b_sims)
testResiduals(hyp3_mod2_DM_b_sims, plot = T)
testQuantiles(hyp3_mod2_DM_b_sims, plot = T) 
testZeroInflation(hyp3_mod2_DM_b_sims, plot = T)
# Model failed

hyp3_mod2_DM_CH_b_sims <- simulateResiduals(hyp3_mod2_DM_CH_b, n = 10000)
plot(hyp3_mod2_DM_CH_b_sims)
testResiduals(hyp3_mod2_DM_CH_b_sims, plot = T)
testQuantiles(hyp3_mod2_DM_CH_b_sims, plot = T) 
testZeroInflation(hyp3_mod2_DM_CH_b_sims, plot = T)
# Model failed

hyp3_mod2_TC_b_sims <- simulateResiduals(hyp3_mod2_TC_b, n = 10000)
plot(hyp3_mod2_TC_b_sims)
testResiduals(hyp3_mod2_TC_b_sims, plot = T)
testQuantiles(hyp3_mod2_TC_b_sims, plot = T) 
testZeroInflation(hyp3_mod2_TC_b_sims, plot = T)
# Only model that doesn't break model assumptions

hyp3_mod2_TC_noWs_b_sims <- simulateResiduals(hyp3_mod2_TC_noWs_b, n = 10000)
plot(hyp3_mod2_TC_noWs_b_sims)
testResiduals(hyp3_mod2_TC_noWs_b_sims, plot = T)
testQuantiles(hyp3_mod2_TC_noWs_b_sims, plot = T) 
testZeroInflation(hyp3_mod2_TC_noWs_b_sims, plot = T)
# Model failed

summary(hyp3_mod1_DM_b)
summary(hyp3_mod1_DM_CH_b)
summary(hyp3_mod1_TC_b)
summary(hyp3_mod1_TC_noWs_b)


summary(hyp3_mod2_DM_b)
summary(hyp3_mod2_DM_CH_b)
summary(hyp3_mod2_TC_b)
summary(hyp3_mod2_TC_noWs_b)




# Rerun GLMMs with best fit with 50 completed data sets --------------------


hyp3_mod2_TC_b_imputdata <- with(data=imputed_fig_data_pmm, glmmTMB(cbind(pollinators, parasites) ~ foundress_number * TC_tmax * TC_prec * TC_ws + asynchrony + flowering_landscape_250_125 * syconium_landscape_250_125  + (1|season) + (1|site/treeID), ziformula=~0, family = binomial(link = "logit")))
hyp3_mod2_TC_b_imputdata_pool_fit <- pool(hyp3_mod2_TC_b_imputdata)
summary(hyp3_mod2_TC_b_imputdata_pool_fit)

# Conclusions --------------------

# Larger foundress number increases the proportion of pollinators produced in a fig relative to parasites.
# Asynchrony is associated with increase production of pollinators
# Increased temperature and increased precipitation are associated with higher pollinator production.
# Interactions among abiotic variables and foudress counts.



# ---------------------------------------------------------------------------- #
# Fig wasp species interaction network inference with Poisson-Lognormal models
# ---------------------------------------------------------------------------- #



# Prepare data
wasp_abund_for_PLN <- prepare_data(counts = wasp_abund, covariates = df_fig_data[,c("site","season")], offset = "TSS")
str(wasp_abund_for_PLN)

# Run network models
wasp_network_offset <- PLNnetwork(Abundance ~ 1 + offset(log(Offset)), data = wasp_abund_for_PLN)

# Select best model
wasp_best_network_offset <- getBestModel(wasp_network_offset)
plot(wasp_best_network_offset)
plot(wasp_best_network_offset, type = "support", output = "corrplot")
# Done



# ---------------------------------------------------------------------------- #
# End of Appendix 6
# ---------------------------------------------------------------------------- #



# Clear working environment 
rm(list = ls())
dev.off()