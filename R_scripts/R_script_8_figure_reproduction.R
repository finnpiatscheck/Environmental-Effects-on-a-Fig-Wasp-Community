# ============================================================================ #
# R script 8: figure reproduction
# ============================================================================ #



# Article title: Landscape-Level Analysis of a Fig-Pollinator-Antagonist Community: Spatial and Temporal Variation in a Fig Wasp Community and its Response to Biotic and Abiotic Factors
# By Finn Piatscheck, Justin van Goor, Derek Houston and John Nason.


# The script allows the reproduction of the following figures, in this order:
# - Figure 4A
# - Supplementary Figure S13
# - Supplementary Figure S14
# - Supplementary Figure S16
# - Supplementary Figure S15
# - Figure 6
# - Supplementary Figure S18
# - Supplementary Figure S17
# - Figure 7
# - Figure 8
# - Figure 9
# - Supplementary Figure S5
# - Supplementary Figure S6
# - Supplementary Figure S7
# - Supplementary Figure S8
# - Supplementary Figure S9
# - Supplementary Figure S10
# - Supplementary Figure S11
# - Supplementary Figure S12
# - Figure 5
# - Figure 10


# Figures S2, S3 and S4 can be reproduced with R script "neighborhood_landscape_analysis.R".
# Figures 1, 2, 3, 4B and S1 were realized with inkscape.


# ---------------------------------------------------------------------------- #
# Load packages and set the working directory
# ---------------------------------------------------------------------------- #



# R version --------------------

# R version used for these analyses: 4.3.1.


# Set working directory --------------------

# To replicate the following analyses, create a folder named "working_directory".
# It should have 4 sub-folders: "raw_data", "generated_data", "tables" and "figures".
# workdir <- "PathToWorkingDirectory/" # Replace here with the path to the directory in which you placed the raw data files.
setwd(workdir)


# Load packages --------------------

# Data manipulation and analysis
library(dplyr) # version 1.1.3
library(tidyr) # version 1.3.0
library(scales) # version 1.3.0
library(reshape2) # version 1.4.4
library(Hmisc) # version 5.1-1
library(DataExplorer) # version 0.8.2
library(sf) # version 1.0-14
library(ggplot2) # version 3.4.4
library(ggridges) # version 0.5.4
library(ggsci) # version 3.0.0
library(cowplot) # version 1.1.1
library(corrplot) # version 0.92
library(egg) # version 0.4.5
library(rnaturalearth) # version 0.3.4
library(rnaturalearthhires) # version 0.2.1
library(ggspatial) # version 1.1.9
library(ggrepel) # version 0.9.4
library(raster) # version 3.6.23
library(PLNmodels) # version 1.0.4
library(igraph) # version 1.5.1
library(ggraph) # version 2.1.0
library(corrplot) # version 0.92
library(automap) # version 1.1-9
library(RColorBrewer) # version 1.1-3
library(viridis) # version 0.6.4



# ------------------------------------------------------------------------------ #
# Import the data, formatting, and plot preparation
# ------------------------------------------------------------------------------ #



# Open df_fig_data.csv and df_tree_data.csv which are formatted files combining fig- and tree-level observations.
# To read the raw data and format them, refer to the R script: R_script_Piatscheck.et.al.2020_Fpet_data_formating.R

df_fig_data <- read.csv(paste(workdir,"generated_data/","df_fig_data.csv", sep = ""))
# df_fig_data is the data set that contains both fig and tree- level observations. 
# However, we measured ecological variables on trees that did not contain E phase syconia.
# These observations are not in this data set.

df_tree_data <- read.csv(paste(workdir,"generated_data/","df_tree_data.csv", sep = ""))
# df_tree_data is a simplified version of the dataset, each row corresponds to one tree observed (at a particular season).
# Fig-level variables are absent here.
# Wasp counts are averages (which can be changed to sums in R script data formating)
# The dataset contains, however, measurements made on trees not containing E phase syconia ("male" phase figs containing wasps), which are absent in the fig-level dataset.


# ____________________________________________________________________________ #
# Format factor variables

df_fig_data$site <- factor(df_fig_data$site, levels = c("158", "172","112","113","95","179","201","96","70", "250"))
df_fig_data$season <- factor(df_fig_data$season, levels = c("F2012", "S2013","F2013","S2014"))
df_tree_data$site <- factor(df_tree_data$site, levels = c("158", "172","112","113","95","179","201","96","70", "250"))
df_tree_data$season <- factor(df_tree_data$season, levels = c("F2012", "S2013","F2013","S2014"))

# Here we do not remove the outliers as they will be visualized in some figures.
# Some will be removed later in the script.


# ____________________________________________________________________________ #
# Set a color palette


# Display colors
display.brewer.all()

# Visualize 14 colors in a pie
# See http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf for more colors.
pie(rep(1,40), col=c("darkcyan", "orange", "deepskyblue", "lightslategray", "firebrick2", 
                     "forestgreen", "blueviolet", "tan4", "yellow","darkgreen", 
                     "red",  "midnightblue", "snow", "mediumaquamarine", "dodgerblue2",
                     "#E31A1C", "green4","#6A3D9A","#FF7F00",  "black",
                     "gold1", "skyblue2","#FB9A99", "palegreen2", "#CAB2D6",
                     "#FDBF6F","gray70", "khaki2","maroon","orchid1",
                     "deeppink1", "blue1", "steelblue4", "darkturquoise", "hotpink",
                     "green1","yellow4","yellow3", "darkorange4","brown"))

# Create my own palettes with colorRampPalette function.

corrPalette <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", 
                                  "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                  "#4393C3", "#2166AC", "#053061"))

waspPalette <- colorRampPalette(c("darkcyan", "orange", "deepskyblue", "firebrick1", "olivedrab1",
                                 "blueviolet", "tan4", "yellow", "mediumaquamarine", "darkgreen", 
                                 "red", "midnightblue", "snow", "lightslategray", "forestgreen",
                                 "gray19"))

pcaPalette <- colorRampPalette(c("black", "lightpink1", "white", "orangered1", "deepskyblue3",
                                 "lemonchiffon3", "lightsalmon", "grey50", "cyan1", "slateblue3",
                                 "saddlebrown", "darkorange1", "springgreen3", "deeppink3", "khaki1",
                                 "seagreen", "mediumorchid1", "dodgerblue4", "greenyellow", "blue2",
                                 "firebrick4", "magenta2", "tomato2", "yellow1", "turquoise1",
                                 "green2", "red2"))


# ____________________________________________________________________________ #
# Labels


# Wasp species labels --------------------

# Partially italicized wasp species labels: 
wasplabels <- c(expression(paste(italic("Pegoscapus"), " sp.")), 
                expression(paste(italic("Idarnes flavicolis"), " sp. 1")),
                expression(paste(italic("Idarnes carme"), " sp. 1")), 
                expression(paste(italic("Idarnes carme"), " sp. 2")), 
                expression(paste(italic("Idarnes"), " males")), 
                expression(paste(italic("Heterandrium"), " sp. 1")), 
                expression(paste(italic("Heterandrium"), " sp. 2")), 
                expression(paste(italic("Ficicola"), " sp.")),
                expression(paste(italic("Physothorax"), " sp.")),
                expression(paste(italic("Sycophila"), " sp.")))

wasplabels_allwasps <- c(expression(paste(italic("Pegoscapus"), " sp.")), 
								expression(paste(italic("Idarnes flavicolis"), " sp. 1")),
								expression(paste(italic("Idarnes carme"), " sp. 1")), 
								expression(paste(italic("Idarnes carme"), " sp. 2")), 
								expression(paste(italic("Idarnes"), " males")), 
								expression(paste(italic("Heterandrium"), " sp. 1")), 
								expression(paste(italic("Heterandrium"), " sp. 2")), 
								expression(paste(italic("Ficicola"), " sp.")),
								expression(paste(italic("Physothorax"), " sp.")),
								expression(paste(italic("Sycophila"), " sp.")),
								expression("All wasp species"))


# Host-related ecological variables --------------------

hostvarlabels <- c("Foundresses", "Tree volume", "Reproduction", "Crop size", "Asynchrony", "Flowering landscape", "Syconium landscape")


# Abiotic variables --------------------

abiovarlabels <- c("DM T.max.", "DM T.min.","DM Prec.","DM V.P.","CHIRPS Prec.","MODIS LST","TC T.max.","TC T.min.","TC Prec.","TC V.P.","TC W.S.")



# Field trips --------------------

fielftriplabels <- c("Fall 2012","Spring 2013","Fall 2013","Spring 2014")


# Function to replace facets labels --------------------

fielftrip_names <- list(
  'F2012'="Fall 2012",
  'S2013'="Spring 2013",
  'F2013'="Fall 2013",
  'S2014'="Spring 2014")

fieldtrip <- function(variable,value){
  return(fielftrip_names[value])
}

# The trick above was found on a forum (I couldn't find the source again), thanks to the person that shared this!



# ---------------------------------------------------------------------------- #
# Figures: Study sites map
# ---------------------------------------------------------------------------- #



# Main study sites coordinates --------------------

Fpet_sites_coord <-  read.csv("D:/Finn/Research/Journal Articles/Manuscripts/2019-Environmental_effects_on_fig_wasp_community/working_directory/raw_data/Fpet_sites_coord.csv")

# Formating
Fpet_sites_coord$pop <- factor(Fpet_sites_coord$pop, levels = c("158","172","112","113","95","179","201","96","70")) # Sites ordered from north (first) to south (last)
Fpet_sites_coord <- Fpet_sites_coord[,c(3,1:2)]

# Extend the extent by several degrees each four cardinal directions
Fpet_sites_extent <- extent(x = c(min(Fpet_sites_coord$longitude) - 3, 
                                  max(Fpet_sites_coord$longitude) + 3,
                                  min(Fpet_sites_coord$latitude) - 1.5,
                                  max(Fpet_sites_coord$latitude) + 1.5))

# Get the country boundaries
mexicoAndUSA <- ne_countries(country = c("mexico","United States of America"), returnclass='sf', scale = "large")
# Actually USA map is not needed for the following dimensions.

# Visualize the map
ggplot(data = mexicoAndUSA) + geom_sf() +
  coord_sf(xlim = Fpet_sites_extent[1:2], ylim = Fpet_sites_extent[3:4], expand = FALSE) +
  geom_point(data = Fpet_sites_coord, aes(x = longitude, y = latitude), size = 3.5, shape = 21, stroke = 1.5, colour = "white", fill = "black") +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) 


# Figure 4A --------------------

ggplot(data = mexicoAndUSA) + geom_sf() +
  coord_sf(xlim = Fpet_sites_extent[1:2], ylim = Fpet_sites_extent[3:4], expand = FALSE) +
  geom_point(data = Fpet_sites_coord, aes(x = longitude, y = latitude), size = 4, shape = 21, stroke = 1.5, colour = "white", fill = "black") +
  geom_text_repel(data = Fpet_sites_coord, aes(x = longitude, y = latitude, label = pop), nudge_x = -0.4, vjust =1, size = 5, fontface = "bold", segment.color = NA)  +
  theme_classic() +
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),
        panel.border = element_rect(fill = NA, size = 0.5), # element_rect(colour = "black", fill=NA, size=1.5)
        text=element_text(size = 16),
        axis.text.x=element_text(colour="black", size = 12, angle = 45, vjust = 0.5),
        axis.text.y=element_text(colour="black", size = 12)) +
  annotate(geom = "text", x = -109.5, y = 30, label = "United \n Mexican States", color = "grey22", size = 4.5) +
  annotate(geom = "text", x = -114, y = 25, label = "Pacific Ocean", fontface = "italic", color = "grey22", size = 4) +
  annotate(geom = "text", x = -111, y = 27, label = "Sea of Cortez", fontface = "italic", color = "grey22", size = 4, angle = -50) +
  annotation_scale(location = "bl", width_hint = 0.5) + 
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.3, "in"), pad_y = unit(0.5, "in"), style = north_arrow_orienteering()) +
  xlab("Longitude") + 
  ylab("Latitude")



# Saved with ggsave(), width = 6 inches, height = 8 inches, dpi = 600.
ggsave(paste0(workdir,"figures/","Fpet_study_sites_map.png",sep=""), plot = last_plot(), device = "png", path = NULL, scale = 1, width = 7, height = 7, units = c("in"), dpi = 600, limitsize = TRUE)

# Clear working environment
rm(Fpet_sites_coord, Fpet_sites_extent, mexicoAndUSA)



# ---------------------------------------------------------------------------- #
# Figures: Data exploration and analysis
# ---------------------------------------------------------------------------- #



# Supplementary Figure S13 --------------------

# Prepare the data
df_fig_reponse_plot <- subset(df_fig_data, select=c(season,site,pollinators,LO1_f,SO1_f,SO2_f,idarnes_m,heterandrium_1,heterandrium_2,ficicola,physothorax,sycophila,wasps_total))
df_fig_reponse_plot_m <- melt(df_fig_reponse_plot, id.vars = c("season", "site"))

ggplot(df_fig_reponse_plot_m, aes(x = value, y = variable, fill = variable)) + 
  geom_density_ridges(alpha = 0.8) + 
  geom_boxplot(alpha = 0.1, outlier.alpha = 0.8, width = 1, colour = "grey") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(15, 10, 10, 10, "pt"),
        plot.title=element_text(size = 20),
        text=element_text(size = 16),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, vjust = 1, margin = margin(3.5,0,0,0, "pt")),
        axis.text.x=element_text(colour="black", size = 12, angle = 0, vjust = 0.5),
        axis.text.y=element_blank(),
        legend.position="right",
        legend.background = element_blank(),
        legend.spacing.y = unit(3, "mm"),
        legend.text=element_text(size=10),
        legend.text.align = 0) +
  # facet_wrap(season ~ .,  ncol=2, strip.position = c("top"), labeller = fieldtrip) + # , scales="free_x"
  xlab("Counts") +
  ylab("Frequency") +
  xlim(c(0,500)) +
  scale_fill_manual(values = rainbow(11), name = "Wasp species", breaks = levels(df_fig_reponse_plot_m$variable), labels = wasplabels_allwasps)
                      
# Saved with ggsave(), width = 8.5 inches, height = 7 inches, dpi = 600.
ggsave(paste0(workdir,"figures/","wasps_density_plot.png",sep=""), plot = last_plot(), device = "png", path = NULL, scale = 1, width = 8.5, height = 7, units = c("in"), dpi = 600, limitsize = TRUE)

# Clear working environment
rm(df_fig_reponse_plot, df_fig_reponse_plot_m)


# Supplementary Figure S14 --------------------

# Prepare the data
df_tree_host_variable_plot <- subset(df_tree_data, select=c(season,site,tree_volume,reproduction,crop_size,asynchrony,flowering_landscape_250_125,syconium_landscape_250_125)) # We chose flowering_landscape_250_50 an syconium_landscape_250_50 based on preliminary analyses and glmm best fits.
df_tree_host_variable_plot_m <- melt(df_tree_host_variable_plot, id.vars = c("season", "site"))

df_foundresses <- subset(df_fig_data, select=c(season,site,foundress_number))
df_foundresses_plot_m <- melt(df_foundresses, id.vars = c("season", "site"))

df_host_variable_plot_m <- rbind(df_foundresses_plot_m, df_tree_host_variable_plot_m)

# Visualize
ggplot(na.omit(df_host_variable_plot_m), aes(x = season, y=value, fill=variable, colour=variable)) +
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(15, 10, 10, 10, "pt"),
        text=element_text(size = 16),
        axis.text.x=element_text( angle = 90, vjust = 0.5),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        strip.text = element_text(size = 12, vjust = 1, margin = margin(3.5,0,0,0, "pt"))) +
  facet_grid(variable ~ site , scales="free_y") +
  xlab("Field season") +
  ylab(NULL) +
  scale_fill_discrete(breaks = c("foundress_number", levels(df_tree_host_variable_plot_m$variable)), labels = hostvarlabels, name = "Host-related\necological variables") +
  scale_colour_discrete(breaks = c("foundress_number", levels(df_tree_host_variable_plot_m$variable)), labels = hostvarlabels, name = "Host-related\necological variables") 
  
# Saved with ggsave(), width = 8.5 inches, height = 7 inches, dpi = 600.
ggsave(paste0(workdir,"figures/","host-related_variables_boxplots.png",sep=""), plot = last_plot(), device = "png", path = NULL, scale = 1, width = 12, height = 8, units = c("in"), dpi = 600, limitsize = TRUE)


# Clear working environment
rm(df_tree_host_variable_plot, df_tree_host_variable_plot_m, df_foundresses, df_foundresses_plot_m ,df_host_variable_plot_m)


# Supplementary Figure S16 --------------------

df_tree_abiotic_variable_plot <- subset(df_tree_data, select=c(season,site,DM_tmax,DM_tmin,DM_prec,DM_vp,CHIRPS_prec,MODIS_LST,TC_tmax,TC_tmin,TC_prec,TC_vp,TC_ws)) # We plot all climate variable although not all are considered in further analyses.
df_tree_abiotic_variable_plot <- df_tree_abiotic_variable_plot[!duplicated(df_tree_abiotic_variable_plot), ]
df_tree_abiotic_variable_plot_m <- melt(df_tree_abiotic_variable_plot, id.vars = c("season", "site"))

# Visualize
ggplot(na.omit(df_tree_abiotic_variable_plot_m), aes(x = season, y=value, fill=variable, colour=variable)) +
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(15, 10, 10, 10, "pt"),
        text=element_text(size = 16),
        axis.text.x=element_text( angle = 90, vjust = 0.5),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        strip.text = element_text(size = 12, vjust = 1, margin = margin(3.5,0,0,0, "pt"))) +
  facet_grid(variable ~ site , scales="free_y") +
  xlab("Field season") +
  ylab(NULL) +
  scale_fill_discrete(breaks = levels(df_tree_abiotic_variable_plot_m$variable), labels = abiovarlabels, name = "Abiotic variables") +
  scale_colour_discrete(breaks = levels(df_tree_abiotic_variable_plot_m$variable), labels = abiovarlabels, name = "Abiotic variables") 

# Saved with ggsave(), width = 8.5 inches, height = 7 inches, dpi = 600.
ggsave(paste0(workdir,"figures/","abiotic_variables_boxplots.png",sep=""), plot = last_plot(), device = "png", path = NULL, scale = 1, width = 12, height = 8, units = c("in"), dpi = 600, limitsize = TRUE)

# Clear working environment
rm(df_tree_abiotic_variable_plot, df_tree_abiotic_variable_plot_m)



# Supplemental Figure S15 --------------------


# Re-create the correlation matrix and tests
df_fig_data_all_corr <- subset(df_fig_data, select=c(foundress_number,crop_size,asynchrony,flowering_landscape_250_125,syconium_landscape_250_125, 
                                                     DM_tmax, DM_tmin,DM_prec,DM_vp,CHIRPS_prec,MODIS_LST,TC_tmax,TC_tmin,TC_prec,TC_vp,TC_ws))
df_fig_data_all_corr <- na.omit(df_fig_data_all_corr) # Remove missing data, which reduces the data set size.
variable_all_corr <- cor(df_fig_data_all_corr, method="spearman") 
colnames(variable_all_corr) <- c("Foundresses", "Crop size", "Asynchrony", "Flowering landscape", "Syconium landscape", "DM T.max.", "DM T.min.","DM Prec.","DM V.P.","CHIRPS Prec.","MODIS LST","TC T.max.","TC T.min.","TC Prec.","TC V.P.","TC W.S.")
rownames(variable_all_corr) <- c("Foundresses", "Crop size", "Asynchrony", "Flowering landscape", "Syconium landscape", "DM T.max.", "DM T.min.","DM Prec.","DM V.P.","CHIRPS Prec.","MODIS LST","TC T.max.","TC T.min.","TC Prec.","TC V.P.","TC W.S.")
diag(variable_all_corr) = NA
variable_all_corr_test_sp <- rcorr(as.matrix(variable_all_corr), type = "spearman") # the function rcorr() resolve the issue above.

# Plot and save with png(), ggsave() with last_plot() won't work here. 
png(width = 10, height=8,file = paste0(workdir,"figures/","ecological_variable_correlation.png", sep=""), units = c("in"), res = 300,  type = "cairo")
corrplot(variable_all_corr, method = "color", type = "upper", 
         tl.col="black", tl.srt=35, order = NULL, diag = TRUE, na.label = "NA",
         addCoef.col = "darkgrey", sig.level = 0.05, insig = "pch", p.mat = variable_all_corr_test_sp$P)
dev.off()



# Clear working environment
rm(df_fig_data_all_corr, df_fig_corr,variable_corr,variable_all_corr_test_sp,variable_all_corr)



# ---------------------------------------------------------------------------- #
# Figures: fig-level response variables visualization
# ---------------------------------------------------------------------------- #



# ____________________________________________________________________________ #
# Wasp proportions per site and season


# Re-format the dataset --------------------

df_plot_proportions <- df_fig_data %>% 
	group_by(season, site) %>%
	summarise(pollinators =sum(pollinators), parasites =sum(parasites))  %>% 
	mutate(freq_pollinators = pollinators / (pollinators + parasites),
				 freq_parasites = parasites / (pollinators + parasites))
df_plot_proportions_red  <- df_plot_proportions[, c(1:2,5:6)]
df_plot_proportions_reshape <- melt(as.data.frame(df_plot_proportions_red), id=c("site", "season"))


# Figure 6 --------------------

ggplot(df_plot_proportions_reshape,aes(x = site, y = value, fill = variable)) + 
	geom_bar(position = "fill",stat = "identity", colour = "black") + # Add " color = "black" but this creates delimitations in the bars
	# or:
	# geom_bar(position = position_fill(), stat = "identity") 
	theme(axis.line.x = element_line(size = 0.5, colour = "black"),
				axis.line.y = element_line(size = 0.5, colour = "black"),
				axis.line = element_line(size=1, colour = "black"),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				panel.border = element_blank(),
				panel.background = element_blank(),
				plot.margin = margin(15, 10, 10, 10, "pt"),
				plot.title=element_text(size = 20),
				text=element_text(size = 16),
				strip.background = element_blank(),
				strip.text = element_text(size = 12, vjust = 2.5, margin = margin(3.5,0,0,0, "pt")),
				axis.text.x=element_text(colour="black", size = 12, angle = 0, vjust = 0.5),
				axis.text.y=element_text(colour="black", size = 12, angle = 90, hjust = 0.5),
				legend.position=c(0.15, 0.68),
				legend.background = element_rect(size=0.5, linetype="solid", colour ="black", fill="white"),
				legend.spacing.y = unit(0, "mm"),
				legend.text=element_text(size=10)) +
	scale_y_continuous(labels = percent_format(), expand = c(0,0)) +
	labs(x = "Sites", y = "Wasp proportions") +
	scale_fill_manual(values = c("#1E90FF", "#E1B378"), name=NULL, # Try with "#FF7F00", "#1E90FF"
										breaks=c("freq_pollinators", "freq_parasites"),
										labels=c("Pollinators", "Non-pollinating fig wasps"))  + 
	facet_wrap( .~ season, ncol=2, scales = "free_x", strip.position = c("top"), labeller = fieldtrip)

# Saved with ggsave(), width = 8.5 inches, height = 7 inches, dpi = 600.
ggsave("wasp_proportions_per_site_and_season.png", plot = last_plot(), device = "png", path = NULL, scale = 1, width = 8.5, height = 7, units = c("in"), dpi = 600, limitsize = TRUE)



# ____________________________________________________________________________ #
# Principal Component Analysis


# Re-run multivariate analysis with RRPP
df_fig_RRPP <- rrpp.data.frame(logY = log(as.matrix(df_fig_data[, c("pollinators", "LO1_f", "SO1_f", "SO2_f", "heterandrium_1", "heterandrium_2", "ficicola", "physothorax", "sycophila")]) + 1), 
															 season = df_fig_data$season,
															 site = df_fig_data$site)
gps.all <- as.factor(paste(df_fig_data$season,df_fig_data$site,sep = "-"))
wasp_abund_fit <- lm.rrpp(logY ~ gps.all, data = df_fig_RRPP, iter = 999, SS.type = "I", print.progress = TRUE)
# The figure might be somewhat different because the random permutations are likely to produce slightly different results.

# Extract the wanted data
pca_PW <- plot(wasp_abund_fit, type = "PC")
pca_plot <- as.data.frame(pca_PW$PC.points)
pca_plot$levels <- gps.all
pca_plot$levels <- factor(pca_plot$levels, levels = c("F2012-158", "F2012-172", "S2013-112", "F2012-113", "F2012-95", "F2012-201","F2012-96", "F2012-70",
																											"S2013-158", "S2013-172", "S2013-201", "S2013-96", "S2013-70",
																											"F2013-172", "F2013113", "F2013-95", "F2013-96", "F2013-70",
																											"S2014-158", "S2014-172", "S2014-112", "S2014113", "S201495", "S2014179", "S2014201", "S201496", "S201470"))

levels(pca_plot$levels) <- c("F2012-s158", "F2012-s172", "S2013-s112", "F2012-s113", "F2012-s95", "F2012-s201","F2012-s96", "F2012-s70",
														 "S2013-s158", "S2013-s172", "S2013-s201", "S2013-s96", "S2013-s70",
														 "F2013-s172", "F2013-s113", "F2013-s95", "F2013-s96", "F2013-s70",
														 "S2014-s158", "S2014-s172", "S2014-s112", "S2014-s113", "S2014-s95", "S2014-s179", "S2014-s201", "S2014-s96", "S2014-s70")

# Count levels for colors
levels(pca_plot$levels)
length(levels(pca_plot$levels)) # 27 colors, we need a 27 color palette which was defined previously.
colourCountPca <- length(levels(pca_plot$levels))


# Supplementary Figure 18 --------------------

ggplot(pca_plot, aes(x=PC1, y=PC2))  + 
	geom_point(aes(fill=factor(levels)), size = 3, shape = 21, stroke = 1, colour = "black") + 
	geom_hline(yintercept=0,linetype=2) + 
	geom_vline(xintercept=0,linetype=2) +
	#scale_fill_discrete(name="Regions") +
	#  stat_ellipse(aes(fill=factor(levels)), geom = "polygon", alpha = 0.2) +
	stat_ellipse(aes(color=factor(levels)), geom = "path", alpha = 1) +
	theme_classic(base_size = 16) +
	theme(plot.margin = margin(15, 10, 10, 10, "pt"),
				axis.text.y=element_text(colour="black", angle = 90),
				legend.position = "right") +
	scale_fill_manual("Site- and season-specific\nwasp community means", values = pcaPalette(colourCountPca)) +
	scale_color_manual(guide = FALSE, values = pcaPalette(colourCountPca)) +
	xlab(pca_PW$plot_args$xlab) +
	ylab(pca_PW$plot.args$ylab)

# Saved with ggsave(), width = 10 inches, height = 7 inches, dpi = 600.
ggsave(paste0(workdir,"figures/","PCA_wasp_community_site_by_season_RRPP.png", sep ="" ), plot = last_plot(), device = "png", path = NULL, scale = 1, width = 10, height = 7, units = c("in"), dpi = 600, limitsize = TRUE)



# ____________________________________________________________________________ #
# Wasp means per site and season

# Re-format the dataset --------------------

df_plot_wasps_means <- df_fig_data %>% 
	group_by(season) %>% 
	summarise(pegoscapus_f = mean(pegoscapus_f), 
						pegoscapus_m = mean(pegoscapus_m),
						LO1=mean(LO1_f), 
						SO1=mean(SO1_f), 
						SO2=mean(SO2_f), 
						idarnes_m=mean(idarnes_m),
						heterandrium_1_f=mean(heterandrium_1_f), 
						heterandrium_1_mw=mean(heterandrium_1_mw),
						heterandrium_1_mu=mean(heterandrium_1_mu),
						heterandrium_2_f=mean(heterandrium_2_f), 
						heterandrium_2_mw=mean(heterandrium_2_mw),
						heterandrium_2_mu=mean(heterandrium_2_mu),
						ficicola_f=mean(ficicola_f), 
						ficicola_m=mean(ficicola_m),
						physothorax_f=mean(physothorax_f),
						physothorax_m=mean(physothorax_m),
						sycophila=mean(sycophila)) %>%
	gather(key = key, value = val, -season)
df_plot_wasps_means$key <- as.factor(df_plot_wasps_means$key)
df_plot_wasps_means$key <- factor(df_plot_wasps_means$key, levels =unique(df_plot_wasps_means$key))
df_plot_wasps_means$sp <- as.factor(c("Pegoscapus", "Pegoscapus","Pegoscapus", "Pegoscapus", "Pegoscapus", "Pegoscapus", "Pegoscapus", "Pegoscapus",
																			"LO1", "LO1", "LO1", "LO1",
																			"SO1", "SO1", "SO1", "SO1", 
																			"SO2", "SO2", "SO2", "SO2", 
																			"Idarnes males", "Idarnes males", "Idarnes males", "Idarnes males",
																			"Heterandrium 1", "Heterandrium 1", "Heterandrium 1", "Heterandrium 1", "Heterandrium 1", "Heterandrium 1", "Heterandrium 1", "Heterandrium 1", "Heterandrium 1", "Heterandrium 1", "Heterandrium 1", "Heterandrium 1", 
																			"Heterandrium 2", "Heterandrium 2", "Heterandrium 2", "Heterandrium 2", "Heterandrium 2", "Heterandrium 2", "Heterandrium 2", "Heterandrium 2", "Heterandrium 2", "Heterandrium 2", "Heterandrium 2", "Heterandrium 2",
																			"Ficicola", "Ficicola", "Ficicola", "Ficicola","Ficicola", "Ficicola","Ficicola", "Ficicola",
																			"Physothorax", "Physothorax", "Physothorax", "Physothorax","Physothorax", "Physothorax","Physothorax", "Physothorax",
																			"Sycophyla", "Sycophyla", "Sycophyla", "Sycophyla"))
df_plot_wasps_means$sp <- factor(df_plot_wasps_means$sp, levels = c("Pegoscapus", "LO1", "SO1", "SO2", "Idarnes males", "Heterandrium 1", "Heterandrium 2", "Ficicola", "Physothorax", "Sycophyla"))

# Prepare some elements that will be added to the plots --------------------

# Count for numbers of color used in the graph
colourCount = length(unique(df_plot_wasps_means$key))

# We want to add some information about number of syconia sampled on the graph
# Total number of syconia sampled
df_plot_F2012 <- subset(df_fig_data, season =="F2012")
length(df_plot_F2012$fig)
# 634
df_plot_S2013 <- subset(df_fig_data, season =="S2013")
length(df_plot_S2013$fig)
# 435
df_plot_F2013 <- subset(df_fig_data, season =="F2013")
length(df_plot_F2013$fig)
# 388
df_plot_S2014 <- subset(df_fig_data, season =="S2014")
length(df_plot_S2014$fig)
# 847
length(df_fig_data$fig)
# Total = 2304

# Total number of wasp sampled
sum(df_plot_F2012$wasps_total)
# 43850
sum(df_plot_S2013$wasps_total)
# 39794
sum(df_plot_F2013$wasps_total)
# 44092
sum(df_plot_S2014$wasps_total)
# 88863

# Add the info above to the plot
facet_text <- c("Fall 2012        \nSyconia: 634 \nWasps: 43850", "Spring 2013  \nSyconia: 435 \nWasps: 39794", "Fall 2013        \nSyconia: 388 \nWasps: 44092", "Spring 2014  \nSyconia: 847 \nWasps: 88863")

# Supplementary Figure 17 --------------------

waspcountplot_means <- ggplot(data=df_plot_wasps_means, aes(x=sp, y=val, fill=key)) +
	geom_bar(colour="black", stat="identity", position = "stack") +
	theme(axis.line.x = element_line(size = 0.5, colour = "black"),
				axis.line.y = element_line(size = 0.5, colour = "black"),
				axis.line = element_line(size=1, colour = "black"),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				panel.border = element_blank(),
				panel.background = element_blank(),
				plot.title=element_text(size = 20),
				text=element_text(size = 16),
				strip.background = element_blank(),
				strip.text = element_text(size = 16),
				axis.text.x=element_text(colour="black", size = 10, angle = 90, vjust = 0.5),
				axis.text.y=element_text(colour="black", size = 12, angle = 90, hjust = 0.5)) +
	labs(x = NULL, y = "Wasp average") +
	facet_wrap( .~ season, ncol=4, strip.position = c("top"), labeller = fieldtrip) +
	scale_fill_manual(values = waspPalette(colourCount)) +
	scale_x_discrete(labels=wasplabels) +
	scale_y_continuous(expand = c(0,0)) +
	theme(legend.position="none") 
waspcountplot_means_tag <- tag_facet(waspcountplot_means, open = "", close = "", tag_pool = facet_text, x = -Inf, y = Inf, fontface = 2) 
# EDIT: tag_facet() removes facet titles, so I added them in "facet_text", fieldtrip() not necessary anymore here (but used later).
waspcountplot_means_tag

# Saved with ggsave(), width = 8.5 inches, height = 5 inches, dpi = 600.
ggsave(paste0(workdir,"figures/","wasp_species_averages.png", sep=""), plot = last_plot(), device = "png", path = NULL, scale = 1, width = 8.5, height = 5, units = c("in"), dpi = 600, limitsize = TRUE)

# Clear memory
rm(df_plot_wasps_means, colourCount, df_plot_F2012, df_plot_S2013, df_plot_F2013, df_plot_S2014, waspcountplot_means,waspcountplot_means_tag )



# ---------------------------------------------------------------------------- #
# Figures: fig- and tree-level predictor variables visualization
# ---------------------------------------------------------------------------- #



# In the following we remove outliers at the tree level as they are previously identified.

# Replace outliers with NAs.
df_tree_data_outlrm <-	df_tree_data
df_tree_data_outlrm[rownames(df_tree_data_outlrm[which(df_tree_data_outlrm$tree_volume > 6000),]),c(which(colnames(df_tree_data_outlrm) == "tree_volume"), which(colnames(df_tree_data_outlrm) == "crop_size"))] <- NA
which(df_tree_data$tree_volume > 6000)
which(df_tree_data_outlrm$tree_volume > 6000) # Outliers removed


# ____________________________________________________________________________ #
# Fig-level ecological variables plots

ggplot(df_fig_data, aes(x = foundress_number)) +
  geom_histogram(stat="count", fill = "#56B4E9", color = "black") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(30, 10, 10, 10, "pt"),
        text=element_text(size = 16),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12, angle = 90, hjust = 0.5)) +
  scale_x_continuous(breaks = seq(0, 10, 1), labels = c("0","1","2","3","4","5","6","7","8","9","10+")) +
  scale_y_continuous(expand = c(0,1)) +
  ylab("Density") +
  xlab("Foundress number")

# Saved with ggsave(), width = 4 inches, height = 4 inches, dpi = 600.
# Part of Figure 7 --------------------

ggsave(paste0(workdir,"figures/", "foundress_counts_histogram.png", sep=""), plot = last_plot(), device = "png", path = NULL, scale = 1, width = 4, height = 4, units = c("in"), dpi = 600, limitsize = TRUE)



# ____________________________________________________________________________ #
# Tree-level ecological variables plots

synchroplot <- ggplot(df_tree_data_outlrm, aes(x = site, y = asynchrony, fill = season)) +
  geom_boxplot() +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(15, 10, 10, 10, "pt"),
        text=element_text(size = 16),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12, angle = 90, hjust = 0.5),
        legend.position=c(0.1, 0.85), 
        legend.background = element_rect(colour='black'),
        legend.title = element_text(colour="black", size = 10),
        legend.text  = element_text(colour="black", size = 10)) +
  scale_y_continuous(name = "Asynchrony index", 
                     breaks = seq(1, 7, 1),
                     limits=c(1, 7)) +
  scale_x_discrete(name = NULL) +
  scale_fill_nejm(name="Sampling season",
                  breaks=c("F2012", "S2013", "F2013", "S2014"),
                  labels=fielftriplabels)

synchroplot
# Saved with ggsave(), width = 8.5 inches, height = 4 inches, dpi = 600.
ggsave(paste0(workdir,"figures/","tree_synchrony_2.png", sep=""), plot = last_plot(), device = "png", path = NULL, scale = 1, width = 8.5, height = 4, units = c("in"), dpi = 600, limitsize = TRUE)


cropsizeplot <- ggplot(df_tree_data_outlrm, aes(x = site, y = crop_size, fill = season)) +
  geom_boxplot() +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(15, 10, 10, 10, "pt"),
        text=element_text(size = 16),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12, angle = 90, hjust = 0.5)) +
  scale_y_continuous(name = "Crop size",
                     breaks = c(0, 1000, 2000),
                     labels = c("Very low", "Medium", "Very high"),
                     limits = c(0, 2000)) +
  scale_x_discrete(name = "Sites") +
  scale_fill_nejm() +
  theme(legend.position="none")

cropsizeplot
# # Saved with ggsave(), width = 8.5 inches, height = 4 inches, dpi = 600.
ggsave(paste0(workdir,"figures/","tree_crop_size.png", sep=""), plot = last_plot(), device = "png", path = NULL, scale = 1, width = 8.5, height = 4, units = c("in"), dpi = 600, limitsize = TRUE)


# Figure 7 --------------------

# # Combine plots
# ggdraw() +
#   draw_plot(synchroplot, x = 0, y = 0.5, width = 1, height = 0.5) +
#   draw_plot(cropsizeplot, x = 0, y = 0, width = 1, height = 0.5) +
#   draw_plot_label(label = c("A", "B"), size = 15,
#                   x = c(0, 0), y = c(1, 0.5))

# # Saved with ggsave(), width = 8.5 inches, height = 8.5 inches, dpi = 600.
# ggsave(paste0(workdir,"figures/","tree_phenology_variables.png", sep=""), plot = last_plot(), device = "png", path = NULL, scale = 1, width = 8.5, height = 8.5, units = c("in"), dpi = 600, limitsize = TRUE)

# Originally I wanted to combine plots on tree reproduction only together with cowplot, I used Inkscape instead to combine Figure 7.

# Clear memory
rm(foundressbar, cropsizeplot, synchroplot)



# ---------------------------------------------------------------------------- #
# Figures: landscape-level predictor variables visualization
# ---------------------------------------------------------------------------- #



# ____________________________________________________________________________ #
# Flowering tree landscape and syconium landscape spatial representation and krigging


# Site 112 --------------------

# Make a grid for the kriging
df_112 <- df_tree_data[which(df_tree_data$site == "112"),c("longitude", "latitude")]
colnames(df_112) <- c("x","y")
coordinates(df_112) <- ~x+y
proj4string(df_112) <- CRS("+init=epsg:4326")
df_112_sp <- spTransform(df_112, CRS("+proj=utm +zone=12"))
df_112_sp@bbox <- as.matrix(extent(df_112_sp)+1000) # extend for a larger representation.
df_112_grid <- makegrid(df_112_sp, n=10000) 
df_112_grid_sp <- spTransform(SpatialPoints(df_112_grid, proj4string = CRS("+proj=utm +zone=12")), CRS("+proj=utm +zone=12"))

# Flowering landscape site 112, spring 2013
fl_land_112_S2013 <- df_tree_data[which(df_tree_data$season == "S2013" & df_tree_data$site == "112"),c("longitude", "latitude", "flowering_landscape_250_125")]
fl_land_112_S2013 <- na.omit(fl_land_112_S2013)
fl_land_112_S2013_sp <- fl_land_112_S2013
colnames(fl_land_112_S2013_sp) <- c("x","y","var")
coordinates(fl_land_112_S2013_sp) <- ~x+y
proj4string(fl_land_112_S2013_sp) <- CRS("+init=epsg:4326")
fl_land_112_S2013_sp <- spTransform(fl_land_112_S2013_sp, CRS("+proj=utm +zone=12"))

# Kriging
fl_land_112_S2013_krig <- autoKrige(var~1, fl_land_112_S2013_sp, df_112_grid_sp)

# Format and plot with ggplot2
fl_land_112_S2013_krig_df <- as.data.frame(fl_land_112_S2013_krig$krige_output)
fl_land_112_S2013_krig_df$var1.pred[fl_land_112_S2013_krig_df$var1.pred<0] <- 0 # We can't have negative values, so we replace them with 0.
fl_land_112_S2013_sf <- st_as_sf(fl_land_112_S2013_sp)

fl_land_112_S2013_plot <- ggplot(data = fl_land_112_S2013_sf) + 
  geom_sf(aes(fill=var), shape = 21, size = 3, stroke = 1) +
  geom_tile(data =as.data.frame(fl_land_112_S2013_krig_df), aes(x = coords.x1, y = coords.x2, fill=var1.pred), alpha = 0.5) +
  theme(legend.position = "none", text=element_text(size = 12),
        # legend.background = element_rect(colour='white'),
        # legend.title = element_text(colour="black", size = 10),
        # legend.text  = element_text(colour="black", size = 10),
        # legend.key.width = unit(0.3,"cm"),
        # legend.key.height = unit(1.5, "cm"),
        axis.text.x=element_text(colour="black", size = 10, angle = 0, vjust = 0.5),  
        axis.text.y=element_text(colour="black", size = 10, angle = 90, hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, size = 0.5),
        panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "white"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(name=NULL, colours = rev(viridis(10)), limits = c(0,4)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-113.10, -113.08), labels = c("","")) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("") +
  ylab("Latitude") +
  ggtitle("Spring 2013")

# Flowering landscape site 112, fall 2013
fl_land_112_F2013 <- df_tree_data[which(df_tree_data$season == "F2013" & df_tree_data$site == "112"),c("longitude", "latitude", "flowering_landscape_250_125")]
fl_land_112_F2013 <- na.omit(fl_land_112_F2013)
fl_land_112_F2013_sp <- fl_land_112_F2013
colnames(fl_land_112_F2013_sp) <- c("x","y","var")
coordinates(fl_land_112_F2013_sp) <- ~x+y
proj4string(fl_land_112_F2013_sp) <- CRS("+init=epsg:4326")
fl_land_112_F2013_sp <- spTransform(fl_land_112_F2013_sp, CRS("+proj=utm +zone=12"))

# Kriging
fl_land_112_F2013_krig <- autoKrige(var~1, fl_land_112_F2013_sp, df_112_grid_sp)

# Format and plot with ggplot2
fl_land_112_F2013_krig_df <- as.data.frame(fl_land_112_F2013_krig$krige_output)
fl_land_112_F2013_krig_df$var1.pred[fl_land_112_F2013_krig_df$var1.pred<0] <- 0 # We can't have negative values, so we replace them with 0.
fl_land_112_F2013_sf <- st_as_sf(fl_land_112_F2013_sp)

fl_land_112_F2013_plot <- ggplot(data = fl_land_112_F2013_sf) + 
  geom_sf(aes(fill=var), shape = 21, size = 3, stroke = 1) +
  geom_tile(data =as.data.frame(fl_land_112_F2013_krig_df), aes(x = coords.x1, y = coords.x2, fill=var1.pred), alpha = 0.5) +
  theme(legend.position = "none", text=element_text(size = 12),
        # legend.background = element_rect(colour='white'),
        # legend.title = element_text(colour="black", size = 10),
        # legend.text  = element_text(colour="black", size = 10),
        # legend.key.width = unit(0.3,"cm"),
        # legend.key.height = unit(1.5, "cm"),
        axis.text.x= element_text(colour="black", size = 12, angle = 0, vjust = 0.5),  
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, size = 0.5),
        panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "white"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(name=NULL, colours = rev(viridis(10)), limits = c(0,3.5)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-113.10, -113.08), labels = c("","")) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("") +
  ylab(NULL) +
  ggtitle("Fall 2013")

# Flowering landscape site 112, spring 2014
fl_land_112_S2014 <- df_tree_data[which(df_tree_data$season == "S2014" & df_tree_data$site == "112"),c("longitude", "latitude", "flowering_landscape_250_125")]
fl_land_112_S2014 <- na.omit(fl_land_112_S2014)
fl_land_112_S2014_sp <- fl_land_112_S2014
colnames(fl_land_112_S2014_sp) <- c("x","y","var")
coordinates(fl_land_112_S2014_sp) <- ~x+y
proj4string(fl_land_112_S2014_sp) <- CRS("+init=epsg:4326")
fl_land_112_S2014_sp <- spTransform(fl_land_112_S2014_sp, CRS("+proj=utm +zone=12"))

# Kriging
fl_land_112_S2014_krig <- autoKrige(var~1, fl_land_112_S2014_sp, df_112_grid_sp)

# Format and plot with ggplot2
fl_land_112_S2014_krig_df <- as.data.frame(fl_land_112_S2014_krig$krige_output)
fl_land_112_S2014_krig_df$var1.pred[fl_land_112_S2014_krig_df$var1.pred<0] <- 0 # We can't have negative values, so we replace them with 0.
fl_land_112_S2014_sf <- st_as_sf(fl_land_112_S2014_sp)

fl_land_112_S2014_plot <- ggplot(data = fl_land_112_S2014_sf) + 
  geom_sf(aes(fill=var), shape = 21, size = 3, stroke = 1) +
  geom_tile(data =as.data.frame(fl_land_112_S2014_krig_df), aes(x = coords.x1, y = coords.x2,  fill=var1.pred), alpha = 0.5) +
  theme(legend.position = "right", text=element_text(size = 12),
        legend.background = element_rect(colour='white'),
        legend.title = element_text(colour="black", size = 10),
        legend.text  = element_text(colour="black", size = 10, angle = -90, hjust = 0.5),
        legend.key.width = unit(0.3,"cm"),
        legend.key.height = unit(1.5, "cm"),
        axis.text.x=element_text(colour="black", size = 12, angle = 0, vjust = 0.5),  
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, size = 0.5),
        panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "white"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(name=NULL, colours = rev(viridis(10)), limits = c(0,3.5)) + 
  scale_x_continuous(expand = c(0, 0), breaks = c(-113.10, -113.08), labels = c("","")) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("") +
  ylab(NULL) +
  ggtitle("Spring 2014")

# Syconium landscape site 112, spring 2013
sy_land_112_S2013 <- df_tree_data[which(df_tree_data$season == "S2013" & df_tree_data$site == "112"),c("longitude", "latitude", "syconium_landscape_250_125")]
sy_land_112_S2013 <- na.omit(sy_land_112_S2013)
sy_land_112_S2013_sp <- sy_land_112_S2013
colnames(sy_land_112_S2013_sp) <- c("x","y","var")
coordinates(sy_land_112_S2013_sp) <- ~x+y
proj4string(sy_land_112_S2013_sp) <- CRS("+init=epsg:4326")
sy_land_112_S2013_sp <- spTransform(sy_land_112_S2013_sp, CRS("+proj=utm +zone=12"))

# Kriging
sy_land_112_S2013_krig <- autoKrige(var~1, sy_land_112_S2013_sp, df_112_grid_sp)

# Format and plot with ggplot2
sy_land_112_S2013_krig_df <- as.data.frame(sy_land_112_S2013_krig$krige_output)
sy_land_112_S2013_krig_df$var1.pred[sy_land_112_S2013_krig_df$var1.pred<0] <- 0 # We can't have negative values, so we replace them with 0.
sy_land_112_S2013_sf <- st_as_sf(sy_land_112_S2013_sp)

sy_land_112_S2013_plot <- ggplot(data = sy_land_112_S2013_sf) + 
  geom_sf(aes(fill=var), shape = 21, size = 3, stroke = 1) +
  geom_tile(data =as.data.frame(sy_land_112_S2013_krig_df), aes(x = coords.x1, y = coords.x2,  fill=var1.pred), alpha = 0.5) +
  theme(legend.position = "none", text=element_text(size = 12),
        # legend.background = element_rect(colour='white'),
        # legend.title = element_text(colour="black", size = 10),
        # legend.text  = element_text(colour="black", size = 10),
        # legend.key.width = unit(0.3,"cm"),
        # legend.key.height = unit(1.5, "cm"),
        axis.text.x=element_text(colour="black", size = 10, angle = 0, vjust = 0.5),  
        axis.text.y=element_text(colour="black", size = 10, angle = 90, hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, size = 0.5),
        panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "white"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(name=NULL, colours = rev(viridis(10)), limits = c(0,810)) + 
  scale_x_continuous(expand = c(0, 0), breaks = c(-113.10, -113.08)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("Spring 2013")

# Syconium landscape site 112, fall 2013
sy_land_112_F2013 <- df_tree_data[which(df_tree_data$season == "F2013" & df_tree_data$site == "112"),c("longitude", "latitude", "syconium_landscape_250_125")]
sy_land_112_F2013 <- na.omit(sy_land_112_F2013)
sy_land_112_F2013_sp <- sy_land_112_F2013
colnames(sy_land_112_F2013_sp) <- c("x","y","var")
coordinates(sy_land_112_F2013_sp) <- ~x+y
proj4string(sy_land_112_F2013_sp) <- CRS("+init=epsg:4326")
sy_land_112_F2013_sp <- spTransform(sy_land_112_F2013_sp, CRS("+proj=utm +zone=12"))

# Kriging
sy_land_112_F2013_krig <- autoKrige(var~1, sy_land_112_F2013_sp, df_112_grid_sp)

# Format and plot with ggplot2
sy_land_112_F2013_krig_df <- as.data.frame(sy_land_112_F2013_krig$krige_output)
sy_land_112_F2013_krig_df$var1.pred[sy_land_112_F2013_krig_df$var1.pred<0] <- 0 # We can't have negative values, so we replace them with 0.
sy_land_112_F2013_sf <- st_as_sf(sy_land_112_F2013_sp)

sy_land_112_F2013_plot <- ggplot(data = sy_land_112_F2013_sf) + 
  geom_sf(aes(fill=var), shape = 21, size = 3, stroke = 1) +
  geom_tile(data =as.data.frame(sy_land_112_F2013_krig_df), aes(x = coords.x1, y = coords.x2, fill=var1.pred), alpha = 0.5) +
  theme(legend.position = "none", text=element_text(size = 12),
        # legend.background = element_rect(colour='white'),
        # legend.title = element_text(colour="black", size = 10),
        # legend.text  = element_text(colour="black", size = 10),
        # legend.key.width = unit(0.3,"cm"),
        # legend.key.height = unit(1.5, "cm"),
        axis.text.x=element_text(colour="black", size = 10, angle = 0, vjust = 0.5),  
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, size = 0.5),
        panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "white"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(name=NULL, colours = rev(viridis(10)), limits = c(0,810)) + 
  scale_x_continuous(expand = c(0, 0), breaks = c(-113.10, -113.08)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Longitude") +
  ylab(NULL) +
  ggtitle("Fall 2013")

# Syconium landscape site 112, spering 2014
sy_land_112_S2014 <- df_tree_data[which(df_tree_data$season == "S2014" & df_tree_data$site == "112"),c("longitude", "latitude", "syconium_landscape_250_125")]
sy_land_112_S2014 <- na.omit(sy_land_112_S2014)
sy_land_112_S2014_sp <- sy_land_112_S2014
colnames(sy_land_112_S2014_sp) <- c("x","y","var")
coordinates(sy_land_112_S2014_sp) <- ~x+y
proj4string(sy_land_112_S2014_sp) <- CRS("+init=epsg:4326")
sy_land_112_S2014_sp <- spTransform(sy_land_112_S2014_sp, CRS("+proj=utm +zone=12"))

# Kriging
sy_land_112_S2014_krig <- autoKrige(var~1, sy_land_112_S2014_sp, df_112_grid_sp)

# Format and plot with ggplot2
sy_land_112_S2014_krig_df <- as.data.frame(sy_land_112_S2014_krig$krige_output)
sy_land_112_S2014_krig_df$var1.pred[sy_land_112_S2014_krig_df$var1.pred<0] <- 0 # We can't have negative values, so we replace them with 0.
sy_land_112_S2014_sf <- st_as_sf(sy_land_112_S2014_sp)

sy_land_112_S2014_plot <- ggplot(data = sy_land_112_S2014_sf) + 
  geom_sf(aes(fill=var), shape = 21, size = 3, stroke = 1) +
  # coord_sf(datum = st_crs(sy_land_112_S2014_sf)) +
  geom_tile(data =as.data.frame(sy_land_112_S2014_krig_df), aes(x = coords.x1, y = coords.x2, fill=var1.pred), alpha = 0.5) +
  theme(legend.position = "right", text=element_text(size = 12),
        legend.background = element_rect(colour='white'),
        legend.title = element_text(colour="black", size = 12, angle = -90),
        legend.text  = element_text(colour="black", size = 10, angle = -90, hjust = 0.5),
        legend.key.width = unit(0.3,"cm"),
        legend.key.height = unit(1.5, "cm"),
        axis.text.x=element_text(colour="black", size = 10, angle = 0, vjust = 0.5),  
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, size = 0.5),
        panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "white"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(name=NULL, colours = rev(viridis(10)), limits = c(0,810)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-113.10, -113.08)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Longitude") +
  ylab(NULL) +
  ggtitle("Spring 2014")

# Figure 8 --------------------

# Combine the plots above
ggdraw() +
  draw_plot(fl_land_112_S2013_plot, x = 0, y = 0.5, width = 0.33, height = 0.5) +
  draw_plot(fl_land_112_F2013_plot, x = 0.33, y = 0.5, width = 0.33, height = 0.5) +
  draw_plot(fl_land_112_S2014_plot, x = 0.66, y = 0.5, width = 0.33, height = 0.5) +
  draw_plot(sy_land_112_S2013_plot, x = 0, y = 0, width = 0.33, height = 0.5) +
  draw_plot(sy_land_112_F2013_plot, x = 0.33, y = 0, width = 0.33, height = 0.5) +
  draw_plot(sy_land_112_S2014_plot, x = 0.66, y = 0, width = 0.33, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C", "D", "E", "F"), size = 15,
                  x = c(0.02, 0.33, 0.64, 0.02, 0.33, 0.64), y = c(1, 1, 1, 0.5, 0.5, 0.5))

# Saved with ggsave(), width = 9 inches, height = 8 inches, dpi = 600.
ggsave(paste0(workdir,"figures/","fl&sy_landscape_site_112.png", sep=""), plot = last_plot(), device = "png", path = NULL, scale = 1, width = 9, height = 8, units = c("in"), dpi = 600, limitsize = TRUE)



# Site 113 --------------------

# Make a grid for the kriging
df_113 <- df_tree_data[which(df_tree_data$site == "113" & df_tree_data$longitude < -112.50),c("longitude", "latitude")]
colnames(df_113) <- c("x","y")
coordinates(df_113) <- ~x+y
proj4string(df_113) <- CRS("+init=epsg:4326")
df_113_sp <- spTransform(df_113, CRS("+proj=utm +zone=12"))
df_113_sp@bbox <- as.matrix(extent(df_113_sp)+200) # extend for a larger representation.
df_113_grid <- makegrid(df_113_sp, n=10000) 
df_113_grid_sp <- spTransform(SpatialPoints(df_113_grid, proj4string = CRS("+proj=utm +zone=12")), CRS("+proj=utm +zone=12"))

# Flowering landscape site 113, spring 2012
fl_land_113_F2012 <- df_tree_data[which(df_tree_data$season == "F2012" & df_tree_data$site == "113" & df_tree_data$longitude < -112.50),c("longitude", "latitude", "flowering_landscape_250_125")]
fl_land_113_F2012 <- na.omit(fl_land_113_F2012)
fl_land_113_F2012_sp <- fl_land_113_F2012
colnames(fl_land_113_F2012_sp) <- c("x","y","var")
coordinates(fl_land_113_F2012_sp) <- ~x+y
proj4string(fl_land_113_F2012_sp) <- CRS("+init=epsg:4326")
fl_land_113_F2012_sp <- spTransform(fl_land_113_F2012_sp, CRS("+proj=utm +zone=12"))

# Kriging
fl_land_113_F2012_krig <- autoKrige(var~1, fl_land_113_F2012_sp, df_113_grid_sp)

# Format and plot with ggplot2
fl_land_113_F2012_krig_df <- as.data.frame(fl_land_113_F2012_krig$krige_output)
fl_land_113_F2012_krig_df$var1.pred[fl_land_113_F2012_krig_df$var1.pred<0] <- 0 # We can't have negative values, so we replace them with 0.
fl_land_113_F2012_sf <- st_as_sf(fl_land_113_F2012_sp)

fl_land_113_F2012_plot <- ggplot(data = fl_land_113_F2012_sf) + 
  geom_sf(aes(fill=var), shape = 21, size = 3, stroke = 1) +
  geom_tile(data =as.data.frame(fl_land_113_F2012_krig_df), aes(x = coords.x1, y = coords.x2, fill=var1.pred), alpha = 0.5) +
  theme(legend.position = "none", text=element_text(size = 12),
        # legend.background = element_rect(colour='white'),
        # legend.title = element_text(colour="black", size = 10),
        # legend.text  = element_text(colour="black", size = 10),
        # legend.key.width = unit(0.3,"cm"),
        # legend.key.height = unit(1.5, "cm"),
        axis.text.x=element_text(colour="black", size = 10, angle = 0, vjust = 0.5),  
        axis.text.y=element_text(colour="black", size = 10, angle = 90, hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, size = 0.5),
        panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "white"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(name=NULL, colours = rev(viridis(10)), limits = c(0,5)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-112.522, -112.518, -112.514), labels = c("","","")) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("") +
  ylab("Latitude") +
  ggtitle("Fall 2012")

# Flowering landscape site 113, spring 2013
fl_land_113_S2013 <- df_tree_data[which(df_tree_data$season == "S2013" & df_tree_data$site == "113" & df_tree_data$longitude < -112.50),c("longitude", "latitude", "flowering_landscape_250_125")]
fl_land_113_S2013 <- na.omit(fl_land_113_S2013)
fl_land_113_S2013_sp <- fl_land_113_S2013
colnames(fl_land_113_S2013_sp) <- c("x","y","var")
coordinates(fl_land_113_S2013_sp) <- ~x+y
proj4string(fl_land_113_S2013_sp) <- CRS("+init=epsg:4326")
fl_land_113_S2013_sp <- spTransform(fl_land_113_S2013_sp, CRS("+proj=utm +zone=12"))

# Kriging
fl_land_113_S2013_krig <- autoKrige(var~1, fl_land_113_S2013_sp, df_113_grid_sp)

# Format and plot with ggplot2
fl_land_113_S2013_krig_df <- as.data.frame(fl_land_113_S2013_krig$krige_output)
fl_land_113_S2013_krig_df$var1.pred[fl_land_113_S2013_krig_df$var1.pred<0] <- 0 # We can't have negative values, so we replace them with 0.
fl_land_113_S2013_sf <- st_as_sf(fl_land_113_S2013_sp)

fl_land_113_S2013_plot <- ggplot(data = fl_land_113_S2013_sf) + 
  geom_sf(aes(fill=var), shape = 21, size = 3, stroke = 1) +
  geom_tile(data =as.data.frame(fl_land_113_S2013_krig_df), aes(x = coords.x1, y = coords.x2, fill=var1.pred), alpha = 0.5) +
  theme(legend.position = "none", text=element_text(size = 12),
        # legend.background = element_rect(colour='white'),
        # legend.title = element_text(colour="black", size = 10),
        # legend.text  = element_text(colour="black", size = 10),
        # legend.key.width = unit(0.3,"cm"),
        # legend.key.height = unit(1.5, "cm"),
        axis.text.x=element_text(colour="black", size = 10, angle = 0, vjust = 0.5),  
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, size = 0.5),
        panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "white"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(name=NULL, colours = rev(viridis(10)), limits = c(0,5)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-112.522, -112.518, -112.514), labels = c("","","")) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("") +
  ylab("Latitude") +
  ggtitle("Spring 2013")

# Flowering landscape site 113, fall 2013
fl_land_113_F2013 <- df_tree_data[which(df_tree_data$season == "F2013" & df_tree_data$site == "113" & df_tree_data$longitude < -112.50),c("longitude", "latitude", "flowering_landscape_250_125")]
fl_land_113_F2013 <- na.omit(fl_land_113_F2013)
fl_land_113_F2013_sp <- fl_land_113_F2013
colnames(fl_land_113_F2013_sp) <- c("x","y","var")
coordinates(fl_land_113_F2013_sp) <- ~x+y
proj4string(fl_land_113_F2013_sp) <- CRS("+init=epsg:4326")
fl_land_113_F2013_sp <- spTransform(fl_land_113_F2013_sp, CRS("+proj=utm +zone=12"))

# Kriging
fl_land_113_F2013_krig <- autoKrige(var~1, fl_land_113_F2013_sp, df_113_grid_sp)

# Format and plot with ggplot2
fl_land_113_F2013_krig_df <- as.data.frame(fl_land_113_F2013_krig$krige_output)
fl_land_113_F2013_krig_df$var1.pred[fl_land_113_F2013_krig_df$var1.pred<0] <- 0 # We can't have negative values, so we replace them with 0.
fl_land_113_F2013_sf <- st_as_sf(fl_land_113_F2013_sp)

fl_land_113_F2013_plot <- ggplot(data = fl_land_113_F2013_sf) + 
  geom_sf(aes(fill=var), shape = 21, size = 3, stroke = 1) +
  geom_tile(data =as.data.frame(fl_land_113_F2013_krig_df), aes(x = coords.x1, y = coords.x2, fill=var1.pred), alpha = 0.5) +
  theme(legend.position = "none", text=element_text(size = 12),
        # legend.background = element_rect(colour='white'),
        # legend.title = element_text(colour="black", size = 10),
        # legend.text  = element_text(colour="black", size = 10),
        # legend.key.width = unit(0.3,"cm"),
        # legend.key.height = unit(1.5, "cm"),
        axis.text.x=element_text(colour="black", size = 10, angle = 0, vjust = 0.5),  
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, size = 0.5),
        panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "white"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(name=NULL, colours = rev(viridis(10)), limits = c(0,5)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-112.522, -112.518, -112.514), labels = c("","","")) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("") +
  ylab("Latitude") +
  ggtitle("Fall 2013")

# Flowering landscape site 113, spring 2014
fl_land_113_S2014 <- df_tree_data[which(df_tree_data$season == "S2014" & df_tree_data$site == "113" & df_tree_data$longitude < -112.50),c("longitude", "latitude", "flowering_landscape_250_125")]
fl_land_113_S2014 <- na.omit(fl_land_113_S2014)
fl_land_113_S2014_sp <- fl_land_113_S2014
colnames(fl_land_113_S2014_sp) <- c("x","y","var")
coordinates(fl_land_113_S2014_sp) <- ~x+y
proj4string(fl_land_113_S2014_sp) <- CRS("+init=epsg:4326")
fl_land_113_S2014_sp <- spTransform(fl_land_113_S2014_sp, CRS("+proj=utm +zone=12"))

# Kriging
fl_land_113_S2014_krig <- autoKrige(var~1, fl_land_113_S2014_sp, df_113_grid_sp)

# Format and plot with ggplot2
fl_land_113_S2014_krig_df <- as.data.frame(fl_land_113_S2014_krig$krige_output)
fl_land_113_S2014_krig_df$var1.pred[fl_land_113_S2014_krig_df$var1.pred<0] <- 0 # We can't have negative values, so we replace them with 0.
fl_land_113_S2014_sf <- st_as_sf(fl_land_113_S2014_sp)

fl_land_113_S2014_plot <- ggplot(data = fl_land_113_S2014_sf) + 
  geom_sf(aes(fill=var), shape = 21, size = 3, stroke = 1) +
  geom_tile(data =as.data.frame(fl_land_113_S2014_krig_df), aes(x = coords.x1, y = coords.x2,  fill=var1.pred), alpha = 0.5) +
  theme(legend.position = "right", text=element_text(size = 12),
        legend.background = element_rect(colour='white'),
        legend.title = element_text(colour="black", size = 10),
        legend.text  = element_text(colour="black", size = 10, angle = -90, hjust = 0.5),
        legend.key.width = unit(0.3,"cm"),
        legend.key.height = unit(1.5, "cm"),
        axis.text.x=element_text(colour="black", size = 10, angle = 0, vjust = 0.5),  
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, size = 0.5),
        panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "white"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(name=NULL, colours = rev(viridis(10)), limits = c(0,5)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-112.522, -112.518, -112.514), labels = c("","","")) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("") +
  ylab("Latitude") +
  ggtitle("Spring 2014")


# Syconium landscape site 113, spring 2012
sy_land_113_F2012 <- df_tree_data[which(df_tree_data$season == "F2012" & df_tree_data$site == "113" & df_tree_data$longitude < -112.50),c("longitude", "latitude", "syconium_landscape_250_125")]
sy_land_113_F2012 <- na.omit(sy_land_113_F2012)
sy_land_113_F2012_sp <- sy_land_113_F2012
colnames(sy_land_113_F2012_sp) <- c("x","y","var")
coordinates(sy_land_113_F2012_sp) <- ~x+y
proj4string(sy_land_113_F2012_sp) <- CRS("+init=epsg:4326")
sy_land_113_F2012_sp <- spTransform(sy_land_113_F2012_sp, CRS("+proj=utm +zone=12"))

# Kriging
sy_land_113_F2012_krig <- autoKrige(var~1, sy_land_113_F2012_sp, df_113_grid_sp)

# Format and plot with ggplot2
sy_land_113_F2012_krig_df <- as.data.frame(sy_land_113_F2012_krig$krige_output)
sy_land_113_F2012_krig_df$var1.pred[sy_land_113_F2012_krig_df$var1.pred<0] <- 0 # We can't have negative values, so we replace them with 0.
sy_land_113_F2012_sf <- st_as_sf(sy_land_113_F2012_sp)

sy_land_113_F2012_plot <- ggplot(data = sy_land_113_F2012_sf) + 
  geom_sf(aes(fill=var), shape = 21, size = 3, stroke = 1) +
  geom_tile(data =as.data.frame(sy_land_113_F2012_krig_df), aes(x = coords.x1, y = coords.x2, fill=var1.pred), alpha = 0.5) +
  theme(legend.position = "none", text=element_text(size = 12),
        # legend.background = element_rect(colour='white'),
        # legend.title = element_text(colour="black", size = 10),
        # legend.text  = element_text(colour="black", size = 10),
        # legend.key.width = unit(0.3,"cm"),
        # legend.key.height = unit(1.5, "cm"),
        axis.text.x=element_text(colour="black", size = 10, angle = 0, vjust = 0.5),  
        axis.text.y=element_text(colour="black", size = 10, angle = 90, hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, size = 0.5),
        panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "white"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(name=NULL, colours = rev(viridis(10)), limits = c(0,602)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-112.522, -112.518, -112.514)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("Fall 2012")

# Syconium landscape site 113, spring 2013
sy_land_113_S2013 <- df_tree_data[which(df_tree_data$season == "S2013" & df_tree_data$site == "113" & df_tree_data$longitude < -112.50),c("longitude", "latitude", "syconium_landscape_250_125")]
sy_land_113_S2013 <- na.omit(sy_land_113_S2013)
sy_land_113_S2013_sp <- sy_land_113_S2013
colnames(sy_land_113_S2013_sp) <- c("x","y","var")
coordinates(sy_land_113_S2013_sp) <- ~x+y
proj4string(sy_land_113_S2013_sp) <- CRS("+init=epsg:4326")
sy_land_113_S2013_sp <- spTransform(sy_land_113_S2013_sp, CRS("+proj=utm +zone=12"))

# Kriging
sy_land_113_S2013_krig <- autoKrige(var~1, sy_land_113_S2013_sp, df_113_grid_sp)

# Format and plot with ggplot2
sy_land_113_S2013_krig_df <- as.data.frame(sy_land_113_S2013_krig$krige_output)
sy_land_113_S2013_krig_df$var1.pred[sy_land_113_S2013_krig_df$var1.pred<0] <- 0 # We can't have negative values, so we replace them with 0.
sy_land_113_S2013_sf <- st_as_sf(sy_land_113_S2013_sp)

sy_land_113_S2013_plot <- ggplot(data = sy_land_113_S2013_sf) + 
  geom_sf(aes(fill=var), shape = 21, size = 3, stroke = 1) +
  geom_tile(data =as.data.frame(sy_land_113_S2013_krig_df), aes(x = coords.x1, y = coords.x2, fill=var1.pred), alpha = 0.5) +
  theme(legend.position = "none", text=element_text(size = 12),
        # legend.background = element_rect(colour='white'),
        # legend.title = element_text(colour="black", size = 10),
        # legend.text  = element_text(colour="black", size = 10),
        # legend.key.width = unit(0.3,"cm"),
        # legend.key.height = unit(1.5, "cm"),
        axis.text.x=element_text(colour="black", size = 10, angle = 0, vjust = 0.5),  
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, size = 0.5),
        panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "white"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(name=NULL, colours = rev(viridis(10)), limits = c(0,602)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-112.522, -112.518, -112.514)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("Spring 2013")

# Syconium landscape site 113, fall 2013
sy_land_113_F2013 <- df_tree_data[which(df_tree_data$season == "F2013" & df_tree_data$site == "113" & df_tree_data$longitude < -112.50),c("longitude", "latitude", "syconium_landscape_250_125")]
sy_land_113_F2013 <- na.omit(sy_land_113_F2013)
sy_land_113_F2013_sp <- sy_land_113_F2013
colnames(sy_land_113_F2013_sp) <- c("x","y","var")
coordinates(sy_land_113_F2013_sp) <- ~x+y
proj4string(sy_land_113_F2013_sp) <- CRS("+init=epsg:4326")
sy_land_113_F2013_sp <- spTransform(sy_land_113_F2013_sp, CRS("+proj=utm +zone=12"))

# Kriging
sy_land_113_F2013_krig <- autoKrige(var~1, sy_land_113_F2013_sp, df_113_grid_sp)

# Format and plot with ggplot2
sy_land_113_F2013_krig_df <- as.data.frame(sy_land_113_F2013_krig$krige_output)
sy_land_113_F2013_krig_df$var1.pred[sy_land_113_F2013_krig_df$var1.pred<0] <- 0 # We can't have negative values, so we replace them with 0.
sy_land_113_F2013_sf <- st_as_sf(sy_land_113_F2013_sp)

sy_land_113_F2013_plot <- ggplot(data = sy_land_113_F2013_sf) + 
  geom_sf(aes(fill=var), shape = 21, size = 3, stroke = 1) +
  geom_tile(data =as.data.frame(sy_land_113_F2013_krig_df), aes(x = coords.x1, y = coords.x2, fill=var1.pred), alpha = 0.5) +
  theme(legend.position = "none", text=element_text(size = 12),
        # legend.background = element_rect(colour='white'),
        # legend.title = element_text(colour="black", size = 10),
        # legend.text  = element_text(colour="black", size = 10),
        # legend.key.width = unit(0.3,"cm"),
        # legend.key.height = unit(2, "cm"),
        axis.text.x=element_text(colour="black", size = 10, angle = 0, vjust = 0.5),  
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, size = 0.5),
        panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "white"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(name=NULL, colours = rev(viridis(10)), limits = c(0,602)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-112.522, -112.518, -112.514)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("Fall 2013")

# Syconium landscape site 113, spring 2014
sy_land_113_S2014 <- df_tree_data[which(df_tree_data$season == "S2014" & df_tree_data$site == "113" & df_tree_data$longitude < -112.50),c("longitude", "latitude", "syconium_landscape_250_125")]
sy_land_113_S2014 <- na.omit(sy_land_113_S2014)
sy_land_113_S2014_sp <- sy_land_113_S2014
colnames(sy_land_113_S2014_sp) <- c("x","y","var")
coordinates(sy_land_113_S2014_sp) <- ~x+y
proj4string(sy_land_113_S2014_sp) <- CRS("+init=epsg:4326")
sy_land_113_S2014_sp <- spTransform(sy_land_113_S2014_sp, CRS("+proj=utm +zone=12"))

# Kriging
sy_land_113_S2014_krig <- autoKrige(var~1, sy_land_113_S2014_sp, df_113_grid_sp)

# Format and plot with ggplot2
sy_land_113_S2014_krig_df <- as.data.frame(sy_land_113_S2014_krig$krige_output)
sy_land_113_S2014_krig_df$var1.pred[sy_land_113_S2014_krig_df$var1.pred<0] <- 0 # We can't have negative values, so we replace them with 0.
sy_land_113_S2014_sf <- st_as_sf(sy_land_113_S2014_sp)

sy_land_113_S2014_plot <- ggplot(data = sy_land_113_S2014_sf) + 
  geom_sf(aes(fill=var), shape = 21, size = 3, stroke = 1) +
  geom_tile(data =as.data.frame(sy_land_113_S2014_krig_df), aes(x = coords.x1, y = coords.x2, fill=var1.pred), alpha = 0.5) +
  theme(legend.position = "right", text=element_text(size = 12),
        legend.background = element_rect(colour='white'),
        legend.title = element_text(colour="black", size = 10),
        legend.text  = element_text(colour="black", size = 10, angle = -90, hjust = 0.5),
        legend.key.width = unit(0.3,"cm"),
        legend.key.height = unit(1.5, "cm"),
        axis.text.x=element_text(colour="black", size = 10, angle = 0, vjust = 0.5),  
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, size = 0.5),
        panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "white"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(name=NULL, colours = rev(viridis(10)), limits = c(0,602)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-112.522, -112.518, -112.514)) + 
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("Spring 2014")


# Figure 9 --------------------

# Combine the plots above
ggdraw() +
  draw_plot(fl_land_113_F2012_plot, x = 0, y = 0.5, width = 0.25, height = 0.5) +
  draw_plot(fl_land_113_S2013_plot, x = 0.25, y = 0.5, width = 0.25, height = 0.5) +
  draw_plot(fl_land_113_F2013_plot, x = 0.50, y = 0.5, width = 0.25, height = 0.5) +
  draw_plot(fl_land_113_S2014_plot, x = 0.75, y = 0.5, width = 0.25, height = 0.5) +
  draw_plot(sy_land_113_F2012_plot, x = 0, y = 0, width = 0.25, height = 0.5) +
  draw_plot(sy_land_113_S2013_plot, x = 0.25, y = 0, width = 0.25, height = 0.5) +
  draw_plot(sy_land_113_F2013_plot, x = 0.50, y = 0, width = 0.25, height = 0.5) +
  draw_plot(sy_land_113_S2014_plot, x = 0.75, y = 0, width = 0.25, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C", "D", "E", "F", "G", "H"), size = 15,
                  x = c(0.02, 0.27, 0.52, 0.75, 0.02, 0.27, 0.52, 0.75), y = c(1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5))

# Saved with ggsave(), width = 8.5 inches, height = 8.5 inches, dpi = 600.
ggsave(paste0(workdir,"figures/","fl&sy_landscape_site_113.png",sep=""), plot = last_plot(), device = "png", path = NULL, scale = 1, width = 16, height = 8, units = c("in"), dpi = 600, limitsize = TRUE)


# Clear working environment --------------------

rm(df_112,df_112_sp,df_112_grid,df_112_grid_sp,df_113,df_113_sp,df_113_grid,df_113_grid_sp)

rm(fl_land_112_S2013,fl_land_112_S2013_sp,fl_land_112_S2013_plot,fl_land_112_S2013_krig,fl_land_112_S2013_krig_df,fl_land_112_S2013_sf,
   fl_land_112_F2013,fl_land_112_F2013_sp,fl_land_112_F2013_plot,fl_land_112_F2013_krig,fl_land_112_F2013_krig_df,fl_land_112_F2013_sf,
   fl_land_112_S2014,fl_land_112_S2014_sp,fl_land_112_S2014_plot,fl_land_112_S2014_krig,fl_land_112_S2014_krig_df,fl_land_112_S2014_sf,
   fl_land_113_F2012,fl_land_113_F2012_sp,fl_land_113_F2012_plot,fl_land_113_F2012_krig,fl_land_113_F2012_krig_df,fl_land_113_F2012_sf,
   fl_land_113_S2013,fl_land_113_S2013_sp,fl_land_113_S2013_plot,fl_land_113_S2013_krig,fl_land_113_S2013_krig_df,fl_land_113_S2013_sf,
   fl_land_113_F2013,fl_land_113_F2013_sp,fl_land_113_F2013_plot,fl_land_113_F2013_krig,fl_land_113_F2013_krig_df,fl_land_113_F2013_sf,
   fl_land_113_S2014,fl_land_113_S2014_sp,fl_land_113_S2014_plot,fl_land_113_S2014_krig,fl_land_113_S2014_krig_df,fl_land_113_S2014_sf)

rm(sy_land_112_S2013,sy_land_112_S2013_sp,sy_land_112_S2013_plot,sy_land_112_S2013_krig,sy_land_112_S2013_krig_df,sy_land_112_S2013_sf,
   sy_land_112_F2013,sy_land_112_F2013_sp,sy_land_112_F2013_plot,sy_land_112_F2013_krig,sy_land_112_F2013_krig_df,sy_land_112_F2013_sf,
   sy_land_112_S2014,sy_land_112_S2014_sp,sy_land_112_S2014_plot,sy_land_112_S2014_krig,sy_land_112_S2014_krig_df,sy_land_112_S2014_sf,
   sy_land_113_F2012,sy_land_113_F2012_sp,sy_land_113_F2012_plot,sy_land_113_F2012_krig,sy_land_113_F2012_krig_df,sy_land_113_F2012_sf,
   sy_land_113_S2013,sy_land_113_S2013_sp,sy_land_113_S2013_plot,sy_land_113_S2013_krig,sy_land_113_S2013_krig_df,sy_land_113_S2013_sf,
   sy_land_113_F2013,sy_land_113_F2013_sp,sy_land_113_F2013_plot,sy_land_113_F2013_krig,sy_land_113_F2013_krig_df,sy_land_113_F2013_sf,
   sy_land_113_S2014,sy_land_113_S2014_sp,sy_land_113_S2014_plot,sy_land_113_S2014_krig,sy_land_113_S2014_krig_df,sy_land_113_S2014_sf)

rm(fl_land_112_S2013_plot,fl_land_112_F2013_plot,fl_land_112_S2014_plot,sy_land_112_S2013_plot,sy_land_112_F2013_plot,sy_land_112_S2014_plot,
   fl_land_113_S2013_plot,fl_land_113_F2013_plot,fl_land_113_S2014_plot,sy_land_113_F2012_plot,sy_land_113_S2013_plot,sy_land_113_F2013_plot,
   sy_land_113_S2014_plot)



# ---------------------------------------------------------------------------- #
# Figures: abiotic factors, daily weather and monthly climate site-specific data
# ---------------------------------------------------------------------------- #



# Colors are:
# "#A3A500": weather stations
# "#F8766D": DAYMET
# "#00B0F6": TerraClimate
# "#E76BF3": CHIRPS
# "#00BF7D": MODIS LST

# Visualize
est_name <- c("weather stations"=1/5,"DAYMET"=1/5,"TerraClimate"=1/5,"CHIRPS"=1/5,"MODIS LST"=1/5)
est_col <- c("#A3A500","#F8766D","#00B0F6","#E76BF3","#00BF7D")
par(lwd=2)
pie(est_name, col = est_col,main = "Corresponding colors", border="white")
par(lwd=2) # Resets lwd size.


# ____________________________________________________________________________ #
# Weather stations data - weather/climate data sets comparison


weather_station_all_estimates <- read.csv("weather_station_all_estimates_2011.csv")
weather_station_all_estimates$weather_station <- factor(weather_station_all_estimates$weather_station,
                                                        levels = c("Chapala","Bahia de los Angeles","San Francisco de la Sierra","Santa Agueda","Ojo de Agua","Loreto DGE","San Ramon","Los Robles","San Bartolo"))

# Format data
weather_station_tmax_melt <- melt(weather_station_all_estimates, id.vars = c("weather_station","date"), measure.vars = , c("tmax","MODIS_lst","DM_tmax","TC_tmax"))
weather_station_tmin_melt <- melt(weather_station_all_estimates, id.vars = c("weather_station","date"), measure.vars = , c("tmin","DM_tmin","TC_tmin"))
weather_station_prec_melt <- melt(weather_station_all_estimates, id.vars = c("weather_station","date"), measure.vars = , c("prcp","CHIRPS_prec","DM_prec","TC_prec"))
weather_station_tmax_melt$date <- as.Date(weather_station_tmax_melt$date)
weather_station_tmin_melt$date <- as.Date(weather_station_tmin_melt$date)
weather_station_prec_melt$date <- as.Date(weather_station_prec_melt$date)


# Options to plot weather/climate data
# geom_line(aes(y = value)) + # Careful with geom_line(), TerraClimate data have missing value, geom_line won't connect the points. We need to subset TC data and omit NA's so dotes are connected
# geom_ribbon(aes(ymin = value - 1, ymax = value + 1, fill = variable), alpha = 0.8) + # I like to add this to geom_line because it highlights a little the line and the "jumps" from different value.
# geom_smooth(aes(y = value), linetype="solid", size=1, se = FALSE) + # Alternatively to geom_line, this display a fit to the data, optional, add:, span = 0.1. I prefer this to represent precipitation, they are hard to see with geom_line().
  

# Supplementary Figure S5 --------------------

ggplot(weather_station_tmax_melt, aes(x = date, colour=variable)) + 
  geom_line(aes(y = value), linetype="solid", size=0.5) + 
  geom_ribbon(aes(ymin = value - 0.5, ymax = value + 0.5, fill = variable), alpha = 0.8) +
  geom_line(data = na.omit(subset(weather_station_tmax_melt, variable =="TC_tmax")), aes(x = date, colour=variable, y = value), size =0.5) + 
  geom_ribbon(data = na.omit(subset(weather_station_tmax_melt, variable =="TC_tmax")), aes(x = date, colour=variable, ymin = value - 0.2, ymax = value + 0.2, fill = variable), alpha = 0.8) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(15, 10, 10, 10, "pt"),
        text=element_text(size = 16),
        axis.text.x=element_text(colour="black", size = 9, angle = 90, vjust = 1), # Change to 0 or 90 depending on facet_wrap
        axis.text.y=element_text(colour="black", size = 9, angle = 0, hjust = 0),
        legend.position = "right", # change depending on facet_wrap ("right" or c(0.15,0.9))
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text( size=10),
        strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        strip.text.x = element_text(size = 9)) +
  scale_y_continuous("Maximum Temperature (°C)", breaks=seq(0,55,10), limits = c(0,55), expand=c(0,0)) +
  scale_x_date("Time", labels = date_format("%Y-%m-%d"), limits = c(as.Date("2011-01-01"), as.Date("2011-12-31")), expand=c(0,0), breaks = date_breaks("2 months")) +
  scale_fill_manual(values=c("#A3A500","#00BF7D","#F8766D","#00B0F6"), name = "Estimates source", breaks = levels(weather_station_tmax_melt$variable), labels = c("Weather station", "MODIS LST", "DAYMET", "TerraClimate")) +
  scale_colour_manual(values=c("#A3A500","#00BF7D","#F8766D","#00B0F6")) +
  facet_wrap(~ weather_station, ncol=3) + # Add/remove for all data or data per stations
  guides(color = FALSE)

# Saved with ggsave(), width = 10 inches, height = 7 inches, dpi = 600.
ggsave(paste0(workdir,"figures/", "weather_station_tmax.png",sep=""), plot = last_plot(), device = "png", path = NULL, scale = 1, width = 10, height = 7, units = c("in"), dpi = 600, limitsize = TRUE)



# Supplementary Figure S6 --------------------

ggplot(weather_station_tmin_melt, aes(x = date, colour=variable)) + 
  geom_line(aes(y = value), linetype="solid", size=0.5) +
  geom_ribbon(aes(ymin = value - 0.5, ymax = value + 0.5, fill = variable), alpha = 0.8) + 
  geom_line(data = na.omit(subset(weather_station_tmin_melt, variable = "TC_tmin")), aes(x = date, colour=variable,y = value), linetype="solid", size=0.5) +
  geom_ribbon(data = na.omit(subset(weather_station_tmin_melt, variable = "TC_tmin")), aes(x = date, colour=variable,ymin = value - 0.2, ymax = value + 0.2, fill = variable), alpha = 0.8) + 
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(15, 10, 10, 10, "pt"),
        text=element_text(size = 16),
        axis.text.x=element_text(colour="black", size = 9, angle = 90, vjust = 1), # Change to 0 or 90 depending on facet_wrap
        axis.text.y=element_text(colour="black", size = 9, angle = 0, hjust = 0),
        legend.position = "right", # change depending on facet_wrap ("right" or c(0.15,0.9))
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text( size=10),
        strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        strip.text.x = element_text(size = 9)) +
  scale_y_continuous("Minimum Temperature (°C)", breaks=seq(-10,35,10), limits = c(-10,35), expand=c(0,0)) +
  scale_x_date("Time", labels = date_format("%Y-%m-%d"), limits = c(as.Date("2011-01-01"), as.Date("2011-12-31")), expand=c(0,0), breaks = date_breaks("2 months")) +
  scale_fill_manual(values=c("#A3A500","#F8766D","#00B0F6"), name = "Estimates source", breaks = levels(weather_station_tmin_melt$variable), labels = c("Weather station", "DAYMET", "TerraClimate")) +
  scale_colour_manual(values=c("#A3A500","#F8766D","#00B0F6")) +
  facet_wrap(~ weather_station, ncol=3) + # Add/remove for all data or data per stations
  guides(color = FALSE)


# Saved with ggsave(), width = 10 inches, height = 7 inches, dpi = 600.
ggsave(paste0(workdir,"figures/", "weather_station_tmin.png", sep=""), plot = last_plot(), device = "png", path = NULL, scale = 1, width = 10, height = 7, units = c("in"), dpi = 600, limitsize = TRUE)



# Supplementary Figure S7 --------------------

ggplot(weather_station_prec_melt, aes(x = date, colour=variable)) + 
  geom_ribbon(aes(ymin = value - 1, ymax = value + 1, fill = variable), alpha = 0.8) + 
  geom_line(aes(y = value), linetype="solid", size=1, alpha=0.7) + 
  geom_bar(data = na.omit(subset(weather_station_prec_melt, variable =="TC_prec")), aes(x = date, y = value), stat="identity", alpha = 0.3, color = "#00B0F6") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(15, 10, 10, 10, "pt"),
        text=element_text(size = 16),
        axis.text.x=element_text(colour="black", size = 9, angle = 90, vjust = 1), # Change to 0 or 90 depending on facet_wrap
        axis.text.y=element_text(colour="black", size = 9, angle = 0, hjust = 0),
        legend.position = "right", # change depending on facet_wrap ("right" or c(0.15,0.9))
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text( size=10),
        strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        strip.text.x = element_text(size = 9)) +
  scale_y_continuous("Precipitation (mm)", breaks=seq(0,65,10), limits = c(0,65), expand=c(0,0)) +
  scale_x_date("Time", labels = date_format("%Y-%m-%d"), limits = c(as.Date("2011-01-01"), as.Date("2011-12-31")), expand=c(0,0), breaks = date_breaks("2 months")) +
  scale_fill_manual(values=c("#A3A500","#E76BF3","#F8766D","#00B0F6"),name = "Estimates source", breaks = levels(weather_station_prec_melt$variable), labels = c("Weather station", "CHIRPS", "DAYMET", "Terra Climate"), drop=FALSE) +
  scale_colour_manual(values=c("#A3A500","#E76BF3","#F8766D","#00B0F6")) +
  facet_wrap(~ weather_station, ncol=3) + # Add/remove for all data or data per stations
  guides(color = FALSE)

# Saved with ggsave(), width = 10 inches, height = 7 inches, dpi = 600.
ggsave(paste0(workdir,"figures/","weather_station_prec.png", sep=""), plot = last_plot(), device = "png", path = NULL, scale = 1, width = 10, height = 7, units = c("in"), dpi = 600, limitsize = TRUE)


# ____________________________________________________________________________ #
# Abiotic variables correlation for 2012-2014


# Import
Fpet_all_estimates <- read.csv(paste0(workdir,"generated_data/","Fpet_all_estimates_2012-2014.csv"))
Fpet_all_estimates$date <- as.Date(Fpet_all_estimates$date)


# Colors for correlation plot
corrPalette <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", 
																	"#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
																	"#4393C3", "#2166AC", "#053061"))

# Compute correlation and plot correlation heatmaps

# Supplementary Figure S8 -------------------- 

daily_abiotic_variable_corr <- cor(na.omit(Fpet_all_estimates[,6:10]), method="spearman") 
colnames(daily_abiotic_variable_corr) <- c("DM Prec.", "DM Tmax", "DM Tmin", "DM VP", "CH. Prec.")
rownames(daily_abiotic_variable_corr) <- c("DM Prec.", "DM Tmax", "DM Tmin", "DM VP", "CH. Prec.")
diag(daily_abiotic_variable_corr) = NA
daily_abiotic_variable_corr_sp <- rcorr(as.matrix(daily_abiotic_variable_corr), type = "spearman") # the function rcorr() resolve the issue above.

png(width = 10, height=8,file = paste0(workdir,"figures/","abiotic_daily_variable_correlation.png", sep=""), units = c("in"), res = 300,  type = "cairo")
corrplot(daily_abiotic_variable_corr, method = "color", type = "upper", col=corrPalette(200), cl.pos = "r",
				 tl.col="black", tl.offset = 0.2, tl.srt = 50 , order = NULL, diag = TRUE, na.label = "NA", 
				 addCoef.col = "darkgrey", sig.level = 0.05, insig = "pch", p.mat = daily_abiotic_variable_corr_sp$P)
dev.off()

# Supplementary Figure S9 --------------------

monthly_abiotic_variable_corr <- cor(na.omit(Fpet_all_estimates[,11:15]), method="spearman") 
colnames(monthly_abiotic_variable_corr) <- c("TC Prec.", "TC Tmax", "TC Tmin", "TC VP", "TC wind")
rownames(monthly_abiotic_variable_corr) <- c("TC Prec.", "TC Tmax", "TC Tmin", "TC VP", "TC wind")
diag(monthly_abiotic_variable_corr) = NA
monthly_abiotic_variable_corr_sp <- rcorr(as.matrix(monthly_abiotic_variable_corr), type = "spearman") # the function rcorr() resolve the issue above.

png(width = 10, height=8,file = paste0(workdir,"figures/","abiotic_monthly_variable_correlation.png", sep=""), units = c("in"), res = 300,  type = "cairo")
corrplot(monthly_abiotic_variable_corr, method = "color", type = "upper", col=corrPalette(200), cl.pos = "r",
				 tl.col="black", tl.offset = 0.2, tl.srt = 50 , order = NULL, diag = TRUE, na.label = "NA", 
				 addCoef.col = "darkgrey", sig.level = 0.05, insig = "pch", p.mat = monthly_abiotic_variable_corr_sp$P)
dev.off()

# Clear memory
rm(daily_abiotic_variable_corr, daily_abiotic_variable_corr_sp, monthly_abiotic_variable_corr, monthly_abiotic_variable_corr_sp)



# ____________________________________________________________________________ #
# Maximum temperature, precipitation and wind speed variations for 2012-2014



# Format
Fpet_all_estimates_T <- Fpet_all_estimates[,c("site","date","latitude","longitude","altitude","DM_tmax","TC_tmax")]
Fpet_all_estimates_T_melt <- melt(Fpet_all_estimates_T, id.vars = c("site","date","latitude","longitude","altitude"))
Fpet_all_estimates_P <- Fpet_all_estimates[,c("site","date","latitude","longitude","altitude","DM_prec", "CHIRPS_prec","TC_prec")]
Fpet_all_estimates_P_melt <- melt(Fpet_all_estimates_P, id.vars = c("site","date","latitude","longitude","altitude"))
Fpet_all_estimates_WS <- Fpet_all_estimates[,c("site","date","latitude","longitude","altitude","TC_ws")]
Fpet_all_estimates_WS_melt <- melt(Fpet_all_estimates_WS, id.vars = c("site","date","latitude","longitude","altitude"))


# Supplementary Figure S10 -------------------- 

ggplot(Fpet_all_estimates_T_melt, aes(x = date, colour=variable)) + 
  geom_line(aes(y = value), linetype="solid", size=0.5) + 
  geom_ribbon(aes(ymin = value - 1, ymax = value + 1, fill = variable), alpha = 0.8) + 
  geom_line(data = na.omit(subset(Fpet_all_estimates_T_melt, variable == "TC_tmax")), aes(x = date, colour=variable, y = value), size =1.2) + 
  geom_ribbon(data = na.omit(Fpet_all_estimates_T_melt), aes(x = date, colour=variable, ymin = value - 0.5, ymax = value + 0.5, fill = variable), alpha = 0.8) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(15, 10, 10, 10, "pt"),
        text=element_text(size = 16),
        axis.text.x=element_text(colour="black", size = 9, angle = 90, vjust = 1), # Change to 0 or 90 depending on facet_wrap
        axis.text.y=element_text(colour="black", size = 9, angle = 0, hjust = 0),
        legend.position = "right", # change depending on facet_wrap ("right" or c(0.15,0.9))
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text( size=10),
        strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        strip.text.x = element_text(size = 9)) +
  scale_y_continuous("Maximum Temperature (°C)", breaks=seq(0,55,10), limits = c(0,55), expand=c(0,0)) +
  scale_x_date("Time", labels = date_format("%Y-%m-%d"), limits = c(as.Date("2012-01-01"), as.Date("2014-12-31")), expand=c(0,0), breaks = date_breaks("4 months")) +
  scale_fill_manual(name = "Estimates source", breaks = levels(Fpet_all_estimates_T_melt$variable), labels = c("DAYMET", "TerraClimate"), values = c("#F8766D","#00B0F6"), drop=FALSE) +
  scale_colour_manual(values = c("#F8766D","#00B0F6")) +
  facet_wrap(~ site, ncol=3) + # Add/remove for all data or data per stations
  guides(color = FALSE)

# Saved with ggsave(), width = 10 inches, height = 7 inches, dpi = 600.
ggsave(paste0(workdir, "figures/", "Tmax_2012-2014.png", sep=""), plot = last_plot(), device = "png", path = NULL, scale = 1, width = 10, height = 7, units = c("in"), dpi = 600, limitsize = TRUE)


# Supplementary Figure S11 -------------------- 

# Plot
ggplot(Fpet_all_estimates_P_melt, aes(x = date, colour=variable)) + 
  geom_ribbon(aes(ymin = value - 1, ymax = value + 1, fill = variable), alpha = 0.8) + 
  geom_line(aes(y = value), linetype="solid", size=1, alpha=0.7) + # For a better fit to the data (but with more data, not spit by stations) add:, span = 0.1 .
  geom_bar(data = na.omit(subset(Fpet_all_estimates_P_melt, variable == "TC_prec")), aes(x = date, y = value), stat="identity", alpha = 0.2, color = "#00994C") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(15, 10, 10, 10, "pt"),
        text=element_text(size = 16),
        axis.text.x=element_text(colour="black", size = 9, angle = 90, vjust = 1), # Change to 0 or 90 depending on facet_wrap
        axis.text.y=element_text(colour="black", size = 9, angle = 0, hjust = 0),
        legend.position = "right", # change depending on facet_wrap ("right" or c(0.15,0.9))
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text( size=10),
        strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        strip.text.x = element_text(size = 9)) +
  scale_y_continuous("Precipitation (mm)", breaks=seq(0,330,50),limits = c(0,330), expand=c(0,0)) +
  scale_x_date("Time", labels = date_format("%Y-%m-%d"), limits = c(as.Date("2012-01-01"), as.Date("2014-12-31")), expand=c(0,0), breaks = date_breaks("4 months")) +
  scale_fill_manual(values=c("#00CCCC","#FF6666","#00994C"), name = "Estimates source", breaks = levels(Fpet_all_estimates_P_melt$variable), labels = c("DAYMET","CHIRPS","TerraClimate"),drop = FALSE) +
  scale_colour_manual(values=c("#00CCCC","#FF6666","#00994C")) +
  facet_wrap(~ site, ncol=3) + # Add/remove for all data or data per stations
  guides(color = FALSE)

# Saved with ggsave(), width = 10 inches, height = 7 inches, dpi = 600.
ggsave(paste0(workdir, "figures/", "Prec_2012-2014.png", sep=""), plot = last_plot(), device = "png", path = NULL, scale = 1, width = 10, height = 7, units = c("in"), dpi = 600, limitsize = TRUE)


# Supplementary Figure S12 -------------------- 

ggplot(na.omit(Fpet_all_estimates_WS_melt), aes(x = date, colour=variable)) + 
  geom_line(aes(y = value),stat="identity") +
  geom_ribbon(aes(ymin = value - 0.05, ymax = value + 0.05, fill = variable), alpha = 0.5) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(15, 10, 10, 10, "pt"),
        text=element_text(size = 16),
        axis.text.x=element_text(colour="black", size = 9, angle = 90, vjust = 1), # Change to 0 or 90 depending on facet_wrap
        axis.text.y=element_text(colour="black", size = 9, angle = 0, hjust = 0),
        legend.position = "right", # change depending on facet_wrap ("right" or c(0.15,0.9))
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text( size=10),
        strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        strip.text.x = element_text(size = 9)) +
  scale_y_continuous("Wind Speed (°C)", breaks=seq(1,4,0.5), limits = c(1,4), expand=c(0,0)) +
  scale_x_date("Time", labels = date_format("%Y-%m-%d"), limits = c(as.Date("2012-01-01"), as.Date("2014-12-31")), expand=c(0,0), breaks = date_breaks("4 months")) +
  scale_fill_manual(name = "Estimates source", breaks = levels(Fpet_all_estimates_WS_melt$variable), labels = c("TerraClimate"), values = c("#00994C")) +
  scale_colour_manual(values = c("#00994C")) +
  facet_wrap(~ site, ncol=3) + # Add/remove for all data or data per stations
  guides(color = FALSE)

# Saved with ggsave(), width = 10 inches, height = 7 inches, dpi = 600.
ggsave(paste0(workdir, "figures/", "WS_2012-2014.png", sep=""), plot = last_plot(), device = "png", path = NULL, scale = 1, width = 10, height = 7, units = c("in"), dpi = 600, limitsize = TRUE)


# ____________________________________________________________________________ #
# Maximum temperature, precipitation on a transect for 2012-2014


# Import data
transect_annual_data <- read.csv("generated_data/transect_annual_data_2012-2014.csv")

# Plots
transect_annual_tmax_2012 <- ggplot(transect_annual_data, aes( x = latitude)) +
  geom_line(aes( y = TC_tmax_2012 ), linetype="solid", colour="black", size=1) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(10, 10, 10, 10, "pt"),
        text=element_text(size = 16),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12, angle = 00, hjust = 0.5)) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"TC_tmax_2012"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"TC_tmax_2012"]+0.02), label = "158", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"TC_tmax_2012"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"TC_tmax_2012"]+0.02), label = "172", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"TC_tmax_2012"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"TC_tmax_2012"]+0.02), label = "112", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"TC_tmax_2012"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"TC_tmax_2012"]+0.02), label = "113", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"TC_tmax_2012"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"TC_tmax_2012"]+0.02), label = "95", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"TC_tmax_2012"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"TC_tmax_2012"]+0.02), label = "179", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"TC_tmax_2012"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"TC_tmax_2012"]+0.02), label = "201", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"TC_tmax_2012"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"TC_tmax_2012"]+0.02), label = "96", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"TC_tmax_2012"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"TC_tmax_2012"]+0.02), label = "70", fontface = "bold", color = "black", size = 4) +
  scale_y_continuous("Annual maximum\n temperature in 2012 (mm)", breaks=seq(0,50,1)) + 
  scale_x_reverse("",breaks=seq(23.5,30,0.5))

transect_annual_tmax_2013 <- ggplot(transect_annual_data, aes( x = latitude)) +
  geom_line(aes( y = TC_tmax_2013 ), linetype="solid", colour="black", size=1) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(10, 10, 10, 10, "pt"),
        text=element_text(size = 16),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12, angle = 00, hjust = 0.5)) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"TC_tmax_2013"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"TC_tmax_2013"]+0.02), label = "158", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"TC_tmax_2013"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"TC_tmax_2013"]+0.02), label = "172", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"TC_tmax_2013"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"TC_tmax_2013"]+0.02), label = "112", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"TC_tmax_2013"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"TC_tmax_2013"]+0.02), label = "113", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"TC_tmax_2013"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"TC_tmax_2013"]+0.02), label = "95", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"TC_tmax_2013"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"TC_tmax_2013"]+0.02), label = "179", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"TC_tmax_2013"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"TC_tmax_2013"]+0.02), label = "201", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"TC_tmax_2013"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"TC_tmax_2013"]+0.02), label = "96", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"TC_tmax_2013"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"TC_tmax_2013"]+0.02), label = "70", fontface = "bold", color = "black", size = 4) +
  scale_y_continuous("Annual maximum\n temperature in 2013 (mm)", breaks=seq(0,50,1)) + 
  scale_x_reverse("",breaks=seq(23.5,30,0.5))

transect_annual_tmax_2014 <- ggplot(transect_annual_data, aes( x = latitude)) +
  geom_line(aes( y = TC_tmax_2014 ), linetype="solid", colour="black", size=1) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(10, 10, 10, 10, "pt"),
        text=element_text(size = 16),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12, angle = 00, hjust = 0.5)) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"TC_tmax_2014"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"TC_tmax_2014"]+0.02), label = "158", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"TC_tmax_2014"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"TC_tmax_2014"]+0.02), label = "172", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"TC_tmax_2014"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"TC_tmax_2014"]+0.02), label = "112", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"TC_tmax_2014"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"TC_tmax_2014"]+0.02), label = "113", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"TC_tmax_2014"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"TC_tmax_2014"]+0.02), label = "95", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"TC_tmax_2014"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"TC_tmax_2014"]+0.02), label = "179", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"TC_tmax_2014"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"TC_tmax_2014"]+0.02), label = "201", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"TC_tmax_2014"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"TC_tmax_2014"]+0.02), label = "96", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"TC_tmax_2014"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"TC_tmax_2014"]+0.02), label = "70", fontface = "bold", color = "black", size = 4) +
  scale_y_continuous("Annual maximum\n temperature in 2014 (mm)", breaks=seq(0,50,1)) + 
  scale_x_reverse("Latitude",breaks=seq(23.5,30,0.5))

transect_annual_prec_2012 <- ggplot(transect_annual_data, aes( x = latitude)) +
  geom_line(aes( y = TC_prec_2012 ), linetype="solid", colour="black", size=1) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(10, 10, 10, 10, "pt"),
        text=element_text(size = 16),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12, angle = 00, hjust = 0.5)) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"TC_prec_2012"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"TC_prec_2012"]+0.25), label = "158", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"TC_prec_2012"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"TC_prec_2012"]+0.25), label = "172", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"TC_prec_2012"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"TC_prec_2012"]+0.25), label = "112", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"TC_prec_2012"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"TC_prec_2012"]+0.25), label = "113", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"TC_prec_2012"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"TC_prec_2012"]+0.25), label = "95", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"TC_prec_2012"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"TC_prec_2012"]+0.25), label = "179", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"TC_prec_2012"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"TC_prec_2012"]+0.25), label = "201", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"TC_prec_2012"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"TC_prec_2012"]+0.25), label = "96", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"TC_prec_2012"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"TC_prec_2012"]+0.25), label = "70", fontface = "bold", color = "black", size = 4) +
  scale_y_continuous("Annual cumulated\n precipitation in 2012 (mm)", breaks=seq(50,1000,50)) + 
  scale_x_reverse("",breaks=seq(23.5,30,0.5))

transect_annual_prec_2013 <- ggplot(transect_annual_data, aes( x = latitude)) +
  geom_line(aes( y = TC_prec_2013 ), linetype="solid", colour="black", size=1) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(10, 10, 10, 10, "pt"),
        text=element_text(size = 16),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12, angle = 00, hjust = 0.5)) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"TC_prec_2013"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"TC_prec_2013"]+0.25), label = "158", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"TC_prec_2013"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"TC_prec_2013"]+0.25), label = "172", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"TC_prec_2013"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"TC_prec_2013"]+0.25), label = "112", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"TC_prec_2013"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"TC_prec_2013"]+0.25), label = "113", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"TC_prec_2013"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"TC_prec_2013"]+0.25), label = "95", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"TC_prec_2013"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"TC_prec_2013"]+0.25), label = "179", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"TC_prec_2013"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"TC_prec_2013"]+0.25), label = "201", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"TC_prec_2013"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"TC_prec_2013"]+0.25), label = "96", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"TC_prec_2013"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"TC_prec_2013"]+0.25), label = "70", fontface = "bold", color = "black", size = 4) +
  scale_y_continuous("Annual cumulated\n precipitation in 2013 (mm)", breaks=seq(50,1000,50)) + 
  scale_x_reverse("",breaks=seq(23.5,30,0.5))

transect_annual_prec_2014 <- ggplot(transect_annual_data, aes( x = latitude)) +
  geom_line(aes( y = TC_prec_2014 ), linetype="solid", colour="black", size=1) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(10, 10, 10, 10, "pt"),
        text=element_text(size = 16),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12, angle = 00, hjust = 0.5)) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"TC_prec_2014"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==158),"TC_prec_2014"]+0.25), label = "158", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"TC_prec_2014"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==172),"TC_prec_2014"]+0.25), label = "172", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"TC_prec_2014"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==112),"TC_prec_2014"]+0.25), label = "112", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"TC_prec_2014"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==113),"TC_prec_2014"]+0.25), label = "113", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"TC_prec_2014"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==95),"TC_prec_2014"]+0.25), label = "95", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"TC_prec_2014"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==179),"TC_prec_2014"]+0.25), label = "179", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"TC_prec_2014"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==201),"TC_prec_2014"]+0.25), label = "201", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"TC_prec_2014"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==96),"TC_prec_2014"]+0.25), label = "96", fontface = "bold", color = "black", size = 4) +
  geom_point(aes(x = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"TC_prec_2014"])), size = 10.5, shape = 21, stroke = 2, colour = "black", fill = "white") +
  annotate(geom = "text", x = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"latitude"]), y = as.vector(transect_annual_data[which(transect_annual_data$pop==70),"TC_prec_2014"]+0.25), label = "70", fontface = "bold", color = "black", size = 4) +
  scale_y_continuous("Annual cumulated\n precipitation in 2014 (mm)", breaks=seq(50,1000,50)) + 
  scale_x_reverse("Latitude",breaks=seq(23.5,30,0.5))


# Figure 5 -------------------- 


# Combine plots
# dev.new(width = 10, height = 10, units = "in", dpi = 600, limitsize = TRUE)
# par(mfrow=c(1,1))

ggdraw() +
  draw_plot(transect_annual_tmax_2012, x = 0, y = 0.66, width = 0.5, height = 0.33) +
  draw_plot(transect_annual_tmax_2013, x = 0, y = 0.33, width = 0.5, height = 0.33) +
  draw_plot(transect_annual_tmax_2014, x = 0, y = 0, width = 0.5, height = 0.33) +
  draw_plot(transect_annual_prec_2012, x = 0.5, y = 0.66, width = 0.5, height = 0.33) +
  draw_plot(transect_annual_prec_2013, x = 0.5, y = 0.33, width = 0.5, height = 0.33) +
  draw_plot(transect_annual_prec_2014, x = 0.5, y = 0, width = 0.5, height = 0.33)

# dev.size(units = "in")
ggsave(paste0(workdir,"figures/","transect_annual_data.png", sep=""), plot = last_plot(), device = "png", path = NULL, scale = 1, width = 12, height = 16, units = c("in"), dpi = 600, limitsize = TRUE)



# ---------------------------------------------------------------------------- #
# Figure: Fig wasp interaction network
# ---------------------------------------------------------------------------- #



# Subset abundances for fig wasp interaction analysis
wasp_abund <- df_fig_analysis[,c("pollinators","LO1_f","SO1_f","SO2_f","heterandrium_1","heterandrium_2","ficicola","physothorax","sycophila")] 
colnames(wasp_abund) <- c("Pegoscapus sp.","Idarnes flavicolis sp.","Idarnes carme sp. 1","Idarnes carme sp. 2","Heterandrium sp. 1","Heterandrium sp. 2","Ficicola sp.","Physothorax sp.","Sycophila sp.")

# Run network models
wasp_network_offset <- PLNnetwork(Abundance ~ 1 + offset(log(Offset)), data = wasp_abund_for_PLN)

# Select best model
wasp_best_network_offset <- getBestModel(wasp_network_offset)

# Create igraph object
wasp_best_network_offset_graph <- plot(wasp_best_network_offset)

# Investigate igraph object
str(wasp_best_network_offset_graph)
vertex_attr(wasp_best_network_offset_graph)
edge_attr(wasp_best_network_offset_graph)


# Figure 10 --------------------

# Make nicer plot with ggraph
wasp_graph <- ggraph(wasp_best_network_offset_graph, layout = "linear", circular = TRUE) + 
	geom_edge_arc(aes(edge_width=width, edge_colour=color, alpha = 0.5)) + # Add:, label = round(width,digits = 2), label_size = 22, for labels; ,check_overlap = TRUE to prevent overlap
	geom_node_point(size = 6) +
	# coord_fixed() +
	theme_graph(background = "white") +
	theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0, "cm")) +
	scale_edge_colour_manual(values = c("#e72a32","#69c8ff"))
wasp_graph

# Add labels
wasp_graph + 	geom_node_text(aes(label=name), size = 5, colour="#007ba7", nudge_x = wasp_graph$data$x * .1, nudge_y = wasp_graph$data$y * .1)
# Labels are hard to work with, they will be added with Inkscape.

# Save
wasp_graph + coord_fixed()
ggsave(paste(workdir,"figures/","wasp_interaction_network_with_numbers.png", sep = ""), plot = last_plot(), device = "png", path = NULL, scale = 1, width = 6, height = 6, units = c("in"), dpi = 600, limitsize = TRUE)

# Partial correlation matrix
wasp_best_network_offset_corrplot_matrix <- as.matrix(wasp_best_network_offset_corrplot)
wasp_best_network_offset_corrplot_matrix
diag(wasp_best_network_offset_corrplot_matrix) <- 0
colnames(wasp_best_network_offset_corrplot_matrix) <- c("Pegoscapus sp.","Idarnes flavicolis sp.","Idarnes carme sp. 1","Idarnes carme sp. 2","Heterandrium sp. 1","Heterandrium sp. 2","Ficicola sp.","Physothorax sp.","Sycophila sp.")
rownames(wasp_best_network_offset_corrplot_matrix) <- c("Pegoscapus sp.","Idarnes flavicolis sp.","Idarnes carme sp. 1","Idarnes carme sp. 2","Heterandrium sp. 1","Heterandrium sp. 2","Ficicola sp.","Physothorax sp.","Sycophila sp.")
corrplot(wasp_best_network_offset_corrplot_matrix, diag = FALSE, type = "upper", tl.col = "black", font = 3,is.corr=FALSE)

# Save
png(paste(workdir,"figures/","wasp_residual_covariance.png", sep = ""), width = 6, height = 6, units = "in", res = 600)
corrplot(wasp_best_network_offset_corrplot_matrix, diag = FALSE, type = "upper", tl.col = "black", font = 3, col.lim = c(-0.3,0.3),is.corr=FALSE)
dev.off()

# Both graphs will be combined with Inkscape.


# ---------------------------------------------------------------------------- #
# End of R script 8
# ---------------------------------------------------------------------------- #



# Clear working environment 
rm(list = ls())
dev.off()
