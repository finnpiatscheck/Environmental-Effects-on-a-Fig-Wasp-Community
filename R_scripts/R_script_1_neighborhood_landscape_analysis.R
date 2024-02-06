# ============================================================================ #
# R script 1: Neighborhood/landscape spatial analysis
# ============================================================================ #



# Article title: Landscape-Level Analysis of a Fig-Pollinator-Antagonist Community: Spatial and Temporal Variation in a Fig Wasp Community and its Response to Biotic and Abiotic Factors
# By Finn Piatscheck, Justin van Goor, Derek Houston and John Nason.


# This R script presents a function computing tree reproductive landscape for each individual tree within a geographic location. 
# More specifically, the function calculates neighborhood indices reflecting the neighboring landscape of reproductive trees or syconia in a given radius with closer trees adding more weight to the variable, and further trees contributing less.
# We followed Nottebrock et al. 2017:
# https://onlinelibrary.wiley.com/doi/abs/10.1111/oik.03438

# This R script also provides code aiming to understand the mechanism of the function used in this study, visualize the output, and choose parameter values.
# All calculations are already done and presented in the output file flowering_and_syconium_landscape.csv that can be reproduced from raw data in the following R script.
# The code to reproduce Supplementary Figures 2, 3 and 4 is provided here.



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

library(geosphere) # Version 1.5-18
library(leaflet) # Version 2.2.1
library(ggplot2) # Version 3.4.3
library(sf) # Version 1.0-14
library(sp) # Version 2.1-1
library(automap) # Version 1.1-9
library(viridis) # Version 0.6.4



# ---------------------------------------------------------------------------- #
# Import the raw data and format flowering variables
# ---------------------------------------------------------------------------- #



# Import the tree data --------------------

df_tree <- read.csv(paste(workdir,"raw_data/","tree_level_variables.csv", sep=""), h=T, na.strings=c("NA",""), stringsAsFactors = FALSE) 


 # Format the data frame and add reproductive variables --------------------

df_tree$fullID <- as.factor(paste(df_tree$season, df_tree$site, df_tree$tree, sep = "-"))
rownames(df_tree) <- df_tree$fullID
df_tree$tree_volume <- (pi/6)*df_tree$tree_diameter*df_tree$tree_diameter*df_tree$tree_height # Note that this formula is equivalent to the one used for fruit_volume: 4/3*pi*d/2*d/2*h/2
df_tree$crop_size <- df_tree$tree_volume*df_tree$reproduction
# Note that these variables will be computed again in a later scripts.


# Remove outliers --------------------

# Preliminary analyses have identified few trees with very high volumes. These are creating outliers in the crop size and syconium landscape variables.
df_tree[which(df_tree$tree_volume > 6000), c("tree_volume","crop_size")] # These values are problematic for downstream analyses and are replaced with NA.
# Trees are: 113-70 and 70-93A (113-73 and 95-02B high too but kept in later analyses because their crop size values aren't problematic)

# We replace these values with missing data (NAs).
df_tree[rownames(df_tree[which(df_tree$tree_volume > 6000),]),c(which(colnames(df_tree) == "tree_volume"), which(colnames(df_tree) == "crop_size"))] <- NA
# Note that these trees will remain considered later but crop size values will similarly be replaced with NAs.

# Add variables if trees were observed with figs releasing wasps (D phase).
# Note: the D phase is short and some figs are maturing a little faster or slower than the others on a synchronous crop, thus we incorporate the E (wasp recently emerged) and C2 (maturing figs near wasp emergence) in it.
# df_tree$flowering_state <- ifelse(df_tree$flowering == "Yes", 1,0) # This if we want to include all flowering trees, whatever flowering stages there are in. This was abandoned in this study in favor for trees with figs releasing (or near releasing) wasps.
df_tree$flowering_overlap <- ifelse(is.na(df_tree$flowering) == TRUE, NA,
                                    ifelse(df_tree$flowering == "No", 0,
                                           ifelse(df_tree$late_interphase > 0 | df_tree$early_male > 0 | df_tree$late_male > 0,
                                                  1,0))) 

# Assign a crop size value to flowering trees observed releasing wasps.
df_tree$neig_overlap_crop_size <-  df_tree$flowering_overlap*df_tree$crop_size

# Split the data set for the 4 collecting trips with variables of interest.
df_F2012 <- subset(df_tree, season == "F2012", select = c(longitude, latitude, flowering_overlap ,neig_overlap_crop_size))
df_S2013 <- subset(df_tree, season == "S2013", select = c(longitude, latitude, flowering_overlap ,neig_overlap_crop_size))
df_F2013 <- subset(df_tree, season == "F2013", select = c(longitude, latitude, flowering_overlap ,neig_overlap_crop_size))
df_S2014 <- subset(df_tree, season == "S2014", select = c(longitude, latitude, flowering_overlap ,neig_overlap_crop_size))



# ---------------------------------------------------------------------------- #
# Visualization of trees on maps
# ---------------------------------------------------------------------------- #



# Format and subset by collecting season --------------------

df_tree_sites <- df_tree[,2:5]
df_tree_sites$ID <- paste(df_tree$site, df_tree$tree, sep = "-")
df_tree_sites <- df_tree_sites[!duplicated(df_tree_sites$ID),]

df_tree_158 <- df_tree_sites[which(df_tree_sites$site == 158),]
df_tree_172 <- df_tree_sites[which(df_tree_sites$site == 172),]
df_tree_112 <- df_tree_sites[which(df_tree_sites$site == 112),]
df_tree_113 <- df_tree_sites[which(df_tree_sites$site == 113),]
df_tree_95 <- df_tree_sites[which(df_tree_sites$site == 95),]
df_tree_179 <- df_tree_sites[which(df_tree_sites$site == 179),]
df_tree_201 <- df_tree_sites[which(df_tree_sites$site == 201),]
df_tree_96 <- df_tree_sites[which(df_tree_sites$site == 96),]
df_tree_70 <- df_tree_sites[which(df_tree_sites$site == 70),]


# Visualize trees geolocation within each site with leaflet --------------------

# Site 158
leaflet(data = df_tree_158) %>% 
  setView(lng = centroid(df_tree_158[,c(4,3)])[1], lat = centroid(df_tree_158[,c(4,3)])[2], zoom = 16) %>% 
  addTiles() %>% 
  addProviderTiles(providers$Esri.WorldImagery) %>% 
  addCircleMarkers(lng = df_tree_158$longitude, lat = df_tree_158$latitude, color = 'blue', weight = 2, radius = 5, popup = ~ID)

# Site 172
leaflet(data = df_tree_172) %>% 
  setView(lng = centroid(df_tree_172[,c(4,3)])[1], lat = centroid(df_tree_172[,c(4,3)])[2], zoom = 16) %>% 
  addTiles() %>% 
  addProviderTiles(providers$Esri.WorldImagery) %>% 
  addCircleMarkers(lng = df_tree_172$longitude, lat = df_tree_172$latitude, color = 'blue', weight = 2, radius = 5, popup = ~ID)

# Site 112
leaflet(data = df_tree_112) %>% 
  setView(lng = centroid(df_tree_112[,c(4,3)])[1], lat = centroid(df_tree_112[,c(4,3)])[2], zoom = 14) %>% 
  addTiles() %>% 
  addProviderTiles(providers$Esri.WorldImagery) %>% 
  addCircleMarkers(lng = df_tree_112$longitude, lat = df_tree_112$latitude, color = 'blue', weight = 2, radius = 5, popup = ~ID)

# Site 113
leaflet(data = df_tree_113) %>% 
  setView(lng = centroid(df_tree_113[,c(4,3)])[1], lat = centroid(df_tree_113[,c(4,3)])[2], zoom = 12) %>% 
  addTiles() %>% 
  addProviderTiles(providers$Esri.WorldImagery) %>% 
  addCircleMarkers(lng = df_tree_113$longitude, lat = df_tree_113$latitude, color = 'blue', weight = 2, radius = 5, popup = ~ID)
# Note that the first trees mapped (01-18) are far from the main site and were not used for this study. They are however present in the raw data (with missing observations) and visible on this map.

# Site 95
leaflet(data = df_tree_95) %>% 
  setView(lng = centroid(df_tree_95[,c(4,3)])[1], lat = centroid(df_tree_95[,c(4,3)])[2], zoom = 15) %>% 
  addTiles() %>% 
  addProviderTiles(providers$Esri.WorldImagery) %>% 
  addCircleMarkers(lng = df_tree_95$longitude, lat = df_tree_95$latitude, color = 'blue', weight = 2, radius = 5, popup = ~ID)

# Site 201
leaflet(data = df_tree_201) %>% 
  setView(lng = centroid(df_tree_201[,c(4,3)])[1], lat = centroid(df_tree_201[,c(4,3)])[2], zoom = 15) %>% 
  addTiles() %>% 
  addProviderTiles(providers$Esri.WorldImagery) %>% 
  addCircleMarkers(lng = df_tree_201$longitude, lat = df_tree_201$latitude, color = 'blue', weight = 2, radius = 5, popup = ~ID)

# Site 96
leaflet(data = df_tree_96) %>% 
  setView(lng = centroid(df_tree_96[,c(4,3)])[1], lat = centroid(df_tree_96[,c(4,3)])[2], zoom = 13) %>% 
  addTiles() %>% 
  addProviderTiles(providers$Esri.WorldImagery) %>% 
  addCircleMarkers(lng = df_tree_96$longitude, lat = df_tree_96$latitude, color = 'blue', weight = 2, radius = 5, popup = ~ID)

# Site 70
leaflet(data = df_tree_70) %>% 
  setView(lng = centroid(df_tree_70[,c(4,3)])[1], lat = centroid(df_tree_70[,c(4,3)])[2], zoom = 14) %>% 
  addTiles() %>% 
  addProviderTiles(providers$Esri.WorldImagery) %>% 
  addCircleMarkers(lng = df_tree_70$longitude, lat = df_tree_70$latitude, color = 'blue', weight = 2, radius = 5, popup = ~ID)



# ---------------------------------------------------------------------------- #
# Function: "neighborVar" for neighborhood variable computation
# ---------------------------------------------------------------------------- #



# ____________________________________________________________________________ #
# Neighbor Analysis function given a radius (in meters): "neighborVar"
# Input is a data frame with
#  - rownames: individual names
#  - longitude
#  - latitude
#  - 1 or more columns including the variables of interest
# Values for "threshold", "alpha" and "vartotest" need to be set by the user

neighborVar  <- function(df, radius, alpha, vartotest){
  require(geosphere)
  n <- nrow(df)
  d.mat <- matrix(0, n, n, byrow = F)
  dimnames(d.mat)<-list(rownames(df), rownames(df))
  for(i in 1:n){
    for(j in 1:n)
      d.mat[i, j] <- distGeo(c(df[i, 1],
                               df[i, 2]), 
                             c(df[j, 1], 
                               df[j, 2]))
  } 
  diag(d.mat) <- NA
  d.tmp<-apply(d.mat, 1, function(x) which(x<=radius))
  output <- c()
  for(i in 1:length(d.tmp)){
    if (length(d.tmp[[i]]) == 0) {
      sub.calc <- 0
    } else {
      sub.l <- as.vector(df[rownames(df)%in%names(d.tmp[[i]]), vartotest])
      sub.d <- as.vector(d.mat[rownames(d.mat)%in%names(d.tmp[i]), colnames(d.mat)%in%names(d.tmp[[i]])])
      sub.calc <- as.numeric(mapply(function(x,y) x*exp(-(y^2/alpha^2)), x=sub.l, y=sub.d))
    }
    if (all(is.na(sub.calc)) == TRUE) {
      output <- c(output, NA)
    } else {
      output <- c(output,sum(sub.calc, na.rm = TRUE))
    }
  }
  output <- as.data.frame(cbind(rownames(df), output))
  output[,2] <- as.numeric(as.character(output[,2]))
  colnames(output) <- c("ID", paste(vartotest, radius, alpha, sep = "_"))
  return(output)
}

# The distance matrix will be calculated every time the function is run. For this study, the function is run 8 times.
# If a large number of neighborhood variables need to be created, it is wiser to calculate the distance matrix first and modify the function accordingly.

# Values to set in the following function:
# radius <- 200 # Change by the user
# alpha <- 100 # Change by the user
# vartotest <- "reproduction" # Change by the user

# Run the function with 
# neighborVar(df, radius, alpha, vartotest) 
# ____________________________________________________________________________ #



# ---------------------------------------------------------------------------- #
# Test the function and analyze the output
# ---------------------------------------------------------------------------- #



# Here we decompose the function to test 
# We also want to visualize the values of x*exp(-(dij^2/alpha^2) given the distance from the focal tree.
# Thus, we use the code below to guide us in choosing the function's parameter values.
# For better understanding we use the variable "flowering_overlap". The values are either 1, 0, or NA.
# Values of 0 will remain 0 at any distance, and values of 1 with decrease with distance depending on the value of alpha.


# The function neighborVar step by step --------------------

df_test <- df_S2014 # We must use one season only.

# Calculate the distance matrix between all trees in meters (can take several minutes). 
# Note that the function neighborVar runs it each time, which is time-consuming for large data set (and not the most elegant coding). 
# To prevent this, I recommend modifying the function and calculating the distance matrix once before running the rest of the function with different radii and alphas.

n <- nrow(df_test)
d.mat <- matrix(0, n, n, byrow = F) # Distance in meters
dimnames(d.mat)<-list(rownames(df_test), rownames(df_test))
for(i in 1:n){
  for(j in 1:n)
    d.mat[i, j] <- distGeo(c(long1=(df_test[i, 1]),
                             lat1=(df_test[i, 2])), 
                           c(long2=(df_test[j, 1]), 
                             lat2=(df_test[j, 2])))
} 

# Remove zeros of the diagonal in the distance matrix.
diag(d.mat) <- NA

# Visualize a preview of the distance matrix.
d.mat[1:10,1:6] 

# Chose values for parameters alpha and radius, and name the variable of interest. These arguments are provided to the function.
radius <- 50
alpha <- 25
vartotest <- "neig_overlap_crop_size" # Change to "neig_overlap_crop_size". 
# Change these values and rerun the following code to create figures with different parameter values.

# Select the trees within the given radius for each focal tree.
d.tmp <- apply(d.mat, 1, function(x) which(x<=radius))
# Check the returned list: a list slice should be named after the focal tree and contains the names and row numbers of the nearby trees within the given radius.
head(d.tmp)
# Some trees will be isolated depending on the radius given.
d.tmp[769] # Here an example of an isolated tree as site 96 that has no neighbor trees within 500m. As a result, its value here is 0.

# Prepare empty outputs.
output_plot <- data.frame() # This is not in the original function but the data will be used for plots.
output <- c()

# Run the loop that calculates x*exp(-(y^2/alpha^2) sums for each focal trees.
for(i in 1:length(d.tmp)){
  if (length(d.tmp[[i]]) == 0) {
    sub.calc <- 0
    sub.calc_plot <- 0
  } else {
    sub.l <- as.vector(df_test[rownames(df_test)%in%names(d.tmp[[i]]), vartotest])
    sub.d <- as.vector(d.mat[rownames(d.mat)%in%names(d.tmp[i]), colnames(d.mat)%in%names(d.tmp[[i]])])
    sub.calc <- as.numeric(mapply(function(x,y) x*exp(-(y^2/alpha^2)), x=sub.l, y=sub.d))
    sub.calc_plot <- mapply(function(x,y) x*exp(-(y^2/alpha^2)), x=sub.l, y=sub.d) # This is not in the original function.
  }
  sub.f <- as.data.frame(cbind(sub.d, sub.calc_plot)) # This is not in the original function.
  output_plot <- rbind(output_plot, sub.f) # This is also not in the original function, but this output will allow to plot the products of x and the negative exponential of the distance.
  if (all(is.na(sub.calc)) == TRUE) {
    output <- c(output, NA)
  } else {
    output <- c(output,sum(sub.calc, na.rm = TRUE))
  }
}

# Format the output.
output <- as.data.frame(cbind(rownames(df_test), output))
output[,2] <- as.numeric(as.character(output[,2]))
colnames(output) <- c("ID", paste(vartotest, radius, alpha, sep = "_"))

# Visualize the output.
head(output, n=50)


# Plot x*exp(-(dij^2/alpha^2)) values relative to the distance from the focal tree --------------------

# This temporary output should have the distance between tree pairs and the neighbor variable value calculated associated with it.
head(output_plot, n=50)

ggplot(output_plot, aes(sub.d, sub.calc_plot)) +
  geom_point() +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 16)) +
  xlab("Distance (m)") +
  ylab("Neighbor flowering lanscape") +
  ggtitle(paste("Radius =", radius,",", "alpha =",alpha))

ggsave(paste(workdir,"figures/","neighbor_flowering_landscape_", radius,"_",alpha ,".png", sep = ""), plot = last_plot(), device = "png", path = NULL, scale = 1, width = 7, height = 7, units = c("in"), dpi = 600, limitsize = TRUE)
# The values behave as expected, their mean decreases with increased distance.
# The code above was repeated with different variables and parameter values to obtain Supplementary Figures 2 and 3.

# Clear some variables from the working environment.
rm(output_plot, sub.calc_plot, sub.l, sub.d, sub.f, sub.calc, i, j, n)


# Test if the function returns correct values --------------------

# We will check if the neighboring flowering landscape values for tree 70 at site 158 returned by the function are correct.
# For this, first, we will visualize the trees within the circle of a given radius around the focal tree 70.
radius <- 250 # Change to create different figures and analyze the outputs.

df_test_158 <- df_test[which(df_test$latitude > 29),] # Select only rows of site 158.
df_test_158$ID <- rownames(df_test_158)
df_test_158_sf <- st_sfc(st_multipoint(cbind(df_test_158$longitude, df_test_158$latitude)), crs = 4326)
df_test_158_sf <- st_cast(df_test_158_sf, "POINT")
df_test_158_sf <- st_transform(df_test_158_sf, "+proj=utm +zone=12")

circle_158_70 <- st_sfc(st_point(c(df_test_158["S2014-158-70",]$longitude, df_test_158["S2014-158-70",]$latitude)), crs = 4326)
circle_158_70 <- st_transform(circle_158_70, "+proj=utm +zone=12")
circle_158_70 <- st_buffer(circle_158_70, radius) # This creates a ring around the focal tree of the given radius.
circle_158_70_longlat <- st_transform(circle_158_70, "+proj=longlat +datum=WGS84")

# Plot with ggplot2
ggplot(df_test_158_sf) + geom_sf() +
  # geom_sf_text(data = df_test_158_sf, aes(label = df_test_158$ID), size = 3, vjust=-1) +
  geom_sf(data = circle_158_70, colour="black", fill=NA) +
  # coord_sf(xlim = c(st_bbox(circle_158_70)[1], st_bbox(circle_158_70)[3]),
  #          ylim = c(st_bbox(circle_158_70)[2], st_bbox(circle_158_70)[4]), expand = TRUE) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 16)) +
  xlab("Latitude") + ylab("Longitude") + ggtitle(paste("Radius =", radius,"m"))

ggsave(paste(workdir,"figures/","site158_tree70_radius_",radius,"m.png", sep = ""), plot = last_plot(), device = "png", path = NULL, scale = 1, width = 7, height = 6, units = c("in"), dpi = 600, limitsize = TRUE)

# The code for the plot above was also run with radii 10m, 50m; 100m, 250m, 500m and 1000m.
# These plots were combined in Supplementary Figure 4.


# Plot with leaflet (the map is interactive, the names of the trees appear if clicked on them).
leaflet(data = df_tree_158) %>% 
  setView(lng = centroid(df_tree_158[,c(4,3)])[1], lat = centroid(df_tree_158[,c(4,3)])[2], zoom = 16) %>% 
  addTiles() %>% 
  addProviderTiles(providers$Esri.WorldImagery) %>% 
  addCircleMarkers(lng = df_tree_158$longitude, lat = df_tree_158$latitude, color = 'blue', weight = 2, radius = 5, popup = ~ID) %>%
  addPolygons(data= circle_158_70_longlat, color = "red", fill=NA)

# The following confirms what we see:
rownames(df_test_158[st_contains(circle_158_70, df_test_158_sf, sparse = FALSE),]) # This function returns the trees within the circle of the indicated radius.
# Check if this is the same in the function temporary output.
names(d.tmp$`S2014-158-70`) # Trees within the neighborhood of 158-70 (158-70 is absent here because its distance 0 to itself was replaced with NA).
setdiff(rownames(df_test_158[st_contains(circle_158_70, df_test_158_sf, sparse = FALSE),]),
        names(d.tmp$`S2014-158-70`))
# Same (158-70 is also present within the circle, as expected).

# The following code tests if the output from the function is similar to the expected result.
# WARNING: the following is for radius = 250 and alpha = 125.
df_test_158[rownames(df_test_158[st_contains(circle_158_70, df_test_158_sf, sparse = FALSE),]),][,c(1,2,4)]
# 3 trees with a reproduction value are within the radius. Others were not flowering (0) or not visited (NA).
# These trees are 158-28, 158-30 and 158-67

# Manually calculate the value that should be associated with tree 158-70 for radius = 250 and alpha = 125. 
x_value1 <- df_test_158["S2014-158-30",vartotest]
x_value2 <- df_test_158["S2014-158-28",vartotest]
x_value3 <- df_test_158["S2014-158-67",vartotest]
# We calculate the distance between these trees to the focal tree with distGeo().
distance_value1 <- distGeo(c(df_test_158["S2014-158-70",]$longitude, df_test_158["S2014-158-70",]$latitude), 
                          c(df_test_158["S2014-158-30",]$longitude, df_test_158["S2014-158-30",]$latitude))
distance_value2 <- distGeo(c(df_test_158["S2014-158-70",]$longitude, df_test_158["S2014-158-70",]$latitude), 
                           c(df_test_158["S2014-158-28",]$longitude, df_test_158["S2014-158-28",]$latitude))
distance_value3 <- distGeo(c(df_test_158["S2014-158-70",]$longitude, df_test_158["S2014-158-70",]$latitude), 
                           c(df_test_158["S2014-158-67",]$longitude, df_test_158["S2014-158-67",]$latitude))
# Check if it is the same as the distance calculated above.
d.mat["S2014-158-70", "S2014-158-30"];distance_value1 # Good
d.mat["S2014-158-70", "S2014-158-28"];distance_value2 # Good
d.mat["S2014-158-70", "S2014-158-67"];distance_value3 # Good
# Now, we calculate ci:
x_value1*exp(-(distance_value1^2/alpha^2)) + x_value2*exp(-(distance_value2^2/alpha^2)) + x_value3*exp(-(distance_value3^2/alpha^2))
# Same as output?
output[98,]
# The values are identical.

# Let's look at the value of the isolated tree 39 at site 96:
# Plot
df_test_96 <- df_test[which(df_test$latitude > 24 | df_test$latitude < 25),] # Select only rows of site 96.
df_test_96$ID <- rownames(df_test_96)
df_test_96_sf <- st_sfc(st_multipoint(cbind(df_test_96$longitude, df_test_96$latitude)), crs = 4326)
df_test_96_sf <- st_cast(df_test_96_sf, "POINT")
df_test_96_sf <- st_transform(df_test_96_sf, "+proj=utm +zone=12")

circle_96_39 <- st_sfc(st_point(c(df_test_96["S2014-96-39",]$longitude, df_test_96["S2014-96-39",]$latitude)), crs = 4326)
circle_96_39 <- st_transform(circle_96_39, "+proj=utm +zone=12")
circle_96_39 <- st_buffer(circle_96_39, radius) # This creates a ring around the focal tree of the given radius.
circle_96_39_longlat <- st_transform(circle_96_39, "+proj=longlat +datum=WGS84")

leaflet(data = df_tree_96) %>% 
  setView(lng = centroid(df_tree_96[,c(4,3)])[1], lat = centroid(df_tree_96[,c(4,3)])[2], zoom = 15) %>% 
  addTiles() %>% 
  addProviderTiles(providers$Esri.WorldImagery) %>% 
  addCircleMarkers(lng = df_tree_96$longitude, lat = df_tree_96$latitude, color = 'blue', weight = 2, radius = 5, popup = ~ID) %>%
  addPolygons(data= circle_96_39_longlat, color = "red", fill=NA)
# This tree is down a river, isolated from the rest of the population.
# Swap south-west on the map to spot the tree surrounded by a red circle.

# Value for tree 96-39 should be 0.
output[769,]

# The above has been tested for different variables and parameter values and always returned the same result as the output.

# Clear working environment
rm(df_test_158, df_test_158_sf, circle_158_70, circle_158_70_longlat, x_value1, x_value2, x_value3, distance_value1, distance_value2, distance_value3,  df_test_96, df_test_96_sf, circle_96_39, circle_96_39_longlat)


# Visualize the ci values on the landscape --------------------

# The following code will plot syconium landscape values for individual trees at site 158 during the spring 2014.
syconium_landscape_158_sp <- as.data.frame(cbind(df_F2012$longitude, df_F2012$latitude, output[,2]))
syconium_landscape_158_sp <- na.omit(syconium_landscape_158_sp[which(syconium_landscape_158_sp[,2] > 29),])
colnames(syconium_landscape_158_sp) <- c("x","y","var")
coordinates(syconium_landscape_158_sp) <- ~x+y
proj4string(syconium_landscape_158_sp) <- CRS("+init=epsg:4326")
syconium_landscape_158_sp <- spTransform(syconium_landscape_158_sp, CRS("+proj=utm +zone=12"))

spplot(syconium_landscape_158_sp,"var",colorkey=TRUE) # spplot of trees and their syconium landscape values

syconium_landscape_158_grid <- makegrid(syconium_landscape_158_sp, n=5000)
syconium_landscape_158_grid <- spTransform(SpatialPoints(syconium_landscape_158_grid, proj4string = CRS(proj4string(syconium_landscape_158_sp))), CRS("+proj=utm +zone=12"))
# plot(syconium_landscape_158_grid); points(syconium_landscape_158_sp,pch=20) # Plot the grid and tree locations

# Next, we interpolate the ci values on the landscape to visualize its variation in space.
# As our goal is not to do a perfect extrapolation but just a simple visualization, we utilize the autoKrige() function for simplification.
syconium_landscape_158_krig <- autoKrige(var~1, syconium_landscape_158_sp, syconium_landscape_158_grid)
plot(syconium_landscape_158_krig)
automapPlot(syconium_landscape_158_krig$krige_output, 
            zcol = "var1.pred", 
            main = "Spring 2014, site 158", 
            sp.layout = list("sp.points", syconium_landscape_158_sp, pch = 3, col = "red", alpha = .5))

# Let's make a nicer plot with ggplot2.
syconium_landscape_158_krig_df <- as.data.frame(syconium_landscape_158_krig$krige_output)
syconium_landscape_158_krig_df$var1.pred[syconium_landscape_158_krig_df$var1.pred<0] <- 0 # We can't have negative values, so we replace them with 0.
syconium_landscape_158_sf <- st_as_sf(syconium_landscape_158_sp)

ggplot(data = syconium_landscape_158_krig_df,aes(x=coords.x1,y=coords.x2)) + 
  geom_tile(aes(fill=var1.pred), alpha = 0.8) +
  geom_point(data = as.data.frame(syconium_landscape_158_sp), aes(x=coords.x1,y=coords.x2), color="blue") +
  theme(legend.position = c(0.15,0.18)) +
  scale_fill_gradientn(colours = rev(inferno(10))) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_equal()

# Plot with longitude and latitude coordinates, sf object projected.
ggplot(data = syconium_landscape_158_sf) + 
  geom_sf(aes(fill=var), shape = 21, size = 3, stroke = 1) +
  geom_tile(data =as.data.frame(syconium_landscape_158_krig_df), aes(x = coords.x1, y = coords.x2, fill=var1.pred), alpha = 0.5) +
  theme(legend.position = "right", text=element_text(size = 12),
        legend.background = element_rect(colour='white'),
        legend.title = element_text(colour="black", size = 10),
        legend.text  = element_text(colour="black", size = 10),
        legend.key.width = unit(0.3,"cm"),
        legend.key.height = unit(2, "cm"),
        axis.text.x=element_text(colour="black", size = 12, angle = 45, vjust = 0.5),  
        axis.text.y=element_text(colour="black", size = 12),
        panel.border = element_rect(fill = NA, size = 0.5),
        panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "white")) +
  scale_fill_gradientn(name=NULL, colours = rev(viridis(10))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(NULL) +
  ylab(NULL) 

# The above can be repeated for different sites and field trips. The code for reproduction figures 10 and 11 is in script "figure_reproduction.R".

# Clear working environment. 
rm(df_test, d.mat, d.tmp, output, df_tree_sites, df_test_158, df_test_158_sf, df_test_96, df_test_96_sf, circle_158_70, circle_158_70_longlat, circle_96_39, circle_96_39_longlat,
   syconium_landscape_158_sp, syconium_landscape_158_grid, syconium_landscape_158_krig, syconium_landscape_158_krig_df, syconium_landscape_158_sf)



# ---------------------------------------------------------------------------- #
# Compute neighborhood variables
# ---------------------------------------------------------------------------- #



# We choose parameter values radiues = 250m and alpha = 125, as the example above. 
# The justification of these values are in Appendix 1.

F2012_flowering_landscape_250_125 <- neighborVar(df_F2012, 250, 125, "flowering_overlap")
S2013_flowering_landscape_250_125 <- neighborVar(df_S2013, 250, 125, "flowering_overlap")
F2013_flowering_landscape_250_125 <- neighborVar(df_F2013, 250, 125, "flowering_overlap")
S2014_flowering_landscape_250_125 <- neighborVar(df_S2014, 250, 125, "flowering_overlap")

F2012_syconium_landscape_250_125 <- neighborVar(df_F2012, 250, 125, "neig_overlap_crop_size")
S2013_syconium_landscape_250_125 <- neighborVar(df_S2013, 250, 125, "neig_overlap_crop_size")
F2013_syconium_landscape_250_125 <- neighborVar(df_F2013, 250, 125, "neig_overlap_crop_size")
S2014_syconium_landscape_250_125 <- neighborVar(df_S2014, 250, 125, "neig_overlap_crop_size")

# Format into 4 variables.
# Each variable value will be associated with a unique ID corresponding to one tree within one site, observed during one season.

flowering_and_syconium_landscape_F2012 <- cbind(F2012_flowering_landscape_250_125,F2012_syconium_landscape_250_125[,2])
flowering_and_syconium_landscape_S2013 <- cbind(S2013_flowering_landscape_250_125,S2013_syconium_landscape_250_125[,2])
flowering_and_syconium_landscape_F2013 <- cbind(F2013_flowering_landscape_250_125,F2013_syconium_landscape_250_125[,2])
flowering_and_syconium_landscape_S2014 <- cbind(S2014_flowering_landscape_250_125,S2014_syconium_landscape_250_125[,2])

colnames(flowering_and_syconium_landscape_F2012) <- gsub("[, 2]", "", colnames(flowering_and_syconium_landscape_F2012), fixed = TRUE)
colnames(flowering_and_syconium_landscape_F2012) <- gsub("F|S|2012_|2013_|2014_", "", colnames(flowering_and_syconium_landscape_F2012))
colnames(flowering_and_syconium_landscape_S2013) <- gsub("[, 2]", "", colnames(flowering_and_syconium_landscape_S2013), fixed = TRUE)
colnames(flowering_and_syconium_landscape_S2013) <- gsub("F|S|2012_|2013_|2014_", "", colnames(flowering_and_syconium_landscape_S2013))
colnames(flowering_and_syconium_landscape_F2013) <- gsub("[, 2]", "", colnames(flowering_and_syconium_landscape_F2013), fixed = TRUE)
colnames(flowering_and_syconium_landscape_F2013) <- gsub("F|S|2012_|2013_|2014_", "", colnames(flowering_and_syconium_landscape_F2013))
colnames(flowering_and_syconium_landscape_S2014) <- gsub("[, 2]", "", colnames(flowering_and_syconium_landscape_S2014), fixed = TRUE)
colnames(flowering_and_syconium_landscape_S2014) <- gsub("F|S|2012_|2013_|2014_", "", colnames(flowering_and_syconium_landscape_S2014))

flowering_and_syconium_landscape <- rbind(flowering_and_syconium_landscape_F2012, flowering_and_syconium_landscape_F2013, flowering_and_syconium_landscape_S2013, flowering_and_syconium_landscape_S2014)
colnames(flowering_and_syconium_landscape) <- c("ID","flowering_landscape_250_125","syconium_landscape_250_125")

# Save data.
write.csv(flowering_and_syconium_landscape, paste(workdir,"generated_data/","flowering_and_syconium_landscape.csv", sep = ""), row.names = FALSE)



# ---------------------------------------------------------------------------- #
# End of R script 1
# ---------------------------------------------------------------------------- #



# Clear working environment.
rm(list = ls())
dev.off()
