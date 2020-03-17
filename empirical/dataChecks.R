# =========================================================================
# FUNCTIONS FOR DATA CHECKING FOR GROUP 3: GEOGRAPHY OF COMMUNITY SYNCHRONY
# =========================================================================

# ---------------------------------------------------------------------------------------------
# Revised by Max Castorani, University of Virginia, castorani@virginia.edu
# Revised on 2018-05-08
#
# Revised by Nina Lany, Michigan State University, lanynina@msu.edu
# Revised on 2018-11-04 ---------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Load or install necessary libraries
for (package in c('tidyverse', 'ggplot2', 'sp', 'rgdal', 'geosphere', 'maps', 'grDevices', 
                  'RColorBrewer', 'reshape2')) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}

# -------------------------------------------------------------------------
# Make list to handle summaries
dat.summary <- list()

# Add number of unique taxa, years, and spatial replicates to data list:
dat.summary$taxa.no <- length(unique(dat$species))

dat.summary$year.min <- min(dat$year)
dat.summary$year.max <- max(dat$year)
dat.summary$year.no <- length(unique(dat$year))

dat.summary$plot.no <- length(unique(dat$plot))
dat.summary$subplot.no <- length(unique(dat$uniqueID))

dat.summary$abund.units <- paste(dat$unitAbund[1], dat$scaleAbund[1], sep = " within ")

# ---------------------------------------------------------------------------------------------------
# Check balanced sampling of species across space and time by inspecting table, and add to data list
dat.summary$subplot.by.year <- print(list(
  ifelse(length(unique(xtabs(~ year + uniqueID, data = dat))) == 1,
         "OK: Balanced number of taxa recorded across space and time.", 
         "ERROR: Unbalanced numbers of observations across space and time, or taxa list not fully propagated across space and time. Inspect matrix below."),
  xtabs(~ year + uniqueID, data = dat)
))

#---------------------------------------------------------------------------------------------------
# ALPHA DIVERSITY (SPECIES RICHNESS) OVER TIME AND SPACE

# Function to examine temporal patterns in observations of the number of species
no.taxa.fun <- function(ex) {
  
  # Number of unique taxa at each site through time
   no.taxa <- ex %>%
    dplyr::filter(abundance > 0) %>%
    dplyr::select(year, species, uniqueID) %>%
    unique() %>%
    mutate(no.taxa = 1) %>%
    group_by(uniqueID, year) %>%
    summarize(no.taxa = sum(no.taxa))
  
   # Summed number of unique taxa among all sites through time
  total.no.taxa <- ex %>%
    dplyr::filter(abundance > 0) %>%
    dplyr::filter(year > 0) %>%
    dplyr::select(year, species) %>%
    unique() %>%
    mutate(no.taxa = 1) %>%
    group_by(year) %>%
    summarize(no.taxa = sum(no.taxa))
  
  return(list("no.taxa" = no.taxa, "total.no.taxa" = total.no.taxa))
}

#---------------------------------------------------------------------------------------------------
# Create color palette for heatmaps
heat.pal.spectral <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

# ---------------------------------------------------------------------------------------------------
# SITE-SPECIFIC AND TOTAL SPECIES ACCUMULATION CURVES

# Make a function that returns the cumulative number of taxa observed for a given set of community data
cuml.taxa.fun <- function(EX){
  taxa.t.list <- list() # Make empty list
  
  # Loop over each year, creating a list that contains the unique taxa found in each year
  for(t in 1:length(unique(EX$year))){
    tmp.dat <- subset(EX, EX$year == t + (min(EX$year) - 1))
    tmp.dat.pres <- subset(tmp.dat, tmp.dat$abundance > 0) 
    taxa.t.list[[t]] <- unique(tmp.dat.pres$species)
  }
  
  # Make cumulative list of taxa through time
  cuml.taxa <- list() # Empty list
  cuml.taxa[[1]] <- taxa.t.list[[1]] # Add the taxa from the first time step 
  
  # Run for-loop to create list of the cumulative taxa, with duplicates
  for(t in 2:length(unique(EX$year))){ 
    cuml.taxa[[t]] <- c(cuml.taxa[[t - 1]], taxa.t.list[[t]])
  }
  
  # Remove duplicates
  cuml.taxa <- lapply(cuml.taxa, function(x){unique(x)})
  
  # Return the number of total unique taxa through time
  cuml.no.taxa <- data.frame("year" = unique(EX$year))
  cuml.no.taxa$no.taxa <- unlist(lapply(cuml.taxa, function(x){length(x)}))
    
  return(cuml.no.taxa)
  }

# Calculate cumulative taxa (i.e., species accumulation) across all sites pooled together
cuml.taxa.all.sites <- cuml.taxa.fun(EX = dat)

# Examine site-level patterns of species accumulation
# First, sort the comm.dat dataframe to make sure its ordered by site
comm.dat <- dat %>% 
  arrange(uniqueID)

# Then, split the data frame, and apply the cuml.taxa.fun() for each site
X <- split(comm.dat, comm.dat$uniqueID)
out <- lapply(X, cuml.taxa.fun)

# For each site in the list, create column within the data frame that adds the static
for (i in 1:length(out)){
  out[[i]]$uniqueID <- names(out[i])
}
# Make the lists a dataframe
output <- do.call("rbind", out)

# Clean up 
row.names(output) <- NULL

cuml.taxa.by.site <- output

rm(output, out, X)

# ---------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------
# SPECIES ACCUMULATION CURVE FUNCTION - OVER SPACE

# Make a function that returns the cumulative number of taxa observed for a given set of community data
cuml.taxa.space.fun <- function(EX){
  taxa.s.list <- list() # Make empty list
  sites <- unique(EX$uniqueID)
  # Loop over each year, creating a list that contains the unique taxa found in each year
  for(t in 1:length(unique(EX$uniqueID))){
    tmp.dat <- subset(EX, EX$uniqueID == sites[t])
    tmp.dat.pres <- subset(tmp.dat, tmp.dat$abundance > 0) 
    taxa.s.list[[t]] <- unique(tmp.dat.pres$species)
  }

  # Make cumulative list of taxa over space
  cuml.taxa.space <- list() # Empty list
  cuml.taxa.space[[1]] <- taxa.s.list[[1]] # Add the taxa from the first sites 
  
  # Run for-loop to create list of the cumulative taxa, with duplicates
  for(t in 2:length(unique(EX$uniqueID))){ 
    cuml.taxa.space[[t]] <- c(cuml.taxa.space[[t - 1]], taxa.s.list[[t]])
  }
  
  # Remove duplicates
  cuml.taxa.space <- lapply(cuml.taxa.space, function(x){unique(x)})
  
  # Return the number of total unique taxa over space
  cuml.no.taxa.space <- data.frame("site" = unique(EX$uniqueID))
  cuml.no.taxa.space$no.taxa <- unlist(lapply(cuml.taxa.space, function(x){length(x)}))
  
  return(cuml.no.taxa.space)
  }







# ---------------------------------------------------------------------------------------------------
# If coordinates exist, do all of the following bits. If not, do nothing
ifelse(exists("coords") == TRUE, 
  { # ---------------------------------------------------------------------------------------------------
    # Format spatial data
    
    # Ensure coordinate data is numeric
    ifelse(is.numeric(coords$latitude) == "TRUE" & is.numeric(coords$longitude) == "TRUE",
           "OK: Latitude and longitude data coded as decimal-degrees", 
           "ERROR: Latitude and longitude data NOT coded decimal-degrees. Please inspect data cleaning script.")
    
    # Project coordinate data from latitude-longitude to a coordinates object
    long.lat <- coords %>% 
      select(longitude, latitude)
    
    coordinates(long.lat) = c("longitude", "latitude") 
    
    crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") # Projection
    proj4string(long.lat) <- crs.geo
    
    # Summarize coordinate data
    dat.summary$bounding.box <- summary(long.lat)
    
    # ---------------------------------------------------------------------------------------------------
    # Calculate spatial summary statistics
    
    # Distance between spatial units (subplots)
    dist.mat <- (distm(long.lat, fun = distVincentyEllipsoid)/1000) # Distance in km
    rownames(dist.mat) <- coords$uniqueID
    colnames(dist.mat) <- coords$uniqueID
    
    # Add distance matrix to data summary list
    dat.summary$dist.mat <- dist.mat
    
    # Maximum and minimum distances among subplots
    dat.summary$dist.min <- min(dist.mat[dist.mat > 0])
    dat.summary$dist.max <- max(dist.mat)
    
    # ---------------------------------------------------------------------------------------------------
    # Plot spatial distribution of subplots
    
    # Create color palette for heatmaps
    heat.pal.spectral <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
    heat.pal.YOR <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))
    
    # Visualize the distance matrix
    dist.df <- melt(dist.mat) %>%
      rename(subplot.1 = Var1,
             subplot.2 = Var2,
             distance.km = value)
    
    # Add distance data frame to data summary list
    dat.summary$dist.df <- dist.df
    
    # Visualize distance matrix
    dist.mat.plot <- ggplot(data = dist.df, aes(x=subplot.1, y=subplot.2, fill=distance.km)) + 
      geom_tile() +
      scale_fill_gradientn(colours = heat.pal.YOR(100), name = "Distance (km)") +
      xlab("Subplot") +
      ylab("Subplot")
    
    dat.summary$dist.mat.plot <- dist.mat.plot
    
    # ---------------------------------------------------------------------------------------------------
    # Visualize subplots in geographic space without any underlying map
    
    dat.summary$map  <- 
      ggplot(data = coords, aes(x=longitude, y=latitude)) + 
      geom_point(fill = "grey", shape = 21, size = 2) +
      xlab("Longitude") +
      ylab("Latitude") +
      theme_bw() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(color = "black"),
            axis.text.y = element_text(color = "black")) +
      coord_fixed(ratio = 1) +  # Put latitude and longitude on the same scale
      scale_x_continuous(expand = c(0.05, 0.05)) + # Add 5% buffer between points and plot area
      scale_y_continuous(expand = c(0.05, 0.05))
    
  },  print("No coordinates available. Skipping spatially-explicit analyses and visualizations"))

