# =========================================================================
# SPATIAL SYNCHRONY FORMATTING AND ANALYSIS FOR SEV DATA
# =========================================================================

# Contributors: Nina Lany, from Max Castorani template
# Created 2018-05-07
rm(list = ls())

# Load or install necessary libraries
for (package in c('tidyverse', 'dplyr', 'ggplot2', 'ecodist', 'abind', 'geosphere', 'rgdal',
                  'maps', 'reshape2', 'codyn', 'igraph', 'vegan', 'devtools')) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}

# Source local version of CommSpatSynch that removed NMDS analysis
source("GeogSynch/Scripts/CommSpatSynch_v3.R")


# -------------------------------------------------------------------------
# Source preliminarily-formatted data
source("Data_cleaning/sev_cleaning.R")

# ---------------------------------------------------------------------------------------------
# Clone specific data sets with a generic name
dat.site <- 'sev'
dat.name <- 'sev_g'
dat.habitat <- 'Black Grama'
dat.domain <- 'terrestrial'

## Site G is black grama-dominated and established in 1999
## Site C is creosote-dominated and established in 1999
## Site B is blue grama dominated and established in 2002

#Select G black-grama
sev<- sev_hierarchies %>%
    filter(habitat == "G") %>%
	mutate_at(vars(c(habitat, project)), as.character) # Ensure that habitat is coded as character

#Select only subplots that were sampled in all years (complete cases):
en <- function(x) {length(unique(x))}
cc <- tapply(sev$year, list(sev$uniqueID, sev$year), en)
cc <- as.data.frame(cc)
cc$uniqueID <- row.names(cc)
cc$Code <- complete.cases(cc)
XX <- subset(cc, Code == "TRUE")
XX$uniqueID
dat <- merge(sev, XX, by = "uniqueID", all.x = F, all.y=T)
dat <- dat[,1:11]
#order by year
dat <- dat[order(dat$uniqueID, dat$year),]

# -----------------------------------------------------
# Compute frequencies of occurrence
tXl<-length(unique(dat$uniqueID))*length(unique(dat$year))
sev_spp<-unique(dat$species)
freq<-NULL
for(spp in sev_spp){
  tmp<-dat[dat$species == spp,]
  freq<-c(freq, sum(tmp$abundance > 0)/tXl)
}

thresh<-0.05
sppmin<-5 #minimum number of species in the community for the analysis to proceed

sppset<-sev_spp[freq > thresh]
if(length(sppset) < sppmin){print("Stop! too few species for analysis")}

dat<-dat[dat$species %in% sppset,]


# -------------------------------------------------------------------------
# Add MCR coordinates manually ##<-these still need to be added if possible################
# coords <- expand.grid(site = unique(dat$site), 
#                       habitat = unique(dat$habitat), 
#                       plot = unique(dat$plot), 
#                       subplot = unique(dat$subplot))
# coords$latitude <- NA
# coords$longitude <- NA
#coords$longitude[coords$plot == "location_1" & coords$habitat == "Fringing" & ]

# -------------------------------------------------------------------------
# Source functions for data checks
source("GeogSynch/Scripts/format_L2_data/geog_synch_data_checks.R")

# -------------------------------------------------------------------------
# What checks are available to view?
names(dat.summary)

# -------------------------------------------------------------------------
# Check number of taxa, nature of time series, nature of spatial units, nature of measurement units
dat.summary$taxa.no

dat.summary$year.min
dat.summary$year.max
dat.summary$year.no

dat.summary$plot.no
dat.summary$subplot.no

dat.summary$abund.units

# -------------------------------------------------------------------------
# Add SEV spatial info to dat.summary because they are in UTM and I can't firgure out how to transform to longlat ##
library(sp)
library(rgdal)
library(maps)

##Note: these coordinates are not geolocated correctly on the Earth, but they should be accurate relative to one another.

#sev_coords <- read.csv(file = "~/Google Drive File Stream/My Drive/LTER_synchrony/Data/Raw_data/Grassland/sev/plot_coordinates/sev_coordinates.csv", stringsAsFactors=F)
sev_coords <-subset(sev_coords, SiteCode=="G")
sev_coords$Plot <-droplevels(sev_coords$Plot)

#there are duplicate entries for most unique quadrats. Perhaps because theexact quadrat sampled shifted at one point over the course of the study? Remove by averaging. 
N <- tapply(sev_coords$Northing, sev_coords$Plot, mean)
E <- tapply(sev_coords$Easting, sev_coords$Plot, mean)
tapply(sev_coords$Northing, sev_coords$Plot, sd) #Looks like the differnce between the two cootds for each quadrat is never greater than 1m.

coords <- cbind(N,E)
sevsp <- SpatialPoints(coords, proj4string = CRS("+init=epsg:6342"))
sevsp <- spTransform(sevsp, CRS("+proj=longlat +ellps=GRS80"))
str(sevsp)

dat.summary$bounding.box <- summary(sevsp)

dist.mat <- (distm(sevsp, fun = distVincentyEllipsoid)/1000) # Distance in km
rownames(dist.mat) <- row.names(N)
colnames(dist.mat) <- row.names(N)
# Add distance matrix to data summary list
dat.summary$dist.mat <- dist.mat

# Maximum and minimum distances among subplots
dat.summary$dist.min <- min(dist.mat[dist.mat > 0])
dat.summary$dist.max <- max(dist.mat)



dat.summary$dist.max
mean(dat.summary$dist.mat)

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
  ggplot(data = as.data.frame(coords), aes(x=E, y=N)) + 
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

# # -------------------------------------------------------------------------
# # Visualize distance matrix of subplots
# dat.summary$dist.mat.plot
# 
# # Visualize map of subplots
# dat.summary$map

# Visualize sampling effort across all subplots
ggplot(data = dat, aes(x = year, y = uniqueID)) +
  geom_raster() +
  theme_bw() +
  xlab("Year") +
  ylab("plot") 

# Visualize subplot-level richness at each time step
no.taxa <- no.taxa.fun(dat) # Result is a list of: (1) no. of taxa at each subplot; (2) no. of taxa at all subplots

# Plot a heatmap of the number of species observed over space (all subplots) and time
ggplot(data = no.taxa$no.taxa, aes(x = year, y = uniqueID, fill = no.taxa)) +
  geom_raster() +
  scale_fill_gradientn(colours = heat.pal.spectral(100), name = "No. of taxa",
                       limits = c(0, max(no.taxa$no.taxa$no.taxa))) + # Make lower bound = 0
  theme_bw() +
  #guides(fill = guide_legend(title = "Number of taxa")) +
  xlab("Year") +
  ylab("Site") +
  theme(aspect.ratio = 1)

# Plot number of taxa through time
ggplot(data=no.taxa$no.taxa, aes(x=year, y=no.taxa)) +
  geom_point(aes(color = uniqueID)) +
  geom_line(aes(color=uniqueID)) +
  geom_point(data=no.taxa$total.no.taxa, aes(x=year, y=no.taxa), color="black", size=3) +
  geom_line(data=no.taxa$total.no.taxa, aes(x=year, y=no.taxa), color="black", size=1) +
  xlab("Year") +
  ylab("Number of taxa observed") +
  guides(color = guide_legend(title = "Site")) +
  ylim(c(0, max(no.taxa$total.no.taxa$no.taxa))) +
  theme_bw() +
  theme(aspect.ratio = 2/3) +
  guides(color = FALSE) +
  scale_x_continuous(breaks = seq(min(cuml.taxa.by.site$year), max(cuml.taxa.by.site$year), by = 2))
ggsave(file = 'GeogSynch/Manuscripts/1_Data_supplemental_methods/sev_g_richness_over_time.pdf', width = 7, height = 4.7, units = "in")# Note that the thick line indicates the total number of taxa among all subplots

# Plot the cumulative number of taxa observed at each subplot, as well as across all subplots together
ggplot(data=cuml.taxa.by.site, aes(x = year, y = no.taxa)) +
  geom_point(aes(color = uniqueID)) +
  geom_line(aes(group = uniqueID, color = uniqueID)) +
  geom_point(data = cuml.taxa.all.sites, aes(x=year, y=no.taxa), size = 3) +
  geom_line(data = cuml.taxa.all.sites, aes(x=year, y=no.taxa), size = 1.5) +
  xlab("Year") +
  ylab("Cumulative number of taxa") +
  guides(color = guide_legend(title = "Site")) +
  ylim(c(0, max(cuml.taxa.all.sites$no.taxa))) +
  theme_bw() +
  theme(aspect.ratio = 2/3) +
  guides(color = FALSE) +
  scale_x_continuous(breaks = seq(min(cuml.taxa.by.site$year), max(cuml.taxa.by.site$year), by = 2))
# Note that the thick line indicates the total number of taxa among all sites
ggsave(file = 'GeogSynch/Manuscripts/1_Data_supplemental_methods/sev_g_sp_acc_curve.pdf', width = 7, height = 4.7, units = "in")

#visualize spp accumulation curve over space and estimate number of species in 'regional species pool'(the asymptote):
no.taxa.space <- cuml.taxa.space.fun(dat)

##### TEST DIFFERENT SPP ACCUM MODELS ###########
#with vegan, fit species accumulation curve:
sites <- as.numeric(rownames(no.taxa.space))
xtmp <- seq(min(sites), max(sites), len=3*length(sites))
## all sites:
S <-no.taxa.space$no.taxa

##The Arrhenius model:(SSarrhenius) is the expression k*area^z. This is the most classical model that can be found in any textbook of ecology (and also in Dengler 2009). Parameter z is the steepness of the species-area curve, and k is the expected number of species in a unit area.
marr <- nls( S ~ SSarrhenius(sites, k, z))
confint(marr) #z = steepness and k = expected number of species

plot(S ~ sites, xlab = "Plots", ylab = "Number of Species", ylim = c(1, max(S)))
lines(xtmp, predict(marr, newdata=data.frame(sites = xtmp)), lwd=2)

## Lomolino: using original names of the parameters (Lomolino 2000):
#The Lomolino model (SSlomolino) is Asym/(1 + slope^log(xmid/area)) (Lomolino 2000, Dengler 2009). Parameter Asym is the asymptotic maximum number of species, slope is the maximum slope of increase of richness, and xmid is the area where half of the maximum richness is achieved.
mlom <- nls(S ~ SSlomolino(sites, Smax, A50, Hill))
# mlom; confint(mlom)
lines(xtmp, predict(mlom, newdata=data.frame(sites=xtmp)), lwd=2, col = 4)

# Michaelis Menten:
mmic <- nls(S ~ SSmicmen(sites, slope, Asym))
lines(xtmp, predict(mmic, newdata = data.frame(sites=xtmp)),lwd =2, col = 5)
mmic; confint(mmic)

## compare models (AIC)
allmods <- list(Arrhenius = marr, Lomolino = mlom, MicMen= mmic)
sapply(allmods, AIC)

# Plot the cumulative number of taxa observed as plots are added, and add the Lomolino line:

pdf(file=paste('GeogSynch/Manuscripts/1_Data_supplemental_methods/', dat.name, '_sp_acc_space.pdf', sep = ""), width = 7, height = 4.7)
plot(sites, no.taxa.space$no.taxa, pch = 19, xaxt="n", bty="l", xlab = "Cumulative number of sites", ylab = "Cumulative number of taxa", cex=1.5, lwd=3, cex.lab=1.5)
axis(side=1, at = sites, labels = seq(1,length(no.taxa.space$site),1))
lines(xtmp, predict(mlom, newdata=data.frame(sites=xtmp)), lwd=2)
dev.off()

coef(mlom)[1]

## Writing L2 and L3 data ##############################################
#make entry (row) for L3 table:
mtdt <- list()
mtdt$dataset <- "sev_g"
mtdt$site <- "sev"
mtdt$initial.year <- dat.summary$year.min
mtdt$study.length <- dat.summary$year.no
mtdt$n.plots <- dat.summary$subplot.no
mtdt$n.taxa <- dat.summary$taxa.no
mtdt$n.taxa.regional <- round(coef(mlom)[1],0)
mtdt$organism <- "terrestrial"
mtdt$taxa <- "black grama"
mtdt$abund.type <- "biomass"
mtdt$abund.units <- dat.summary$abund.units
mtdt$extent <-  dat.summary$dist.max#/1000 #TEMPORARY UNTIL WE HEAR BACK FROM SCOTT
mtdt$interplot.dist <- mean(dat.summary$dist.mat[lower.tri(dat.summary$dist.mat, diag=F)])#/1000 #TEMPORARY UNTIL WE HEAR BACK FROM SCOTT
mtdt <- data.frame(mtdt)
#write metadata
write.csv(mtdt, file = "GeogSynch/Scripts/make_site_table/site_summaries/sev_g_metadata.csv", row.names=F)

#write L2 data
write.csv(dat, file = "GeogSynch/L2_data/sev_g.csv", row.names=F)



# -------------------------------------------------------------------------
# Run analyses

str(dat)
# Convert from long to wide and back to long to be sure that we have fully propagated taxa
dat_wide <- spread(dat, key = species, value = abundance, fill = 0)
dat <- gather(dat_wide, key = species, value = abundance, -year, -uniqueID, -site, -habitat, -project, -plot, -subplot, -unitAbund, -scaleAbund)



sev_data_array <- array(NA, dim=c(length(unique(dat$species)), # No. of species
                                  length(unique(dat$uniqueID)), # No. of unique subplots
                                  length(unique(dat$year))))  # No. of years

for(spp in unique(dat$species)){
  
  mydat <- dat %>%
    filter(species == spp)
    
    mydat <- spread(mydat, key = year, value = abundance, fill=0) 
  
    mydat <- mydat %>%
      select(-(uniqueID:species))
  
  sev_data_array[unique(dat$species) == spp, , ] <- as.matrix(mydat)
}


## Running analyses #######################################################


sev_results<-CommSpatSynch(inarray=sev_data_array, patch.het=NULL)
sev_results$n.spp <- dim(sev_data_array)[1]

sev_g_results<-sev_results

save(sev_g_results, file="GeogSynch/Scripts/format_L2_data/Output/sev_g_results.RData")
