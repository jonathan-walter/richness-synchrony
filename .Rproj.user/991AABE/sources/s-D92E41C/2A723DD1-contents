# =============================================================================================
# SPATIAL SYNCHRONY FORMATTING AND ANALYSIS FOR MCR DATA
# =============================================================================================
# ---------------------------------------------------------------------------------------------
# Revised by Max Castorani, University of Virginia, castorani@virginia.edu
# Revised on 2018-05-09
# ---------------------------------------------------------------------------------------------
rm(list = ls())
# ---------------------------------------------------------------------------------------------
## Data manipulation packages
# ---------------------------------------------------------------------------------------------
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



# ---------------------------------------------------------------------------------------------
# Load datasets
# ---------------------------------------------------------------------------------------------
# Source preliminarily-formatted data
source("Data_cleaning/mcr_cleaning.R")

# ---------------------------------------------------------------------------------------------
# Clone specific data sets with a generic name
dat.site <- 'mcr'
dat.name <- 'mcr_coral_algae_fringing'
dat.habitat <- 'Fringing'
dat.domain <- 'marine'

dat <- mcr; rm(mcr)

coords <- mcr_coordinates; rm(mcr_coordinates)

# ---------------------------------------------------------------------------------------------
# Subset data based on habitat or project

# First, find the years for which we have all data projects
# min.yr.algae <- min(dat[dat$project == "algae",]$year)
# min.yr.coral <- min(dat[dat$project == "coral",]$year)
# max.yr.algae <- max(dat[dat$project == "algae",]$year)
# max.yr.coral <- max(dat[dat$project == "coral",]$year)
# min.yr <- max(c(min.yr.algae, min.yr.coral))
# max.yr <- min(c(max.yr.algae, max.yr.coral))
# 
min.yr=2006
max.yr=2015
dat <- dat %>%
  dplyr::filter(year >= min.yr & year <= max.yr) %>%
  droplevels(.)
# 
# rm(min.yr.algae, min.yr.coral, max.yr.algae, max.yr.coral, min.yr, max.yr)

# Second, subset by habitat type
dat <- dat %>%
  dplyr::filter(habitat == dat.habitat) %>%
  droplevels(.)

coords <- coords %>%
  dplyr::filter(habitat == dat.habitat) %>%
  droplevels(.)


# ---------------------------------------------------------------------------------------------
# Compute frequencies of occurrence and exclude very rare species
# ---------------------------------------------------------------------------------------------

tXl <- length(unique(dat$uniqueID))*length(unique(dat$year)) # How many combinations of years and spatial locations?
no_spp <- unique(dat$species) # How many taxa?
freq <- NULL

for(spp in no_spp){  # Break up the dataset into separate datasets for each taxon
  tmp <- dat[dat$species == spp,]
  freq <- c(freq, sum(tmp$abundance > 0)/tXl)  # Then, calculate the frequency of occurrence for all taxa relative to the number of sampling time-location combinations
}

thresh <- 0.05 # Minimum proportion of sampling events that a taxon must be present in the data
sppmin <- 5 # Minimum number of species in the community for the analysis to proceed

sppset <- no_spp[freq > thresh]
if(length(sppset) < sppmin){print("Stop! Too few species for analysis")}

dat <- dat[dat$species %in% sppset, ]


# ---------------------------------------------------------------------------------------------
# For later reumanplatz package, create a data array based on taxon abundance at each space-time sampling point
# ---------------------------------------------------------------------------------------------

# First spread the abundance of each taxon over years
dat.spread <- dat %>%
  group_by(species) %>%
  spread(key = year, value = abundance, fill = 0) %>%
  dplyr::select(site, habitat, plot, subplot, uniqueID, species, unitAbund, scaleAbund, everything(), -guild)

# Create empty array for each taxon at each time point and location
data_array <- array(NA, dim=c(length(unique(dat$species)),  # No. of taxa
                              length(unique(dat$uniqueID)), # No. of unique subplots
                              length(unique(dat$year))))    # No. of years

# Fill array with the biomass of each taxon
for(spp in unique(dat.spread$species)){
  data_array[unique(dat.spread$species) == spp, , ] <- dat.spread %>%
    dplyr::filter(species == spp) %>%
    ungroup() %>%
    dplyr::select(-c(site:scaleAbund)) %>%
    as.matrix(.)
}

# Remove temporary data objects
rm(dat.spread, spp, tmp)


# ---------------------------------------------------------------------------------------------
# Check data
# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# Source functions for data checks
source("GeogSynch/Scripts/format_L2_data/geog_synch_data_checks.R")

# ---------------------------------------------------------------------------------------------
# What checks are available to view?
names(dat.summary)

# ---------------------------------------------------------------------------------------------
# Check number of taxa, nature of time series, nature of spatial units, nature of measurement units
dat.summary$taxa.no

dat.summary$year.min
dat.summary$year.max
dat.summary$year.no

dat.summary$plot.no
dat.summary$subplot.no

dat.summary$abund.units


# ---------------------------------------------------------------------------------------------
# Visualize distance matrix of subplots
dat.summary$dist.mat.plot

# Visualize map of subplots
dat.summary$map
max(dat.summary$dist.mat)

# Visualize sampling effort across all subplots
ggplot(data = dat, aes(x = year, y = uniqueID)) +
  geom_point(size = 2) +
  theme_bw() +
  xlab("Year") +
  ylab("plot") 

# Visualize subplot-level richness at each time step
no.taxa <- no.taxa.fun(dat) # Result is a list of: (1) no. of taxa at each subplot; (2) no. of taxa at all subplots

# Join descriptive information to data frame on species richness
no.taxa$no.taxa <-  left_join(no.taxa$no.taxa, 
                              dat[dat$year == min(dat$year), c('site', 'habitat', 'plot', 'subplot', 'uniqueID')], 
                                  by = 'uniqueID') %>%
  dplyr::select(uniqueID, site, habitat, plot, subplot, year, no.taxa)

# Plot a heatmap of the number of species observed over space (all subplots) and time
ggplot(data = no.taxa$no.taxa, aes(x = year, y = uniqueID, fill = no.taxa)) +
  facet_wrap(~ habitat, scales = 'free') +
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
ggsave(file = paste('GeogSynch/Manuscripts/1_Data_supplemental_methods/', dat.name, '_richness_over_time.pdf', sep= ""), width = 7, height = 4.7, units = "in")# Note that the thick line indicates the total number of taxa among all subplots

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
ggsave(file = paste('GeogSynch/Manuscripts/1_Data_supplemental_methods/', dat.name, '_sp_acc_curve.pdf', sep = ""), width = 7, height = 4.7, units = "in")


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
#mlom; confint(mlom)
lines(xtmp, predict(mlom, newdata=data.frame(sites=xtmp)), lwd=2, col = 4)

## Michaelis Menten:
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

# ---------------------------------------------------------------------------------------------
## Writing L2 and L3 data ##############################################
#make entry (row) for L3 table:
mtdt <- list()
mtdt$dataset <- dat.name
mtdt$site <- dat.site
mtdt$initial.year <- dat.summary$year.min
mtdt$study.length <- dat.summary$year.no
mtdt$n.plots <- dat.summary$subplot.no
mtdt$n.taxa <- dat.summary$taxa.no
mtdt$n.taxa.regional <- round(coef(mlom)[1],0)
mtdt$organism <- dat.domain
mtdt$taxa <- "coral & algae"
mtdt$abund.type <- "percent cover"
mtdt$abund.units <- dat.summary$abund.units
mtdt$extent <- dat.summary$dist.max
mtdt$interplot.dist <- mean(dat.summary$dist.mat[lower.tri(dat.summary$dist.mat, diag=F)])
mtdt <- data.frame(mtdt)
#write metadata
write.csv(mtdt, file = paste("GeogSynch/Scripts/make_site_table/site_summaries/", dat.name, "_metadata.csv", sep =""), row.names=F)

#write L2 data
write.csv(dat, file = "GeogSynch/L2_data/mcr_coral_algae_fringing.csv", row.names=F)
# ---------------------------------------------------------------------------------------------
# Run synchrony analyses 
# ---------------------------------------------------------------------------------------------

no_habitats <- c(rep(1, dim(data_array)[2]))

mcr_coral_algae_fringing_results <- CommSpatSynch(inarray = data_array)
mcr_coral_algae_fringing_results$n.spp <- dim(data_array)[1]

save(mcr_coral_algae_fringing_results, file = paste("GeogSynch/Scripts/format_L2_data/Output/", dat.name, "_results.RData", sep =""))
