# =============================================================================================
# SPATIAL SYNCHRONY FORMATTING AND ANALYSIS FOR MCR DATA
# =============================================================================================
# ---------------------------------------------------------------------------------------------
# Revised by Max Castorani, University of Virginia, castorani@virginia.edu
# Revised on 2018-05-09
# ---------------------------------------------------------------------------------------------
rm(list=ls())
# ---------------------------------------------------------------------------------------------
# Source combined data from EDI
# ---------------------------------------------------------------------------------------------
# Package ID: edi.437.1 Cataloging System:https://pasta.edirepository.org.
# Data set title: Benthic algae and sessile invertebrate survey data from LTER and other reef ecosystems.
# Data set creator:  Lauren Hallett - University of Oregon
# Data set creator:  Max Castorani - University of Virginia
# Data set creator:  Nina Lany - US Forest Service
# Data set creator:  Daniel Reed - SBC LTER 
# Data set creator:  Peter Edmunds - USVI, MCR LTER 
# Data set creator:  Robert Carpenter - MCR LTER 
# Contact:  Nina Lany - US Forest Service - Nina.Lany@usda.gov
# Stylesheet v2.7 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl2  <- "https://pasta-s.lternet.edu/package/data/eml/edi/437/1/d9c5d73d11fc6a16ec41139b63b27751" 
infile2 <- tempfile()
download.file(inUrl2,infile2,method="curl")


marine <-read.csv(infile2,header=F 
                  ,skip=1
                  ,sep=","  
                  ,quot='"' 
                  , col.names=c(
                    "year",     
                    "site",     
                    "habitat",     
                    "plot",     
                    "subplot",     
                    "uniqueID",     
                    "guild",     
                    "species",     
                    "abundance",     
                    "unitAbund",     
                    "scaleAbund"    ), check.names=TRUE, stringsAsFactors=FALSE)

rm(infile2,inUrl2)

# Source analysis routine script
source(here("CommSpatSynch_v3.R"))

# ---------------------------------------------------------------------------------------------
# select data
# ---------------------------------------------------------------------------------------------
# Source preliminarily-formatted data
dat<-marine[marine$site=="mcr",]
dat<-dat[dat$habitat=="Fringing",]

min.yr=2006
max.yr=2015
dat <- dat %>%
  dplyr::filter(year >= min.yr & year <= max.yr) %>%
  droplevels(.)

# ---------------------------------------------------------------------------------------------
# Clone specific data sets with a generic name
dat.site <- 'mcr'
dat.name <- 'mcr_coral_algae_fringing'
dat.habitat <- 'Fringing'
dat.domain <- 'marine'

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
#  create a data array based on taxon abundance at each space-time sampling point
# ---------------------------------------------------------------------------------------------

data_array<-array(0, dim=c(length(unique(dat$species)),length(unique(dat$uniqueID)),
                           length(unique(dat$year))))

for(spp in 1:length(unique(dat$species))){
  for(loc in 1:length(unique(dat$uniqueID))){
    for(yr in 1:length(unique(dat$year))){
      if(any(dat$species==unique(dat$species)[spp] & dat$uniqueID==unique(dat$uniqueID)[loc]
             & dat$year==unique(dat$year)[yr])){
        data_array[spp,loc,yr]<-dat$abundance[dat$species==unique(dat$species)[spp] & dat$uniqueID==unique(dat$uniqueID)[loc]
                                              & dat$year==unique(dat$year)[yr]]
      }
    }
  }
}


# ---------------------------------------------------------------------------------------------
# Check data
# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# Source functions for data checks
source(here("empirical/dataChecks.R"))

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
# Note that the thick line indicates the total number of taxa among all subplots
ggsave(file = here("figures/supplement",paste(dat.name, '_richness_over_time.pdf', sep= "")), width = 7, height = 4.7, units = "in")

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
ggsave(file = here("figures/supplement",paste(dat.name, '_sp_acc_curve.pdf', sep = "")), width = 7, height = 4.7, units = "in")

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

pdf(file=here("figures/supplement",paste(dat.name, '_sp_acc_space.pdf', sep = "")), width = 7, height = 4.7)
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
mtdt$extent <- 15.674
mtdt$interplot.dist <- 9.659
mtdt <- data.frame(mtdt)
#write metadata
write.csv(mtdt, file = here("empirical/site_summary_table",paste(dat.name, "_metadata.csv", sep ="")), row.names=F)

#write L2 data
#write.csv(dat, file = "GeogSynch/L2_data/mcr_coral_algae_fringing.csv", row.names=F)
# ---------------------------------------------------------------------------------------------
# Run synchrony analyses 
# ---------------------------------------------------------------------------------------------

mcr_coral_algae_fringing_results <- CommSpatSynch(inarray = data_array, do.Modularity=F, do.Mantel=F)
mcr_coral_algae_fringing_results$n.spp <- dim(data_array)[1]

saveRDS(mcr_coral_algae_fringing_results, file = here("empirical/analyses_by_site/output",paste(dat.name, "_results.rds", sep ="")))
