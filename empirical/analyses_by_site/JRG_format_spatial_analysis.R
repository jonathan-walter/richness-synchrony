# =============================================================================================
# SPATIAL SYNCHRONY FORMATTING AND ANALYSIS FOR JASPER RIDGE DATA
# =============================================================================================
# ---------------------------------------------------------------------------------------------
# by Jon Walter, University of Virginia, jaw3es@virginia.edu
# Created 2019-11-21
# ---------------------------------------------------------------------------------------------
rm(list=ls())
# ---------------------------------------------------------------------------------------------
# Source combined data from EDI
# ---------------------------------------------------------------------------------------------
# Package ID: edi.436.1 Cataloging System:https://pasta.edirepository.org.
# Data set title: Plant survey data from LTER and other US grasslands.
# Data set creator:  Lauren Hallett - University of Oregon
# Data set creator:  Max Castorani -  University of Virginia
# Data set creator:  Nina Lany - US Forest Service
# Data set creator:  Peter Adler - Utah State University 
# Data set creator:  Scott Collins - SEV LTER 
# Data set creator:  Kay Gross - KBS LTER 
# Data set creator:  Richard Hobbs - Jasper Ridge Biological Preserve 
# Data set creator:  Esteban  Muldavin - JRN LTER, SEV LTER 
# Data set creator:  Jennifer Rudgers - JRN LTER 
# Data set creator:  David Tilman - CDR LTER 
# Contact:  Nina Lany - US Forest Service - Nina.Lany@usda.gov
# Stylesheet v2.7 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta-s.lternet.edu/package/data/eml/edi/436/1/d9c5d73d11fc6a16ec41139b63b27751" 
infile1 <- tempfile()
download.file(inUrl1,infile1,method="curl")


grasslands <-read.csv(infile1,header=F 
                      ,skip=1
                      ,sep=","  
                      ,quot='"' 
                      , col.names=c(
                        "year",     
                        "site",     
                        "habitat",     
                        "project",     
                        "plot",     
                        "subplot",     
                        "uniqueID",     
                        "species",     
                        "unitAbund",     
                        "scaleAbund",     
                        "abundance"    ), check.names=TRUE, stringsAsFactors=F)

rm(infile1,inUrl1)

# Source analysis routine script
source(here("CommSpatSynch_v3.R"))
## Data cleaning and formatting #######################################################

dat<-grasslands[grasslands$site=="jrg",]

# rm(list=ls())
# library(dplyr)
# library(vegan)
# 
# source("Data_cleaning/jrg_cleaning.R")
# source("GeogSynch/Scripts/CommSpatSynch_v3.R")
# 
# 
# jrg_crosssite <- jrg_crosssite %>%
#   spread(key = species, value = abundance, fill = 0) %>%
#   gather(key = species, value = abundance, AGHE:VUMI)
# 
# jrg.coords <- jrg.coords %>%
#   arrange(uniqueID)
# coords<-as.data.frame(jrg.coords[,2:3]) %>%
# #coords <- coords %>%
# rename(latitude = coordY, longitude = coordX)
# 
# jrg.soildepth <-jrg.soildepth %>%
#   arrange(uniqueID)
# 
# jrg_habitats <- rep(1, 12)

# -----------------------------------------------------
# Compute frequencies of occurrence
tXl<-length(unique(dat$uniqueID))*length(unique(dat$year))
jrg_spp<-unique(dat$species)
freq<-NULL
for(spp in jrg_spp){
  tmp<-dat[dat$species == spp,]
  freq<-c(freq, sum(tmp$abundance > 0)/tXl)
}

# ------------------------------------------------------------------------
# Remove un-abundant species
thresh<-0.05
sppmin<-5 #minimum number of species in the community for the analysis to proceed
  
sppset<-jrg_spp[freq > thresh]
if(length(sppset)<sppmin){print("stop! too few species to proceed")}

dat<-dat[dat$species %in% sppset,]

# ------------------------------------------------------------------------
# Create data array
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

# -------------------------------------------------------------------------
# Source functions for data checks
source(here("empirical/dataChecks.R"))

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

# -------------------------------------------------------------------------
# Visualize distance matrix of subplots
dat.summary$dist.mat.plot

# Visualize map of subplots
dat.summary$map

# Visualize sampling effort across all subplots
ggplot(data = dat, aes(x = year, y = uniqueID)) +
  geom_point(size = 3) +
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
# Note that the thick line indicates the total number of taxa among all subplots
ggsave(file=here("figures/supplement/jrg_richness_over_time.pdf"), width=7, height=4.7, units="in")

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
ggsave(file=here("figures/supplement/jrg_sp_acc_curve.pdf"), width = 7, height = 4.7, units = "in")

#############################
#visualize spp accumulation curve over space and estimate number of species in 'regional species pool'(the asymptote):
no.taxa.space <- cuml.taxa.space.fun(dat)

##### TEST DIFFERENT SPP ACCUM MODELS ###########
#with vegan, fit species accumulation curve:
sites <- as.numeric(as.factor(no.taxa.space$site))
xtmp <- seq(min(sites), max(sites), len=3*length(sites))
## all sites:
S <-no.taxa.space$no.taxa

##The Arrhenius model:(SSarrhenius) is the expression k*area^z. This is the most classical model that can be found in any textbook of ecology (and also in Dengler 2009). Parameter z is the steepness of the species-area curve, and k is the expected number of species in a unit area.
marr <- nls( S ~ SSarrhenius(sites, k, z))
confint(marr) #z = steepness and k = expected number of species

plot(S ~ sites, xlab = "Time (years)", ylab = "Number of Species", ylim = c(1, max(S)))
lines(xtmp, predict(marr, newdata=data.frame(sites = xtmp)), lwd=2)

## Lomolino: using original names of the parameters (Lomolino 2000):
#The Lomolino model (SSlomolino) is Asym/(1 + slope^log(xmid/area)) (Lomolino 2000, Dengler 2009). Parameter Asym is the asymptotic maximum number of species, slope is the maximum slope of increase of richness, and xmid is the area where half of the maximum richness is achieved.
# mlom <- nls(S ~ SSlomolino(sites, Smax, A50, Hill))
# mlom; confint(mlom)
# lines(xtmp, predict(mlom, newdata=data.frame(sites=xtmp)), lwd=2, col = 4)

## Michaelis Menten:
mmic <- nls(S ~ SSmicmen(sites, slope, Asym))
lines(xtmp, predict(mmic, newdata = data.frame(sites=xtmp)),lwd =2, col = 5)
mmic; confint(mmic)

## compare models (AIC)
allmods <- list(Arrhenius = marr, MicMen= mmic)
sapply(allmods, AIC)

# Plot the cumulative number of taxa observed as plots are added, and add the MicMen line:
pdf(file=here("figures/supplement/jrg_sp_acc_space.pdf"), width = 7, height = 4.7)
plot(as.numeric(no.taxa.space$site), no.taxa.space$no.taxa, pch = 19, xaxt="n", bty="l", xlab = "Cumulative number of sites", ylab = "Cumulative number of taxa", cex=1.5, lwd=3, cex.lab=1.5, ylim = c(0,25))
axis(side=1, at = as.numeric(no.taxa.space$site), labels = seq(1,length(no.taxa.space$site),1))
lines(xtmp, predict(mmic, newdata=data.frame(sites=xtmp)), lwd=2)
dev.off()

coef(mmic)[1]


## Writing L2 and L3 data ##############################################
#make entry (row) for L3 table:
mtdt <- list()
mtdt$dataset <- "jrg"
mtdt$site <- "jrg"
mtdt$initial.year <- dat.summary$year.min
mtdt$study.length <- dat.summary$year.no
mtdt$n.plots <- dat.summary$subplot.no
mtdt$n.taxa <- dat.summary$taxa.no
mtdt$n.taxa.regional <- round(coef(mmic)[1],0)
mtdt$organism <- "terrestrial"
mtdt$taxa <- "grass & forbs"
mtdt$abund.type <- "percent cover"
mtdt$abund.units <- dat.summary$abund.units
mtdt$extent <- 0.0278
mtdt$interplot.dist <-0.0143
mtdt <- data.frame(mtdt)
#write metadata
write.csv(mtdt, file=here("empirical/site_summary_table/jrg_metadata.csv"), row.names=F)
#write L2 data
#write.csv(dat, file = "GeogSynch/L2_data/jrg.csv", row.names=F)



## Running analyses #######################################################
jrg_results<-CommSpatSynch(inarray=data_array, do.Modularity=FALSE, do.Mantel=FALSE)
jrg_results$n.spp <- dim(data_array)[1]

#write to output folder
write_rds(jrg_results, file=here("empirica/analyses_by_site/output/jrg_results.rds"))

