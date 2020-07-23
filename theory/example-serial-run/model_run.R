# This script runs 10 simulations of the theoretical model on your computer 
# without using parallel processing.
# This script calls:
# 1. example10run.csv, which has parameters for 10 runs
#    these 10 parameters are a subset of those included in "example-cluster-run"
# 2. theory_functions.R, which has the theoretical model functions and analyses
#    note that the analyses are a subset of those in CommSpatSynch_v3.R, but paired 
#    down for computational efficiency.

# Set the working directory and load necessary data and libraries
library(here)
setwd(here("theory/example-serial-run"))

# load functions for theory model and necessary pacakges for analysis 
source("theory_functions.R") 
library(codyn) 
library(ecodist) 
library(igraph) 
library(vegan) 
library(wsyn) 

# read in parameters
params <- read.csv("example10run.csv")

# create output objects
num <- dim(params)[1]

AvgPlotRich <- rep(NA,num)
Evenness <- rep(NA,num)
Turnover <- rep(NA,num)
Jaccard <- rep(NA,num)
CVTotBioMass <- rep(NA,num)
LoreauSynch <- rep(NA,num)
VarRatio <- rep(NA,num)
rSppMean <- rep(NA,num)
rTotBiomass <- rep(NA,num)
rRichness <- rep(NA,num)
rEvenness <- rep(NA,num)
rTurnover <- rep(NA,num)
sd.rSppMean <- rep(NA,num)
sd.rTotBiomass <- rep(NA,num)
sd.rRichness <- rep(NA,num)
sd.rEvenness <- rep(NA,num)
sd.rTurnover <- rep(NA,num)

# Run the model for each run
for (i in 1:num) {
  
  patch_het <- params$patch_het[i]
  a <- params$a[i]
  species <- params$species[i]
  patches <- params$patches[i]
  
  r_temp <- params$r_temp[i]
  r <- runif(species, .5-r_temp, .5+r_temp)
  beta_temp <- params$beta_temp[i]
  beta <- matrix(runif(species*species, 0, beta_temp), nrow=species, ncol=species)
  sigma_env_temp <- params$sigma_env_temp[i]
  sigma_env <- rnorm(species, 0, sigma_env_temp)
  disp_temp <- params$disp_temp[i]
  disp_rate <- rep(disp_temp, species)
  
  # set parameters that are constant across runs
  time <- 100
  b <- (1-a^2)^.5
  K <- round(rlnorm(species, meanlog=3, sdlog=1))
  K <- ifelse(K < 2, 2, K)
  sigma_dem <- runif(species, 0, .75)
  
  climate <- array(NA, dim=c(time))
  climate[1] <- 0
  
  # set starting conditions equal to carrying capacities
  N <- array(rep(NA, species*time*patches), dim=c(species, patches, time))
  N_pre <- array(rep(NA, species*time*patches), dim=c(species, patches, time))
  N[,,1] <- K
  
  # set demographic and environmental timeseries
  dem <- array(rnorm((time-1)*patches*species, 0, 1), dim=c(patches, species, (time-1)))
  
  # determine climate for each patch type
  for (t in 1:(time-1)) {
    climate[t+1] <- a*climate[t]+b*rnorm(1,0,1)
  }
  
  # determine how much variation there is in the climate timeseries between patches
  env <- matrix(NA, nrow=patches, ncol=(time-1))
  for (p in 1: patches) {
    env[p,] <- rnorm((time-1), climate, patch_het)
  }
  
  # run abundance model
  for (t in 1:(time-1)) {
    pre_disp <- N[,,t]
    # dispersal; temp update N with the pre-growth population sizes
    N_pre[,,t] <- run.global.dispersal(pre_disp, disp_rate, patches, species)
    # patch model
    post_growth <- run.patch.model(N_pre, r, K, beta, sigma_env, sigma_dem, species, patches, t, env, dem)
    # update N array
    N[,,t+1] <- post_growth
  }
  
  # remove species that went extinct
  temp <- 0
  counter <- 0
  for (s in 1:species) {
    if(all(N[s,,50:100]==0)) {
      temp <- c(temp, s)
      counter <- counter + 1
    }
  }
  
  if(counter > 0) {
    temp <- temp[2:(counter+1)]
    N_new <- N[-temp,,50:100]
  } else {
    N_new <- N[,,50:100]
  }
  
  # run synchrony analyses on abundance model output
  if (sum(N_new > 0)) {
    results <- CommSpatSynch(inarray = N_new)
  } else {
    comm.vars<-rep(NA, 7)
    names(comm.vars)<-c("AvgPlotRich","Evenness","Turnover","Jaccard","CVTotBiomass","LoreauSynch","VarRatio")
    synch.vars<-rep(NA, 10)
    names(synch.vars)<-c("rSppMean","rTotBiomass","rRichness","rEvenness","rTurnover",
                         "sd.rSppMean","sd.rTotBiomass","sd.rRichness","sd.rEvenness","sd.rTurnover")
     results <- list(comm.vars=comm.vars, synch.vars=synch.vars)
  }
  
  # save results
  AvgPlotRich[i] <- results$comm.vars$AvgPlotRich
  Evenness[i] <- results$comm.vars$Evenness
  Turnover[i] <- results$comm.vars$Turnover
  Jaccard[i] <- results$comm.vars$Jaccard
  CVTotBioMass[i] <- results$comm.vars$CVTotBiomass
  LoreauSynch[i] <- results$comm.vars$LoreauSynch
  VarRatio[i] <- results$comm.vars$VarRatio
  rSppMean[i] <- results$synch.vars[1]
  rTotBiomass[i] <- results$synch.vars[2]
  rRichness[i] <- results$synch.vars[3]
  rEvenness[i] <- results$synch.vars[4]
  rTurnover[i] <- results$synch.vars[5]
  sd.rSppMean[i] <- results$synch.vars[6]
  sd.rTotBiomass[i] <- results$synch.vars[7]
  sd.rRichness[i] <- results$synch.vars[8]
  sd.rEvenness[i] <- results$synch.vars[9]
  sd.rTurnover[i] <- results$synch.vars[10]
  
  print(i)
}

# create dataframe of results
final <- data.frame(AvgPlotRich,
                Evenness,
                Turnover,
                Jaccard,
                CVTotBioMass,
                LoreauSynch,
                VarRatio,
                rSppMean,
                rTotBiomass,
                rRichness,
                rEvenness,
                rTurnover,
                sd.rSppMean,
                sd.rTotBiomass,
                sd.rRichness,
                sd.rEvenness,
                sd.rTurnover, row.names=NULL)

write.csv(final, file="example10run_results.csv", row.names=FALSE)
save(final, file = "example10run_results.rdata")
