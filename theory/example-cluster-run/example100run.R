# This script runs 100 simulations of the theoretical model on a cluster using Rmpi.
# This script calls:
# 1. example100run.csv, which has parameters for 100 runs
# 2. theory_functions.R, which has the theoretical model functions and analyses
#    note that the analyses are a subset of those in CommSpatSynch_v3.R, but paired 
#    down for computational efficiency.

# This script is called by example100run.sh, which runs the model using the 
# University of Minnesota Super Computer. It will need to be altered for your cluster.

# Set the number of processors and number of simulations to be run
nProc <- 24      # this is processor*nodes
NumSims <- 100   # match number of rows in params.csv

# Set the working directory and load necessary data and libraries
setwd("~/GeoSynchSimMesabi/")
library(parallel)
library(Rmpi)

params <- read.csv("example100run.csv")

# Simulation function to pass to to the worker nodes. The function performs
#    a single simulation, and will then be repeated enough times to produce
#    the desired number of simulations.
SimFunc <- function(i){
  
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
  sigma_dem <- runif(species, 0, .75)
  
  climate <- array(NA, dim=c(time))
  climate[1] <- 0
  
  # set starting conditions equal to carrying capacities
  N <- array(rep(NA, species*time*patches), dim=c(species, patches, time))
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
    # patch model
    pre_disp <- run.patch.model(N, r, K, beta, sigma_env, sigma_dem, species, patches, t, env, dem)
    # dispersal
    post_disp <- run.global.dispersal(pre_disp, disp_rate, patches, species)
    # update N array
    N[,,t+1] <- post_disp
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
  
  return(results)
}

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

ObjectsToImport <- c("params")
clusterExport(cl, ObjectsToImport)

# Run any commands necessary on the processors before running the simulation. 
temp2 <- clusterEvalQ(cl, setwd("~/GeoSynchSimMesabi/") )
temp2 <- clusterEvalQ(cl, source("theory_functions.R") )
temp2 <- clusterEvalQ(cl, library(codyn) )
temp2 <- clusterEvalQ(cl, library(ecodist) )
temp2 <- clusterEvalQ(cl, library(igraph) )
temp2 <- clusterEvalQ(cl, library(vegan) )
temp2 <- clusterEvalQ(cl, library(wsyn) )

# Run the simulations on the cluster. 
Sims <- clusterApply(cl, x = 1:NumSims, fun = SimFunc)

# or save as an r dataframe
results <- Sims

save(results, file = "example100run_results.rdata")


