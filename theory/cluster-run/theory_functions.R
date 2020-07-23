# Lauren Shoemaker
# Last updated June 26, 2019

# ----------------------------------------------------------------------------------------------------------------------
# model functions
# N = population size
# r = density independent growth rate
# K = carrying capacity
# beta = competition coefficient
# sigma_env = effect of environmental stochasticity
# sigma_dem = effect of demographic stochasticity
# species = number of species
# patches = number of patches
# t = timepoint
# env = environmetal stochasticity
# dem = demographic stochasticity
run.patch.model <- function(N, r, K, beta, sigma_env, sigma_dem, species, patches, t, env, dem) {
  N_temp <- matrix(rep(NA, patches*species), nrow=species, ncol=patches)
  for (p in 1:patches) {
    for (s in 1:species) {
      if(N[s,p,t] > 0) {
        N_temp[s,p] <- N[s,p,t]*exp(r[s]*(1-N[s,p,t]/K[s]-sum(beta[s,]*N[,p,t]/K) + beta[s,s]*N[s,p,t]/K[s]) + 
                                      sigma_env[s]*env[p,t] + sigma_dem[s]*dem[p,s,t]/sqrt(N[s,p,t]))
        if(N_temp[s,p] <0) {
          N_temp[s,p] <- 0
        }
      }else {
        N_temp[s,p] <- 0
      }
    }
  }
  
  N_temp <- round(N_temp)
  
  return(N_temp)
}

# pre_disp = pre-dispersal population size in each patch
# disp_rate = dispersal rate
# patches = number of patches
# species = number of species
run.global.dispersal <- function(pre_disp, disp_rate, patches, species){
  
  # keep track of who stays and who goes
  stay <- matrix(rep(NA, species*patches), nrow=species, ncol=patches)
  where <- matrix(rep(0, species*patches), nrow=species, ncol=patches)
  
  for (p in 1:patches) {
    # determine who disperses
    to_disperse <- rbinom(species, pre_disp[,p], disp_rate)
    stay[,p] <- pre_disp[,p]-to_disperse
    for (s in 1:species) {
      if(to_disperse[s] > 0) {
        vals <- seq(1:patches)
        vals <- vals[-p]
        temp_where <- sample(vals, to_disperse[s], replace=TRUE)
        for (x in 1:patches) {
          where[s,x] <- sum(temp_where==x) + where[s,x]
        }
      }
    }
  }
  total <- stay + where
  return(total)
}

# ----------------------------------------------------------------------------------------------------------------------
# analysis function

CommSpatSynch<-function(inarray, clnlev=2, corrmeth="spearman", patch.het=NULL, do.SynchSig=F, 
                        do.Modularity=F, do.Mantel=F){
  
  ## Data Cleaning Function ########################################
  clean<-function(indata, clnlev){
    
    if(!clnlev %in% c(1,2)){stop("clnlev must be 1 or 2")}
    
    out<-matrix(NA, nrow(indata), ncol(indata))
    
    y<-indata#[rowSums(indata)!=0,]
    x<-1:ncol(indata)
    
    for(mm in 1:nrow(y))
      
      if(clnlev>=1){
        for(mm in 1:nrow(y)){
          if(sum(y[mm,], na.rm=T)==0){
            out[mm,]<-y[mm,]
          }
          else{
            out[mm,]<-residuals(lm(y[mm,]~x, na.action=na.exclude))
          }
        }
      }
    if(clnlev==2){
      for(mm in 1:nrow(out)){
        if(sum(out[mm,], na.rm=T)==0){next}
        out[mm,]<-out[mm,]/sd(out[mm,])
      }
    }
    return(out)
  }
  
  if(length(unique(patch.het))==1){patch.het<-NULL}
  
  # Get dims of data array
  nspp<-dim(inarray)[1]
  nlocs<-dim(inarray)[2]
  tmax<-dim(inarray)[3]
  
  yrs=1:tmax
  
  ## Initialize outputs ################################################################
  comm.vars<-rep(NA, 7)
  names(comm.vars)<-c("AvgPlotRich","Evenness","Turnover","Jaccard","CVTotBiomass","LoreauSynch","VarRatio")
  synch.vars<-rep(NA, 10)
  names(synch.vars)<-c("rSppMean","rTotBiomass","rRichness","rEvenness","rTurnover",
                       "sd.rSppMean","sd.rTotBiomass","sd.rRichness","sd.rEvenness","sd.rTurnover")
  synch.sig<-rep(NA, 5)
  names(synch.sig)<-c("p.rSppMean","p.rTotBiomass","p.rRichness","p.rEvenness","p.rTurnover")
  module.info<-rep(NA, 10)
  names(module.info)<-c("nMod.SppMean","nMod.TotBiomass","nMod.Richness","nMod.Evenness","nMod.Turnover",
                        "Mod.SppMean","Mod.TotBiomass","Mod.Richness","Mod.Evenness","Mod.Turnover")
  mantel.results<-matrix(NA, nrow=15, ncol=2)
  colnames(mantel.results)<-c("coeff","pval")
  rownames(mantel.results)<-c("rTotBiomassXrSppMean","rTotBiomassXrRichness","rTotBiomassXrEvenness","rTotBiomassXrTurnover","rTotBiomassXJaccard",
                              "rSppMeanXrRichness","rSppMeanXrEvenness","rSppMeanXrTurnover", "rSppMeanXJaccard",
                              "rRichnessXrEvenness","rRichnessXrTurnover","rRichnessXJaccard",
                              "rEvennessXrTurnover","rEvennessXJaccard",
                              "rTurnoverXJaccard")
  module.correspondence<-rep(NA, 5)
  names(module.correspondence)<-c("rSppMean","rTotBiomass","rRichness","rEvenness","rTurnover")
  
  ## Start Chugging! ########################################################################
  
  ## Matrices of Community Totals --  use this to calculate community-wide richness and evenness
  CommSums<-apply(inarray, 1, sum, na.rm=T)
  CommAvgs<-apply(inarray, 1, mean, na.rm=T)
  sppXtime<-apply(inarray, c(1,3), sum,  na.rm=T)
  sppXplot<-apply(inarray, c(2,1), mean,  na.rm=T)
  
  TotBiomassXtime<-apply(inarray,3,sum,na.rm=T)
  comm.vars[5]<-sd(TotBiomassXtime)/mean(TotBiomassXtime)
  
  ## Loreau synchrony
  indices<-which(sppXtime==sppXtime, arr.ind=T)
  data_long<-data.frame(time = indices[,2], species = indices[,1], abundance=c(sppXtime))
  comm.vars[6]<-synchrony(data_long, time.var="time", species.var="species", abundance.var="abundance")
  
  ## Variance Ratio
  comm.vars[7]<-variance_ratio(data_long, time.var="time", species.var="species", abundance.var="abundance", bootnumber=1)[4]
  
  ## Beta diversity (Jaccard similarity)
  jaccsim<-1-full(distance(sppXplot>0, method="jaccard"))
  comm.vars[4]<-mean(jaccsim[lower.tri(jaccsim)], na.rm=T)
  
  ## Species average synchrony
  synchrony_array<-array(NA,dim=c(nspp,nlocs,nlocs))
  for(spp in 1:nspp){
    synchrony_array[spp,,]<-suppressWarnings(cor(t(clean(inarray[spp,,],clnlev)),method=corrmeth,use="pairwise.complete.obs"))
  }
  
  syn.spp.mean<-apply(synchrony_array, c(2,3), mean, na.rm=T)
  syn.spp.mean[is.nan(syn.spp.mean)]<-0 ##Some sites never share a species--set correlation here to 0
  synch.vars[1]<-mean(syn.spp.mean[lower.tri(syn.spp.mean)])
  synch.vars[6]<-sd(syn.spp.mean[lower.tri(syn.spp.mean)])
  
  ## Total biomass
  total.biomass<-apply(inarray, c(2,3), sum, na.rm=T)
  syn.total.biomass<-suppressWarnings(cor(t(clean(total.biomass,clnlev)),method=corrmeth))
  syn.total.biomass[is.na(syn.total.biomass)]<-0
  synch.vars[2]<-mean(syn.total.biomass[lower.tri(syn.total.biomass)])
  synch.vars[7]<-sd(syn.total.biomass[lower.tri(syn.total.biomass)])
  
  ## Richness
  richness<-apply(inarray>0, c(2,3), sum, na.rm=T)
  syn.richness<-suppressWarnings(cor(t(clean(richness,clnlev)),method=corrmeth))
  syn.richness[is.na(syn.richness)]<-0
  synch.vars[3]<-mean(syn.richness[lower.tri(syn.richness)])
  synch.vars[8]<-sd(syn.richness[lower.tri(syn.richness)])
  comm.vars[1]<-mean(richness, na.rm=T)
  
  ## Evenness
  fn.evenness<-function(x,totP=NULL,S=NULL){
    x<-x[!is.na(x) & x > 0]
    if(is.null(totP)){totP = sum(x)}
    if(is.null(S)){S=length(x)}
    Pi = x/totP
    H = -sum(Pi * log(Pi))
    return(H/log(S))
  }
  
  evenness<-apply(inarray, c(2,3), fn.evenness, S=nspp)
  comm.vars[2]<-fn.evenness(CommAvgs)
  syn.evenness<-suppressWarnings(cor(t(clean(evenness,clnlev)),method=corrmeth))
  syn.evenness[is.na(syn.evenness)]<-0
  synch.vars[4]<-mean(syn.evenness[lower.tri(syn.evenness)])
  synch.vars[9]<-sd(syn.evenness[lower.tri(syn.evenness)])
  
  ## Turnover rate
  # turnrate.vec<-rep(NA, times=tmax-1)
  # for(tt in 1:(tmax-1)){
  #   com.t1<-sppXtime[,tt]>0
  #   com.t2<-sppXtime[,tt+1]>0
  #   turnrate.vec[tt]<-sum(com.t1!=com.t2)/(sum(com.t1)+sum(com.t2))
  # }
  # comm.vars[3]<-mean(turnrate.vec)
  
  turnrate.mat<-matrix(NA, nrow=nlocs, ncol=tmax-1)
  for(loc in 1:nlocs){
    for(tt in 1:(tmax-1)){
      com.t1<-inarray[,loc,tt]>0
      com.t2<-inarray[,loc,tt+1]>0
      turnrate.mat[loc,tt]<-sum(com.t1!=com.t2)/(sum(com.t1)+sum(com.t2))
    }
  }
  turnrate.mat[is.na(turnrate.mat)]<-0
  comm.vars[3]<-mean(turnrate.mat, na.rm=T) ##<---- this is average plot-level temporal turnover rates
  
  syn.turnover<-suppressWarnings(cor(t(clean(turnrate.mat,clnlev)),method=corrmeth))
  syn.turnover[is.na(syn.turnover)]<-0
  synch.vars[5]<-mean(syn.turnover[lower.tri(syn.turnover)], na.rm=T)
  synch.vars[10]<-sd(syn.turnover[lower.tri(syn.turnover)], na.rm=T)
  
  return(list(comm.vars=comm.vars, synch.vars=synch.vars)) 
}