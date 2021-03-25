#' Perform exploratory community synchrony analyses on a data array
#' 
#' UPDATES COMPLETED 5/31/2019, TESTED on JRG DATA
#' 
#' @param inarray a species x locations x time data array
#' @param clnlev cleaning level options. If clnlev==1, detrend the time series. If clnlev==2, detrend and rescale,
#' @param corrmeth correlation method, one of "spearman","pearson", or "kendall"
#' @param patch.het a factor vector specifying the patch type for each location. If NULL (the default) this is ignored.
#' @param do.SynchSig T of F, significance test spatial synchrony.
#' @param do.Modularity T or F, do modularity analyses.
#' @param do.Mantel T or F, do Mantel test analyses.
#' 
#' @return \code{CommSpatSynch}
#' 
#' @note Detailed methods given in a working document provided to the LTER Synchrony working group.
#' As of 6 July 2018, NMDS is nolonger computed and placeholders have been removed from output.
#' 
#' @author Jonathan Walter, \email{jaw3es@@virginia.edu}
#'
#' @references Walter, J. A., et al. (2017) The geography of spatial synchrony. Ecology Letters. doi: 10.1111/ele.12782
#'  
#' @export

CommSpatSynch<-function(inarray, clnlev=2, corrmeth="spearman", patch.het=NULL, do.SynchSig=T, 
                        do.Modularity=T, do.Mantel=T){
  
  library(codyn)
  library(ecodist)
  library(vegan)
  library(wsyn)
  library(betapart)
  
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
  comm.vars<-rep(NA, 9)
  names(comm.vars)<-c("AvgPlotRich","Evenness","Turnover","Jaccard","Jacc.tu","Jacc.ne","CVTotBiomass","LoreauSynch","VarRatio")
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
  sppXtime<-apply(inarray, c(1,3), mean,  na.rm=T)
  sppXplot<-apply(inarray, c(2,1), mean,  na.rm=T)
  
  TotBiomassXtime<-apply(inarray,3,sum,na.rm=T)
  comm.vars[7]<-sd(TotBiomassXtime)/mean(TotBiomassXtime)
  
  ## Loreau synchrony
  indices<-which(sppXtime==sppXtime, arr.ind=T)
  data_long<-data.frame(time = indices[,2], species = indices[,1], abundance=c(sppXtime))
  comm.vars[8]<-synchrony(data_long, time.var="time", species.var="species", abundance.var="abundance")
  
  ## Variance Ratio
  comm.vars[9]<-variance_ratio(data_long, time.var="time", species.var="species", abundance.var="abundance", bootnumber=1)[4]
  
  ## Beta diversity (Jaccard similarity)
  btcore<-betapart.core(ifelse(sppXplot>0,1,0))
  jaccdiss.part<-beta.pair(btcore, index.family = "jaccard")
  jaccsim<-full(1-jaccdiss.part$beta.jac)
  jacctu<-full(1-jaccdiss.part$beta.jtu)
  jaccne<-full(1-jaccdiss.part$beta.jne)
  comm.vars[4]<-mean(jaccsim[lower.tri(jaccsim)], na.rm=T)
  comm.vars[5]<-mean(jacctu[lower.tri(jacctu)], na.rm=T)
  comm.vars[6]<-mean(jaccne[lower.tri(jaccne)], na.rm=T)
  
  
  # jaccsim<-1-full(distance(sppXplot>0, method="jaccard"))
  # comm.vars[4]<-mean(jaccsim[lower.tri(jaccsim)], na.rm=T)
  
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
  
  ## Sig test community synchronies ########################################

  if(do.SynchSig){
    sigtest<-function(emp, surr){
      percentile<-ecdf(surr)
      return(1-percentile(emp))
    }
    surrsynch<-function(surrogs){
      surrsynchmats<-lapply(lapply(surrogs,t), cor, method=corrmeth)
      out<-rep(NA, length(surrogs))
      for(ii in 1:length(surrogs)){
        out[ii]<-mean(surrsynchmats[[ii]][lower.tri(surrsynchmats[[ii]])], na.rm=T)
      }
      return(out)
    }
    ##Average population synchrony
    surrcommarray<-array(NA, dim=c(dim(inarray),1000))
    for(spp in 1:dim(inarray)[1]){
      surrcommarray[spp,,,]<-simplify2array(wsyn::surrog(clean(inarray[spp,,],clnlev), 
                                                         1000, "aaft", syncpres=F),higher=T)
    }
    surrogsynch<-rep(NA, 1000)
    for(rep in 1:1000){
      surrcomm<-surrcommarray[,,,rep]
      surrsyncharray<-array(NA,dim=c(nspp,nlocs,nlocs))
      for(spp in 1:nspp){
        surrsyncharray[spp,,]<-suppressWarnings(cor(t(surrcomm[spp,,]),method=corrmeth,use="pairwise.complete.obs"))
      }
      surrsynchmean<-apply(surrsyncharray, c(2,3), mean, na.rm=T)
      surrsynchmean[is.nan(surrsynchmean)]<-0 ##Some sites never share a species--set correlation here to 0
      surrogsynch[rep]<-mean(surrsynchmean[lower.tri(surrsynchmean)],na.rm=T)
    }
    synch.sig[1]<-sigtest(emp=mean(syn.spp.mean[lower.tri(syn.spp.mean)]), surr=surrogsynch)
    rm(surrsyncharray)
    #Synchrony in total biomass
    surr.TotBiomass<-wsyn::surrog(cleandat(total.biomass,yrs,clev=3)$cdat, 1000, "aaft", syncpres=F)
    synch.sig[2]<-sigtest(emp=mean(syn.total.biomass[lower.tri(syn.total.biomass)]), surr=surrsynch(surr.TotBiomass))
    #Synchrony in richness
    surr.Richness<-wsyn::surrog(cleandat(richness,yrs,clev=3)$cdat, 1000, "aaft", syncpres=F)
    synch.sig[3]<-sigtest(emp=mean(syn.richness[lower.tri(syn.richness)]), surr=surrsynch(surr.Richness))
    #Synchrony in evenness
    surr.Evenness<-wsyn::surrog(cleandat(evenness,yrs,clev=3)$cdat, 1000, "aaft", syncpres=F)
    synch.sig[4]<-sigtest(emp=mean(syn.evenness[lower.tri(syn.evenness)]), surr=surrsynch(surr.Evenness))
    #Synchrony in turnover rate
    surr.Turnover<-wsyn::surrog(cleandat(turnrate.mat,1:ncol(turnrate.mat),clev=3)$cdat, 1000, "aaft", syncpres=F)
    synch.sig[5]<-sigtest(emp=mean(syn.turnover[lower.tri(syn.turnover)]), surr=surrsynch(surr.Evenness))
  }
  
  ## Mantel Tests ##########################################################
  #TODO: renovate these codes so output goes to the right places

  if(do.Mantel){
    mantel.results[1,]<-ecodist::mantel(as.dist(syn.total.biomass)~as.dist(syn.spp.mean))[c(1,4)]
    mantel.results[2,]<-ecodist::mantel(as.dist(syn.total.biomass)~as.dist(syn.richness))[c(1,4)]
    mantel.results[3,]<-ecodist::mantel(as.dist(syn.total.biomass)~as.dist(syn.evenness))[c(1,4)]
    mantel.results[4,]<-ecodist::mantel(as.dist(syn.total.biomass)~as.dist(syn.turnover))[c(1,4)]
    mantel.results[5,]<-ecodist::mantel(as.dist(syn.total.biomass)~as.dist(jaccsim))[c(1,4)]
    mantel.results[6,]<-ecodist::mantel(as.dist(syn.spp.mean)~as.dist(syn.richness))[c(1,4)]
    mantel.results[7,]<-ecodist::mantel(as.dist(syn.spp.mean)~as.dist(syn.evenness))[c(1,4)]
    mantel.results[8,]<-ecodist::mantel(as.dist(syn.spp.mean)~as.dist(syn.turnover))[c(1,4)]
    mantel.results[9,]<-ecodist::mantel(as.dist(syn.spp.mean)~as.dist(jaccsim))[c(1,4)]
    mantel.results[10,]<-ecodist::mantel(as.dist(syn.richness)~ as.dist(syn.evenness))[c(1,4)]
    mantel.results[11,]<-ecodist::mantel(as.dist(syn.richness) ~ as.dist(syn.turnover))[c(1,4)]
    mantel.results[12,]<-ecodist::mantel(as.dist(syn.richness) ~ as.dist(jaccsim))[c(1,4)]
    mantel.results[13,]<-ecodist::mantel(as.dist(syn.evenness) ~ as.dist(syn.turnover))[c(1,4)]
    mantel.results[14,]<-ecodist::mantel(as.dist(syn.evenness) ~ as.dist(jaccsim))[c(1,4)]
    mantel.results[15,]<-ecodist::mantel(as.dist(syn.turnover) ~ as.dist(jaccsim))[c(1,4)]
  }

  if(do.Modularity){
  ## modularity
    diag(syn.spp.mean)<-0
    mod.syn.spp.mean<-cluseigen(syn.spp.mean)
    diag(syn.total.biomass)<-0
    mod.total.biomass<-cluseigen(syn.total.biomass)
    diag(syn.richness)<-0
    mod.syn.richness<-cluseigen(syn.richness)
    diag(syn.evenness)<-0
    mod.syn.evenness<-cluseigen(syn.evenness)
    diag(syn.turnover)<-0
    mod.syn.turnover<-cluseigen(syn.turnover)
    
    module.info[1]<-max(unlist(mod.syn.spp.mean))
    module.info[2]<-max(unlist(mod.total.biomass))
    module.info[3]<-max(unlist(mod.syn.richness))
    module.info[4]<-max(unlist(mod.syn.evenness))
    module.info[5]<-max(unlist(mod.syn.turnover))
    module.info[6]<-wsyn::modularity(syn.spp.mean, mod.syn.spp.mean[[length(mod.syn.spp.mean)]])
    module.info[7]<-wsyn::modularity(syn.total.biomass, mod.total.biomass[[length(mod.total.biomass)]])
    module.info[8]<-wsyn::modularity(syn.richness, mod.syn.richness[[length(mod.syn.richness)]])
    module.info[9]<-wsyn::modularity(syn.evenness, mod.syn.evenness[[length(mod.syn.evenness)]])
    module.info[10]<-wsyn::modularity(syn.turnover, mod.syn.turnover[[length(mod.syn.turnover)]])
    
    ##Test of association between habitat patch classes and module membership
    if(!is.null(patch.het)){ 
      if(!max(unlist(mod.syn.spp.mean))==1){module.correspondence[1]<-chisq.test(table(as.factor(patch.het), as.factor(mod.syn.spp.mean$membership)))$p.value}
      if(!max(unlist(mod.total.biomass))==1){module.correspondence[2]<-chisq.test(table(as.factor(patch.het), as.factor(mod.total.biomass$membership)))$p.value}
      if(!max(unlist(mod.syn.richness))==1){module.correspondence[3]<-chisq.test(table(as.factor(patch.het), as.factor(mod.syn.richness$membership)))$p.value}
      if(!max(unlist(mod.syn.evenness))==1){module.correspondence[4]<-chisq.test(table(as.factor(patch.het), as.factor(mod.syn.evenness$membership)))$p.value}
      if(!max(unlist(mod.syn.evenness))==1){module.correspondence[6]<-chisq.test(table(as.factor(patch.het), as.factor(mod.syn.turnover$membership)))$p.value}
    }
  }
  
  return(list(comm.vars=comm.vars, synch.vars=synch.vars, module.info=module.info, mantel.results=mantel.results, module.correspondence=module.correspondence,
              synch.sig=synch.sig)) 
}