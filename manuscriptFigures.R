## This script produces manuscript-quality figures synthesizing empirical and theoretical results.

rm(list=ls())

library(here)
library(RColorBrewer)

#load empirical data
emp.dat<-read.csv(here("empirical/analysisvars_table.csv"), stringsAsFactors = F)

#load theory data -- this CSV was assembled offline from cluster output. We are providing the
#results, and simulation model function R script; see manuscript for a description of the
#distributions parameter values were drawn from.
#thry.dat<-read.csv(here("theory/pre_disp_theory_output.csv"))

# local dispersal
thry.dat<-read.csv(here("theory/newv2_local_params2500_results.csv"))

## ---------------------------------------------------------------------------------------
## Figure 1: Is there spatial synchrony in species richness?

emp.fig1<-emp.dat[,colnames(emp.dat) %in% c("dataset","site","organism","rRichness","sd.rRichness")]
emp.fig1<-emp.fig1[order(emp.fig1$organism, decreasing=T),]
emp.fig1$displayname<-c("HAY","JRG","JRN-BASN","JRN-IBPE","JRN-SUMM","KNZ-UP","KNZ-LOW","SEV-B",
                        "SEV-C","SEV-G","DRT","LOK","MAU","MCR-BACK","MCR-FRNG","MCR-OUT","MDK","SBC","UPK","USVI")

bb<-rgb(0,114,178,255,maxColorValue=255)
bg<-rgb(0,148,115,255,maxColorValue=255)

pdf(here("figures/main/Fig1_richsynch_combined.pdf"), width=6.5, height=3.25)

layout(matrix(c(1,2),nrow=1),widths=c(0.45,0.55))

par(mar=c(1.5,3.7,1.5,1.1), oma=c(1.3,0,0,0), tcl=-0.4, mgp=c(3,0.5,0))

hist(thry.dat$rRichness, breaks=seq(-.1,1,by=.1), freq=TRUE, xlim=c(-2/20,1), 
          xlab="", main="", xaxt="n", yaxt="n", ylab="", col="lightgrey")
axis(side=1, labels = T, at=c(0, .25, .5, .75, 1), pos=0/20)                             
axis(side=2, at=seq(0,700,100), labels=seq(0,700,100)/2500, las=1, pos=-2/20)
mtext("Frequency", side=2, line=2.3)
mtext("Richness Synchrony", side=1, line=0, outer=T)
mtext("Theoretical", line=0.25, side=3)
mtext("A)",at=-.1)

par(mar=c(1.5,4.1,1.5,1.5))

barplot(height=emp.fig1$rRichness, names.arg=emp.fig1$displayname, las=2, horiz=T, xlab="",
        col=c(rep(bg,10),rep(bb,10)), xlim=c(0,1),offset=0, xaxt="n", cex.names=0.7)
axis(side=1, labels=T, at=c(0,0.25,.5,0.75,1),pos=0)
mtext("Empirical", side=3, line=0.25)
legend("topright",legend=rev(c("Grassland","Marine")),fill=rev(c(bg,bb)), cex=0.9, bty="n")
mtext("C)",at=0)

par(fig=c(0.17,0.4,0.4,0.85), mar=c(2.1,2.1,1,0.2), mgp=c(0.8,0.1,0), cex=0.6, tcl=-0.2, new=T)
hist(emp.dat$rRichness, breaks=seq(-.1,1,by=0.1), freq=TRUE, xlim=c(-.1,1), main="", xlab="", ylab="",
     xaxt="n",yaxt="n", col="lightgrey")
axis(side=1, labels = T, at=c(0, .25, .5, .75, 1), pos=0/20)                             
axis(side=2, at=0:4, labels=c("0","","2","","4"), las=1, pos=-2/20, mgp=c(0.8,0.3,0))
mtext("Empirical",cex=0.6)
mtext("B)",at=-0.1,cex=0.6)

dev.off()

##---------------------------------------------------------------------------------
## Figure 2: Theoretical parameter sensitivity

thry.std<-thry.dat[,colnames(thry.dat) %in% 
                     c("rRichness","AvgPlotRich","Evenness","Jaccard","Turnover","Npatches","growrt","competition"
                       ,"env_var","disp","patch_het","autocorr","jacc.tu","jacc.ne")]
thry.std$AvgPlotRich<-scale(thry.std$AvgPlotRich)
thry.std$Evenness<-scale(thry.std$Evenness)
thry.std$Jaccard<-scale(thry.std$Jaccard)
thry.std$Turnover<-scale(thry.std$Turnover)
thry.std$Npatches<-scale(thry.std$Npatches)
thry.std$growrt<-scale(thry.std$growrt)
thry.std$competition<-scale(thry.std$competition)
thry.std$env_var<-scale(thry.std$env_var)
thry.std$disp<-scale(thry.std$disp)
thry.std$patch_het<-scale(thry.std$patch_het)
thry.std$autocorr<-scale(thry.std$autocorr)
thry.std$jacc.ne<-scale(thry.std$jacc.ne)
thry.std$jacc.tu<-scale(thry.std$jacc.tu)

param.effects<-lm(rRichness~env_var + autocorr + Npatches + patch_het + growrt + competition + disp, data=thry.std)
summary(param.effects)
# param.effects2<-lm(rRichness~env_var + autocorr + Npatches + patch_het + growrt + competition + disp + env_var:patch_het, data=thry.std)
# summary(param.effects2)

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

#pal=brewer.pal(3,"Accent")

pal=c(rgb(0,0,0,255,maxColorValue=255),rgb(230,159,0,255,maxColorValue=255),rgb(86,180,233,255,maxColorValue=255))


pdf(here("figures/main/Fig2_thry_effects.pdf"), width=3.25, height=3.25)

par(mar=c(6.1,3.5,0.7,0.7), tcl=-0.4, mgp=c(1,0.5,0))

barplot(param.effects$coefficients[2:8], ylim=c(-.02, .15), names.arg=c("Temporal \n Variation", "Temporal \n Autocorrelation", 
                                                                           "# Patches", "Patch \n Heterogeneity",
                                                                           "Growth \n Rate", "Competitive \n Strength",
                                                                           "Dispersal"), 
        cex.names = .8, las=2, cex.axis = .8, col=pal[c(1,1,2,2,3,3,3)])
abline(h=0)
error.bar(seq(from=0.7,by=1.2,length.out = 7),param.effects$coefficients[-1],summary(param.effects)$coefficients[-1,2],length=0.05)
mtext("Effect size",2,line=2.2)
legend("topright",fill=pal,legend=c("Abiotic temporal","Abiotic spatial","Demographic"), bty="n")

dev.off()

##--------------------------------------------------------------------
## Figure 3: Empirical effects

emp.std<-emp.dat[,colnames(emp.dat) %in% c("rRichness","organism","extent","n.taxa.regional","Evenness","Jaccard",
                                           "Jacc.tu","Jacc.ne","Turnover","AvgPlotRich","n.plots")]
emp.std$extent<-scale(emp.std$extent)
emp.std$n.taxa.regional<-scale(emp.std$n.taxa.regional)
emp.std$Evenness<-scale(emp.std$Evenness)
emp.std$Jaccard<-scale(emp.std$Jaccard)
emp.std$Jacc.tu<-scale(emp.std$Jacc.tu)
emp.std$Jacc.ne<-scale(emp.std$Jacc.ne)
emp.std$Turnover<-scale(emp.std$Turnover)
emp.std$AvgPlotRich<-scale(emp.std$AvgPlotRich)
emp.std$n.plots<-scale(emp.std$n.plots)

emp.effects<-lm(rRichness~organism+extent+AvgPlotRich+Evenness+Jaccard+Turnover, data=emp.std)
summary(emp.effects)

emp.effects2<-lm(rRichness~organism+extent+AvgPlotRich+Evenness+Jacc.tu+Jacc.ne+Turnover, data=emp.std)
summary(emp.effects2)

thry.effects2<-lm(rRichness~AvgPlotRich+Evenness+Jaccard+Turnover, data=thry.std)
summary(thry.effects2)

thry.effects3<-lm(rRichness~AvgPlotRich+Evenness+jacc.tu+jacc.ne+Turnover, data=thry.std)
summary(thry.effects3)


pdf(here("figures/main/Fig3_emp_effects.pdf"), width=3.25, height=3.25)

par(mar=c(4.1,3.5,0.7,0.7), tcl=-0.4, mgp=c(1,0.5,0))

barplot(rbind(emp.effects$coefficients[2:7], c(NA,NA,thry.effects2$coefficients[-1])),
        ylim=c(-0.11, 0.32), names.arg=c("Biome", "Extent", "Richness", "Evenness", "Beta \n diversity", "Turnover"),
        cex.names = .8, las=2, cex.axis = .8, beside=T,
        col=rep(c("grey40","grey80"),6))
abline(h=0)
error.bar(seq(from=1.5,by=3,length.out = 7),emp.effects$coefficients[-1],summary(emp.effects)$coefficients[-1,2],length=0.03)
error.bar(seq(from=8.5,by=3,length.out = 5), thry.effects2$coefficients[-1],summary(thry.effects2)$coefficients[-1,2],length=0.03)
mtext("Effect size",2,line=2.2)
legend("topleft",fill=c("grey40","grey80"),legend=c("Empirical","Theoretical"),bty="n", cex=0.8, inset=c(0.075,0))

dev.off()

# pdf(here("figures/main/Fig3_emp_effects.pdf"), width=3.25, height=3.25)
# 
# par(mar=c(4.1,3.5,0.7,0.7), tcl=-0.4, mgp=c(1,0.5,0))
# 
# barplot(rbind(emp.effects2$coefficients[-1], c(NA,NA,thry.effects3$coefficients[-1])), 
#         ylim=c(-0.11, 0.35), names.arg=c("Biome", "Extent", "Richness", "Evenness", "Beta-repl.","Beta-nest.", "Turnover"), 
#         cex.names = .8, las=2, cex.axis = .8, beside=T,
#         col=rep(c("grey40","grey80"),7))
# abline(h=0)
# error.bar(seq(from=1.5,by=3,length.out = 8),emp.effects2$coefficients[-1],summary(emp.effects2)$coefficients[-1,2],length=0.03)
# error.bar(seq(from=8.5,by=3,length.out = 6), thry.effects3$coefficients[-1],summary(thry.effects3)$coefficients[-1,2],length=0.03)
# mtext("Effect size",2,line=2.2)
# legend("topleft",fill=c("grey40","grey80"),legend=c("Empirical","Theoretical"),bty="n", cex=0.8, inset=c(0.075,0))
# 
# dev.off()

##--------------------------------------------------------------------------
## Figure 4: community cv ~ richness synchrony and avg richness

ptcol<-c(bb,bg,bg,bg,bg,bg,bg,bg,bb,bb,bb,bb,bb,bb,bb,bg,bg,bg,bb,bb)

fit3<-lm(I(-1*thry.dat$cv)~thry.dat$rRichness)
cor.test(thry.dat$cv,thry.dat$rRichness)
summary(fit3)

fit4<-lm(I(-1*emp.dat$CVTotBiomass)~emp.dat$rRichness)
cor.test(-1*emp.dat$CVTotBiomass,emp.dat$rRichness)
summary(fit4)

cor.test(-1*thry.dat$cv,thry.dat$AvgPlotRich)
fit5<-lm(I(-1*thry.dat$cv)~thry.dat$AvgPlotRich)
summary(fit5)

fit6<-lm(I(-1*emp.dat$CVTotBiomass)~emp.dat$AvgPlotRich)
cor.test(-1*emp.dat$CVTotBiomass,emp.dat$AvgPlotRich)
summary(fit6)

#### 

pdf(here("figures/main/Fig4_cv_richsynch.pdf"), width=6.5, height=6.5)

par(mar=c(3.1,3.1,1.75,0.5), tcl=-0.4, mgp=c(1.75,0.5,0), mfcol=c(2,2), oma=c(0,0,0,0.6))

plot(thry.dat$rRichness, -1*thry.dat$cv, xlab="Richness synchrony", ylab="Stability", xlim=c(0,1),
     pch=20, col="grey")
abline(fit3, lwd=2)
#mtext("Theoretical",3,line=0.2)
mtext(expression(paste(italic(R)^2,"=0.22, ",hat(beta),"=-0.29")),3,line=0.1,cex=0.9)
#mtext(expression(paste(italic(R2),"=0.48, ",italic(beta),"=0.28")),3,line=-1.3,cex=0.9)
mtext("A)",at=0.01,line=0.3)

plot(emp.dat$rRichness, -1*emp.dat$CVTotBiomass, xlab="Richness synchrony", ylab="Stability", xlim=c(0,1), pch=20, col=ptcol, cex=1.75)
abline(fit4,lwd=2)
#mtext("Empirical",3,line=0.2)
mtext(expression(paste(italic(R)^2,"=0.42, ",hat(beta),"=-0.58")),3,line=0.1,cex=0.9)
#mtext(expression(paste(italic(r),"=0.65, ",italic(p),"=0.002, ",italic(beta),"=0.58")),3,line=-1.3,cex=0.9)
mtext("C)",at=0.01,line=0.3)

plot(thry.dat$AvgPlotRich, -1*thry.dat$cv, xlab="Richness", ylab="Stability", pch=20, col="grey")
abline(fit5,lwd=2)
mtext(expression(paste(italic(R)^2,"=0.08, ",hat(beta),"=0.004")),3,line=0.1,cex=0.9)
#mtext(expression(paste(italic(R2),"=0.02, ",italic(beta),"=-0.0006")),3,line=-1.3,cex=0.9)
mtext("B)",at=2,line=0.3)

plot(emp.dat$AvgPlotRich, -1*emp.dat$CVTotBiomass, xlab="Richness", ylab="Stability", pch=20, col=ptcol,cex=1.75)
abline(fit6,lwd=2)
mtext(expression(paste(italic(R)^2,"=0.13, ",hat(beta),"=0.02")),3,line=0.1,cex=0.9)
#mtext(expression(paste(italic(R),"=-0.37, ",italic(p),"=0.11, ",italic(beta),"=-0.02")),3,line=-1.3,cex=0.9)
legend("bottomright",legend=rev(c("Grassland","Marine")),pch=20,col=rev(c(bg,bb)), inset=c(0.03,0.1), cex=1.1)
mtext("D)",at=4.5,line=0.3)

dev.off()

# fit3a<-lm(thry.dat$cv~scale(thry.dat$rRichness))
# fit5a<-lm(thry.dat$cv~scale(thry.dat$AvgPlotRich))
# fit4a<-lm(emp.dat$CVTotBiomass~scale(emp.dat$rRichness))
# fit6a<-lm(emp.dat$CVTotBiomass~scale(emp.dat$AvgPlotRich))

summary(lm(I(-1*thry.dat$cv)~thry.dat$rRichness))
summary(lm(I(-1*thry.dat$cv)~thry.dat$AvgPlotRich))
summary(lm(I(-1*thry.dat$cv)~thry.dat$Evenness))
summary(lm(I(-1*thry.dat$cv)~thry.dat$Turnover))
summary(lm(I(-1*thry.dat$cv)~thry.dat$Jaccard))


summary(lm(I(-1*emp.dat$CVTotBiomass)~emp.dat$rRichness))
summary(lm(I(-1*emp.dat$CVTotBiomass)~emp.dat$AvgPlotRich))
summary(lm(I(-1*emp.dat$CVTotBiomass)~emp.dat$Evenness))
summary(lm(I(-1*emp.dat$CVTotBiomass)~emp.dat$Turnover))
summary(lm(I(-1*emp.dat$CVTotBiomass)~emp.dat$Jaccard))

          
## use a bootstrapping procedure to check confidence interval on R^2

r2_boot<-NULL
for(ii in 1:1000){
  
  samp<-sample(1:20, 20, replace=TRUE)
  cv.ii <- -1*emp.dat$CVTotBiomass[samp]
  turnover.ii <- emp.dat$Turnover[samp]
  
  r2_boot<-c(r2_boot, summary(lm(cv.ii~turnover.ii))$r.squared)
  
}

quantile(r2_boot, c(0.025,0.5,0.975))