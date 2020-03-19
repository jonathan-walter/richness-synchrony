# Do lake morphology, nutrients, or climate intensification explain variability in HAB trends?

# NOTE: run bloom metric trends first to get trend result .csv files 

rm(list=ls())

if (!require(here)) install.packages("dplyr")
if (!require(here)) install.packages("ncf")
if (!require(here)) install.packages("nlme")
if (!require(here)) install.packages("car")
if (!require(here)) install.packages("here")

library(dplyr)
library(ncf)
library(nlme)
library(car)
library(here)

covars<-read.csv(here("lake_covariates.csv"), stringsAsFactors=F)
res.avg.chla<-read.csv(here("results_AvgChl_vs_Year_gls.csv"), stringsAsFactors=F)
res.p95.chla<-read.csv(here("results_p95Chl_vs_Year_gls.csv"), stringsAsFactors=F)
res.tIMug<-read.csv(here("results_tIMug_vs_Year_gls.csv"), stringsAsFactors=F)

#Investigate missingness in covariates ----------------------------------------------------------
sum(is.na(covars$dr)) #27
sum(is.na(covars$maxdepth)) #37
sum(is.na(covars$lake_area_ha)) #24
sum(is.na(covars$tn.trend)) #141
sum(is.na(covars$tp.trend)) #130
sum(is.na(covars$np)) #145
sum(is.na(covars$wetland_pct)) #14
sum(is.na(covars$agri_pct)) #14
sum(is.na(covars$devel_pct)) #14

#Scale numeric covariates ------------------------------------------------------------------------

covars$lake_area_ha<-scale(log10(covars$lake_area_ha))
covars$maxdepth<-scale(log10(covars$maxdepth))
covars$dr<-scale(log10(covars$dr))
covars$tn.trend<-scale(covars$tn.trend)
covars$tp.trend<-scale(covars$tp.trend)
covars$np<-scale(covars$np)
covars$trend.pptp95<-scale(covars$trend.pptp95)
covars$trend.tmaxmean<-scale(covars$trend.tmaxmean)
covars$agri_pct<-scale(covars$agri_pct)
covars$devel_pct<-scale(covars$devel_pct)
covars$wetland_pct<-scale(covars$wetland_pct)


# Average chlorophyll-a ---------------------------------------------------------------------------

dat.avg.chla<-inner_join(covars, res.avg.chla)
dat.avg.chla<-dat.avg.chla[!dat.avg.chla$recovery.flag,]
#table(dat.avg.chla$class)

lm.avg.chla<-lm(b1 ~ lake_area_ha + maxdepth + dr + tn.trend + tp.trend + np + trend.pptp95 + 
                  trend.tmaxmean + agri_pct + devel_pct + wetland_pct, data=dat.avg.chla, na.action="na.exclude")
# there was spatial autocorrelation so looking at the GLS version
lm.avg.chla2<-gls(b1 ~ lake_area_ha + maxdepth + dr + tn.trend + tp.trend + np + trend.pptp95
                  + trend.tmaxmean  + agri_pct + devel_pct + wetland_pct, data=dat.avg.chla,
                 correlation=corExp(form=~nhd_long + nhd_lat), na.action="na.exclude")
summary(lm.avg.chla)
summary(lm.avg.chla2)
vif(lm.avg.chla)
hist(lm.avg.chla$residuals)
rescor.avg.chla<-spline.correlog(x=dat.avg.chla$nhd_long, y=dat.avg.chla$nhd_lat, z=residuals(lm.avg.chla), latlon = T, na.rm=T)
plot(rescor.avg.chla)

# 95th percentile chlorophyll-a -------------------------------------------------------------------

dat.p95.chla<-inner_join(covars, res.p95.chla)
dat.p95.chla<-dat.p95.chla[!dat.p95.chla$recovery.flag,]

lm.p95.chla<-lm(b1 ~ lake_area_ha + maxdepth + dr + tn.trend + tp.trend + np + trend.pptp95 + 
                  trend.tmaxmean + agri_pct + devel_pct + wetland_pct, data=dat.p95.chla, na.action="na.exclude")
# there was spatial autocorrelation so looking at the GLS version
lm.p95.chla2<-gls(b1 ~ lake_area_ha + maxdepth + dr + tn.trend + tp.trend + np + trend.pptp95 + 
                    trend.tmaxmean +trend.tmaxmean + agri_pct + devel_pct + wetland_pct, data=dat.p95.chla,
                  correlation=corExp(form=~nhd_long + nhd_lat),na.action="na.exclude")
summary(lm.p95.chla)
summary(lm.p95.chla2)
vif(lm.p95.chla)
hist(lm.p95.chla$residuals)
rescor.p95.chla<-spline.correlog(x=dat.p95.chla$nhd_long, y=dat.p95.chla$nhd_lat, z=residuals(lm.p95.chla), latlon=T, na.rm=T)
plot(rescor.p95.chla)

# Proportion of Observations Impaired -------------------------------------------------------------

dat.tIMug<-right_join(covars, res.tIMug)
dat.tIMug<-dat.tIMug[!dat.tIMug$recovery.flag,]

lm.tIMug<-lm(b1 ~ lake_area_ha + maxdepth + dr + tn.trend + tp.trend + np +
                  trend.pptp95 + trend.tmaxmean + agri_pct + devel_pct + wetland_pct, data=dat.tIMug)
summary(lm.tIMug)
vif(lm.tIMug)
hist(lm.tIMug$residuals)
rescor.tIMug<-spline.correlog(x=dat.tIMug$nhd_long, y=dat.tIMug$nhd_lat, z=residuals(lm.tIMug))
plot(rescor.tIMug)

# Make a barplot of coefficients -----------------------------------------------------------------
coeffmat<-rbind(lm.avg.chla2$coefficients,
                lm.p95.chla2$coefficients,
                lm.tIMug$coefficients)

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

se<-c(summary(lm.avg.chla2)$tTable[-1,2],
      summary(lm.p95.chla2)$tTable[-1,2],
      summary(lm.tIMug)$coefficients[-1,2])


pdf(here("covariates_barplot.pdf"), width=6.5, height=4)

#Colors for plot
durationcol <- rgb(102,102,102, max = 255, alpha = 200)
severitycol <- rgb(16,78,139, max = 255, alpha = 200)
intensecol <- rgb(0,139,139, max = 255, alpha = 200)

par(mar=c(6.1,4.1,1.1,1.1))

bp<-barplot(coeffmat[,-1], beside=T, col=c(intensecol, severitycol, durationcol),
            names.arg=c("Lake Area","Lake Depth","CA:SA","TN trend","TP trend","N:P","Precip trend","Temp trend","% Agriculture","% Developed","% Wetland"),
            ylab="Standardized regression coefficient", ylim=c(-0.05,0.08),las=2)

#add error bars
error.bar(x=c(t(bp)), y=c(t(coeffmat[,-1])),upper=2*se,length=0.01, col="black")

legend("top", pch=15, col=c(intensecol, severitycol, durationcol),
       legend=c("Intensity","Severity","Duration"), ncol=3, bty="n", pt.cex=2)
abline(0,0)
dev.off()
