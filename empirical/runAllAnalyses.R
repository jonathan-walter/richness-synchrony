## This script organizes for reproducibility all analyses for the richness synchrony manuscript
## by Walter and colleagues. 

rm(list = ls())

# Load or install necessary libraries
for (package in c('tidyverse', 'dplyr', 'ggplot2', 'ecodist', 'abind', 'geosphere', 'rgdal',
                  'maps', 'reshape2', 'codyn', 'igraph', 'vegan', 'devtools','here','wsyn')) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}
rm(package)
# Source analysis routine script
source(here("CommSpatSynch_v3.R"))

## ------------------------------------------------------------------------------------------------
## Read in all-sites data files from EDI

# GRASSLANDS --------------------------------------------------------------------------------------
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

# MARINE ------------------------------------------------------------------------------------------
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

## Data and Site-by-Site Results Generation -------------------------------------------------------
source(here("empirical/analyses_by_site/DRT_coral_format_spatial_analysis.R"))
source(here("empirical/analyses_by_site/HAY_format_spatial_analysis.R"))
source(here("empirical/analyses_by_site/JRG_format_spatial_analysis.R"))
source(here("empirical/analyses_by_site/JRN_BASN_format_spatial_analysis.R"))
source(here("empirical/analyses_by_site/JRN_IBPE_format_spatial_analysis.R"))
source(here("empirical/analyses_by_site/JRN_SUMM_format_spatial_analysis.R"))
source(here("empirical/analyses_by_site/KNZ__lowland-t_format_spatial_analysis.R"))
source(here("empirical/analyses_by_site/KNZ__upland-f_format_spatial_analysis.R"))
source(here("empirical/analyses_by_site/LOK_coral_format_spatial_analysis.R"))
source(here("empirical/analyses_by_site/MAU_coral_format_spatial_analysis.R"))
source(here("empirical/analyses_by_site/MCR_coral_algae_backreef_format_spatial_analysis.R"))
source(here("empirical/analyses_by_site/MCR_coral_algae_fringing_format_spatial_analysis.R"))
source(here("empirical/analyses_by_site/MCR_coral_algae_outer10_format_spatial_analysis.R"))
source(here("empirical/analyses_by_site/MDK_coral_format_spatial_analysis.R"))
source(here("empirical/analyses_by_site/SBC_format_spatial_analysis.R"))
source(here("empirical/analyses_by_site/SEV_B_format_spatial_analysis.R"))
source(here("empirical/analyses_by_site/SEV_C_format_spatial_analysis.R"))
source(here("empirical/analyses_by_site/SEV_G_format_spatial_analysis.R"))
source(here("empirical/analyses_by_site/UPK_coral_format_spatial_analysis.R"))
source(here("empirical/analyses_by_site/USVI_coral_format_spatial_analysis.R"))

## Make table of site information -----------------------------------------------------------------
source("./GeogSynch/Scripts/make_site_table/make_site_table.R")

## Collate results of synchrony significance tests ------------------------------------------------
source("./GeogSynch/Scripts/empirical_analyses/synch_sig_summary.R")

## Make 