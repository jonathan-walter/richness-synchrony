## This script organizes for reproducibility all analyses for the richness synchrony manuscript
## by Walter and colleagues. 

rm(list = ls())

# Load or install necessary libraries
for (package in c('tidyverse', 'dplyr', 'ggplot2', 'ecodist', 'abind', 'geosphere', 'rgdal',
                  'maps', 'reshape2', 'codyn', 'igraph', 'vegan', 'devtools','here','wsyn','betapart')) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}
rm(package)

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
source(here("empirical/site_summary_table/makeSiteTable.R"))

## Make manuscript figures ------------------------------------------------------------------------
source(here("manuscriptFigures.R"))
