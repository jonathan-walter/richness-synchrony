# ---------------------------------------- #
# Make site table - 9 May 2018     #
# ------------------------------------ #

# Clear environment
rm(list = ls())

library(here)

file.list <- list.files(here("empirical/site_summary_table"), pattern = ".csv")

metadata_table <- read.csv(paste(here("empirical/site_summary_table/"), file.list[1], sep=""))
for (file in file.list[2:length(file.list)]){
	 a <- read.csv(paste(here("empirical/site_summary_table/"), file, sep=""))
	 metadata_table <- rbind(metadata_table,a)
}

## ---------------------------------------
## Assemble a table of results from CommSpatSynch

file.list <- list.files(here("empirical/analyses_by_site/output"), pattern = ".rds")

dataset.name <- simplify2array(strsplit(file.list,"_results"))[1,]
results_by_site<-list()
for(i in 1:length(file.list)){	
	results_by_site[[paste0(dataset.name[i])]]<-readRDS(paste0(here("empirical/analyses_by_site/output/"),file.list[i]))
}

response_vars <- NULL
for(i in 1:length(dataset.name)){
  temp <- c(unlist(eval(parse(text = paste0("results_by_site$",dataset.name[i],"$comm.vars")))),
            eval(parse(text = paste0("results_by_site$",dataset.name[i],"$synch.vars"))),
            eval(parse(text = paste0("results_by_site$",dataset.name[i],"$synch.sig"))))
  
  response_vars <- rbind(response_vars, temp)
}

response_vars <- as.data.frame(response_vars)
response_vars$dataset<-dataset.name

response_vars <-response_vars[,colnames(response_vars) %in% c("AvgPlotRich","Evenness","Turnover","Jaccard","Jacc.tu","Jacc.ne","CVTotBiomass","LoreauSynch",
                                                             "VarRatio","rRichness","sd.rRichness","p.rRichness", "rSppMean","dataset")]

## ---------------------------------------
## Merge environmental variables and response variables into the table ##

analysisvars_table <- merge(metadata_table, response_vars, by.x = "dataset", by.y="dataset", all = T)

write.csv(analysisvars_table, file = paste0(here("empirical"),"/analysisvars_table.csv"), row.names=F)

#make table in LaTeX format for manuscript:
#install.packages('xtable')
library('xtable')

table1 <- metadata_table[,c("dataset", "initial.year", "study.length", "n.plots", "extent", "taxa", "n.taxa","abund.type","abund.units")]
print(xtable(table1))

table2 <- response_vars[,c("dataset","rRichness","p.rRichness")]
print(xtable(table2,digits=3))
