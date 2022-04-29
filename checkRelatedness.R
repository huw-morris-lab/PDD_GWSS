#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

my_packages <- c("ggplot2", "plinkQC")                       # Specify packages needed
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]    
if(length(not_installed)) install.packages(not_installed)                               # Install uninstalled packages

library(plinkQC)
library(ggplot2)

qcdir <- args[1]          #path/to/directory/with/QC/results containing name.imiss and name.genome results as returned by plink –missing and plink –genome
name <- "IBD_check"

relatednessQC <- evaluate_check_relatedness(
  qcdir,
  name,
  highIBDTh = 0.125,      #Threshold for acceptable proportion of IBD between pair of individuals
  imissTh = 0.03,         #Threshold for acceptable missing genotype rate in any individual
  interactive = FALSE,
  verbose = FALSE
)

list2env(relatednessQC,envir = .GlobalEnv) #to separate the objects

write.table(failIDs,"IBD_DUPLICATES.txt",quote=FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(fail_highIBD[,c(1:5,10)],"IBD_pairs.txt",quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

pdf("Relatedness.pdf")
print(p_IBD)
dev.off() 

q("no")

