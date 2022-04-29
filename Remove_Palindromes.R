
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(dplyr) 

#Import bim file 
bim = read.table(args[1], header=F)

#Get indices of A/T and G/C SNPs
w = which((bim$V5=="A" & bim$V6=="T") |
(bim$V5=="T" & bim$V6=="A") |
(bim$V5=="C" & bim$V6=="G") |
(bim$V5=="G" & bim$V6=="C"))

#Extract A/T and G/C SNPs
at.cg.snps = bim[w,]

#Export
write.table(at.cg.snps$V2,"at-cg-snps.txt", row.names = F, col.names = F, quote = F)
q("no")
