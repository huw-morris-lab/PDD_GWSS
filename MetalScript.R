#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#Script to generate files ready for meta-analysis using METAL (https://genome.sph.umich.edu/wiki/METAL_Documentation)
#args[1] = input CPH results file
#args[2] = input frequency file
#args[3] = output file

#Load libraries
library(dplyr)
library(data.table)
library(tidyr)

#Import data
data <- fread(args[1], header=TRUE)
freq <- fread(args[2])

freq <- freq %>%
	select(-CHR, -NCHROBS) %>%
	rename(A1_freq = A1,
		A2_freq = A2)

#Sort by p value
data <- data %>% arrange(Pvalue)

#Split SNP ID into chr, position and alleles
data_split <- data %>%
	separate(SNP, into = c("chr", "bp", "REF", "ALT", "A1"))
  
#Check instances where A1 is not the same as ALT 
#A1 is the effect allele in the survival model
data_split  <- data_split %>%
  mutate(effect = ifelse(A1 == ALT, ALT,
                          ifelse(A1 == REF, REF, NA)),
         noneffect = ifelse(A1 == ALT, REF,
                          ifelse(A1 == REF, ALT, NA)))
                          
#Remove indels
data_split_noindels <- data_split %>%
	filter(!is.na(effect)) %>%
	filter(!is.na(noneffect))
  
#Merge frequencies file with results file
data_split_noindels_freq <- data_split_noindels %>%
	mutate(SNP = paste(chr,bp,REF,ALT, sep = ":")) %>%
	inner_join(freq, by = "SNP")
  
#Remove any A1 allele mismatches
data_split_noindels_freq <- data_split_noindels_freq %>%
	filter(A1==A1_freq)
  
#Export for METAL
export_METAL_hg19 <- data_split_noindels_freq %>%
	mutate(SNP_new = paste(chr, bp, sep = ":")) %>%
	select(SNP_new, effect, noneffect, Coeff, se, Pvalue, N, MAF) %>%
	rename(SNP = SNP_new,
	effect_allele = effect,
		noneffect_allele = noneffect,
		beta = Coeff)

fwrite(export_METAL_hg19, file=args[3] , quote = F, sep = "\t", row.names = F)

q("no")


