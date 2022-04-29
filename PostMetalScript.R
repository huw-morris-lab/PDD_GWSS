#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#args[1] = file with results from METAL

#Load libraries
library(data.table)
library(dplyr)
library(tidyr)

#Read in meta analysis results
data <- fread(args[1])

#Filter out SNPs with total sample size < N total
length_data <- max(data$TotalSampleSize)
data_filtered <- data %>%
	filter(TotalSampleSize >= length_data)

#Sort by p value
data_filtered_sorted <- data_filtered %>%
	arrange(`P-value`)

#Filter out SNPs with HetPVal < 0.05 (Cochran's Q-test for heterogeneity) and HetISq > 80
data_filtered_sorted_het <-data_filtered_sorted %>%
	filter(HetPVal > 0.05) %>%
	filter(HetISq < 80)

#Remove variants with MAF variability > 15%
data_filtered_sorted_het_MAF <- data_filtered_sorted_het %>%
	mutate(MAF_variability = MaxFreq - MinFreq) %>%
	filter(MAF_variability <= 0.15)

#Calculate HR + 95CIs 
data_filtered_sorted_het_MAF$HR <- exp(data_filtered_sorted_het_MAF$Effect) #calculate HR
data_filtered_sorted_het_MAF$lower95 <- exp(data_filtered_sorted_het_MAF$Effect - 1.96*data_filtered_sorted_het_MAF$StdErr) #calculate lower 95% CI
data_filtered_sorted_het_MAF$upper95 <- exp(data_filtered_sorted_het_MAF$Effect + 1.96*data_filtered_sorted_het_MAF$StdErr) #calculate upper 95% CI

#Create columns CHR and BP
data_filtered_sorted_het_MAF <- separate(data = data_filtered_sorted_het_MAF, col = MarkerName, into = c("CHR", "BP"), sep = "\\:")
#Create column rsID with CHR:BP
data_filtered_sorted_het_MAF$rsID <- paste(data_filtered_sorted_het_MAF$CHR,data_filtered_sorted_het_MAF$BP, sep = ":")

#Export for FUMA
#Allele1=non-effect allele, Allele2=effect allele 
export_FUMA <- data_filtered_sorted_het_MAF %>%
	select(CHR,BP,rsID,Allele1,Allele2,`P-value`,Effect,StdErr,HR,lower95,upper95,TotalSampleSize) %>%
	rename(pval = `P-value`,
       non_effect = Allele1,
       effect = Allele2)
		
fwrite(export_FUMA, "metaanalysis_dementia.txt", quote = F, row.names = F, col.names = T, sep = "\t")

q("no")

