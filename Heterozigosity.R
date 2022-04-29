#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(tidyverse)

#Import heterozygosity file (.het)
sample_het <- read.table(args[1], header = TRUE)

#Import call rate file (.miss)
sample_callrate <- read.table(args[2], header = TRUE)

#Calculate proportion of heterozygosity
sample_het <- sample_het %>% 
  mutate(het = (N.NM. - O.HOM.)/N.NM.)

#Calculate mean and SD of heterozygosity
summary <- sample_het %>% 
  summarise(mean_het = mean(het),
            sd_het = sd(het))

mean_het <- summary[1,1]
sd_het <- summary[1,2]

#Write list of samples who are >2 SDs from mean of heterozygosity
sample_het <- sample_het %>% 
  mutate(remove_het = ifelse(het > 2*sd_het + mean_het, "remove",
                             ifelse(het < mean_het - 2*sd_het, "remove", "keep")))

#Merge 
sample_stats <- sample_het %>% 
  left_join(sample_callrate, by = c("FID", "IID"))

#Calculate genotyping rate
sample_stats <- sample_stats %>% 
  mutate(callrate = 1 - F_MISS)

#Plot 
pdf("Heterozigosity_rate vs Missingness.pdf")
ggplot(data = sample_stats, mapping = aes(x = het, y = callrate, color = remove_het)) +
  geom_point() +
  theme_bw()
dev.off()

#Write list of samples to remove - if heterozygosity outliers or if call rate is <98%
samples_to_remove <- sample_stats %>% 
  filter(remove_het == "remove" | callrate < 0.98) %>% 
  select(FID, IID)

write.table(samples_to_remove, "samples_to_remove.txt",
            quote=F, col.names = F, row.names = F)

q("no")
