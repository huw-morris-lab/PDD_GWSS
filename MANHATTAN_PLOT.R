#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#Based on https://github.com/danielroelfs/danielroelfs.com/tree/main/docs/blog/how-i-create-manhattan-plots-using-ggplot

#Load libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(ggtext)
library(ggrepel)

# Load results file from meta-analysis containing the following columns: CHR, BP, SNP, P
gwas_data <- fread(args[1], header = T)

# Order positions
gwas_data_ordered <- gwas_data[
  with(gwas_data, order(CHR, BP)), ]

# Select only SNPs with P>0.05
gwas_data_ordered_trimed <- gwas_data_ordered %>% 
  filter(-log10(pval)>1)

data_cum <- gwas_data_ordered_trimed %>% 
  group_by(CHR) %>% 
  summarise(max_bp = as.numeric(max(BP))) %>% 
  mutate(bp_add = cumsum(max_bp)-max_bp) %>% 
  select(-max_bp) %>%
  left_join(gwas_data_ordered, .,by = c("CHR"="CHR")) %>% 
  arrange(CHR,BP) %>%
  mutate(bp_cum = BP + bp_add)

# Prepare X axis
axis_set <- data_cum %>% 
  group_by(CHR) %>% 
  summarize(center = mean(bp_cum))

# Prepare Y axis
ylim <- data_cum %>% 
  filter(pval == min(pval)) %>% 
  mutate(ylim = abs(floor(log10(pval))) + 1) %>% 
  pull(ylim)

sig <- 5e-8

# Create object snpsOfInterest to highlight top SNPs
snpsOfInterest <- c("2:142000271","5:1415068","6:7320911","19:45411941"))

data_cum <- data_cum %>% 
 mutate( is_highlight=ifelse(rsID %in% snpsOfInterest, "yes", "no")) %>%
  mutate( is_annotate=-log10(pval) > -log10(sig), "yes", "no") 

#Nominate nearest gene of top SNPs
data_cum$GENE <- ifelse(data_cum$SNP == "2:142000271", "LRP1B", 
                    ifelse(data_cum$SNP == "5:1415068", "SLC6A3", 
                      ifelse(data_cum$SNP == "6:7320911", "SSR1",
                        ifelse(data_cum$SNP == "19:45411941", "APOE", NA))))
                        
# Make the plot
manhplot <- ggplot(data_cum, aes(x=bp_cum, y=-log10(P))) +
    # Show all points
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("#E69F00","#56B4E9"), 22 )) +
    # Add significance line
    geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
    # custom X and Y axis:
    scale_x_continuous( label = axis_set$CHR, breaks= axis_set$center ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +     # remove space between plot area and x axis
    # Add highlighted points
    geom_point(data=subset(data_cum, is_highlight=="yes"), color="#b30000", size=1.3) +
    # Add label using ggrepel to avoid overlapping
    geom_text_repel(data=subset(data_cum, is_annotate=="yes"), aes(label=GENE), size=2.5, min.segment.length = 0) +
    # Edit axis labels
    labs(x = "Chromosome", 
       y = "-log<sub>10</sub>(p)") + 
    # Custom the theme:
    theme_classic() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_markdown(),
      axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
    )
ggsave(args[2], manhplot, device = "jpeg", width = 10, height = 8, units = "cm", scale=2, dpi="print")

q("no")
