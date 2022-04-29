#!/usr/bin/env Rscript

#Load libraries
library(data.table)
library(dplyr)

#Import data
  data_eigen=fread("cohortAndHapmap_PCA.eigenvec")
  names(data_eigen)[1] <- "FID"
  fam=fread("cohortAndHapmap.fam")
  names(fam)[1] <- "FID"
  names(fam)[2] <- "IID"
  merge_eigenvec.fam <- merge(data_eigen, fam, by=c("FID","IID")) 
  redux <- merge_eigenvec.fam[,c(1,2,16,3:12)]

#Select EUR (4=CEU & 12=TSI) from HapMap  
  EUR <- redux %>% filter(V6=="4" | V6=="12") 

#Calculate mean +/- SD for PC1:10 
  mean.pc1 <- mean(EUR$PC1) 
  mean.pc2 <- mean(EUR$PC2)
  mean.pc3 <- mean(EUR$PC3) 
  mean.pc4 <- mean(EUR$PC4)
  mean.pc5 <- mean(EUR$PC5) 
  mean.pc6 <- mean(EUR$PC6) 
  mean.pc7 <- mean(EUR$PC7)
  mean.pc8 <- mean(EUR$PC8) 
  mean.pc9 <- mean(EUR$PC9)
  mean.pc10 <- mean(EUR$PC10) 
  sd.pc1 <- sd(EUR$PC1)
  sd.pc2 <- sd(EUR$PC2)
  sd.pc3 <- sd(EUR$PC3)
  sd.pc4 <- sd(EUR$PC4)
  sd.pc5 <- sd(EUR$PC5)
  sd.pc6 <- sd(EUR$PC6)
  sd.pc7 <- sd(EUR$PC7)
  sd.pc8 <- sd(EUR$PC8)
  sd.pc9 <- sd(EUR$PC9)
  sd.pc10 <- sd(EUR$PC10)

#Calculate 6SD from mean for PC1:10
  x.sd.pos.pc1 <- mean.pc1 + 6*sd.pc1
  x.sd.neg.pc1 <- mean.pc1 - 6*sd.pc1
  x.sd.pos.pc2 <- mean.pc2 + 6*sd.pc2
  x.sd.neg.pc2 <- mean.pc2 - 6*sd.pc2
  x.sd.pos.pc3 <- mean.pc3 + 6*sd.pc3
  x.sd.neg.pc3 <- mean.pc3 - 6*sd.pc3
  x.sd.pos.pc4 <- mean.pc4 + 6*sd.pc4
  x.sd.neg.pc4 <- mean.pc4 - 6*sd.pc4
  x.sd.pos.pc5 <- mean.pc5 + 6*sd.pc5
  x.sd.neg.pc5 <- mean.pc5 - 6*sd.pc5
  x.sd.pos.pc6 <- mean.pc6 + 6*sd.pc6
  x.sd.neg.pc6 <- mean.pc6 - 6*sd.pc6
  x.sd.pos.pc7 <- mean.pc7 + 6*sd.pc7
  x.sd.neg.pc7 <- mean.pc7 - 6*sd.pc7
  x.sd.pos.pc8 <- mean.pc8 + 6*sd.pc8
  x.sd.neg.pc8 <- mean.pc8 - 6*sd.pc8
  x.sd.pos.pc9 <- mean.pc9 + 6*sd.pc9
  x.sd.neg.pc9 <- mean.pc9 - 6*sd.pc9
  x.sd.pos.pc10 <- mean.pc10 + 6*sd.pc10
  x.sd.neg.pc10 <- mean.pc10 - 6*sd.pc10

#Subset individuals within 6SD of mean for PC1:10 
  redux_filtered <- subset(redux, redux$PC1< x.sd.pos.pc1 & redux$PC1> x.sd.neg.pc1  
                         & redux$PC2< x.sd.pos.pc2 & redux$PC2> x.sd.neg.pc2
                         & redux$PC3< x.sd.pos.pc3 & redux$PC3> x.sd.neg.pc3
                         & redux$PC4< x.sd.pos.pc4 & redux$PC4> x.sd.neg.pc4
                         & redux$PC5< x.sd.pos.pc5 & redux$PC5> x.sd.neg.pc5
                         & redux$PC6< x.sd.pos.pc6 & redux$PC6> x.sd.neg.pc6
                         & redux$PC7< x.sd.pos.pc7 & redux$PC7> x.sd.neg.pc7
                         & redux$PC8< x.sd.pos.pc8 & redux$PC8> x.sd.neg.pc8
                         & redux$PC9< x.sd.pos.pc9 & redux$PC9> x.sd.neg.pc9
                         & redux$PC10< x.sd.pos.pc10 & redux$PC10> x.sd.neg.pc10)

  redux_filtered <- redux_filtered %>% filter (V6=="-9" | V6=="1" | V6=="2") #Select only study samples

#Export list of EUR to keep
  write.table(redux_filtered [,c(1,2)], "Individuals_to_keep_PCA.txt", sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE) #This creates a text file that keeps only columns 1 (individual ID) and 2 (family ID)

  q("no")
