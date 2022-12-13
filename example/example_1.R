library("ExpSigfinder")
library(tidyverse)

df <- read.table("pig_a_mutation_matrix.tsv", sep = "\t", header = T, as.is = T) 
df <- df %>%
  mutate(bg_profile=rowMeans(df[,c("VehCtrlRep3c.1","VehCtrlRep4c.1","VehCtrlRep2c.1")]))

bg_mean <- sum(df$bg_profile)
bg_min <- min(colSums(df[,c("VehCtrlRep3c.1","VehCtrlRep4c.1","VehCtrlRep2c.1")]))
# High
Wrap_KOSig(df,"bg_profile",c("EMShighRep2c.1","EMShighRep3c.1","EMShighRep4c.1"),100, sum(df$bg_profile),2,"EMShigh")
sig <- read.table("EMShigh.txt", sep = "\t", header = T, as.is = T)
plotCountbasis(sig,1,6,9,paste0("MShigh",".pdf"))
plotPercentagebasis(sig,1,6,9,paste0("MShigh","_percentage.pdf"))

# Med
Wrap_KOSig(df,"bg_profile",c("EMSmedRep4c.1","EMSmedRep2c.1","EMSmedRep3c.1"),100, sum(df$bg_profile),2,"EMSmed")
sig <- read.table("EMSmed.txt", sep = "\t", header = T, as.is = T)
plotCountbasis(sig,1,6,9,paste0("EMSmed",".pdf"))
plotPercentagebasis(sig,1,6,9,paste0("EMSmed","_percentage.pdf"))

# Low
Wrap_KOSig(df,"bg_profile",c("EMSlowRep2c.1","EMSlowRep3c.1","EMSlowRep4c.1"),100, 400,2.5,"EMSlow")
sig <- read.table("EMSlow.txt", sep = "\t", header = T, as.is = T)
plotCountbasis(sig,1,6,9,paste0("EMSlow",".pdf"))
plotPercentagebasis(sig,1,6,9,paste0("EMSlow","_percentage.pdf"))
