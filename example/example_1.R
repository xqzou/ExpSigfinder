library("ExpSigfinder")
library(tidyverse)

df <- read.table("example.tsv", sep = "\t", header = T, as.is = T) 
df <- df %>%
  mutate(bg_profile=rowMeans(df[,c("treatment.1","treatment.2","treatment.3")]))

# calculate the average number of mutation burden in controls
bg_mean <- sum(df$bg_profile)


Wrap_KOSig(df,"bg_profile",c("treatment.1","treatment.2","treatment.3"),100, bg_mean,2,"treatment")
sig <- read.table("treatment.txt", sep = "\t", header = T, as.is = T)
plotCountbasis(sig,1,6,9,paste0("treatment",".pdf"))
plotPercentagebasis(sig,1,6,9,paste0("treatment","_percentage.pdf"))