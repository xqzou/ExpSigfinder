# ExpSigfinder
Find new mutational signatures from experiments

## Installation

```{r, eval = FALSE}
# Install the released version from Github
git clone https://github.com/xqzou/ExpSigfinder.git
cd ExpSigfinder
R CMD INSTALL .
```


## Example
Apply MMRDetect to 26 breast cancers

```{r, eval = FALSE}
# import the packages 
library("ExpSigfinder")
library(tidyverse)

# calculate the average number of mutation burden in controls
bg_mean <- sum(df$bg_profile)

# High
Wrap_KOSig(df,"bg_profile",c("EMShighRep2c.1","EMShighRep3c.1","EMShighRep4c.1"),100, bg_mean,2,"EMShigh")
sig <- read.table("EMShigh.txt", sep = "\t", header = T, as.is = T)
plotCountbasis(sig,1,6,9,paste0("MShigh",".pdf"))
plotPercentagebasis(sig,1,6,9,paste0("MShigh","_percentage.pdf"))

# Med
Wrap_KOSig(df,"bg_profile",c("EMSmedRep4c.1","EMSmedRep2c.1","EMSmedRep3c.1"),100, bg_mean*0.8,2,"EMSmed")
sig <- read.table("EMSmed.txt", sep = "\t", header = T, as.is = T)
plotCountbasis(sig,1,6,9,paste0("EMSmed",".pdf"))
plotPercentagebasis(sig,1,6,9,paste0("EMSmed","_percentage.pdf"))

# Low
Wrap_KOSig(df,"bg_profile",c("EMSlowRep2c.1","EMSlowRep3c.1","EMSlowRep4c.1"),100, bg_mean*0.8,3,"EMSlow")
sig <- read.table("EMSlow.txt", sep = "\t", header = T, as.is = T)
plotCountbasis(sig,1,6,9,paste0("EMSlow",".pdf"))
plotPercentagebasis(sig,1,6,9,paste0("EMSlow","_percentage.pdf"))
```
