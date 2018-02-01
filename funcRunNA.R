library(partykit)
# setwd(
#     "C:/Users/Afrooz/Google Drive/education/Computational Science/matching/GmatchRepo"
# )
setwd(
  "C:/Users/Afrooz/Google Drive/education/Computational Science/matching/iterOptMatchNA/iterOptMatchNA"
)
source("growNA.R")
source("growtreeNA.R")
source("partitionNA.R")
source("ordinalizeNA.R")
source("splitruleNA.R")
source("rfConstNA.R")
source("smdNA.R")
source("x2DistNA.R")
source("x2DistanceNA.R")
source("optMatchNA.R")
source("summaryGmatchNA.R")
source("iterOptNA.R")

#### Data ####

# nonOutlierLowADOS = read.csv("lowADOSNonOutlierBefore.csv", header = T)
# nonOutlierLowADOSNAVar= as.data.frame(lapply(nonOutlierLowADOS[c(-2,-1)],function(x) x[sample(c(TRUE,NA),prob=c(0.85,0.10), size=length(x), replace=TRUE)]))  
# nonOutlierLowADOS=data.frame(nonOutlierLowADOS[,c(1,2)],nonOutlierLowADOSNAVar)
data = read.csv("Abide_Combined_Full_Demographics.csv",
                na.strings = c("", "NA"))
selecteddata = data[, c(
    "SUB_ID",
    "DX_GROUP",
    "RMSD",
    "AGE_AT_SCAN",
    "ï..SITE_ID",
    "FIQ",
    "SEX",
    "HANDEDNESS_CATEGORY",
    "PERCENT_GOODTP",
    "EYE_STATUS_AT_SCAN"
)]
selecteddata = subset(
    selecteddata,
    AGE_AT_SCAN <= 18 & AGE_AT_SCAN >= 7 & PERCENT_GOODTP >= .80 & 
        EYE_STATUS_AT_SCAN == 1 & HANDEDNESS_CATEGORY<=3 & RMSD <0.16
) 

print(sapply(selecteddata, is.factor))
catVars = c("DX_GROUP","SEX", "HANDEDNESS_CATEGORY", "EYE_STATUS_AT_SCAN" )
selecteddata[,catVars] <- lapply(selecteddata[,catVars] , factor)

### Formula ####
#type names correctly#form <-group ~ RMSD.PRE.CENSORING + Age + WASI.NVIQ + Gender + Handedness 
form <-DX_GROUP ~ RMSD + AGE_AT_SCAN + FIQ +SEX + HANDEDNESS_CATEGORY 


### Debug function ####
#debug(x2Dist)
#debug(rfConst)
debug(Gmatch)

Gmatch (selecteddata, form, 10,"opt-one-to-one","chi")
#Gmatch (nonOutlierLowADOS, form, 1000,"1To3Dist")#1:3
#Gmatch (nonOutlierLowADOS, form, 1000,"propensity" )
#Gmatch (Gmatch (nonOutlierLowADOSNA, form, 1000,"propensity" ), form, 1000,"propensity" )
#Gmatch (nonOutlierLowADOS, form, 1000, "propensity")#to get propensity score
##Gmatch (nonOutlierLowADOS, form, 1000,"opt-coars-exact-rev","chi")
#Gmatch (nonOutlierLowADOS, form, 10,"opt-coars-exact-rev")
### Statistics of the dataset ####


mean(subset(optData, group==0)$Age)
sd(subset(optData, group==0)$Age)
range(subset(optData, group==0)$Age)

mean(subset(optData, group==0)$RMSD.PRE.censoring)
sd(subset(optData, group==0)$RMSD.PRE.censoring)
range(subset(optData, group==0)$RMSD.PRE.censoring)

mean(subset(optData, group==0)$WASI.NVIQ)
sd(subset(optData, group==0)$WASI.NVIQ)
range(subset(optData, group==0)$WASI.NVIQ)

table(subset(optData, group==0)$Handedness)
table(subset(dataNewPro, group==1)$Gender)
