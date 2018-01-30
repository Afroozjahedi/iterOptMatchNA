#install.packages("partykit")
library(partykit)
setwd("/Users/mo/Documents/Thesis/Abide/Excel_Final")
source("grow.R")
source("growtree.R")
source("partition.R")
source("ordinalize.R")
source("splitrule.R")
source("rfConst.R")
source("smd.R")
source("x2Dist.R")
source("x2Distance.R")
#source("0-1Distance.R")
#install.packages("optmatch")
library(optmatch)
source("optMatch.R")
source("summaryGmatch.R")
source("iterOpt.R")

#### Data ####

data = read.csv("Abide_Combined_Full_Demographics.csv", na.strings = c("","NA"))
selecteddata = data[,c("DX_GROUP","RMSD","AGE_AT_SCAN","SITE_ID","FIQ","SEX","HANDEDNESS_CATEGORY","PERCENT_GOODTP","EYE_STATUS_AT_SCAN")]#x,"HANDEDNESS_CATEGORY","EYE_STATUS")] #1:20 is detemining rows and the c,"" is picking out coloums
selecteddata = subset(selecteddata,AGE_AT_SCAN <= 18 & AGE_AT_SCAN >= 7 & PERCENT_GOODTP >=.80 &EYE_STATUS_AT_SCAN==1) #picks out all RMSD values below the given  number
#which(data$RMSD<=.16&data$AGE_AT_SCAN<=18&data$AGE_AT_SCAN>=7&data$PERCENT_GOODTP>=.80)#shows me which participant was dropped because of RMSD values
#selecteddata$HANDEDNESS_CATEGORY = ifelse (selecteddata$HANDEDNESS_CATEGORY == "R",1,2)#


### Formula ####
form <-DX_GROUP ~ RMSD + AGE_AT_SCAN + PIQ + SEX + HANDEDNESS_CATEGORY + SITE_ID + PERCENT_GOODTP + EYE_STATUS #changes these to actual names in file

### Debug function ####
debug(Gmatch) #shows where the problem is coming from

Gmatch (selecteddata, form, 1000,"opt-one-to-one","chi")
##Gmatch (nonOutlierLowADOS, form, 1000,"opt-coars-exact-rev","chi")
#Gmatch (nonOutlierLowADOSNA, form, 10,"opt-coars-exact-rev")

