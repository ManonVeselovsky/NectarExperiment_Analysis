### Script to organize and process raw data from 2022 monarch experiment
### Written by D. Rivest and M. Veselovsky
### Last updated 2023-08-14

# clear the R environment
rm(list=ls())

#### Set the working directory to the location of the project file "NectarExperiment_Analysis"
#setwd("NectarExperiment_Analysis/")
setwd("C:/Users/Manon/Documents/MSc Biology/Master's Project/Statistics/NectarExperiment_Analysis")
project_directory = getwd()
# Load data
TreatmentData<-read.csv("data/TreatmentData.csv")
treatmentdat = TreatmentData
ForewingMeasurements = read.csv("data/ForewingMeasurements.csv")
FloralMeasurements = read.csv("data/FloralMeasurements.csv")
FatData = read.csv("data/FatExtractions.csv")
temp_data = read.csv("data/TemperatureData.csv")
library(dplyr)
library(tidyr) #to separate date from time in tempdat

############ RAW WEIGHTS ##############################

# Start with a dataframe that only contains the columns of ID, Trial Day, and Weight
experimentdat<-subset(TreatmentData, select=c("ID", "TrialDay", "Weight"))

# Keep only rows for weights measured on day 0 (i.e., original weight)
experimentdat<-experimentdat[which(experimentdat$TrialDay=="0"),]

#rename the Weight column
colnames(experimentdat)[colnames(experimentdat)=="Weight"] <- "RawWeight_day0"

#remove TrialDay Column
experimentdat$TrialDay<-NULL

# add column for day 1

foo1<-subset(TreatmentData, select=c("ID", "TrialDay", "Weight"))
foo1<-foo1[which(foo1$TrialDay=="1"),]
colnames(foo1)[colnames(foo1)=="Weight"] <- "RawWeight_day1"
experimentdat<-merge(experimentdat, foo1[,c("ID","RawWeight_day1")], by.x=c("ID"), all=TRUE)

#add column for raw weights for day 5
foo1<-subset(TreatmentData, select=c("ID", "TrialDay", "Weight"))
foo1<-foo1[which(foo1$TrialDay=="5"),]
colnames(foo1)[colnames(foo1)=="Weight"] <- "RawWeight_day5"
experimentdat<-merge(experimentdat, foo1[,c("ID","RawWeight_day5")], by.x=c("ID"), all=TRUE)

#add column for raw weights for day 7
foo1<-subset(TreatmentData, select=c("ID", "TrialDay", "Weight"))
foo1<-foo1[which(foo1$TrialDay=="7"),]
colnames(foo1)[colnames(foo1)=="Weight"] <- "RawWeight_day7"
experimentdat<-merge(experimentdat, foo1[,c("ID","RawWeight_day7")], by.x=c("ID"), all=TRUE)

#add column for raw weights for day 10
foo1<-subset(TreatmentData, select=c("ID", "TrialDay", "Weight"))
foo1<-foo1[which(foo1$TrialDay=="10"),]
colnames(foo1)[colnames(foo1)=="Weight"] <- "RawWeight_day10"
experimentdat<-merge(experimentdat, foo1[,c("ID","RawWeight_day10")], by.x=c("ID"), all=TRUE)

############ ADD COLUMNS FOR TREATMENT DATA ######################
#Add additional variables from the experimental treatments (Couple, cohort, plant species, experimental location, number of butterflies in the
# enclosure, enclosure colour, sex of butterfly)
foo2<-subset(TreatmentData, select=c("ID", "TrialDay", "Couple", "Cohort", "Plant", "FlowerHeads","ExpLoc", "NumButterflies","EnclCol", "Sex", "EmergDate"))
foo2<-foo2[which(foo2$TrialDay=="0"),]
experimentdat<-merge(experimentdat, foo2[,c("ID", "Couple", "Cohort", "Plant","ExpLoc","FlowerHeads", "NumButterflies","EnclCol", "Sex", "EmergDate")], by.x=c("ID"), all=TRUE)


experimentdat_83 = experimentdat
############ FOREWING MEASUREMENTS ###################

forewingdat = subset(ForewingMeasurements, select=c("ID", "ForewingLength", "ForewingDamage"))
experimentdat<-merge(experimentdat, forewingdat[,c("ID", "ForewingLength", "ForewingDamage")], by.x=c("ID"), all=TRUE)
treatmentdat<-merge(treatmentdat, forewingdat[,c("ID", "ForewingLength", "ForewingDamage")], by.x=c("ID"), all=TRUE)

# for butterflies that were weighed while dead (see notes on raw datasheet), remove those values
experimentdat = subset(experimentdat, experimentdat$ID != "VH23")

# # remove individuals that do not have a day 0 weight
# experimentdat = subset(experimentdat, experimentdat$ID != "VH24")
# experimentdat = subset(experimentdat, experimentdat$ID != "VH30")

# # remove individuals that do not have a day 5 weight
# experimentdat = subset(experimentdat, experimentdat$ID != "Unknown1")

# remove individuals that do not have a day 7 weight
experimentdat = subset(experimentdat, experimentdat$ID != "Jacob8")
experimentdat = subset(experimentdat, experimentdat$ID != "TWA74")
experimentdat = subset(experimentdat, experimentdat$ID != "Unknown1")
experimentdat = subset(experimentdat, experimentdat$ID != "VH24")


# remove individuals that do not have the number of flowering heads recorded
experimentdat = subset(experimentdat, experimentdat$ID != "VH23")

expdat_78 = experimentdat
############ RAW WEIGHT GAIN

#Raw weight gain day 1 calculation (Day 1 weight minus day 0 weight)
experimentdat$RawWeightGain_day1 <- experimentdat$RawWeight_day1-experimentdat$RawWeight_day0

#repeat for day 5
experimentdat$RawWeightGain_day5 <- experimentdat$RawWeight_day5-experimentdat$RawWeight_day0

#repeat for day 7
experimentdat$RawWeightGain_day7 <- experimentdat$RawWeight_day7-experimentdat$RawWeight_day0

#repeat for day 10
experimentdat$RawWeightGain_day10 <- experimentdat$RawWeight_day10-experimentdat$RawWeight_day0

############ RELATIVE WEIGHT GAIN #############
#add a column for relative weight gain on day 0 (i.e., 0!), which is a proportional change
experimentdat$RelWeightGain_day0 <- 0

# repeat for day 1
experimentdat$RelWeightGain_day1<-experimentdat$RawWeightGain_day1/experimentdat$RawWeight_day0


#repeat for day 5
experimentdat$RelWeightGain_day5<-experimentdat$RawWeightGain_day5/experimentdat$RawWeight_day0

#repeat for day 7
experimentdat$RelWeightGain_day7<-experimentdat$RawWeightGain_day7/experimentdat$RawWeight_day0

#repeat for day 10
experimentdat$RelWeightGain_day10<-experimentdat$RawWeightGain_day10/experimentdat$RawWeight_day0


############ ADD FLORAL AREA DATA ################

floraldat = FloralMeasurements

floraldat$Radius = floraldat$Diameter/2

# Calculate floral resource surface area for a single floral unit with height (as bottomless cylinders)
cylFlowers = floraldat[which(floraldat$AreaType=="CircHeight"),]

#Calculate floral area for a single floral unit with height (bottomless cylinder) (2*pi*r*h + pi*r^2)
cylFlowers$SurfaceArea = 2*pi*cylFlowers$Radius*cylFlowers$Height + pi*cylFlowers$Radius^2


# Calculate floral resource surface area a single floral unit without height (simple circles)
circFlowers = floraldat[which(floraldat$AreaType=="Circular"),]

# Circular flower surface area (pi r^2)
circFlowers$SurfaceArea = 2*pi*circFlowers$Radius^2

# Combine cylindrical and simple circle surface area flowers into single temp database
foo5 = rbind(cylFlowers, circFlowers)
#Reduce the number of columns to only those needed (plant sp., location, flower area type, surface area)
foo6 = subset(foo5, select=c("Species", "ExpLoc","AreaType", "SurfaceArea"))

# get the mean values of surface area by species in new temp database:
foo6 = foo6 %>% group_by(Species,ExpLoc,AreaType) %>% summarise_each(funs(mean))

# Rename plant species column from "Species" to "Plant" to match up with my treatment data
colnames(foo6)[colnames(foo6) == "Species"] ="Plant"


# Merge average surface area into temp database of master sheet
foo7 <- merge(experimentdat,foo6[,c("SurfaceArea","ExpLoc","Plant")],by=c("Plant","ExpLoc"))

#Calculate the total SA by multiplying the number of flowering heads by the average SA of a flowering head for that species
foo7$TotalSA = foo7$FlowerHeads*foo7$SurfaceArea


foo8 = merge(treatmentdat,foo6[,c("SurfaceArea","ExpLoc","Plant")],by=c("Plant","ExpLoc"))
foo8$TotalSA = foo8$FlowerHeads*foo8$SurfaceArea

# Merge the TotaSA column into the processed database
experimentdat <- merge(experimentdat,foo7[,c("TotalSA","ID")],by=c("ID"),all.x = TRUE)
treatmentdat = foo8

################### FAT DATA #################
fatdat = FatData
fatdat <- fatdat %>% 
  rename("DryMass" = "Dry.Mass","DryLeanMass"="Dry.Lean.Mass","WaterMass"="Water.Mass")
experimentdat <- merge(experimentdat,fatdat[,c("DryMass","DryLeanMass", "RelDryFat", "DryFatMass","WaterMass","ID")],by.x=c("ID"),all.x=TRUE)


###### Re-order the data so that plants are in order of most visited to least visited
# #Subset out the individual plants - order to be solalt, buddav, symeri, ********CHECK eutmac, echpur, rudhir, helhel (most to least visited)
# soldat = subset(experimentdat, Plant == "SolAlt")
# buddat = subset(experimentdat, Plant =="BudDav")
# symdat = subset(experimentdat, Plant =="SymEri")
# eutdat = subset(experimentdat, Plant =="EutMac")
# echdat = subset(experimentdat, Plant =="EchPur")
# ruddat = subset(experimentdat, Plant =="RudHir")
# heldat = subset(experimentdat, Plant =="HelHel")
# 
# foo9 = rbind(soldat, buddat,symdat,eutdat, echdat,ruddat,heldat)
# experimentdat_o = foo9
# 
# soldat = subset(treatmentdat, Plant == "SolAlt")
# buddat = subset(treatmentdat, Plant =="BudDav")
# symdat = subset(treatmentdat, Plant =="SymEri")
# eutdat = subset(treatmentdat, Plant =="EutMac")
# echdat = subset(treatmentdat, Plant =="EchPur")
# ruddat = subset(treatmentdat, Plant =="RudHir")
# heldat = subset(treatmentdat, Plant =="HelHel")
# 
# foo10 = rbind(soldat, buddat,symdat,eutdat, echdat,ruddat,heldat)
# treatmentdat_o = foo10
# 
# 

#Rename plants - add column alphabetical order prefix in the order of most visited to least visited
newdat = treatmentdat %>%
  mutate(AlphPlant = case_when(
    Plant %in% c("SolAlt") ~ "1_SolAlt",
    Plant %in% c("BudDav") ~ "2_BudDav",
    Plant %in% c("SymEri") ~ "3_SymEri",
    Plant %in% c("EchPur") ~ "4_EchPur",
    Plant %in% c("RudHir") ~ "5_RudHir",
    Plant %in% c("HelHel") ~ "6_HelHel",
    Plant %in% c("EutMac") ~ "7_EutMac",
    TRUE ~ Plant
  ))

treatmentdat=newdat

newdat2 = experimentdat %>%
  mutate(AlphPlant = case_when(
    Plant %in% c("SolAlt") ~ "1_SolAlt",
    Plant %in% c("BudDav") ~ "2_BudDav",
    Plant %in% c("SymEri") ~ "3_SymEri",
    Plant %in% c("EchPur") ~ "4_EchPur",
    Plant %in% c("RudHir") ~ "5_RudHir",
    Plant %in% c("HelHel") ~ "6_HelHel",
    Plant %in% c("EutMac") ~ "7_EutMac",
    TRUE ~ Plant
  ))

experimentdat=newdat2

######### Temperature data - remove unnecessary columns, rename columns #############
summary(temp_data)
foo9<-subset(temp_data, select=c(1:7))

names(foo9)[5] = "GMT"
names(foo9)[6] = "Temperature"
names(foo9)[7] = "Light.Lux"

# separate date and time into different columns
summary(foo9$Date.Time.GMT)
library(stringr)
foo9[c('Date', 'Time','AM.PM')] <- str_split_fixed(foo9$Date.Time.GMT, ' ')

tempdat = foo9
############ FINAL CLEAN-UP #####################


summarydat = experimentdat
TreatmentData_p = treatmentdat

############ SAVE THE PROCESSED FILE
setwd("processed/")
write.csv(summarydat, file= "summarydat.csv")
write.csv(TreatmentData_p,file="TreatmentData_p.csv")
write.csv(tempdat, file="tempdat.csv")
setwd(project_directory)



