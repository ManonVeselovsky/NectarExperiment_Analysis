### Script to organize and process raw data from 2022 monarch experiment
### Written by D. Rivest and M. Veselovsky
### Last updated 2023-01-12

# clear the R environment
rm(list=ls())

#### Set the working directory to the location of the project file "NectarExperiment_Analysis"
#setwd("NectarExperiment_Analysis/")
setwd("C:/Users/Manon/Documents/MSc Biology/Master's Project/Statistics/NectarExperiment_Analysis")
project_directory = getwd()
# Load data
TreatmentData<-read.csv("data/TreatmentData.csv")
ForewingMeasurements = read.csv("data/ForewingMeasurements.csv")
FloralMeasurements = read.csv("data/FloralMeasurements.csv")
library(dplyr)

############ RAW WEIGHTS ##############################
# Start with a dataframe including only the raw weight at day 0
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


############ FOREWING MEASUREMENTS ####################3

forewingdat = subset(ForewingMeasurements, select=c("ID", "ForewingLength", "ForewingDamage"))
experimentdat<-merge(experimentdat, forewingdat[,c("ID", "ForewingLength", "ForewingDamage")], by.x=c("ID"), all=TRUE)

# for butterflies that were weighed while dead, remove those values (i.e., VH23)
experimentdat = subset(experimentdat, experimentdat$ID != "VH23")

# remove individuals that do not have a day 0 weight
experimentdat = subset(experimentdat, experimentdat$ID != "VH24")
experimentdat = subset(experimentdat, experimentdat$ID != "VH30")

# remove individuals that do not have a day 5 weight
experimentdat = subset(experimentdat, experimentdat$ID != "Unknown1")

# remove individuals that do not have a day 7 weight
experimentdat = subset(experimentdat, experimentdat$ID != "Jacob8")
experimentdat = subset(experimentdat, experimentdat$ID != "TWA74")
experimentdat = subset(experimentdat, experimentdat$ID != "Unknown1")
experimentdat = subset(experimentdat, experimentdat$ID != "VH24")



############ RAW WEIGHT GAIN
#add a column for raw weight gain on day 0 (i.e., 0!), which is just the difference
experimentdat$RawWeightGain_day0 <- 0

#Repeat for day1
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

#calculate total SA and add into experimentdat
experimentdat$TotalSA = foo7$FlowerHeads*foo7$SurfaceArea


############ FINAL CLEAN-UP #####################


summarydat = experimentdat


############ SAVE THE PROCESSED FILE
setwd("processed/")
write.csv(summarydat, file= "summarydat.csv")
setwd(project_directory)

