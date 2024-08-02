### Statistical analysis script
### Written by M. Veselovsky
### Last modified 2024-04-24

# clear the R environment
rm(list=ls())

#setwd("NectarExperiment_Analysis/")

#Load data and libraries
summarydat<-read.csv("processed/summarydat.csv")
TreatmentData<-read.csv("processed/TreatmentData_p.csv")#raw treatment data with all trial days
data = summarydat

library(lmerTest)
library(car)
library(effects)
library(ggplot2)
library(emmeans)
library(multcomp)
library(sjPlot)
library(sjmisc)
library(ggeffects)
library(ggbeeswarm)
library(performance) #for check_model() assumptions function



# create separate working databases for greenhouse plants, goldenrod (for GH-Field comparison), and field plants
GH_data = subset(summarydat, summarydat$ExpLoc == "GH") # data from greenhouse
GH_F_data = summarydat[which(summarydat$Plant=="SolAlt"),] #all data on goldenrod, for GH-Field comparison
Field_data = subset(summarydat,summarydat$ExpLoc == "Field") #data from field
SolAlt_F_data = subset(Field_data,Field_data$Plant == "SolAlt") #Subset SolAlt in field to check for enclosure effects

######## 0. Exploration of trial days, Descriptive plots

trialday_data = subset(TreatmentData, TreatmentData$TrialDay != "10")
trialday_data = subset(trialday_data, trialday_data$TrialDay != "3")
trialday_data = subset(trialday_data, trialday_data$TrialDay != "2")
trialday_data = subset(trialday_data, trialday_data$TrialDay != "8")
trialday_data = subset(trialday_data, trialday_data$TrialDay != "11")

scatterplot(Weight~TrialDay,data=TreatmentData,ylab="Weight (g)",xlab="Trial day")
scatterplot(Weight~EmergDate,data=TreatmentData,ylab="Weight (g)",xlab="Julian date")
scatterplot(DryFatMass~EmergDate +(1|Plant),data=GH_data,ylab="% body fat",xlab="Body mass - dry (g)",smooth=FALSE)


boxplot(Weight~TrialDay,data=trialday_data)
day_factor = trialday_data
day_factor$TrialDay = as.factor(day_factor$TrialDay)
trialday_lmer = lmer(Weight~TrialDay + ExpLoc+TrialDay:ExpLoc+ (1|ID),data=day_factor)
Anova(trialday_lmer)
plot(allEffects(trialday_lmer),ylim=c(0.3,0.65))


#######Plot the range of start dates for the different plant treatment trials
# basic plot
p=ggplot(trialday_data, aes(EmergDate, as.factor(Plant):ExpLoc,color=Plant,show.legend=FALSE)) + 
  geom_line() + geom_point()+ 
  labs(y = "Plant", x = "Start of trial (Julian date)")
#remove legend 
p_nolegend = p + theme(legend.position="none")
#Make theme of the plot black and white (apart from plant color) and remove unnecessary elements
p_nolegend + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(),
                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


########## 1. GREEHOUSE DATA ANALYSIS ##############
################## 1.1 Normality check ###############
opar <- par(mfrow = c(1, 1))

#Check for normality of the butterflies raw weights on day 7 of their respective trials
weight_distr <- ggplot(GH_data, aes(x = RawWeight_day7))
weight_distr <- weight_distr +
  # add data density smooth
  geom_density() +
  # add rug (bars at the bottom of the plot)
  geom_rug() +
  # add black semitransparent histogram
  geom_histogram(aes(y = ..density..),
                 color = "black",
                 alpha = 0.3) +
  # add normal curve in red, with mean and sd from fklength
  stat_function(fun = dnorm,
                args = list(
                  mean = mean(GH_data$RawWeight_day7),
                  sd = sd(GH_data$RawWeight_day7)
                ),
                color = "red")

#display graph
weight_distr
qqnorm(GH_data$RawWeight_day7)
qqline(GH_data$RawWeight_day7)


#Check for normality of fat weight in butterflies after trial finished
fat_distr <- ggplot(GH_data, aes(x = DryFatMass))
fat_distr <- fat_distr +
  # add data density smooth
  geom_density() +
  # add rug (bars at the bottom of the plot)
  geom_rug() +
  # add black semitransparent histogram
  geom_histogram(aes(y = ..density..),
                 color = "black",
                 alpha = 0.3) +
  # add normal curve in red, with mean and sd from DryFatMass
  stat_function(fun = dnorm,
                args = list(
                  mean = mean(GH_data$DryFatMass),
                  sd = sd(GH_data$DryFatMass)
                ),
                color = "red")

#display graph
fat_distr
qqnorm(GH_data$DryFatMass)
qqline(GH_data$DryFatMass)


# Formal test for normality (Shapiro)
shapiro.test(GH_data$RawWeight_day7) #weights normal
shapiro.test(GH_data$DryFatMass) #fat non-normal

# boxplots of raw weight on day 7 by plant species in the greenhouse
myplot<-ggplot(data=GH_data, aes(x=Plant, y=RawWeight_day7),na.action=na.exclude)
myplot+geom_boxplot(notch=FALSE)

#########plots for sex differences
myplot<-ggplot(data=GH_data, aes(x=Sex, y=DryFatMass),na.action=na.exclude)
myplot+geom_boxplot(notch=FALSE)


myplot<-ggplot(data=GH_data, aes(x=Sex, y=WaterMass),na.action=na.exclude)
myplot+geom_boxplot(notch=FALSE)

myplot<-ggplot(data=GH_data, aes(x=Sex, y=DryLeanMass),na.action=na.exclude)
myplot+geom_boxplot(notch=FALSE)

myplot<-ggplot(data=GH_data, aes(x=Sex, y=ForewingLength),na.action=na.exclude)
myplot+geom_boxplot(notch=FALSE)


#################### 1.2 MODEL BUILDING - LMER Weight ###########################
gh_TreatmentData = subset(TreatmentData, TreatmentData$ExpLoc == "GH")

overall_weight = lmer(Weight~Plant+OrigWeight+Sex+TrialDay+EnclCol+(1|ID),data=gh_TreatmentData)
#removed dateweighed (collinearity with trialday) and cohort (i don't think it is relevant)
weight_effects=allEffects(overall_weight)


check_model(overall_weight)
Anova(overall_weight)
plot(weight_effects)
plot(weight_effects$Plant)


# gh_TreatmentData$EmergDate_c = gh_TreatmentData$EmergDate-mean(gh_TreatmentData$EmergDate)
# gh_TreatmentData$TrialDay_f = as.factor(gh_TreatmentData$TrialDay)

# overall_weight_sc = lmer(Weight~Plant + EmergDate+Sex+TrialDay+EnclCol+Cohort+(1|ID),data=gh_TreatmentData)
# weight_effects = allEffects(overall_weight_sc)
# 
# Anova(overall_weight_sc, type="3")
plot(weight_effects)

effects_weight = allEffects(overall_weight)
plot(effects_weight)
 summary(overall_weight)
check_model(overall_weight)

emmeans(overall_weight, list(pairwise~Plant), adjust="tukey")
# Check correlation of starting weight and forewing length (to see if I need both in the model)
# The a-priori expectation is that they will be highly correlated, in which case I think forewing length is
# the ideal variable to retain as it does not change with adult age or time of day (feeding could change weight slightly).
# Weights could fluctuate depending on the time they were taken
# since emergence (metabolism as they have not eaten, and drying time for the fluid they expel
# as they eclose from their chrysalises)

# fwl_day0w_lm = lm(RawWeight_day0 ~ ForewingLength, data=GH_data)
# Anova(fwl_day0w_lm, type=3)# Significantly related, so I will keep forewing length and not use weight
# plot(allEffects(fwl_day0w_lm))
# summary(fwl_day0w_lm) 

# # Estimate the change in raw weight on day 7 of trial based on predictor variables
# weight_lm = lm(RawWeight_day7~ Plant + Sex + ForewingLength + EmergDate + TotalSA, data=GH_data)
# 
# 
# #Test for collinearity among predictor variables (GVIF<4 for continuous acceptable,
# # GVIF^(1/(2*df)) < 2 acceptable for categorical)
# vif(weight_lm) #all good
# opar <- par(mfrow = c(2, 2))
# 
# plot(weight_lm)
# plot(weight_lm, resid(., scaled=TRUE) ~ fitted(.), abline = 0,pch=16,xlab="Fitted values",ylab="Standardised residuals")
# Anova(weight_lm, type=3)
# summary(weight_lm)
# 
# effects_info = allEffects(weight_lm)
# effects_info
# plot(effects_info)
# Anova(weight_lm,type=3)

############ 1.3 Multiple comparisons test for Weights
# 
# emmeans(weight_lm, list(pairwise~Plant), adjust="tukey")

######## 2. COMPARISON OF GH TO FIELD DATA (USING GOLDENROD) ##############
#Create a database of field and gh goldenrod to compare experiments
GH_F_data = summarydat[which(summarydat$Plant=="SolAlt"),]
summary(GH_F_data)

#Check for normality on raw day 7 weights
hist(GH_F_data$RawWeight_day7)


### linear model to compare GH method to field method (ExpLoc)
gh_f_lm = lm(RawWeight_day7~ Sex + ForewingLength + ExpLoc, data=GH_F_data)
Anova(gh_f_lm, type=3)

plot(allEffects(gh_f_lm))
summary(gh_f_lm)


################ 3. PLANT SURFACE AREA DIFFERENCES BETWEEN SPECIES ############

########### 3.1 SA - GREENHOUSE PLANT COMPARISONS
## Histogram of surface area by plant species
myplot<-ggplot(data=GH_data, aes(x=Plant, y=log(TotalSA)),na.action=na.exclude)
myplot+geom_boxplot(notch=TRUE)

#### Compare surface area of the greenhouse plants
## scatter plot by plant species
ggplot(data = GH_data) +
  aes(y = Plant, x = log(TotalSA)) +
  geom_beeswarm() +
  coord_flip()
boxplot(TotalSA~Plant,data=GH_data)


########### 3.2 SA - FIELD SPECIES
ggplot(data = Field_data) +
  aes(y = Plant, x = TotalSA) +
  geom_beeswarm() +
  coord_flip()

ggplot(data = Field_data) +
  aes(y = Plant, x = log(TotalSA)) +
  geom_beeswarm() +
  coord_flip()

myplot<-ggplot(data=Field_data, aes(x=Plant, y=TotalSA),na.action=na.exclude)
myplot+geom_boxplot(notch=TRUE)
# 
# SA_field_lm = lm(TotalSA~Plant,data=Field_data)
# Anova(SA_field_lm)
# plot(allEffects(SA_field_lm))


###################### 4 FIELD PLANTS & COMPARISONS #####################

# Because greenhouse SolAlt (Tall goldenrod) had a different effect from field SolAlt,
# I will specifically compare the two field species (EutMac/Joe-pye weed & SolAlt)
# as they cannot be pooled with the greenhouse data
summary(Field_data)


##################### 4.1 CHECK FOR FIELD ENCLOSURE EFFECTS ####################
# Check for enclosure effects with just SolAlt data (SolAlt had two different types
# of enclosures in the trial)
encl_lm = lm(RawWeight_day7~ EnclCol, data=SolAlt_F_data)

Anova(encl_lm) #EnclCol least significant, remove and compare
summary(encl_lm)

#################### 4.2 COMPARE FIELD SPECIES ########################
# I had a limited number of field plots (only 3 plots of JPW)
# so the first set of field trials had 2 females per enclosure for EutMac,
# and a mix of double females, single males, and single females for SolAlt.
# For the second cohort, I put additional double females and males for both species
# and added a single female and single male for EutMac. Unfortunately, I did not get
# more male replicates for EutMac so I will exclude males from the field analysis
females_data = subset(Field_data,Field_data$Sex == "F")

#### build model using females, without EnclCol
females_data$NumButterflies = as.factor(females_data$NumButterflies)
field_lm = lm(RawWeight_day7~ Plant + ForewingLength + NumButterflies + TotalSA, data=females_data)
Anova(field_lm, type=3)
summary(field_lm)
plot(allEffects(field_lm))

field_fat = lm(DryFatMass ~ Plant + ForewingLength + NumButterflies + TotalSA, data=females_data)

Anova(field_fat)
summary(field_fat)

############### 5. BODY COMPOSITION ANALYSES ##########################

############### 5.1 Exploration ##################################
opar <- par(mfrow = c(1, 1))
boxplot(WaterMass~Plant,data=GH_data,ylab="Body water content (g)")
water_lm = lm(WaterMass~Plant,data=GH_data)
summary(water_lm)

ggplot(data, aes(x = Plant, y = WaterMass, color = ExpLoc)) +  # ggplot function
  geom_point()

scatterplot(DryFatMass~RawWeight_day7,data=data, ylab="Body fat (g)",xlab="Body mass - wet (g)")
fat_lm =lm(DryFatMass~RawWeight_day7,data=data)
summary(fat_lm)

ggplot(data, aes(x = RawWeight_day7, y = log(DryFatMass), color = ExpLoc)) +  # ggplot function
  geom_point()+
  geom_smooth(method=lm , se=TRUE) +

scatterplot(RelDryFat~RawWeight_day7,data=data,ylab="% body fat",xlab="Body mass - dry (g)")
scatterplot(RelDryFat~DryMass,data=data,ylab="% body fat",xlab="Body mass - dry (g)")

ggplot(data, aes(x = Plant, y = RelDryFat, color = ExpLoc)) +  # ggplot function
  geom_point()+
  geom_boxplot(notch=FALSE)


relfat_lm = lm(RelDryFat~RawWeight_day7+Plant,data=data)
summary(relfat_lm)

scatterplot(DryMass~RawWeight_day7,data=data,ylab="Dry body mass (g)",xlab="Body mass - wet (g)")
drymass_lm = lm(DryMass~RawWeight_day7+Plant,data=data)
summary(drymass_lm)

scatterplot(DryFatMass~DryLeanMass,data=data,ylab="Body fat (g)",xlab="Lean mass (g)")

scatterplot(DryMass~EmergDate,data=data)
scatterplot(DryFatMass~EmergDate,data=data)


########################## 5.2 GREENHOUSE FAT MODELS ####################################
gh_fat_lm = lm(DryFatMass~ Plant + Sex + ForewingLength + EmergDate + TotalSA, data=GH_data)
gh_fat_lm_log = lm(log(DryFatMass)~ Plant + Sex + ForewingLength + EmergDate + TotalSA, data=GH_data)

#Test for collinearity among predictor variables (GVIF<4 for continuous acceptable,
# GVIF^(1/(2*df)) < 2 acceptable for categorical)
vif(gh_fat_lm)

opar <- par(mfrow = c(2, 2))
plot(gh_fat_lm)
plot(gh_fat_lm_log)

check_model(gh_fat_lm)
check_model(gh_fat_lm_log)

# Anova(gh_fat_lm)
Anova(gh_fat_lm_log)
summary(gh_fat_lm_log)

# pred_gh_fat_lm=allEffects(gh_fat_lm)
pred_gh_fat_log=allEffects(gh_fat_lm_log)


# plot(pred_gh_fat_lm)
plot(pred_gh_fat_log)
# Anova(gh_fat_lm_log)

predict_plant = predictorEffect("Plant",gh_fat_lm_log)

str(predict_plant)
predict_plant$fit = exp(predict_plant$fit)
predict_plant$upper=exp(predict_plant$upper)
predict_plant$lower=exp(predict_plant$lower)
plot(predict_plant,ylab ="Fat mass (g)")
plot(allEffects(gh_fat_lm_log))

Anova(gh_fat_lm_log)
emmeans(gh_fat_lm_log, list(pairwise~Plant), adjust="tukey")


################### GLM Using Gamma distribution for Fat #######################
# gh_fat_lm_gm = glm(DryFatMass ~ Sex + Plant + ForewingLength + EmergDate + TotalSA, data=GH_data,family="Gamma")
# opar <- par(mfrow = c(2, 2))
# plot(gh_fat_lm_gm)
# 
# Anova(gh_fat_lm_gm)
# summary(gh_fat_lm_gm)
# pred_gh_fat_gm=allEffects(gh_fat_lm_gm)
# plot(pred_gh_fat_gm)
# 
# str(pred_gh_fat_gm$Plant)
# str(pred_gh_fat_gm$Plant$response)
# str(pred_gh_fat_gm$fit)
# str(pred_gh_fat_gm)
# 
# summary(gh_fat_lm_gm)
# str(pred_gh_fat_gm)
# plot(pred_gh_fat_gm)
# check_model(gh_fat_lm_gm)
# 
# plot(pred_gh_fat_gm)
# plot(pred_gh_fat_log)
# emmeans(gh_fat_lm_gm, list(pairwise~Plant), adjust="tukey")

# ## Plot of effects using ggplot glm
# ggplot(GH_data, aes(EmergDate, DryFatMass)) + 
#   geom_smooth(method = "glm", formula = DryFatMass~ Plant + Sex + ForewingLength + EmergDate + TotalSA, colour = "black",
#               linetype = 2, fill = "gray80", alpha = 0.2,
#               method.args = list(family = gamma)) +
#   geom_rug(sides = "b") +
#   theme_bw() +
#   labs(y = "Fat weight (g)", x = "Plant") +
#   theme(text = element_text(size = 16),
#         plot.margin = margin(50, 50, 50, 50),
#         axis.title.x = element_text(vjust = -8),
#         axis.title.y = element_text(vjust = 10),
#         plot.title = element_text(vjust = 8))

########## relative fat content (%)
gh_relfat_lm = lm(RelDryFat~ Plant + EmergDate + TotalSA, data=GH_data)
gh_relfat_log = lm(log(RelDryFat)~ Plant + EmergDate + TotalSA, data=GH_data)
vif(gh_relfat_lm)
Anova(gh_relfat_lm)

summary(gh_relfat_lm)
check_model(gh_relfat_lm)
check_model(gh_relfat_log)
emmeans(gh_relfat_lm, list(pairwise~Plant), adjust="tukey")



######################## 5.3 FIELD FAT ANALYSIS #######################
opar <- par(mfrow = c(2, 2))
field_fat_lm = lm(DryFatMass~ EmergDate + TotalSA + Plant, data=Field_data)
Anova(field_fat_lm)
plot(allEffects(field_fat_lm))
vif(field_fat_lm)

opar <- par(mfrow = c(2, 2))
plot(field_fat_lm)

field_fat_glm = glm(DryFatMass~ EmergDate, data=Field_data,family="Gamma")
Anova(field_fat_glm)
check_model(field_fat_glm)

plot(allEffects(field_fat_glm))



#################### 6. LEAN MASS (PROTEIN) ANALYSIS ###############

lean_distr <- ggplot(GH_data, aes(x = DryLeanMass))
lean_distr <- lean_distr +
  # add data density smooth
  geom_density() +
  # add rug (bars at the bottom of the plot)
  geom_rug() +
  # add black semitransparent histogram
  geom_histogram(aes(y = ..density..),
                 color = "black",
                 alpha = 0.3) +
  # add normal curve in red, with mean and sd from fklength
  stat_function(fun = dnorm,
                args = list(
                  mean = mean(GH_data$DryLeanMass),
                  sd = sd(GH_data$DryLeanMass)
                ),
                color = "red")

#display graph
lean_distr
#generate the QQ plot
qqnorm(GH_data$DryLeanMass)
qqline(GH_data$DryLeanMass)

# 
# ggplot(data, aes(x = RawWeight_day7, y = DryLeanMass, color = ExpLoc)) +  # ggplot function
#   geom_point()+
#   geom_smooth(method=lm , color="red", se=TRUE)

gh_prot_lm = lm(DryLeanMass ~ Sex + Plant + ForewingLength + EmergDate, data=GH_data)
opar <- par(mfrow = c(2, 2))
plot(gh_prot_lm)

check_model(gh_prot_lm)
Anova(gh_prot_lm)

pred_prot_lm = allEffects(gh_prot_lm)
plot(pred_prot_lm)
emmeans(gh_prot_lm, list(pairwise~Plant), adjust="tukey")


##################### 7. WATER ############################

water_distr <- ggplot(GH_data, aes(x = WaterMass))
water_distr <- water_distr +
  # add data density smooth
  geom_density() +
  # add rug (bars at the bottom of the plot)
  geom_rug() +
  # add black semitransparent histogram
  geom_histogram(aes(y = ..density..),
                 color = "black",
                 alpha = 0.3) +
  # add normal curve in red, with mean and sd from fklength
  stat_function(fun = dnorm,
                args = list(
                  mean = mean(GH_data$WaterMass),
                  sd = sd(GH_data$WaterMass)
                ),
                color = "red")

#display graph
water_distr
qqnorm(GH_data$DryLeanMass)
qqline(GH_data$DryLeanMass)

gh_water_lm = lm(WaterMass ~ Sex + Plant + ForewingLength + EmergDate, data=GH_data)

opar <- par(mfrow = c(2, 2))

plot(gh_water_lm)
check_model(gh_water_lm)
water_effects = allEffects(gh_water_lm)

plot(water_effects)
Anova(gh_water_lm)

emmeans(gh_water_lm, list(pairwise~Plant), adjust="tukey")


