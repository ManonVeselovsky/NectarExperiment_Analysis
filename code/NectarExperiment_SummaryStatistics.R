### Statistical analysis script
### Written by M. Veselovsky
### Last modified 2023-01-13

# clear the R environment
rm(list=ls())

#setwd("NectarExperiment_Analysis/")

#Load data and libraries
summarydat<-read.csv("processed/summarydat.csv")
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
  # add normal curve in red, with mean and sd from fklength
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



#################### 1.2 MODEL BUILDING - LMER ###########################

# Check correlation of starting weight and forewing length (to see if I need both in the model)
# The a-priori expectation is that they will be highly correlated, in which case I think forewing length is
# the ideal variable to retain as it does not change with adult age or time of day (feeding could change weight slightly).
# Weights could fluctuate depending on the time they were taken
# since emergence (metabolism as they have not eaten, and drying time for the fluid they expel
# as they eclose from their chrysalises)

fwl_day0w_lm = lm(RawWeight_day0 ~ ForewingLength, data=GH_data)
Anova(fwl_day0w_lm, type=3)# Significantly related, so I will keep forewing length and not use weight
plot(allEffects(fwl_day0w_lm))
summary(fwl_day0w_lm) 

# Estimate the change in raw weight on day 7 of trial based on predictor variables
weight_lm = lm(RawWeight_day7~ Plant + Sex + ForewingLength + EmergDate + TotalSA, data=GH_data)


#Test for collinearity among predictor variables (GVIF<4 for continuous acceptable,
# GVIF^(1/(2*df)) < 2 acceptable for categorical)
vif(weight_lm) #all good
opar <- par(mfrow = c(2, 2))

plot(weight_lm)
plot(weight_lm, resid(., scaled=TRUE) ~ fitted(.), abline = 0,pch=16,xlab="Fitted values",ylab="Standardised residuals")
Anova(weight_lm, type=3)
summary(weight_lm)

effects_info = allEffects(weight_lm)
effects_info
plot(effects_info)


############ 1.3 Multiple comparisons test

emmeans(weight_lm, list(pairwise~Plant), adjust="tukey")

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

########### 3.1 SA - GREENHOUSE COMPARISON
## Histogram of surface area by plant species
myplot<-ggplot(data=GH_data, aes(x=Plant, y=log(TotalSA)),na.action=na.exclude)
myplot+geom_boxplot(notch=TRUE)
stripchart(TotalSA ~ Plant, data=GH_data,              # Data
                      method = "jitter", # Random noise
                      pch = 19,          # Pch symbols
                      col = 4,           # Color of the symbol
                      add = TRUE) 

#### Compare surface area of the greenhouse plants
## scatter plot by plant species
ggplot(data = GH_data) +
  aes(y = Plant, x = log(TotalSA)) +
  geom_beeswarm() +
  coord_flip()
boxplot(TotalSA~Plant,data=GH_data)
#Model 
SA_comparison_GH = lm(TotalSA ~ Plant, data=GH_data)
summary(SA_comparison_GH)

#model effects
plot(allEffects(SA_comparison_GH), ylab="Floral Surface Area")
Anova(SA_comparison_GH)


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

SA_field_lm = lm(TotalSA~Plant,data=Field_data)
Anova(SA_field_lm)
plot(allEffects(SA_field_lm))


###################### 4 FIELD PLANTS & COMPARISONS #####################

# Because greenhouse SolAlt (Tall goldenrod) had a different effect from field SolAlt,
# I will specifically compare the two field species (EutMac/Joe-pye weed & SolAlt)
# as they cannot be pooled with the greenhouse data
summary(Field_data)


##################### 4.1 CHECK FOR ENCLOSURE EFFECTS ####################
# Check for enclosure effects with just SolAlt data (SolAlt had two different types
# of enclosures in the trial)
encl_lm = lm(RawWeight_day7~ Sex + ForewingLength + EmergDate + EnclCol, data=SolAlt_F_data)

Anova(encl_lm) #EnclCol least significant, remove and compare
encl_reduced = lm(RawWeight_day7~ Sex + ForewingLength + EmergDate, data=SolAlt_F_data)
anova(encl_reduced, encl_lm) # more complex (with EnclCol) is not significantly better, leave out EnclCol

Anova(encl_reduced) #EmergDate least significant, remove and compare
encl_reduced2 = lm(RawWeight_day7~ Sex + ForewingLength, data=SolAlt_F_data)
anova(encl_reduced2, encl_reduced) # more complex (with EmergDate) is not significantly better
Anova(encl_reduced2)

#Sex least significant, remove and compare
encl_reduced3 = lm(RawWeight_day7~ ForewingLength, data=SolAlt_F_data)
anova(encl_reduced3, encl_reduced2) # more complex (with EmergDate) is not significantly better


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
summary(females_data)

field_fat = lm(RawWeight_day7~ Plant + ForewingLength + NumButterflies + TotalSA, data=females_data)

## remove plant, least significant term, and compare
#field_reduced = lm(RawWeight_day7~ ForewingLength + NumButterflies, data=females_data)

# anova(field_reduced, field_lm) #more complex is not significantly better
# Anova(field_reduced) #remove numButterflies
# 
# field_reduced2 = lm(RawWeight_day7~ ForewingLength, data=females_data)
# anova(field_reduced2,field_reduced)
# Anova(field_reduced2)



############### 5. FAT EXTRACTION ANALYSES ##########################



############### 5.1 Exploration ##################################
opar <- par(mfrow = c(1, 1))
boxplot(WaterMass~Plant,data=data,ylab="Body water content (g)")
water_lm = lm(WaterMass~Plant,data=data)
summary(water_lm)

ggplot(data, aes(x = Plant, y = WaterMass, color = ExpLoc)) +  # ggplot function
  geom_point()

scatterplot(DryFatMass~RawWeight_day7,data=data, ylab="Body fat (g)",xlab="Body mass - wet (g)")
fat_lm =lm(DryFatMass~RawWeight_day7,data=data)
summary(fat_lm)

ggplot(data, aes(x = RawWeight_day7, y = log(DryFatMass), color = ExpLoc)) +  # ggplot function
  geom_point()+
  geom_smooth(method=lm , color="red", se=TRUE) +

scatterplot(RelDryFat~RawWeight_day7,data=data,ylab="% body fat",xlab="Body mass - dry (g)")
scatterplot(RelDryFat~DryMass,data=data,ylab="% body fat",xlab="Body mass - dry (g)")

ggplot(data, aes(x = Plant, y = RelDryFat, color = ExpLoc)) +  # ggplot function
  geom_point()


relfat_lm = lm(RelDryFat~RawWeight_day7,data=data)
summary(relfat_lm)

scatterplot(DryMass~RawWeight_day7,data=data,ylab="Dry body mass (g)",xlab="Body mass - wet (g)")
drymass_lm = lm(DryMass~RawWeight_day7,data=data)
summary(drymass_lm)

scatterplot(DryFatMass~DryLeanMass,data=data,ylab="Body fat (g)",xlab="Lean mass (g)")

scatterplot(DryMass~EmergDate,data=data)


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
# Anova(gh_fat_lm_log)
# 
# summary(gh_fat_lm_log)
# 
# pred_gh_fat_lm=allEffects(gh_fat_lm)
# pred_gh_fat_log=allEffects(gh_fat_lm_log)


# plot(pred_gh_fat_lm)
# plot(pred_gh_fat_log)
# Anova(gh_fat_lm_log)

predict_plant = predictorEffect("Plant",gh_fat_lm_log)

str(predict_plant)
predict_plant$fit = exp(predict_plant$fit)
predict_plant$upper=exp(predict_plant$upper)
predict_plant$lower=exp(predict_plant$lower)
plot(predict_plant,ylab ="Fat mass (g)")
plot(allEffects(gh_fat_lm_log))

emmeans(gh_fat_lm_log, list(pairwise~Plant), adjust="tukey")


################### GLM Using Gamma distribution for Fat #######################
gh_fat_lm_gm = glm(DryFatMass ~ Sex + Plant + ForewingLength + EmergDate + TotalSA, data=GH_data,family="Gamma")
opar <- par(mfrow = c(2, 2))
plot(gh_fat_lm_gm)

Anova(gh_fat_lm_gm)
summary(gh_fat_lm_gm)
pred_gh_fat_gm=allEffects(gh_fat_lm_gm)
plot(pred_gh_fat_gm)

str(pred_gh_fat_gm$Plant)
str(pred_gh_fat_gm$Plant$response)
str(pred_gh_fat_gm$fit)
str(pred_gh_fat_gm)

summary(gh_fat_lm_gm)
str(pred_gh_fat_gm)
plot(pred_gh_fat_gm)
check_model(gh_fat_lm_gm)

plot(pred_gh_fat_gm)
emmeans(gh_fat_lm_gm, list(pairwise~Plant), adjust="tukey")

## Plot of effects using ggplot glm
ggplot(GH_data, aes(EmergDate, DryFatMass)) + 
  geom_smooth(method = "glm", formula = DryFatMass~ Plant + Sex + ForewingLength + EmergDate + TotalSA, colour = "black",
              linetype = 2, fill = "gray80", alpha = 0.2,
              method.args = list(family = gamma)) +
  geom_rug(sides = "b") +
  theme_bw() +
  labs(y = "Fat weight (g)", x = "Plant") +
  theme(text = element_text(size = 16),
        plot.margin = margin(50, 50, 50, 50),
        axis.title.x = element_text(vjust = -8),
        axis.title.y = element_text(vjust = 10),
        plot.title = element_text(vjust = 8))






########## relative fat content (%)
gh_relfat_lm = lm(RelDryFat~ Plant + EmergDate + TotalSA, data=GH_data)
vif(gh_relfat_lm)
Anova(gh_relfat_lm)

emmeans(gh_relfat_lm, list(pairwise~Plant), adjust="tukey")



######################## 5.3 FIELD FAT ANALYSIS #######################
opar <- par(mfrow = c(2, 2))
field_fat_lm = lm(DryFatMass~ EmergDate + TotalSA, data=Field_data)
Anova(field_fat_lm)
plot(allEffects(field_fat_lm))
vif(field_fat_lm)

opar <- par(mfrow = c(2, 2))
plot(field_fat_lm)

field_fat_glm = glm(DryFatMass~ EmergDate, data=Field_data,family="Gamma")
Anova(field_fat_glm)

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




##################### WATER ############################

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


