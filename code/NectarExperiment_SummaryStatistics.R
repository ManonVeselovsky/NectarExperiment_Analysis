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
require(ggplot2)
library(emmeans)
library(multcomp)


# create separate working databases for greenhouse plants, goldenrod (for GH-Field comparison), and field plants
GH_data = subset(summarydat, summarydat$ExpLoc == "GH")
GH_F_data = summarydat[which(summarydat$Plant=="SolAlt"),] #all data on goldenrod, for GH-Field comparison
Field_data = subset(summarydat,summarydat$ExpLoc == "Field")
SolAlt_F_data = subset(Field_data,Field_data$Plant == "SolAlt") #Subset SolAlt in field to check for enclosure effects


########## GREEHOUSE DATA ANALYSIS ##############
################## Normality check ###############
opar <- par(mfrow = c(1, 1))

#Check for normality of raw weight on day 7
mygraph <- ggplot(GH_data, aes(x = RawWeight_day7))
mygraph <- mygraph +
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
mygraph
qqnorm(GH_data$RawWeight_day7)
qqline(GH_data$RawWeight_day7)

# Shapiro formal test for normality
shapiro.test(GH_data$RawWeight_day7)



#################### MODEL BUILDING ####################################

#check collinearity of starting weight and forewing length (to see if I need both)
fwl_day0w_lm = lm(RawWeight_day0 ~ ForewingLength, data=GH_data)
Anova(fwl_day0w_lm, type=3)
summary(fwl_day0w_lm)

# Estimate the change in raw weight on day 7 of trial based on predictor variables
lmer.model = lmer(RawWeight_day7~ Plant + Sex + ForewingLength + EmergDate + Cohort + EnclCol + (1|Couple), data=GH_data)


#Test for collinearity among predictor variables (GVIF<4 for continuous acceptable,
# GVIF^(1/(2*df)) < 2 acceptable for categorical)
vif(lmer.model) #all good
opar <- par(mfrow = c(2, 2))
Anova(lmer.model, type=3)
plot(lmer.model)


## Reduce model to remove insignificant terms, starting with cohort
reduced_1 = lmer(RawWeight_day7~ Plant + Sex + ForewingLength + EmergDate + EnclCol + (1|Couple), data=GH_data)

#test to see if more complex model is significantly better at describing the variation
anova(reduced_1, lmer.model) #not significantly better, leave cohort out
Anova(reduced_1, type=3)

# Remove next least significant term (EmergDate)
reduced_2 = lmer(RawWeight_day7~ Plant + Sex + ForewingLength + EnclCol + (1|Couple), data=GH_data)

#test to see if more complex is significantly better
anova(reduced_2, reduced_1) #not significantly better, keep the simpler model
Anova(reduced_2, type=3)

# remove next least significant (EnclCol)
reduced_3 = lmer(RawWeight_day7~ Plant + Sex + ForewingLength + (1|Couple), data=GH_data)
#test to see if reduced is significantly better
anova(reduced_3, reduced_2)

#more complex is significantly better, keep reduced_2
bestmodel = reduced_2


#dd <- dredge(lmer.model)
# 
# # get models within 4 units of AICc from the best model
# top.models.1 <- get.models(dd, subset = delta < 4)
# avgmodel1<-model.avg(top.models.1) # compute average parameters
# summary(avgmodel1) #display averaged model
# confint(avgmodel1) #display CI for averaged coefficients
# model.avg(object = top.models.1)



############ Multiple comparisons test
m1<-bestmodel
library(multcomp)
#WHICH groups are different from each other
g<-glht(m1, mcp(Plant="Tukey")); confint(g)
#OR
library(emmeans)
emmeans(m1, list(pairwise~Plant), adjust="tukey")

plot(allEffects(m1))

# diagnostic plots
plot(m1)
qqnorm(resid(m1))
qqline(resid(m1))
require("lattice")
qqmath(m1,id=0.05)

######## COMPARISON OF GH TO FIELD DATA (USING GOLDENROD) ##############
#Create a database of field and gh goldenrod to compare experiments
GH_F_data = summarydat[which(summarydat$Plant=="SolAlt"),]
summary(GH_F_data)

#Check for normality on raw day 7 weights
hist(GH_F_data$RawWeight_day7)

#linear mixed effects model to compare GH to field
#gh_f_lmer = lmer(RawWeight_day7~ Sex + ForewingLength + ExpLoc + (1|Couple), data=GH_F_data)
#Anova(gh_f_lmer, type=3)
#plot(allEffects(gh_f_lmer))

### linear model to compare GH method to field method (ExpLoc)
gh_f_lm = lm(RawWeight_day7~ Sex + ForewingLength + ExpLoc, data=GH_F_data)
Anova(gh_f_lm, type=3)

plot(allEffects(gh_f_lm))
summary(gh_f_lm)


################ CHECK PLANT SURFACE AREA
SA_comparison = lm(TotalSA ~ Plant + ExpLoc + Plant:ExpLoc, data=data)
summary(SA_comparison)
# emmeans(SA_comparison, list(pairwise~Plant), adjust="tukey")
# TukeyHSD(SA_comparison)
plot(allEffects(SA_comparison))
Anova(SA_comparison)

SA_reduced = lm(TotalSA ~ Plant + ExpLoc, data=data)

anova(SA_reduced, SA_comparison)

SA_reduced2 = lm(TotalSA ~ Plant, data=data)

Anova(SA_reduced2)
anova(SA_reduced2, SA_reduced)

Anova(SA_reduced2)
summary(SA_reduced2)

########## Redo but with lm

lm.model = lm(RawWeight_day7~ Plant + Sex + ForewingLength + EmergDate + Cohort + EnclCol, data=GH_data)

#Test for collinearity among predictor variables (GVIF<4 for continuous acceptable,
# GVIF^(1/(2*df)) < 2 acceptable for categorical)
vif(lm.model) #remove emerg date, CORRELATED WITH COHORT

lm.model = lm(RawWeight_day7~ Plant + Sex + ForewingLength + Cohort + EnclCol, data=GH_data)
vif(lm.model)
opar <- par(mfrow = c(2, 2))
Anova(lm.model, type=3)
plot(lm.model)

## Reduce model to remove insignificant terms, starting with cohort
reduced_1 = lm(RawWeight_day7~ Plant + Sex + ForewingLength + EnclCol, data=GH_data)

#test to see if more complex model is significantly better at describing the variation
anova(reduced_1, lm.model) #not significantly better, leave cohort out
Anova(reduced_1, type=3)

# Remove next least significant term (EnclCol)
reduced_2 = lm(RawWeight_day7~ Plant + Sex + ForewingLength, data=GH_data)

#test to see if more complex is significantly better
anova(reduced_2, reduced_1) #significantly better, keep the more complex model (reduced 1)
bestmodel = reduced_1
Anova(bestmodel, type=3)


############ Multiple comparisons test for lm ###################
m1<-bestmodel
library(multcomp)
#WHICH groups are different from each other
g<-glht(m1, mcp(Plant="Tukey")); confint(g)
#OR
library(emmeans)
emmeans(m1, list(pairwise~Plant), adjust="tukey")

plot(allEffects(m1))

# diagnostic plots
plot(m1)
require("lattice")
qqnorm(resid(m1))
qqline(resid(m1))
summary(bestmodel)
Anova(bestmodel, type=3)



###################### Compare field plants ##############

summary(Field_data)

#subset for just females in the field (used doubles in JPW enclosures and only females for those, only 1 male on JPW field)
females_data = subset(Field_data,Field_data$Sex == "F")


#### Check for enclosure effects with just SolAlt data
encl_lm = lm(RawWeight_day7~ Sex + ForewingLength + EnclCol, data=SolAlt_F_data)

Anova(encl_lm) #Encl col least significant, remove and compare
encl_reduced = lm(RawWeight_day7~ Sex + ForewingLength, data=SolAlt_F_data)
anova(encl_reduced, encl_lm) # more complex (with EnclCol) is not significantly better, can pool across enclcol

#### build model using females, without EnclCol
field_lm = lm(RawWeight_day7~ Plant + ForewingLength + NumButterflies, data=females_data)
Anova(field_lm)
summary(field_lm)
plot(allEffects(field_lm))

## remove plant, least significant term, and compare
field_reduced = lm(RawWeight_day7~ ForewingLength + NumButterflies, data=females_data)

anova(field_reduced, field_lm) #more complex is not significantly better
Anova(field_reduced) #remove numButterflies

field_reduced2 = lm(RawWeight_day7~ ForewingLength, data=females_data)
anova(field_reduced2,field_reduced)
Anova(field_reduced2)
