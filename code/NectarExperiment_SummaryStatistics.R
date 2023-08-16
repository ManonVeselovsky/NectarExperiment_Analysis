### Statistical analysis script
### Written by M. Veselovsky
### Last modified 2023-01-13

# clear the R environment
rm(list=ls())
usethis::git_sitrep()
usethis::use_git_config(
  user.name = "ManonCV",
  user.email = "mvese092@uottawa.ca"
)
usethis::create_github_token()
credentials::set_github_pat("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")
usethis::use_git()
usethis::use_github()
usethis::git_vaccinate()


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
library(sjPlot)
library(sjmisc)
library(ggplot2)
library(ggeffects)


# create separate working databases for greenhouse plants, goldenrod (for GH-Field comparison), and field plants
GH_data = subset(summarydat, summarydat$ExpLoc == "GH") # data from greenhouse
GH_F_data = summarydat[which(summarydat$Plant=="SolAlt"),] #all data on goldenrod, for GH-Field comparison
Field_data = subset(summarydat,summarydat$ExpLoc == "Field") #data from field
SolAlt_F_data = subset(Field_data,Field_data$Plant == "SolAlt") #Subset SolAlt in field to check for enclosure effects


########## 1. GREEHOUSE DATA ANALYSIS ##############
################## 1.1 Normality check ###############
opar <- par(mfrow = c(1, 1))

#Check for normality of the butterflies raw weights on day 7 of their respective trials
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

# Formal test for normality (Shapiro)
shapiro.test(GH_data$RawWeight_day7)

# boxplots of raw weight on day 7 by plant species in the greenhouse
myplot<-ggplot(data=GH_data, aes(x=Plant, y=RawWeight_day7),na.action=na.exclude)
myplot+geom_boxplot(notch=TRUE)




#################### 1.2 MODEL BUILDING - LMER ###########################

# Check collinearity of starting weight and forewing length (to see if I need both in the model)
# The a-priori expectation is that they will be correlated, with forewing length the ideal variable
# as it does not change with adult age. Weights could fluctuate depending on the time they were taken
# since emergence (metabolism as they have not eaten, and drying time for the fluid they expel
# as they eclose from their chrysalises)

fwl_day0w_lm = lm(RawWeight_day0 ~ ForewingLength, data=GH_data)
Anova(fwl_day0w_lm, type=3)# Significantly related, so I will keep forewing length and not use weight
plot(allEffects(fwl_day0w_lm))
summary(fwl_day0w_lm) 

# Estimate the change in raw weight on day 7 of trial based on predictor variables
lmer.model = lm(RawWeight_day7~ Plant + Sex + ForewingLength
                  + Cohort + EnclCol + TotalSA, data=GH_data)

#Test for collinearity among predictor variables (GVIF<4 for continuous acceptable,
# GVIF^(1/(2*df)) < 2 acceptable for categorical)
vif(lmer.model) #all good
Anova(lmer.model, type=3)
summary(lmer.model)
plot(lmer.model)
plot(lmer.model, resid(., scaled=TRUE) ~ fitted(.), abline = 0,pch=16,xlab="Fitted values",ylab="Standardised residuals")

effects_info = allEffects(lmer.model)
effects_info
plot(effects_info)
#plot(lmer.model, resid(., scaled=TRUE) ~ fitted(.), abline = 0,pch=16,col=GH_data$Couple,xlab="Fitted values",ylab="Standardised residuals")

## Reduce model to remove insignificant terms, starting with cohort
# reduced_1 = lmer(RawWeight_day7~ Plant + Sex + ForewingLength + EmergDate + EnclCol + (1|Couple), data=GH_data)
# 
# #test to see if more complex model is significantly better at describing the variation
# anova(reduced_1, lmer.model) #not significantly better, leave cohort out
# Anova(reduced_1, type=3)
# 
# # Remove next least significant term (EmergDate)
# reduced_2 = lmer(RawWeight_day7~ Plant + Sex + ForewingLength + EnclCol + (1|Couple), data=GH_data)
# 
# #test to see if more complex is significantly better
# anova(reduced_2, reduced_1) #not significantly better, keep the simpler model
# Anova(reduced_2, type=3)
# 
# # remove next least significant (EnclCol)
# reduced_3 = lmer(RawWeight_day7~ Plant + Sex + ForewingLength + (1|Couple), data=GH_data)
# #test to see if reduced is significantly better
# anova(reduced_3, reduced_2)

#more complex is significantly better, keep reduced_2
# bestmodel = reduced_2
bestmodel = lmer.model
Anova(bestmodel, type=3)
summary(bestmodel)
plot(allEffects(bestmodel))
#library(MuMIn)
#dd <- dredge(lmer.model)
 
# # get models within 4 units of AICc from the best model
# top.models.1 <- get.models(dd, subset = delta < 4)
# avgmodel1<-model.avg(top.models.1) # compute average parameters
# summary(avgmodel1) #display averaged model
# confint(avgmodel1) #display CI for averaged coefficients
# model.avg(object = top.models.1)



############ 1.3 Multiple comparisons test
m1<-bestmodel
Anova(bestmodel)
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

######## 2. COMPARISON OF GH TO FIELD DATA (USING GOLDENROD) ##############
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


################ 3. PLANT SURFACE AREA DIFFERENCES BETWEEN SPECIES ############

########### 3.1 GREENHOUSE COMPARISON
## Histogram of surface area by plant species
myplot<-ggplot(data=GH_data, aes(x=Plant, y=log(TotalSA)),na.action=na.exclude)
myplot+geom_boxplot(notch=TRUE)

#Compare surface area of the greenhouse plants
SA_comparison_GH = lm(TotalSA ~ Plant, data=GH_data)

summary(SA_comparison_GH)
plot(allEffects(SA_comparison_GH), ylab="Floral Surface Area")
Anova(SA_comparison_GH)


########### 3.2 FIELD COMPARISON
myplot<-ggplot(data=Field_data, aes(x=Plant, y=TotalSA),na.action=na.exclude)
myplot+geom_boxplot(notch=TRUE)

##################### 4. MODEL BUILDING WITH LM #################################
# Individuals that were collected in the field did not have a known parent couple.
# These individuals were excluded from my LMER as I had the grouping "Couple" -->
# Repeating the analysis here with a simple lm WITHOUT the "couple" variable
cohort_lm = lm(RawWeight_day7 ~ Cohort, data=GH_data)
summary(cohort_lm)
plot(allEffects(cohort_lm))
plot(cohort_lm)

lm.model = lm(RawWeight_day7~ Plant + Sex + ForewingLength + EmergDate + TotalSA, data=GH_data)

#Test for collinearity among predictor variables (GVIF<4 for continuous acceptable,
# GVIF^(1/(2*df)) < 2 acceptable for categorical)
vif(lm.model) #remove emerg date, CORRELATED WITH COHORT

lm.model = lm(RawWeight_day7~ Plant + Sex + ForewingLength + EmergDate + TotalSA, data=GH_data)
vif(lm.model)
opar <- par(mfrow = c(2, 2))
Anova(lm.model, type=3)
library(knitr)
kable(lm.model, digits = 3)
summary(lm.model)
plot(lm.model)

lm_effects = allEffects(lm.model)
lm_effects
plot(lm_effects, ylab = "Day 7 Weight (g)")
plot(allEffects(m1)$Plant, )
summary(lm.model)
lm_effects
## Reduce model to remove insignificant terms, starting with cohort
# reduced_1 = lm(RawWeight_day7~ Plant + Sex + ForewingLength + EnclCol, data=GH_data)
# 
# #test to see if more complex model is significantly better at describing the variation
# anova(reduced_1, lm.model) #not significantly better, leave cohort out
# Anova(reduced_1, type=3)

# Remove next least significant term (EnclCol)
# reduced_2 = lm(RawWeight_day7~ Plant + Sex + ForewingLength, data=GH_data)
# 
# #test to see if more complex is significantly better
# anova(reduced_2, reduced_1) #significantly better, keep the more complex model (reduced 1)
# bestmodel = reduced_1
bestmodel=lm.model
Anova(bestmodel, type=3)
summary(bestmodel)

############ 4.1 MULTIPLE COMPARISONS TEST FOR LM ###################

# Check for significant differences between plant species and their effect on day7 raw weight
# using a Tukey adjustment for multiple comparisons/hypothesis tests
m1<-lm.model
library(multcomp)
#WHICH groups are different from each other
g<-glht(m1, mcp(Plant="Tukey")); confint(g)
#OR
library(emmeans)
emmeans(m1, list(pairwise~Plant), adjust="tukey")
effects_lm.model = allEffects(m1)
plot(allEffects(m1))
plot(allEffects(m1)$Plant, ylab = "Weight (g)", xlab = "Plant")
str(effects_lm.model)
ggplot(GH_data, aes(x, y)) +
  geom_smooth(method = "lm", formula = RawWeight_day7 ~ Plant + Sex + ForewingLength + EmergDate + EnclCol, colour = "black",
              linetype = 2, fill = "gray80", alpha = 0.2) +
  geom_rug(sides = "b") +
  theme_bw() +
  labs(y = "Weight (g)", x = "Plant",
       title = "Plant effect plot") +
  theme(text = element_text(size = 16),
        plot.margin = margin(50, 50, 50, 50),
        axis.title.x = element_text(vjust = -8),
        axis.title.y = element_text(vjust = 10),
        plot.title = element_text(vjust = 8))

# diagnostic plots
plot(m1)
require("lattice")
qqnorm(resid(m1))
qqline(resid(m1))
summary(bestmodel)
plot_model(bestmodel, type="pred", terms = c("Plant"),ylab="Weight (g)",)
Anova(bestmodel, type=3)
ggpredict(lm, terms = c("Plant","RawWeight_day7"))



###################### 5. FIELD PLANTS & COMPARISONS #####################

# Because greenhouse SolAlt (Tall goldenrod) had a different effect from field SolAlt,
# I will specifically compare the two field species (EutMac/Joe-pye weed & SolAlt)
# as they cannot be pooled with the greenhouse data
summary(Field_data)


##################### 5.1 CHECK FOR ENCLOSURE EFFECTS ####################
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


#################### 5.2 COMPARE FIELD SPECIES ########################
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

## remove plant, least significant term, and compare
#field_reduced = lm(RawWeight_day7~ ForewingLength + NumButterflies, data=females_data)

# anova(field_reduced, field_lm) #more complex is not significantly better
# Anova(field_reduced) #remove numButterflies
# 
# field_reduced2 = lm(RawWeight_day7~ ForewingLength, data=females_data)
# anova(field_reduced2,field_reduced)
# Anova(field_reduced2)

