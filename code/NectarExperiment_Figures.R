### Statistical analysis script
### Written by M. Veselovsky
### Last modified 2024-04-24
rm(list=ls())

library(lmerTest)
library(sjPlot) #For model summary figures
library(sjmisc) #for model tables
library(sjlabelled) #for model tables
library(ggplot2)
library(car)
library(lattice) #For creating plots with multiple panels vertically, common x values
library(plotrix) #to get mean and standard errors of columns
library(effects)
library(ggeffects)
library(emmeans)
library(tidyverse) #re-order categorical variables in effects plots
library(dplyr)

#raw treatment data with all trial days
TreatmentData<-read.csv("processed/TreatmentData_p.csv")
summarydat<-read.csv("processed/summarydat.csv")
tempdat = read.csv("processed/tempdat.csv")

# create separate working databases for greenhouse plants, goldenrod (for GH-Field comparison), and field plants
GH_data = subset(summarydat, summarydat$ExpLoc == "GH") # data from greenhouse
GH_F_data = summarydat[which(summarydat$Plant=="SolAlt"),] #all data on goldenrod, for GH-Field comparison
Field_data = subset(summarydat,summarydat$ExpLoc == "Field") #data from field
SolAlt_F_data = subset(Field_data,Field_data$Plant == "SolAlt") #Subset SolAlt in field to check for enclosure effects
gh_TreatmentData = subset(TreatmentData, TreatmentData$ExpLoc == "GH")
gh_field_weight = TreatmentData[which(TreatmentData$Plant=="SolAlt"),]

## Create a theme to use on all ggpredict plots
my_plot_theme = theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size = 14),
        axis.title=element_text(size=18),
        line=element_line(linewidth=0.75)
        )

################ FINAL MODELS FOR EACH VARIABLE ######################

######### WEIGHT EFFECTS 
#remove cohort
overall_weight = lmer(Weight~AlphPlant+ForewingLength+DateWeighed+Sex+TrialDay+EnclCol+(1|ID),data=gh_TreatmentData)
table_weight = lmer(Weight~Plant+OrigWeight+DateWeighed+Sex+TrialDay+EnclCol+(1|ID),data=gh_TreatmentData)

check_model(overall_weight)

#remove dateweighed - collinearity
overall_weight = lmer(Weight~AlphPlant+OrigWeight+Sex+TrialDay+EnclCol+(1|ID),data=gh_TreatmentData)
table_weight = lmer(Weight~Plant+OrigWeight+Sex+TrialDay+EnclCol+(1|ID),data=gh_TreatmentData)
check_model(overall_weight)



summary(overall_weight)
Anova(table_weight)
effects_weight = allEffects(overall_weight)
plot(ggpredict(overall_weight))
Anova(table_weight)

plot(effects_weight$AlphPlant,ylab="Weight (g)")
effects_weight_gg = ggpredict(overall_weight, ci.lvl = 0.95)
weight_plot = plot(effects_weight_gg$AlphPlant) +
  labs(y = "Fat mass (g)",x = "Plant",title="") +
  my_plot_theme

weight_plot

overallplot
# overallplot + scale_x_discrete(
#   labels = c(
#     "1_SolAlt" = "Goldenrod",
#     "2_BudDav" = "Butterfly bush",
#     "3_SymEri" = "Heath Aster",
#     "4_EchPur" = "Coneflower",
#     "5_RudHir" = "Black-eyed Susan",
#     "6_HelHel" = "Ox-eye")
#   )

emmeans(table_weight,list(pairwise~Plant), adjust="tukey")
tab_model(table_weight)

###########FAT PLOTS
gh_fat_lm_log = lm(log(DryFatMass)~ AlphPlant + Sex + ForewingLength + EmergDate + TotalSA, data=GH_data)
tab_fat = lm(log(DryFatMass)~ Plant + Sex + ForewingLength + EmergDate + TotalSA, data=GH_data)

simfat = lm(log(DryFatMass)~ AlphPlant, data=GH_data)
Anova(simfat)
plot(ggpredict(simfat))
summary(simfat)

# plot(allEffects(tab_fat)) + my_plot_theme
effects_fat = ggpredict(gh_fat_lm_log, ci.lvl = 0.95)

fatplot = plot(effects_fat$AlphPlant) +
  labs(y = "Fat (g)",x = "Plant",title="") +
  my_plot_theme

#Plots
fatplot
tab_model(tab_fat)
emmeans(tab_fat,list(pairwise~Plant), adjust="tukey")
Anova(gh_fat_lm_log)

############LEAN MASS PLOTS
gh_prot_lm = lm(DryLeanMass ~ AlphPlant + Sex + ForewingLength + EmergDate + TotalSA, data=GH_data)
tab_lean = lm(DryLeanMass ~ Plant + Sex + ForewingLength + EmergDate + TotalSA, data=GH_data)
effects_prot = allEffects(gh_prot_lm)
plot(effects_prot)
effects_lean = ggpredict(gh_prot_lm, ci.lvl = 0.95)

leanplot = plot(effects_lean$AlphPlant) +
  labs(y = "Lean mass (g)",x = "Plant",title="") +
  my_plot_theme

leanplot

Anova(gh_prot_lm)
# plot(effects_prot,ylab="Lean mass (g)")
# plot(effects_prot$AlphPlant,ylab = "Lean mass(g)")
tab_model(tab_lean)
emmeans(tab_lean,list(pairwise~Plant),adjust="tukey")

################# WATER PLOTS
gh_water_lm = lm(WaterMass ~ Plant + Sex + ForewingLength + EmergDate + TotalSA, data=GH_data)
effects_prot = allEffects(gh_water_lm)
Anova(gh_water_lm)
effects_water = allEffects(gh_water_lm)
plot(effects_water)

############### GH VS FIELD

ghf_weight = lmer(Weight~ ExpLoc+TrialDay+Sex+(1|ID), data=gh_field_weight)
library(performance)
check_model(ghf_weight)
tab_model(ghf_weight)
effects_ghf = ggpredict(ghf_weight)
plot(effects_ghf$ExpLoc) + my_plot_theme + labs(y = "Weight (g)",x = "Location",title="")
Anova(ghf_weight)


ghf_fat = lm(DryFatMass~ExpLoc + Sex, data=GH_F_data)
check_model(ghf_fat)
Anova(ghf_fat)
ghf_fat_pred = ggpredict(ghf_fat)
ghf_plot = plot(ghf_fat_pred)

ghf_plot = plot(ghf_fat_pred$ExpLoc) +
  labs(y = "Fat (g)",x = "Location",title="") +
  my_plot_theme

Anova(ghf_fat)
ghf_plot
################# DESCRIPTIVE SCATTERPLOTS ###########

scatterplot(Weight~TrialDay,data=TreatmentData,ylab="Weight (g)",xlab="Trial day")
scatterplot(Weight~EmergDate,data=TreatmentData,ylab="Weight (g)",xlab="Julian date")
scatterplot(DryFatMass~EmergDate +(1|Plant),data=GH_data,ylab="Body fat (g)",xlab="Trial start date",smooth=FALSE,grid=FALSE)

boxplot(Weight~TrialDay,data=trialday_data)
day_factor = trialday_data
day_factor$TrialDay = as.factor(day_factor$TrialDay)
trialday_lmer = lmer(Weight~TrialDay + ExpLoc+TrialDay:ExpLoc+ (1|ID),data=day_factor)
Anova(trialday_lmer)
plot(allEffects(trialday_lmer),ylim=c(0.3,0.65))


############### Males vs females ##################
boxplot(RelWeightGain_day7*100~Sex,data=GH_data, ylab="% Weight change",xlab="Sex",ylim=c(-60,0))


data_msd = GH_data %>%
  group_by(Sex) %>%
  summarise_at(vars(RelWeightGain_day7),
               list(mean = mean,
                    sd = sd))
data_msd

ggplot(GH_data, aes(Sex, mean)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = mean - sd, 
                    ymax = mean + std.error(mean)),
                width = 0.1, linewidth = 0.1) +
  theme_bw() +
  labs(x = "Sex",
       y = "% Weight change")

################### TRIAL DATE RANGES ######################

# basic plot
p=ggplot(summarydat, aes(EmergDate, as.factor(Plant),color=ExpLoc,show.legend=TRUE)) + 
  geom_point()+ geom_jitter(width=0,height=0.2)+
  labs(y = "Plant", x = "Start of trial (Julian date)")
#remove legend 
p_nolegend = p + theme(legend.position="none")
#Make theme of the plot black and white (apart from plant color) and remove unnecessary elements
p_nolegend + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(),
                   panel.grixd.minor = element_blank(), axis.line = element_line(colour = "black"))

mean(GH_data$RelWeightGain_day7)
std.error(GH_data$RelWeightGain_day7)

mean(Field_data$RelWeightGain_day7)
std.error(Field_data$RelWeightGain_day7)

######################## Temperature data ###########################
fieldtemp = subset(tempdat,tempdat$Location == "Fletcher")
ghtemp = subset(tempdat,tempdat$Location =="Greenhouse")
summary(ghtemp)


########### GH EFFECTS #####################################################

####################### Fat effects ###########
plot(allEffects(gh_fat_lm_log), ylab = "Log(Fat mass (g))")
emmeans(gh_fat_lm_log,list(pairwise~Plant), adjust="tukey")
Anova(gh_fat_lm_log)
tab_model(gh_fat_lm_log)

###################### Protein effects ############
plot(allEffects(gh_prot_lm))
Anova(gh_prot_lm)
emmeans(gh_prot_lm,list(pairwise~Plant), adjust="tukey")
tab_model(gh_prot_lm)
summary(gh_prot_lm)

###################### Water effects##############
plot(allEffects(gh_water_lm))
Anova(gh_water_lm)
summary(gh_water_lm)
emmeans(gh_water_lm,list(pairwise~Plant), adjust="tukey")
tab_model(gh_water_lm)


############ Table of all physiological effects #########
tab_model(gh_fat_lm_log, gh_prot_lm, gh_water_lm, auto.label=FALSE)

##################### Overall weight ###############



#weight
overall_weight = lmer(Weight~Plant+OrigWeight+DateWeighed+Sex+TrialDay+EnclCol+Cohort+(1|ID),data=gh_TreatmentData)
effects_weight_gg = ggpredict(overall_weight, ci.lvl = 0.95)
plot(effects_weight_gg$Plant) + labs(y = "Weight (g)") + my_plot_theme

ordered_data = gh_TreatmentData$Plant == "SolAlt" | gh_TreatmentData$Plant == "BudDav" | gh_TreatmentData$Plant == "SymEri" | gh_TreatmentData$Plant=="EchPur"| gh_TreatmentData$Plant=="RudHir"| gh_TreatmentData$Plant=="HelHel"

ordered_data$Plant = factor(gh_TreatmentData$Plant, levels=c("SolAlt", "BudDav", "SymEri", "EchPur", "RudHir", "HelHel"))

overall_weight_ordered = lmer(Weight~Plant+OrigWeight+DateWeighed+Sex+TrialDay+EnclCol+Cohort+(1|ID),data=ordered_data)

effects_ordered = ggpredict(overall_weight_ordered,ci.lvl=0.95)

plot(effects_ordered)


glimpse(effects_weight_gg)

effects_weight_gg %>%
  arrange(Plant) %>%
  mutate(Plant = factor(Plant, levels=c("SolAlt", "BudDav", "SymEri", "EchPur", "RudHir", "HelHel"))) %>%
  plot( aes(x=Plant, y=Plant$predicted)) +
  geom_segment( aes(xend=Plaant, yend=0)) +
  geom_point( size=4, color="orange") +
  theme_bw() +
  xlab("")

