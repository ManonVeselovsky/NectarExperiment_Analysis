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
library(performance)
library(multcomp) #for adding multcomp tukey to summary
library(patchwork) #for creating a multi-panel plot

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

GH_data$Fat_mg = GH_data$DryFatMass*1000

desired_order <- c("SymEri", "BudDav", "EchPur", "RudHir", "SolAlt", "HelHel")

# Reorder the Plant factor in the original data
GH_data$Plant <- factor(GH_data$Plant, levels = desired_order)
gh_TreatmentData$Plant = factor(gh_TreatmentData$Plant, levels = desired_order)


## Create a theme to use on all ggpredict plots
my_plot_theme = theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),  # Standardize axis line thickness
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        line = element_line(linewidth = 0.5)  # Standardize all other line elements
  )

################ FINAL MODELS FOR EACH VARIABLE ######################

######### WEIGHT EFFECTS
#Remove TrialDay = 0 because I already have an original weight column which is the TrialDay = 0 weight
alt_data = gh_TreatmentData %>% 
  filter(TrialDay %in% c("1", "5","7", "10", "11"))
overall_weight = lmer(Weight~Plant+TotalSA+ForewingLength+OrigWeight+Sex+TrialDay+DateWeighed+EnclCol+(1|ID),data=alt_data)

#Create a table of the model output in a word doc
tab_model(overall_weight, 
          file = "tables/Weight_Model_Results.doc",   # Export to a Word document
          show.p = TRUE)                   # Add significance stars


check_model(overall_weight) #Assumptions look good, weird error bars with TrialDay VIF
vif(overall_weight) #I don't see a problem with TrialDay here so I will leave it

Anova(overall_weight) #Plant, FWL, Start weight, sex, and trial day all significant

# Tukey-adjusted pairwise comparisons for Plant groups
tukey_weight <- emmeans(overall_weight, pairwise ~ Plant, adjust = "tukey")

# Extracting the compact letter display (CLD) to get group letters
cld_weight <- cld(tukey_weight$emmeans, Letters = letters)

# Extract predictions specifically for the Plant variable
effects_weight_plant <- ggpredict(overall_weight, terms = "Plant", ci.lvl = 0.95)

# Change the order for plant levels to the order from most visited to least visited
effects_weight_plant$x <- factor(effects_weight_plant$x)



# Ensure Plant is treated as a factor in the effects_weight_plant data
effects_weight_plant$x <- as.factor(effects_weight_plant$x)
# Reorder the Plant levels in cld_weight to match effects_weight_plant$x
cld_weight$Plant <- factor(cld_weight$Plant, levels = levels(effects_weight_plant$x))

# Ensure the Tukey-adjusted letters are ordered according to Plant levels
cld_weight <- cld_weight[order(cld_weight$Plant),]

# Now recreate the labels_df with the correctly ordered Plant factor
labels_df <- data.frame(
  Plant = effects_weight_plant$x,   # Correct factor levels from effects_weight_plant
  Label = cld_weight$.group,        # Use reordered Tukey letters
  y = effects_weight_plant$conf.high + 0.01  # Position just above the upper confidence interval
)

# Create plot using the specific Plant predictions
weight_plot <- ggplot(effects_weight_plant, show_residuals=TRUE, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(y = "Wet weight (g)", x = "", title = "") +
  my_plot_theme +
  geom_text(data = labels_df, aes(x = Plant, y = y, label = Label), vjust = -0.5) + # Add significance letters
  
  # 1. Change x-axis labels
  scale_x_discrete(labels = c("SymEri" = "S. ericoides", 
                              "BudDav" = "B. davidii", 
                              "EchPur" = "E. purpurea", 
                              "RudHir" = "R. hirta", 
                              "SolAlt" = "S. altissima", 
                              "HelHel" = "H. helianthoides")) +
  
  # 2. Angle the x-axis labels
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12),  # Increase font size for x-axis labels
        axis.text.y = element_text(size = 12),                          # Increase font size for y-axis labels
        axis.title.x = element_text(size = 12),                        # Increase font size for x-axis title
        axis.title.y = element_text(size = 12),                        # Increase font size for y-axis title
        panel.background = element_blank(),                             # Remove background color
        plot.background = element_blank(),                               # Remove plot background color
        panel.border = element_blank(),                                  # Remove panel border
        axis.line = element_line(color = "black"))                  # Keep axis lines

# Display plot with updated x-axis labels and correct significance letters
weight_plot

tukey_p_values <- summary(tukey_weight$emmeans)$p.value
tukey_p_values_df <- data.frame(Plant = names(tukey_p_values), Tukey_p_value = tukey_p_values)

# Display the model output with tab_model
tab_model(overall_weight)


###########FAT PLOTS

# Fit the linear model
overall_fat <- lm(log(DryFatMass) ~ Plant + Sex + ForewingLength + EmergDate + TotalSA, data = GH_data)

#Create a table of the model output in a word doc
tab_model(overall_fat, 
          file = "tables/Fat_Results.doc",   # Export to a Word document
          show.p = TRUE)                   # Add significance stars


# Tukey-adjusted pairwise comparisons for Plant groups
tukey_fat <- emmeans(overall_fat, pairwise ~ Plant, adjust = "tukey")

# Extracting the compact letter display (CLD) to get group letters
cld_fat <- cld(tukey_fat$emmeans, Letters = letters)

# Extract predictions specifically for the Plant variable
effects_fat_plant <- ggpredict(overall_fat, terms = "Plant", ci.lvl = 0.95)

# Convert ggpredict object to a data frame
effects_fat_plant_df <- as.data.frame(effects_fat_plant)

# Round predicted values and confidence intervals to 4 decimal places
effects_fat_plant <- effects_fat_plant %>%
  mutate(predicted = round(predicted, 4),
         conf.low = round(conf.low, 4),
         conf.high = round(conf.high, 4))

print(effects_fat_plant_df)

# Create a data frame for the plot labels, positioning the labels slightly above the upper confidence interval
fat_df <- data.frame(
  Plant = factor(cld_fat$Plant),  # Ensure matching factor levels
  Label = cld_fat$.group,       # Tukey-adjusted significant difference letters
  y = effects_fat_plant$conf.high + 0.005  # Position just above the upper confidence interval
)


# Adjust the tukey label positions for specific plants
fat_df$y[fat_df$Plant == "EchPur"] <- 0.033
fat_df$y[fat_df$Plant == "RudHir"] <- 0.01
fat_df$y[fat_df$Plant == "BudDav"] <- 0.0035
fat_df$y[fat_df$Plant == "SolAlt"] <- 0.024

# Also make sure the Plant column in fat_df is ordered correctly
fat_df$Plant <- factor(fat_df$Plant)

fat_plot <- ggplot(effects_fat_plant, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_jitter(data = GH_data, aes(x = Plant, y = DryFatMass), width = 0.1, height = 0, size = 1.25, alpha=.25) +  # Raw data points
  geom_text(data = fat_df, aes(x = Plant, y = y, label = Label), vjust = -0.5) +
  scale_x_discrete(labels = c("SymEri" = "S. ericoides", 
                              "BudDav" = "B. davidii", 
                              "EchPur" = "E. purpurea", 
                              "RudHir" = "R. hirta", 
                              "SolAlt" = "S. altissima", 
                              "HelHel" = "H. helianthoides")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12),  # Increase font size for x-axis labels
        axis.text.y = element_text(size = 12),                          # Increase font size for y-axis labels
        axis.title.x = element_text(size = 12),                        # Increase font size for x-axis title
        axis.title.y = element_text(size = 12),                        # Increase font size for y-axis title
        panel.background = element_blank(),                             # Remove background color
        plot.background = element_blank(),                               # Remove plot background color
        panel.border = element_blank(),                                  # Remove panel border
        axis.line = element_line(color = "black")) +                   # Keep axis lines
  labs(x = NULL, y = "Fat mass (g)")  # Remove x-axis title and set y-axis title

print(fat_plot)

tab_model(overall_fat)
emmeans(overall_fat,list(pairwise~Plant), adjust="tukey")
Anova(overall_fat)
summary(overall_fat)

############LEAN MASS PLOTS
overall_lean = lm(DryLeanMass ~ Plant + Sex + ForewingLength + EmergDate + TotalSA, data=GH_data)

#Create a table of the model output in a word doc
tab_model(overall_lean, 
          file = "tables/Lean_Results.doc",   # Export to a Word document
          show.p = TRUE)                   # Add significance stars


effects_lean = ggpredict(overall_lean, ci.lvl = 0.95)
Anova(overall_lean)

# Tukey-adjusted pairwise comparisons for Plant groups
tukey_lean <- emmeans(overall_lean, pairwise ~ Plant, adjust = "tukey")

# Extracting the compact letter display (CLD) to get group letters
cld_lean <- cld(tukey_lean$emmeans, Letters = letters)

# Extract predictions specifically for the Plant variable
effects_lean_plant <- ggpredict(overall_lean, terms = "Plant", ci.lvl = 0.95)


# Ensure Plant is treated as a factor in the effects_lean_plant data
effects_lean_plant$x <- as.factor(effects_lean_plant$x)

# Create a data frame for the plot labels, using the upper confidence interval for positioning
lean_df <- data.frame(
  Plant = effects_lean_plant$x,  # Ensure matching factor levels with the predictions
  Label = cld_lean$.group,       # Tukey-adjusted significant difference letters
  y = effects_lean_plant$conf.high + 0.01  # Position just above the upper confidence interval
)

# Change the order for plant levels to the order from most visited to least visited
effects_lean_plant$x <- factor(effects_lean_plant$x)

# Also make sure the Plant column in lean_df is ordered correctly
lean_df$Plant <- factor(lean_df$Plant)

# Create plot using the specific Plant predictions
lean_plot <- ggplot(effects_lean_plant, show_residuals=TRUE, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_jitter(data = GH_data, aes(x = Plant, y = DryLeanMass), width = 0.1, height = 0, size = 1.25, alpha=.25) +  # Raw data points
  labs(y = "Lean mass (g)", x = "", title = "") +
  my_plot_theme +
  # geom_text(data = lean_df, aes(x = Plant, y = y, label = Label), vjust = -0.5) + # Add significance letters
  # 
  # 1. Change x-axis labels
  scale_x_discrete(labels = c("SymEri" = "S. ericoides", 
                              "BudDav" = "B. davidii", 
                              "EchPur" = "E. purpurea", 
                              "RudHir" = "R. hirta", 
                              "SolAlt" = "S. altissima", 
                              "HelHel" = "H. helianthoides")) +
  
  # 2. Angle the x-axis labels
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12),  # Increase font size for x-axis labels
        axis.text.y = element_text(size = 12),                          # Increase font size for y-axis labels
        axis.title.x = element_text(size = 12),                        # Increase font size for x-axis title
        axis.title.y = element_text(size = 12),                        # Increase font size for y-axis title
        panel.background = element_blank(),                             # Remove background color
        plot.background = element_blank(),                               # Remove plot background color
        panel.border = element_blank(),                                  # Remove panel border
        axis.line = element_line(color = "black"))

# Display plot with updated x-axis labels
lean_plot


################# WATER PLOTS
# Fit the linear model for WaterMass
gh_water_lm <- lm(WaterMass ~ Plant + Sex + ForewingLength + EmergDate + TotalSA, data = GH_data)
vif(gh_water_lm)

# Check model assumptions and perform ANOVA
Anova(gh_water_lm)
check_model(gh_water_lm)

# Tukey-adjusted pairwise comparisons for Plant groups (even if not significant)
tukey_water <- emmeans(gh_water_lm, pairwise ~ Plant, adjust = "tukey")

# Extracting the compact letter display (CLD) to get group letters
cld_water <- cld(tukey_water$emmeans, Letters = letters)

# Extract predictions specifically for the Plant variable
effects_water_plant <- ggpredict(gh_water_lm, terms = "Plant", ci.lvl = 0.95)

# Convert ggpredict object to a data frame
effects_water_plant_df <- as.data.frame(effects_water_plant)

# Round predicted values and confidence intervals
effects_water_plant <- effects_water_plant %>%
  mutate(predicted = round(predicted, 4),
         conf.low = round(conf.low, 4),
         conf.high = round(conf.high, 4))

# # Create a data frame for the plot labels, positioning the labels slightly above the upper confidence interval
# water_df <- data.frame(
#   Plant = factor(cld_water$Plant),  # Ensure matching factor levels
#   Label = cld_water$.group,         # Tukey-adjusted significant difference letters
#   y = effects_water_plant$conf.high + 0.005  # Position just above the upper confidence interval
# )

# Create the water mass plot
water_plot <- ggplot(effects_water_plant, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_jitter(data = GH_data, aes(x = Plant, y = WaterMass), 
              width = 0.1, height = 0, size = 1.5, alpha = 0.25) +  # Raw data points colored by ID
  # geom_text(data = water_df, aes(x = Plant, y = y, label = Label), vjust = -0.5) +
  scale_x_discrete(labels = c("SymEri" = "S. ericoides", 
                              "BudDav" = "B. davidii", 
                              "EchPur" = "E. purpurea", 
                              "RudHir" = "R. hirta", 
                              "SolAlt" = "S. altissima", 
                              "HelHel" = "H. helianthoides")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12),  # Increase font size for x-axis labels
        axis.text.y = element_text(size = 12),                          # Increase font size for y-axis labels
        axis.title.x = element_text(size = 12),                        # Increase font size for x-axis title
        axis.title.y = element_text(size = 12),                        # Increase font size for y-axis title
        panel.background = element_blank(),                             # Remove background color
        plot.background = element_blank(),                               # Remove plot background color
        panel.border = element_blank(),                                  # Remove panel border
        axis.line = element_line(color = "black")) +                   # Keep axis lines
  labs(x = NULL, y = "Water mass (g)") +  # Remove x-axis title and set y-axis title
  scale_color_viridis_d()  # Color scale for ID

# Print the water plot
print(water_plot)

tab_water = lm(WaterMass ~ Plant + Sex + ForewingLength + EmergDate + TotalSA, data=GH_data)
tab_model(tab_water)
Anova(tab_water)


multi_panel_plot <- (weight_plot + fat_plot) / (lean_plot + water_plot) + 
  plot_annotation(tag_levels = 'a')&
  theme(plot.tag.position = c(0.05, 1))

# Display the combined plot
multi_panel_plot

ggsave("plots/comboplot_weight+BodyComp.png", multi_panel_plot, width = 27, height = 25, units = "cm", dpi = 300,)

############### GH VS FIELD
ghf_weight = lmer(Weight~ExpLoc+TotalSA+OrigWeight+ForewingLength+Sex+TrialDay+(1|ID),data=gh_field_weight)
summary(ghf_weight)
check_model(ghf_weight)
tab_model(ghf_weight)
effects_ghf = ggpredict(ghf_weight)
weight_loc = plot(effects_ghf$ExpLoc, rawdata=TRUE,jitter=c(0.35,0), alpha=.40) +
  labs(y = "Weight (g)",x = "Location",title="") 
weight_loc
Anova(ghf_weight)
summary(ghf_weight)
summary(gh_field_weight)
ghf_fat = lm(DryFatMass~ExpLoc+ForewingLength, data=GH_F_data)
Anova(ghf_fat)
check_model(ghf_fat)
ghf_fat_pred = ggpredict(ghf_fat)
ghf_plot = plot(ghf_fat_pred)

ghf_plot = plot(ghf_fat_pred$ExpLoc) +
  labs(y = "Fat (g)",x = "Location",title="") 

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
                   panel.grixd.minor = element_blank(), axis.line = element_line(color = "black"))

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

