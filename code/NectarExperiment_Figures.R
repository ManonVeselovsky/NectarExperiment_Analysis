### Statistical analysis script
### Written by M. Veselovsky
### Last modified 2024-09-23
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

desired_order <- c("SymEri", "SolAlt", "BudDav", "EchPur", "RudHir", "HelHel")

# Rename "GH" to "Greenhouse" in the ExpLoc variable
GH_F_data$ExpLoc[GH_F_data$ExpLoc == "GH"] <- "Greenhouse"
gh_field_weight$ExpLoc[gh_field_weight$ExpLoc == "GH"] <- "Greenhouse"

# Convert ExpLoc to a factor with proper levels
GH_F_data$ExpLoc <- factor(GH_F_data$ExpLoc, levels = c("Field", "Greenhouse"))

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

##############################################################################################################################
############################ FINAL MODELS FOR EACH VARIABLE ##################################################

########## 0. General statistics

forewing_range <- range(GH_data$ForewingLength, na.rm = TRUE) #get the range of forewing lengths
forewing_range
mean_forewing = mean(summarydat$ForewingLength) #get the mean forewing length
mean_forewing

forewing_males <- mean(summarydat$ForewingLength[summarydat$Sex == "M"], na.rm = TRUE)
forewing_males

forewing_females <- mean(summarydat$ForewingLength[summarydat$Sex == "F"], na.rm = TRUE)
forewing_females
############################################################################################################
############################ 1. Greenhouse Plants ########################################################
######### 1.1 WEIGHT EFFECTS
mean_weight <- mean(gh_TreatmentData$Weight[gh_TreatmentData$TrialDay != 0], na.rm = TRUE)
mean_weight

# Calculate standard error
std_error_weight_day7 <- sd(gh_TreatmentData$Weight[gh_TreatmentData$TrialDay == 7], na.rm = TRUE) / 
  sqrt(sum(!is.na(gh_TreatmentData$Weight[gh_TreatmentData$TrialDay == 7])))


# Display results
mean_weight
std_error_weight_day7


#Remove TrialDay = 0 because I already have an original weight (OrigWeight) column which is the TrialDay = 0 weight
alt_data = gh_TreatmentData %>% 
  filter(TrialDay %in% c("1", "5","7", "10", "11"))
overall_weight = lmer(Weight~Plant+TotalSA+OrigWeight+Sex+ForewingLength+TrialDay+DateWeighed+EnclCol+(1|ID),data=alt_data)

#Create a table of the model output in a word doc
tab_model(overall_weight, 
          file = "tables/Weight_Model_Results.doc",   # Export to a Word document
          show.p = TRUE, # Add significance stars
          show.stat = TRUE)                   

summary(overall_weight)

check_model(overall_weight) #Assumptions look good, weird error bars with TrialDay VIF
#Test for collinearity among predictor variables (GVIF<4 for continuous acceptable,
# GVIF^(1/(2*df)) < 2 acceptable for categorical)
vif(overall_weight) #no problem here

Anova(overall_weight) #Plant, FWL, Start weight, Sex, and Trial day all significant
summary(overall_weight)

# Tukey-adjusted pairwise comparisons for Plant groups
tukey_weight <- emmeans(overall_weight, pairwise ~ Plant, adjust = "tukey")

tukey_weight
# Extracting the compact letter display (CLD) to get group letters
cld_weight <- cld(tukey_weight$emmeans, Letters = letters)

# Extract predictions specifically for the Plant variable
effects_weight_plant <- ggpredict(overall_weight, terms = "Plant", ci.lvl = 0.95)

# Ensure Plant is treated as a factor in the effects_weight_plant data
effects_weight_plant$x <- as.factor(effects_weight_plant$x)
# # Reorder the Plant levels in cld_weight to match effects_weight_plant$x
# cld_weight$Plant <- factor(cld_weight$Plant, levels = levels(effects_weight_plant$x))

# Ensure the Tukey-adjusted letters are ordered according to Plant levels
cld_weight <- cld_weight[order(cld_weight$Plant),]

# Now recreate the labels_df with the correctly ordered Plant factor
labels_df <- data.frame(
  Plant = effects_weight_plant$x,   # Correct factor levels from effects_weight_plant
  Label = cld_weight$.group,        # Use reordered Tukey letters
  y = effects_weight_plant$conf.high + 0.01  # Position just above the upper confidence interval
)

# Create plot using the specific Plant predictions with adjusted legend and angled x-axis labels
weight_plot <- ggplot(effects_weight_plant, show_residuals=TRUE, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_jitter(data = alt_data, aes(x = Plant, y = Weight, color = as.factor(TrialDay)), 
              width = 0.1, height = 0, alpha = 0.6) +  # Add raw data points colored by TrialDay
  labs(y = "Wet weight (g)", x = "", title = "", color = "Trial Day") +  # Label legend as "Trial Day"
  my_plot_theme +
  geom_text(data = labels_df, aes(x = Plant, y = y, label = Label), vjust = -0.5) +  # Add significance letters
  scale_x_discrete(expand = c(0.1, 0.1)) +  # Reduce space between groups
  
  # 1. Change x-axis labels
  scale_x_discrete(labels = c("SymEri" = "S. ericoides", 
                              "BudDav" = "B. davidii", 
                              "EchPur" = "E. purpurea", 
                              "RudHir" = "R. hirta", 
                              "SolAlt" = "S. altissima", 
                              "HelHel" = "H. helianthoides")) +
  
  # 2. Angle x-axis labels
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12),  # Angle the x-axis labels
        axis.text.y = element_text(size = 12),                          # Increase font size for y-axis labels
        axis.title.x = element_text(size = 12),                        # Increase font size for x-axis title
        axis.title.y = element_text(size = 12),                        # Increase font size for y-axis title
        panel.background = element_blank(),                             # Remove background color
        plot.background = element_blank(),                               # Remove plot background color
        panel.border = element_blank(),                                  # Remove panel border
        axis.line = element_line(color = "black"),                       # Keep axis lines
        legend.position = c(0.97, 0.88),  # Adjust legend position to avoid overlap with data
        legend.background = element_blank(),                            # Remove background color for legend
        legend.box.background = element_blank())                        # Remove border around the legend box

# Display plot with updated x-axis labels, adjusted legend, and significance letters
weight_plot


summary(overall_weight)

tukey_p_values <- summary(tukey_weight$emmeans)$p.value
tukey_p_values_df <- data.frame(Plant = names(tukey_p_values), Tukey_p_value = tukey_p_values)

# Display the model output with tab_model
tab_model(overall_weight)




#####################################################################################################
############################ 1.2 FAT EFFECTS ########################################################

#Get the mean fat mass of butterflies
mean_fat <- mean(GH_data$DryFatMass, na.rm = TRUE)

# Calculate standard error
se_fat <- sd(GH_data$DryFatMass, na.rm = TRUE) / 
  sqrt(sum(!is.na(GH_data$DryFatMass)))

mean_fat
se_fat

# Display results
mean_weight_day7
std_error_weight_day7


# Fit the linear model
overall_fat <- lm(log(DryFatMass) ~ Plant + Sex + ForewingLength + EmergDate + TotalSA, data = GH_data)

e_date = lm(DryFatMass~Plant + EmergDate,GH_data)
Anova(e_date)

fat_alt = lm(DryFatMass~Plant + Sex + ForewingLength + EmergDate, data = GH_data)

Anova(fat_alt)

#Create a table of the model output in a word doc
tab_model(overall_fat, 
          file = "tables/Fat_Results.doc",   # Export to a Word document
          show.p = TRUE,                   # Add significance stars
          show.stat = TRUE
)

# Tukey-adjusted pairwise comparisons for Plant groups
tukey_fat <- emmeans(overall_fat, pairwise ~ Plant, adjust = "tukey")

tukey_fat
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

# Reorder the Plant levels in cld_fat to match effects_fat_plant$x
cld_fat$Plant <- factor(cld_fat$Plant, levels = levels(effects_fat_plant$x))

# Reorder cld_fat based on the Plant levels
cld_fat <- cld_fat[order(cld_fat$Plant),]

# Now recreate the labels_df with the correctly ordered Plant factor and Tukey letters
fat_df <- data.frame(
  Plant = effects_fat_plant$x,   # Correct factor levels from effects_fat_plant
  Label = cld_fat$.group,        # Use reordered Tukey letters
  y = effects_fat_plant$conf.high + 0.001  # Position just above the upper confidence interval
)


fat_plot <- ggplot(effects_fat_plant, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_jitter(data = GH_data, aes(x = Plant, y = DryFatMass), width = 0.095, height = 0, size = 1.25, alpha=.25) +  # Raw data points
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



#####################################################################################################
############################ 1.3 LEAN MASS EFFECTS ########################################################
overall_lean = lm(DryLeanMass ~ Plant + Sex + ForewingLength + EmergDate + TotalSA, data=GH_data)

#Create a table of the model output in a word doc
tab_model(overall_lean, 
          file = "tables/Lean_Results.doc",   # Export to a Word document
          show.p = TRUE, # Add significance stars
          show.stat = TRUE)                   


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
  geom_jitter(data = GH_data, aes(x = Plant, y = DryLeanMass), width = 0.09, height = 0, size = 1.25, alpha=.25) +  # Raw data points
  labs(y = "Lean mass (g)", x = "", title = "") +
  my_plot_theme +
  # geom_text(data = lean_df, aes(x = Plant, y = y, label = Label), vjust = -0.5) + # Add significance letters
  # 
  # 1. Change x-axis labels
  scale_x_discrete(labels = c("SymEri" = "Symphyotrichum ericoides", 
                              "BudDav" = "Buddleja davidii", 
                              "EchPur" = "Echinacea purpurea", 
                              "RudHir" = "Rudbeckia hirta", 
                              "SolAlt" = "Solidago altissima", 
                              "HelHel" = "Heliopsis helianthoides")) +
  
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
Anova(overall_lean)

# Round predicted values and confidence intervals to 4 decimal places
effects_fat_plant <- effects_fat_plant %>%
  mutate(predicted = round(predicted, 4),
         conf.low = round(conf.low, 4),
         conf.high = round(conf.high, 4))

print(effects_lean_plant)



#####################################################################################################
############################ 1.4 WATER EFFECTS ########################################################
# Fit the linear model for WaterMass
overall_water <- lm(WaterMass ~ Plant + Sex + ForewingLength + EmergDate + TotalSA, data = GH_data)
vif(overall_water)

# Check model assumptions and perform ANOVA
Anova(overall_water)
check_model(overall_water)

# Tukey-adjusted pairwise comparisons for Plant groups (even if not significant)
tukey_water <- emmeans(overall_water, pairwise ~ Plant, adjust = "tukey")

# Extracting the compact letter display (CLD) to get group letters
cld_water <- cld(tukey_water$emmeans, Letters = letters)

# Extract predictions specifically for the Plant variable
effects_water_plant <- ggpredict(overall_water, terms = "Plant", ci.lvl = 0.95)

# Convert ggpredict object to a data frame
effects_water_plant_df <- as.data.frame(effects_water_plant)

# Round predicted values and confidence intervals
effects_water_plant <- effects_water_plant %>%
  mutate(predicted = round(predicted, 4),
         conf.low = round(conf.low, 4),
         conf.high = round(conf.high, 4))

# # Create a data frame for the plot labels, positioning the labels slightly above the upper confidence interval
water_df <- data.frame(
  Plant = factor(cld_water$Plant),  # Ensure matching factor levels
  Label = cld_water$.group,         # Tukey-adjusted significant difference letters
  y = effects_water_plant$conf.high + 0.005  # Position just above the upper confidence interval
)

# Create the water mass plot
water_plot <- ggplot(effects_water_plant, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_jitter(data = GH_data, aes(x = Plant, y = WaterMass), 
              width = 0.095, height = 0, size = 1.5, alpha = 0.25) +  # Raw data points
  # geom_text(data = water_df, aes(x = Plant, y = y, label = Label), vjust = -0.5) +
  scale_x_discrete(labels = c("SymEri" = "Symphyotrichum ericoides", 
                              "BudDav" = "Buddleja davidii", 
                              "EchPur" = "Echinacea purpurea", 
                              "RudHir" = "Rudbeckia hirta", 
                              "SolAlt" = "Solidago altissima", 
                              "HelHel" = "Heliopsis helianthoides")) +
  theme(axis.text.x = element_text(angle = 15, hjust = 0.5, size = 12),  # Increase font size for x-axis labels
        axis.text.y = element_text(size = 12),                          # Increase font size for y-axis labels
        axis.title.x = element_text(size = 12),                        # Increase font size for x-axis title
        axis.title.y = element_text(size = 12),                        # Increase font size for y-axis title
        panel.background = element_blank(),                             # Remove background color
        plot.background = element_blank(),                               # Remove plot background color
        panel.border = element_blank(),                                  # Remove panel border
        axis.line = element_line(color = "black")) +                   # Keep axis lines
  labs(x = NULL, y = "Water mass (g)") +  # Remove x-axis title and set y-axis title
  scale_color_viridis_d()  # Color scale for ID
  


# Print the water ploteffects_water_plant
print(water_plot)

tab_model(overall_water, 
          file = "tables/Water_Results.doc",   # Export to a Word document
          show.p = TRUE, # Add significance stars
          show.stat = TRUE) 
Anova(overall_water)
summary(overall_water)

plot(allEffects(overall_water))




# Remove x-axis labels from the upper panels and legend from weight_plot
weight_panel <- weight_plot + theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  legend.position = "none")  # Remove the legend

fat_panel <- fat_plot + theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank()
  )

# Create the multi-panel figure
gh_multi_panel_plot <- (weight_panel + fat_panel) / (lean_plot + water_plot) + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag.position = c(0.05, 1))  # Adjust the position of the plot tags

# Display the multi-panel plot
gh_multi_panel_plot

#save the multi-panel plot
ggsave("plots/comboplot_weight+BodyComp.png", gh_multi_panel_plot, width = 27, height = 25, units = "cm", dpi = 300,)



#####################################################################################################
############################## 2. GH VS FIELD ########################################################
############################ 2.1 WEIGHT ########################################################

ghf_weight = lmer(Weight~ExpLoc+TotalSA+OrigWeight+ForewingLength+Sex+TrialDay+(1|ID),data=gh_field_weight)

# Extract the model summary
summary_model <- summary(ghf_weight)

# Extract the p-value for ExpLoc directly from the fixed effects
p_value_exploc <- summary_model$coefficients["ExpLocGreenhouse", "Pr(>|t|)"]

# Create a label for significance based on the p-value
significance_label <- ifelse(p_value_exploc < 0.05, "*", "ns")

# Extract predictions for ExpLoc (without Tukey adjustment)
effects_weight_exploc <- ggpredict(ghf_weight, terms = "ExpLoc", ci.lvl = 0.95)

# Ensure ExpLoc is treated as a factor
effects_weight_exploc$x <- factor(effects_weight_exploc$x)

# Create a label dataframe for positioning significance label
labels_df_exploc <- data.frame(
  ExpLoc = effects_weight_exploc$x,
  Label = significance_label,
  y = effects_weight_exploc$conf.high + 0.01  # Position just above upper confidence interval
)

# Create the plot for ExpLoc with adjusted legend and no border
exploc_weight_plot <- ggplot(effects_weight_exploc, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_jitter(data = gh_field_weight, aes(x = ExpLoc, y = Weight, color = as.factor(TrialDay)), 
              0.095, height = 0, alpha = 0.6) +  # Add raw data points colored by TrialDay
  labs(y = "Weight (g)", x = "", title = "", color = "Trial Day") +  # Label legend as "Trial Day"
  my_plot_theme +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),  # Center x-axis labels
        axis.text.y = element_text(size = 12),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = c(0.85, 0.8),  # Move the legend slightly to the right
        legend.background = element_blank(),
        legend.box.background = element_blank())  # Remove the black border around the legend box

# Display the plot
exploc_weight_plot






summary(ghf_weight)
check_model(ghf_weight)
tab_model(ghf_weight)
effects_ghf = ggpredict(ghf_weight)
weight_loc = plot(effects_ghf$ExpLoc, rawdata=TRUE,jitter=c(0.35,0), alpha=.40) +
  labs(y = "Weight (g)",x = "",title="") 
weight_loc
Anova(ghf_weight)
summary(ghf_weight)
summary(gh_field_weight)


# Extract predictions for ExpLoc from ghf_weight
effects_ghf_weight = ggpredict(ghf_weight, terms = "ExpLoc", ci.lvl = 0.95)

# Ensure ExpLoc is treated as a factor
effects_ghf_weight$x <- as.factor(effects_ghf_weight$x)

# Create plot for Weight by Location with updated ExpLoc labels
weight_loc_plot <- ggplot(effects_ghf_weight, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(y = "Weight (g)", x = "", title = "") +
  # geom_jitter(data = gh_field_weight, aes(x = ExpLoc, y = Weight), 
  #             width = 0.35, alpha = 0.4, size = 1.5) + # Add raw data points with jitter
  theme(axis.text.x = element_text(size = 14),  # Increase x-axis label size
        # axis.text.x = element_text(angle = 30, hjust = 0.5, size = 12),  # Angle x-axis labels
        axis.text.y = element_text(size = 14),                          # Y-axis labels
        axis.title.x = element_text(size = 14),                         # X-axis title
        axis.title.y = element_text(size = 14),                         # Y-axis title
        panel.background = element_blank(),                             # Remove background
        plot.background = element_blank(),                              # Remove plot background
        panel.border = element_blank(),                                 # Remove panel border
        axis.line = element_line(color = "black"))                      # Keep axis lines

# Display the plot
weight_loc_plot


#####################################################################################################
############################ 2.2 GH-Field Fat comparison ############################

# Fit the model
ghf_fat <- lm(DryFatMass ~ ExpLoc + ForewingLength, data = GH_F_data)

summary(GH_F_data)
summary(ghf_fat)
# Check assumptions and model summary
check_model(ghf_fat)
summary(ghf_fat)


# Tukey-adjusted pairwise comparisons for ExpLoc
tukey_fat <- emmeans(ghf_fat, pairwise ~ ExpLoc, adjust = "tukey")
cld_fat <- cld(tukey_fat$emmeans, Letters = letters)


# Generate predictions for ExpLoc variable
gh_fat_pred_exploc <- ggpredict(ghf_fat, terms = "ExpLoc")

# Convert predictions to a data frame for plotting
effects_ghf_fat_df <- as.data.frame(gh_fat_pred_exploc)

# Create a data frame for significance letters with adjusted position
labels_df <- data.frame(
  ExpLoc = cld_fat$ExpLoc,
  Label = cld_fat$.group,
  y = c(
    max(effects_ghf_fat_df$conf.high[effects_ghf_fat_df$x == "Greenhouse"] + 0.002),  # Position for Greenhouse
    max(effects_ghf_fat_df$conf.high[effects_ghf_fat_df$x == "Field"] + 0.002)  # Position for Field
  )
)

# Create the plot for Fat by Location with significance letters
fat_loc_plot <- ggplot(effects_ghf_fat_df, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(y = "Fat (g)", x = "", title = "") +
  geom_jitter(data = GH_F_data, aes(x = ExpLoc, y = DryFatMass), 
              width = 0.35, alpha = 0.4, size = 1.5) +  # Add raw data points with jitter
  geom_text(data = labels_df, aes(x = ExpLoc, y = y, label = Label), vjust = -0.5, size = 5) +  # Add significance letters
  theme(axis.text.x = element_text(hjust = 0.5, size = 16),  # Increase x-axis label size
        axis.text.y = element_text(size = 12),                          # Y-axis labels
        axis.title.x = element_text(size = 12),                         # X-axis title
        axis.title.y = element_text(size = 12),                         # Y-axis title
        panel.background = element_blank(),                             # Remove background
        plot.background = element_blank(),                              # Remove plot background
        panel.border = element_blank(),                                 # Remove panel border
        axis.line = element_line(color = "black"))                      # Keep axis lines

# Display the plot
fat_loc_plot


gh_f_multi_panel_plot <- (weight_loc_plot + fat_loc_plot) + 
  plot_annotation(tag_levels = 'a')&
  theme(plot.tag.position = c(0.05, 1))

# Display the combined plot
gh_f_multi_panel_plot

ggsave("plots/comboplot_GH_Field.png", gh_f_multi_panel_plot, width = 27, height = 25, units = "cm", dpi = 300,)


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
p
#remove legend 
p_nolegend = p + theme(legend.position="none")
p_nolegend
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


