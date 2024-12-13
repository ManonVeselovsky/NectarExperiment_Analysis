### Statistical analysis script
### Written by M. Veselovsky
### Last modified 2024-12-10
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
library(officer)
library(knitr)
library(ggeffects)
library(emmeans)
library(broom)
library(tidyverse) #re-order categorical variables in effects plots
library(dplyr)
library(performance)
library(multcomp) #for adding multcomp tukey to summary
library(patchwork) #for creating a multi-panel plot
library(lubridate)
#raw treatment data with all trial days
TreatmentData<-read.csv("processed/TreatmentData_p.csv")
summarydat<-read.csv("processed/summarydat.csv")
tempdat = read.csv("processed/tempdat.csv")

# create separate working databases for greenhouse plants, goldenrod (for GH-Field comparison), and field plants
GH_data = subset(summarydat, summarydat$ExpLoc == "GH") # data from greenhouse
GH_F_data = summarydat[which(summarydat$Plant=="SolAlt"),] #all data on goldenrod, for GH-Field comparison
Field_data = subset(summarydat,summarydat$ExpLoc == "Field") #data from field for fat
SolAlt_F_data = subset(Field_data,Field_data$Plant == "SolAlt") #Subset SolAlt in field to check for enclosure effects
gh_TreatmentData = subset(TreatmentData, TreatmentData$ExpLoc == "GH")
field_TreatmentData = subset(TreatmentData, TreatmentData$ExpLoc == "Field")
gh_field_weight = TreatmentData[which(TreatmentData$Plant=="SolAlt"),]

desired_order <- c("SymEri", "SolAlt", "BudDav", "EchPur", "RudHir", "HelHel")


# Check for duplicates
GH_data %>%
  filter(ID == "MVD47")  # Filter rows for ID MVD47

# Remove the first occurrence of MVD47 as the second one is the correct entry
GH_data <- GH_data %>%
  group_by(ID) %>%
  filter(!(ID == "MVD47" & row_number() == 1)) %>%
  ungroup()


# Reorder the Plant factor in the original data
GH_data$Plant <- factor(GH_data$Plant, levels = desired_order)
gh_TreatmentData$Plant = factor(gh_TreatmentData$Plant, levels = desired_order)

# Rename "GH" to "Greenhouse" in the ExpLoc variable
GH_F_data$ExpLoc[GH_F_data$ExpLoc == "GH"] <- "Greenhouse"
gh_field_weight$ExpLoc[gh_field_weight$ExpLoc == "GH"] <- "Greenhouse"

# Convert ExpLoc to a factor with proper levels
GH_F_data$ExpLoc <- factor(GH_F_data$ExpLoc, levels = c("Field", "Greenhouse"))

Field_data$Plant[Field_data$Plant == "EutMac"] = "E. maculatum"
Field_data$Plant[Field_data$Plant == "SolAlt"] = "S. altissima"

field_TreatmentData$Plant[field_TreatmentData$Plant == "EutMac"] = "E. maculatum"
field_TreatmentData$Plant[field_TreatmentData$Plant == "SolAlt"] = "S. altissima"

## Create a theme to use for plots
my_plot_theme = theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),  # Standardize axis line thickness
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
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



butterfly_count_interact_weight <- gh_TreatmentData %>%
  group_by(Plant) %>%
  summarise(TotalButterflies = n_distinct(ID))

print(butterfly_count_interact_weight)

# Assuming gh_TreatmentData has 'Plant' and 'DateWeighed' columns (as Julian date)
plant_plot = ggplot(gh_TreatmentData, aes(x = DateWeighed, y = Plant, color = Plant)) +
  geom_jitter(width = 0, height = 0.2, alpha = 0.6) +  # Adjust jitter to avoid overlap
  labs(
    title = "",
    x = "Julian Date of Weight Measurement",
    y = "Plant Species"
  ) +
  scale_x_continuous(
    breaks = seq(min(gh_TreatmentData$DateWeighed), max(gh_TreatmentData$DateWeighed), by = 5),  # Set tick marks every 5 Julian days
    labels = scales::number_format()  # Keeps Julian format
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    panel.grid = element_blank(),  # Remove gridlines
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.ticks = element_line(color = "black")  # Add tick marks
  )
ggsave("plots/s3_PlantOverlap.jpeg", plant_plot, width = 27, height = 25, units = "cm", dpi = 300,)


############################################################################################################
############################ 1. Greenhouse Plants ########################################################
######### 1.1 WEIGHT EFFECTS
mean_weight <- mean(gh_TreatmentData$Weight[gh_TreatmentData$TrialDay != 0], na.rm = TRUE)
mean_weight
count(gh_TreatmentData)

# Calculate standard error
std_error_weight_day7 <- sd(gh_TreatmentData$Weight[gh_TreatmentData$TrialDay == 7], na.rm = TRUE) / 
  sqrt(sum(!is.na(gh_TreatmentData$Weight[gh_TreatmentData$TrialDay == 7])))


# Display results
mean_weight
std_error_weight_day7

# Create the plot
ggplot(gh_TreatmentData, aes(x = TrialDay, y = Weight, group = ID, color = as.factor(ID))) +
  geom_line(alpha = 0.7, size = 1.2) + # Individual butterfly trajectories
  geom_smooth(aes(group = 1), method = "loess", se = TRUE, color = "black", size = 1.5) + # Overall trend
  labs(
    title = "Individual Butterfly Weight Trajectories with Overall Trend",
    x = "Trial Day",
    y = "Butterfly Weight (mg)",
    color = "Butterfly ID"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "gray", size = 0.5),
    panel.grid.minor = element_blank()
  )

# overall_weight = lmer(Weight~Plant+TotalSA+ForewingLength+OrigWeight+Sex+TrialDay+DateWeighed+EnclCol+(1|ID),data=gh_TreatmentData)
 
interact_weight = lmer(Weight~Plant*TrialDay+TotalSA+ForewingLength+Sex+EnclCol+(1|ID),data=gh_TreatmentData)

# Extract the data used in the model
model_data_interact_weight <- interact_weight@frame

# Count the number of unique butterflies tested on each plant
butterfly_count_by_plant <- model_data_interact_weight %>%
  group_by(Plant) %>%
  summarise(TotalButterflies = n_distinct(ID))

print(butterfly_count_by_plant)


check_model(interact_weight) #Trialday*Plant interaction shows multicollinearity but this is inherent with an interaction term
# Use the model without interaction terms to check multicollinearity
collinearity_check = lmer(Weight~Plant + TrialDay+TotalSA+ForewingLength+Sex+EnclCol+(1|ID),data=gh_TreatmentData)
check_model(collinearity_check) #no issues with collinearity or other 

# Step 1: Get emmeans for the interaction at TrialDay = 7
interaction_emmeans <- emmeans(interact_weight, ~ Plant * TrialDay, 
                               at = list(TrialDay = 7))  # Fix TrialDay at 7
tukey_weight = emmeans(interact_weight, ~ Plant * TrialDay, 
                       at = list(TrialDay = 7))  # Fix TrialDay at 7

tukey_weight

# Step 2: Perform pairwise comparisons for Plant (Tukey-adjusted)
plant_comparisons <- contrast(interaction_emmeans, 
                              method = "pairwise", 
                              simple = "each",  # Compare within each TrialDay
                              combine = FALSE,  # Only focus on Plant comparisons
                              adjust = "tukey")  # Force Tukey adjustment

# View pairwise results
print(plant_comparisons)

# Step 3: Generate compact letter display for Plant at TrialDay = 7
plant_cld <- cld(interaction_emmeans, 
                 by = "TrialDay",   # Generate letters per TrialDay
                 adjust = "tukey",  # It is forcing a sidak adjustment
                 #but the CLD for sidak matches the data from the tukey adjustment so I am running with this
                 Letters = letters)

# View the compact letter display
print(plant_cld)

# Step 4: Prepare the emmeans data frame for plotting
interaction_emmeans_df <- as.data.frame(interaction_emmeans) %>%
  filter(TrialDay == 7) %>%  # Only keep TrialDay = 7
  mutate(.group = plant_cld$.group)  # Add significance letters

# Step 5: Plot estimated means and raw data
interaction_plot <- ggplot(interaction_emmeans_df, aes(x = Plant, y = emmean)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(width = 0.5)) +
  # Add raw data points for TrialDay = 7
  geom_jitter(data = subset(gh_TreatmentData, TrialDay == 7), 
              aes(x = Plant, y = Weight), width = 0.095, height = 0, size = 1.25, alpha=.25) +
  # Add significance letters above error bars
  # Add significance letters above the upper confidence intervals
  geom_text(aes(label = .group, y = upper.CL + 0.02), 
            position = position_dodge(width = 0.5), size = 5) +
  labs(x = "", y = "Weight (g)", title = "") +
  theme_minimal() +
  my_plot_theme +
  scale_x_discrete(labels = c("SymEri" = "S. ericoides", 
                              "BudDav" = "B. davidii", 
                              "EchPur" = "E. purpurea", 
                              "RudHir" = "R. hirta", 
                              "SolAlt" = "S. altissima", 
                              "HelHel" = "H. helianthoides"))

print(interaction_plot)


Anova(interact_weight, type = 3)
summary(interact_weight)

#Create a table of the model output in a word doc
tab_model(interact_weight, 
          file = "tables/Weight_Model_Results.doc",   # Export to a Word document
          show.p = TRUE, # Add significance stars
          show.stat = TRUE,
          show.se = TRUE,
          show.df = TRUE)                   

plot(ggpredict(interact_weight))

#####################################################################################################
############################ 1.2 FAT EFFECTS ########################################################

# Fit the linear model
overall_fat <- lm(log(DryFatMass) ~ Plant + Sex + ForewingLength +RawWeight_day0 + TotalSA, data = GH_data)
normal_fat = lm(DryFatMass ~ Plant + Sex + ForewingLength +RawWeight_day0 + TotalSA, data = GH_data)
plot(allEffects(normal_fat))

check_model(overall_fat)
# Get the model data from the overall_fat model
model_data_overall_fat <- model.frame(overall_fat)

# Count the number of butterflies tested on each plant (by counting rows for each Plant)
butterfly_count_by_plant_overall_fat <- model_data_overall_fat %>%
  group_by(Plant) %>%
  summarise(TotalButterflies = n())
plot(allEffects(overall_fat))
print(butterfly_count_by_plant_overall_fat)

Anova(overall_fat)
#Create a table of the model output in a word doc
tab_model(overall_fat, 
          file = "tables/Fat_Results.doc",   # Export to a Word document
          show.p = TRUE,                   
          show.stat = TRUE,
          show.se = TRUE,
          show.df = TRUE
          )

# Tukey-adjusted pairwise comparisons for Plant groups
tukey_fat <- emmeans(overall_fat, pairwise ~ Plant, adjust = "tukey")

tukey_contrasts <- tukey_fat$contrasts #Extract the contrast letters from the mult-comp

# tidy_fat <- tidy(tukey_contrasts)
# tidy_fat <- tidy_fat %>%
#   mutate(across(where(is.numeric), round, 3))


# Create a Word document and add the table as a proper Word table
my_doc <- read_docx() %>%
  body_add_par("Tukey-adjusted pairwise comparisons for Plant groups") %>%
  body_add_table(value = tidy_fat) %>%  # Add the dataframe directly as a table
  print(target = "tables/Tukey_Contrasts_fat.docx")  # Export the document to the target path

tukey_weight

tukey_fat
# Extracting the compact letter display (CLD) to get group letters
cld_fat <- cld(tukey_fat$emmeans, Letters = letters)

# Extract predictions specifically for the Plant variable
effects_fat_plant <- ggpredict(overall_fat, terms = "Plant", ci.lvl = 0.95)

plot(ggpredict(overall_fat, terms = "TotalSA",ci.lvl=0.95)+ggplot(gh))

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
  geom_text(data = fat_df, aes(x = Plant, y = y, label = Label), vjust = -0.5, size = 5) +
  scale_x_discrete(labels = c("SymEri" = "S. ericoides", 
                              "BudDav" = "B. davidii", 
                              "EchPur" = "E. purpurea", 
                              "RudHir" = "R. hirta", 
                              "SolAlt" = "S. altissima", 
                              "HelHel" = "H. helianthoides")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 13),  # Increase font size for x-axis labels
        axis.text.y = element_text(size = 13),                          # Increase font size for y-axis labels
        axis.title.x = element_text(size = 13),                        # Increase font size for x-axis title
        axis.title.y = element_text(size = 13),                        # Increase font size for y-axis title
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
check_model(overall_fat)
vif(overall_fat)
sum(overall_fat$model$Sex == "M")
sum(overall_fat$model$Sex == "F")


#####################################################################################################
############################ 1.3 LEAN MASS EFFECTS ########################################################
overall_lean = lm(DryLeanMass ~ Plant + Sex + ForewingLength + RawWeight_day0 + TotalSA, data=GH_data)

#Create a table of the model output in a word doc
tab_model(overall_lean, 
          file = "tables/Lean_Results.doc",   # Export to a Word document
          show.p = TRUE, # Add significance stars
          show.stat = TRUE,
          show.se = TRUE,
          show.df = TRUE)                   

check_model(overall_lean)

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

# Define the correct plant order
plant_levels <- c("SymEri", "BudDav", "EchPur", "RudHir", "SolAlt", "HelHel")

# Ensure Plant is treated as a factor in the effects_lean_plant data
effects_lean_plant$x <- factor(effects_lean_plant$x, levels = c("SymEri", "SolAlt", "BudDav", "EchPur", "RudHir", "HelHel"))
cld_lean$Plant <- factor(cld_lean$Plant, levels = c("SymEri", "SolAlt", "BudDav", "EchPur", "RudHir", "HelHel"))

# Merge cld_lean with effects_lean_plant to ensure letters align with predictions
lean_df <- merge(
  data.frame(Plant = effects_lean_plant$x, y = effects_lean_plant$conf.high + 0.01),
  data.frame(Plant = cld_lean$Plant, Label = cld_lean$.group),
  by = "Plant"
)

# Ensure GH_data has the correct levels for comparison
GH_data$Plant <- factor(GH_data$Plant, levels = c("SymEri", "SolAlt", "BudDav", "EchPur", "RudHir", "HelHel"))

# Create the plot with consistent theming and style
lean_plot <- ggplot(effects_lean_plant, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_jitter(data = GH_data, aes(x = Plant, y = DryLeanMass), width = 0.09, height = 0, size = 1.25, alpha = .25) +  # Raw data points
  labs(y = "Lean mass (g)", x = "", title = "") +
  my_plot_theme +
  # Add the Tukey significance letters correctly aligned with the factor levels
  geom_text(data = lean_df, aes(x = Plant, y = y, label = Label), size = 5, vjust = -0.5) +
  
  # Ensure the x-axis has the desired plant name order
  scale_x_discrete(labels = c(
    "SymEri" = "Symphyotrichum ericoides", 
    "SolAlt" = "Solidago altissima", 
    "BudDav" = "Buddleja davidii", 
    "EchPur" = "Echinacea purpurea", 
    "RudHir" = "Rudbeckia hirta", 
    "HelHel" = "Heliopsis helianthoides"
  )) +
  
  # Adjust theme for consistent look
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 20, unit = "mm"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 13),  # Adjust x-axis text
        axis.text.y = element_text(size = 13),                          # Adjust y-axis text
        axis.title.x = element_text(size = 13),                        # Adjust x-axis title
        axis.title.y = element_text(size = 13),                        # Adjust y-axis title
        panel.background = element_blank(),                             # Remove panel background
        plot.background = element_blank(),                               # Remove plot background
        panel.border = element_blank(),                                  # Remove borders
        axis.line = element_line(color = "black"))




# Display the plot
lean_plot
Anova(overall_lean)

# Tukey-adjusted pairwise comparisons for Plant groups (even if not significant)
tukey_lean <- emmeans(overall_lean, pairwise ~ Plant, adjust = "tukey")
tukey_lean

# Round predicted values and confidence intervals to 4 decimal places
effects_fat_plant <- effects_fat_plant %>%
  mutate(predicted = round(predicted, 4),
         conf.low = round(conf.low, 4),
         conf.high = round(conf.high, 4))

print(effects_lean_plant)



#####################################################################################################
############################ 1.4 WATER EFFECTS ########################################################
# Fit the linear model for WaterMass
overall_water <- lm(WaterMass ~ Plant + Sex + ForewingLength + RawWeight_day0 + TotalSA, data = GH_data)
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
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 13),  # Increase font size for x-axis labels
        axis.text.y = element_text(size = 13),                          # Increase font size for y-axis labels
        axis.title.x = element_text(size = 13),                        # Increase font size for x-axis title
        axis.title.y = element_text(size = 13),                        # Increase font size for y-axis title
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
          show.stat = TRUE,
          show.se = TRUE,
          show.df = TRUE
          ) 
Anova(overall_water)
summary(overall_water)

plot(allEffects(overall_water))



interact_weight
# Remove x-axis labels from the upper panels and legend from weight_plot
weight_panel <- interaction_plot + theme(
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

gh_multi_panel_plot = gh_multi_panel_plot + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 20))  # Increase the left margin (l = 20)


# Display the multi-panel plot
gh_multi_panel_plot 

#save the multi-panel plot
ggsave("plots/fig3_comboplot_weight+BodyComp.png", gh_multi_panel_plot, width = 27, height = 25, units = "cm", dpi = 300,)



#####################################################################################################
############################## 2. GH VS FIELD ########################################################
############################ 2.1 WEIGHT ########################################################

ghf_weight = lmer(Weight~ExpLoc*TrialDay+ForewingLength+Sex+(1|ID),data=gh_field_weight)

check_model(ghf_weight)
Anova(ghf_weight, type = 3)
# Extract the model summary
summary_model <- summary(ghf_weight)
summary_model

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
ghf_weight_summary = summary(ghf_weight)

# Extract the p-value for ExpLoc (ExpLocGreenhouse as per the correct coefficient name)
p_weight_exploc <- ghf_weight_summary$coefficients["ExpLocGreenhouse", "Pr(>|t|)"]

# Format the p-value: if it's less than 0.001, display "<0.001"; if greater than 0.05, display nothing
if (p_weight_exploc < 0.001) {
  p_weight_label <- "***"
} else if (p_weight_exploc < 0.01) {
  p_weight_label <- "**"
} else if (p_weight_exploc < 0.05) {
  p_weight_label = "*"
} else {
  p_weight_label <- ""  # No label if p > 0.05
} 


exploc_weight_plot <- ggplot(effects_weight_exploc, show_residuals = TRUE, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_jitter(data = gh_field_weight, aes(x = ExpLoc, y = Weight), 
              width = 0.1, height = 0, alpha = 0.25) +  # Raw data points colored by TrialDay
  labs(y = "Wet weight (g)", x = "") +  
  my_plot_theme +
  
  # Add the p-value label, now adjusted for position
  annotate("text", x = 2.5, y = 0.72, 
           label = p_weight_label, hjust = 1.5, vjust = 1.5, size = 7, color = "black", fontface = "bold") +
  
  # Adjust x-axis labels
  scale_x_discrete(labels = c("GH" = "Greenhouse", "Field" = "Field")) +
  
  theme(axis.text.x = element_text(size = 13, hjust = 0.5),  # Center the x-axis labels
        axis.text.y = element_text(size = 13),                # Increase font size for y-axis labels
        axis.title.x = element_text(size = 13),               # Increase font size for x-axis title
        axis.title.y = element_text(size = 13),               # Increase font size for y-axis title
        panel.background = element_blank(),                   # Remove background color
        plot.background = element_blank(),                    # Remove plot background color
        panel.border = element_blank(),                       # Remove panel border
        axis.line = element_line(color = "black"),            # Keep axis lines
        legend.position = c(0.97, 0.88),                      # Adjust legend position to avoid overlap
        legend.background = element_blank(),                  # Remove background color for legend
        legend.box.background = element_blank())              # Remove border around the legend box


# Display the plot
exploc_weight_plot

summary(ghf_weight)
check_model(ghf_weight)
tab_model(ghf_weight, 
          file = "tables/ghf_Weight.doc",   # Export to a Word document
          show.p = TRUE, # Add significance stars
          show.stat = TRUE,
          show.se = TRUE,
          show.df = TRUE
) 


#####################################################################################################
############################ 2.2 GH-Field Fat comparison ############################

# Fit the model
ghf_fat <- lm(DryFatMass ~ ExpLoc + Sex+RawWeight_day0, data = GH_F_data)

check_model(ghf_fat)
vif(ghf_fat)

# Get the summary of the model
ghf_fat_summary <- summary(ghf_fat)

anova_field_fat <- Anova(field_fat)
p_field_fat <- anova_field_fat["Plant", "Pr(>F)"]  # Adjust depending on output format
significance_fat <- ifelse(p_field_fat < 0.001, "***",
                           ifelse(p_field_fat < 0.01, "**",
                                  ifelse(p_plant_fat < 0.05, "*", "")))

# Extract the predicted effects as a data frame
effects_ghf_fat_df <- as.data.frame(ggpredict(ghf_fat, terms = "ExpLoc"))

check_model(ghf_fat)
Anova(ghf_fat)

# Extract the p-value for ExpLoc (ExpLocGreenhouse as per the correct coefficient name)
p_fat_exploc <- ghf_fat_summary$coefficients["ExpLocGreenhouse", "Pr(>|t|)"]

# Format the p-value: if it's less than 0.001, display "<0.001"; if greater than 0.05, display nothing
if (p_fat_exploc < 0.001) {
  p_fat_label <- "***"
} else if (p_fat_exploc < 0.01) {
  p_fat_label <- "**"
} else if (p_fat_exploc < 0.05) {
  p_fat_label = "*"
} else {
  p_fat_label <- ""  # No label if p > 0.05
} 

# Create the plot for DryFatMass by ExpLoc
fat_loc_plot <- ggplot(effects_ghf_fat_df, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(y = "Fat (g)", x = "", title = "") +
  
  # Add raw data points with jitter
  geom_jitter(data = GH_F_data, aes(x = ExpLoc, y = DryFatMass), 
              width = 0.1, alpha = 0.25, size = 1.5) +
  
  # Add the p-value label, now adjusted for position
  annotate("text", x = 2.25, y = 0.06, 
           label = p_fat_label, hjust = 1.5, vjust = 1.5, size = 7, color = "black", fontface = "bold") +
  
  # Adjust theme and formatting
  theme(axis.text.x = element_text(hjust = 0.5, size = 13),  # X-axis label size
        axis.text.y = element_text(size = 13),               # Y-axis label size
        axis.title.x = element_text(size = 13),              # X-axis title size
        axis.title.y = element_text(size = 13),              # Y-axis title size
        panel.background = element_blank(),                  # Remove background
        plot.background = element_blank(),                   # Remove plot background
        panel.border = element_blank(),                      # Remove panel border
        axis.line = element_line(color = "black"))           # Keep axis lines

# Display the plot
fat_loc_plot
summary(ghf_fat)
check_model(ghf_fat)

tab_model(ghf_fat, 
          file = "tables/ghf_Fat.doc",   # Export to a Word document
          show.p = TRUE, # Add significance stars
          show.stat = TRUE,
          show.se = TRUE,
          show.df = TRUE
) 

############################ 2.3 GH-Field Lean mass comparison #############################
# Fit the model
ghf_lean <- lm(DryLeanMass ~ ExpLoc + Sex+RawWeight_day0, data = GH_F_data)

# Calculate VIF for the model
vif_values <- vif(ghf_lean)

vif_values

# Check if the VIF output is a matrix/data frame or a vector
if (is.matrix(vif_values) || is.data.frame(vif_values)) {
  # For multiple degrees of freedom, adjust each GVIF^(1/(2*Df))
  gvif_adjusted <- vif_values[, "GVIF"]^(1 / (2 * vif_values[, "Df"]))
} else {
  # If VIF is a simple vector, no adjustment needed for Df (assumes Df = 1 for each)
  gvif_adjusted <- vif_values^(1 / 2)
}

# Display the adjusted GVIF values
gvif_adjusted

# Extract the predicted effects as a data frame
effects_ghf_lean_df <- as.data.frame(ggpredict(ghf_lean, terms = "ExpLoc"))

check_model(ghf_lean)
Anova(ghf_lean)

ghf_lean_summary = summary(ghf_lean)

# Extract the p-lean for ExpLoc (ExpLocGreenhouse as per the correct coefficient name)
p_lean_exploc <- ghf_lean_summary$coefficients["ExpLocGreenhouse", "Pr(>|t|)"]

if (p_lean_exploc < 0.001) {
  p_lean_label <- "***"
} else if (p_lean_exploc < 0.01) {
  p_lean_label <- "**"
} else if (p_lean_exploc < 0.05) {
  p_lean_label = "*"
} else {
  p_lean_label <- ""  # No label if p > 0.05
} 

# Create the plot for DryLeanMass by ExpLoc
lean_loc_plot <- ggplot(effects_ghf_lean_df, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(y = "Lean mass (g)", x = "", title = "") +
  
  # Add raw data points with jitter
  geom_jitter(data = GH_F_data, aes(x = ExpLoc, y = DryLeanMass), 
              width = 0.1, alpha = 0.25, size = 1.5) +
  
  # Add the p-lean label, now adjusted for position
  annotate("text", x = 2.25, y = 0.15, 
           label = p_lean_label, hjust = 1.5, vjust = 1.5, size = 7, color = "black", fontface = "bold") +
  
  # Adjust theme and formatting
  theme(axis.text.x = element_text(hjust = 0.5, size = 13),  # X-axis label size
        axis.text.y = element_text(size = 13),               # Y-axis label size
        axis.title.x = element_text(size = 13),              # X-axis title size
        axis.title.y = element_text(size = 13),              # Y-axis title size
        panel.background = element_blank(),                  # Remove background
        plot.background = element_blank(),                   # Remove plot background
        panel.border = element_blank(),                      # Remove panel border
        axis.line = element_line(color = "black"))           # Keep axis lines

# Display the plot
lean_loc_plot

summary(ghf_lean)
check_model(ghf_lean)

tab_model(ghf_lean, 
          file = "tables/ghf_Lean.doc",   # Export to a Word document
          show.p = TRUE, # Add significance stars
          show.stat = TRUE,
          show.se = TRUE,
          show.df = TRUE
) 



############### 2.4 Greenhouse-Field Water Mass ################################

# Fit the model
ghf_water <- lm(WaterMass ~ ExpLoc+Sex+RawWeight_day0, data = GH_F_data)

# Get the summary of the model
ghf_water_summary <- summary(ghf_water)

# Extract the predicted effects as a data frame
effects_ghf_water_df <- as.data.frame(ggpredict(ghf_water, terms = "ExpLoc"))

check_model(ghf_water)
Anova(ghf_water)

# Create the plot for WaterMass by ExpLoc
water_loc_plot <- ggplot(effects_ghf_water_df, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(y = "Water (g)", x = "", title = "") +
  
  # Add raw data points with jitter
  geom_jitter(data = GH_F_data, aes(x = ExpLoc, y = WaterMass), 
              width = 0.1, alpha = 0.25, size = 1.5) +
  
  # Adjust theme and formatting
  theme(axis.text.x = element_text(hjust = 0.5, size = 13),  # X-axis label size
        axis.text.y = element_text(size = 13),               # Y-axis label size
        axis.title.x = element_text(size = 13),              # X-axis title size
        axis.title.y = element_text(size = 13),              # Y-axis title size
        panel.background = element_blank(),                  # Remove background
        plot.background = element_blank(),                   # Remove plot background
        panel.border = element_blank(),                      # Remove panel border
        axis.line = element_line(color = "black"))           # Keep axis lines

# Display the plot
water_loc_plot

summary(ghf_water)
tab_model(ghf_water, 
          file = "tables/ghf_Water.doc",   # Export to a Word document
          show.p = TRUE, # Add significance stars
          show.stat = TRUE,
          show.se = TRUE,
          show.df = TRUE
) 


gh_f_multi_panel_plot <- (exploc_weight_plot + fat_loc_plot) / (lean_loc_plot + water_loc_plot) + 
  plot_annotation(tag_levels = 'a')&
  theme(plot.tag.position = c(0.05, 1))

# Display the combined plot
gh_f_multi_panel_plot

ggsave("plots/comboplot_GH_Field.png", gh_f_multi_panel_plot, width = 27, height = 20, units = "cm", dpi = 300,)

###############################################################################################
############################ 3. FIELD SolAlt & EutMac COMPARISON##############################

field_fem_weight = subset(field_TreatmentData,field_TreatmentData$Sex == "F") #remove males as there were not many
field_fem_weight = subset(field_fem_weight,field_fem_weight$EnclCol != "FBrown") #Remove brown enclosures

field_fem_comp = subset(Field_data,Field_data$Sex == "F")
field_fem_comp = subset(field_fem_comp,field_fem_comp$EnclCol != "FBrown")

# Add EnclID to the fat dataset by matching the ID from weight data
field_fem_comp$EnclID <- field_fem_weight$EnclID[match(field_fem_comp$ID, field_fem_weight$ID)]
summary(field_fem_comp)

sum(Field_data$Sex == "M")
sum(Field_data$Sex == "F")
sum(field_fem_comp$Plant == "E. maculatum")
sum(field_fem_comp$Plant == "S. altissima")
sum(Field_data$EnclCol =="FBrown")


sum(field_fem_comp$EnclCol =="FBrown")
sum(field_fem_comp$EnclCol =="FBlack")

############################ 3.1 FIELD WEIGHT ###############################

field_weight = lmer(Weight ~ Plant*TrialDay+TotalSA+(1|EnclID/ID), data = field_fem_weight)
check_model(field_weight)



# Get predicted values for Plant effect (fixing random effects)
predicted_field_weight <- ggpredict(field_weight, terms = "Plant")

# Calculate residuals from the model
field_fem_weight$residuals_field_weight <- residuals(field_weight)

# Account for random effects: generate conditional residuals
predicted_values <- predict(field_weight, re.form = ~(1 | EnclID/ID), allow.new.levels = TRUE)
field_fem_weight$conditional_residuals <- field_fem_weight$Weight - predicted_values

# Calculate the mean predicted value for each Plant level
mean_predicted_values <- predicted_field_weight$predicted

# Recalculate residuals relative to the predicted mean for each factor level of Plant
field_fem_weight$centered_residuals <- field_fem_weight$residuals_field_weight + 
  mean_predicted_values[match(field_fem_weight$Plant, predicted_field_weight$x)]

# Create the plot for the residuals with the same level of grey as your raw data points
field_weight_plot <- ggplot(predicted_field_weight, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(y = "Wet weight (g)", x = "", title = "") +
  
  # Add residual points with jitter, using the same grey level for styling as in your previous example
  geom_jitter(data = field_fem_weight, aes(x = Plant, y = centered_residuals), 
              width = 0.1, alpha = 0.25, size = 1.5, color = "grey30") +  # Use grey30 for residual points
  
  # Add significance annotation
  annotate("text", x = length(unique(predicted_field_weight$x)) + 0.3, 
           y = max(predicted_field_weight$predicted) + 0.1, 
           label = significance, size = 6) +
  
  # Adjust theme and formatting to preserve your initial design
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 13),  # X-axis label size
    axis.text.y = element_text(size = 13),               # Y-axis label size
    axis.title.x = element_text(size = 13),              # X-axis title size
    axis.title.y = element_text(size = 13),              # Y-axis title size
    panel.background = element_blank(),                  # Remove background
    plot.background = element_blank(),                   # Remove plot background
    panel.border = element_blank(),                      # Remove panel border
    axis.line = element_line(color = "black")            # Keep axis lines
  )

# Display the plot
field_weight_plot


anova_field_weight <- Anova(field_weight, type = "III")  # Specify type = "III"
anova_field_weight

# Extract fixed-effects predictions
field_fem_weight <- field_fem_weight %>%
  mutate(
    # Predicted values without random effects
    fixed_effects_pred = predict(field_weight, re.form = NA),
    # Partial residuals for TotalSA
    partial_residuals_totalsa = residuals(field_weight) + fixed_effects_pred
  )

# Select the variables we need for plotting (TotalSA and residuals)
partial_residual_data <- field_fem_weight %>%
  select(TotalSA, partial_residuals_totalsa)

# Get predictions for TotalSA from ggeffects
predicted_field_weight <- ggeffects::ggpredict(field_weight, terms = "TotalSA") %>%
  rename(x = x, predicted = predicted, conf.low = conf.low, conf.high = conf.high)

# Add partial residuals
partial_residual_data <- field_fem_weight %>%
  select(TotalSA, partial_residuals_totalsa)








p_field_plant <- anova_field_weight["Plant", "Pr(>Chisq)"]  # Adjust depending on output format
significance <- ifelse(p_field_plant < 0.001, "***",
                       ifelse(p_field_plant < 0.01, "**",
                              ifelse(p_field_plant < 0.05, "*", "")))

summary(field_weight)
plot(allEffects(field_weight))

tab_model(field_weight)
# Predicted effects for the model
predicted_field_weight <- ggpredict(field_weight, terms = "Plant")

# Create the plot for DryFatMass by ExpLoc
field_weight_plot_rawdata <- ggplot(predicted_field_weight, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(y = "Wet weight (g)", x = "", title = "") +
  
  # Add raw data points with jitter
  geom_jitter(data = field_fem_weight, aes(x = Plant, y = Weight), 
              width = 0.1, alpha = 0.25, size = 1.5) +
  
  # Add significance annotation
  annotate("text", x = length(unique(predicted_field_weight$x)) + 0.3, 
           y = max(predicted_field_weight$predicted) + 0.1, 
           label = significance, size = 6) +
  
  # Adjust theme and formatting
  theme(
        axis.text.x = element_text(hjust = 0.5, size = 13),  # X-axis label size
        axis.text.y = element_text(size = 13),               # Y-axis label size
        axis.title.x = element_text(size = 13),              # X-axis title size
        axis.title.y = element_text(size = 13),              # Y-axis title size
        panel.background = element_blank(),                  # Remove background
        plot.background = element_blank(),                   # Remove plot background
        panel.border = element_blank(),                      # Remove panel border
        axis.line = element_line(color = "black"))           # Keep axis lines


# Print the plot
print(field_weight_plot_rawdata)

############################# 3.2 FIELD FAT ##########################
# Fit the model
field_fat <- lmer(DryFatMass ~ Plant+TotalSA+RawWeight_day0+(1|EnclID), data = field_fem_comp)
check_model(field_fat)

Anova(field_fat)

# Predict effects for the model
predicted_field_fat <- ggpredict(field_fat, terms = "Plant")

tab_model(field_fat)
plot(allEffects(field_fat))

# Create the plot for DryFatMass by ExpLoc
fat_field_plot_raw <- ggplot(predicted_field_fat, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(y = "Fat (g)", x = "", title = "") +
  
  # Add raw data points with jitter
  geom_jitter(data = field_fem_comp, aes(x = Plant, y = DryFatMass), 
              width = 0.1, alpha = 0.25, size = 1.5) +
  
  
  # Adjust theme and formatting
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 13),  # X-axis label size
    axis.text.y = element_text(size = 13),               # Y-axis label size
    axis.title.x = element_text(size = 13),              # X-axis title size
    axis.title.y = element_text(size = 13),              # Y-axis title size
    panel.background = element_blank(),                  # Remove background
    plot.background = element_blank(),                   # Remove plot background
    panel.border = element_blank(),                      # Remove panel border
    axis.line = element_line(color = "black"))           # Keep axis lines

fat_field_plot_raw




# Step 1: Get residuals for the model
field_fem_comp$residuals_fat <- residuals(field_fat)

# Step 2: Get predicted values for the model
predicted_field_fat <- ggpredict(field_fat, terms = "Plant")

# Step 3: Extract the random effects for EnclID
random_effects <- ranef(field_fat)$EnclID

# Step 4: Make sure random_effects is in the correct format for matching
random_effects <- as.data.frame(random_effects)

# Step 5: Add the corresponding random effect for each EnclID into the data
field_fem_comp$random_effect <- random_effects[match(field_fem_comp$EnclID, rownames(random_effects)), 1]

# Step 6: Center the residuals by adding the mean predicted value and the random effect
field_fem_comp$centered_residuals_fat <- field_fem_comp$residuals_fat + mean(predicted_field_fat$predicted) + field_fem_comp$random_effect

# Step 7: Plot the residuals for the field fat model
field_fat_residual_plot <- ggplot(predicted_field_fat, aes(x = x, y = predicted)) +
  geom_point(size = 3) +  # Predicted points
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(y = "Dry fat mass (g)", x = "", title = "") +
  
  # Add residual data points (from field_fem_comp) with jitter
  geom_jitter(data = field_fem_comp, aes(x = Plant, y = centered_residuals_fat), 
              width = 0.1, alpha = 0.25, size = 1.5, color = "grey30") +  # Residuals in grey
  
  # Formatting theme as per original plot
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 13),  # X-axis label size
    axis.text.y = element_text(size = 13),               # Y-axis label size
    axis.title.x = element_text(size = 13),              # X-axis title size
    axis.title.y = element_text(size = 13),              # Y-axis title size
    panel.background = element_blank(),                  # Remove background
    plot.background = element_blank(),                   # Remove plot background
    panel.border = element_blank(),                      # Remove panel border
    axis.line = element_line(color = "black")            # Keep axis lines
  )

# Step 8: Display the residual plot for the 'fat' model
field_fat_residual_plot


############################# 3.3 Field lean #####################################

# Fit the model
field_lean <- lmer(DryLeanMass ~ Plant+TotalSA+ForewingLength + (1|EnclID), data = field_fem_comp)
check_model(field_lean)

Anova(field_lean)

# Predict effects for the model
predicted_field_lean <- ggpredict(field_lean, terms = "Plant")

# Create the plot for DryleanMass by ExpLoc
lean_field_plot_raw <- ggplot(predicted_field_lean, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(y = "Lean (g)", x = "", title = "") +
  
  # Add raw data points with jitter
  geom_jitter(data = field_fem_comp, aes(x = Plant, y = DryLeanMass), 
              width = 0.1, alpha = 0.25, size = 1.5) +
  
  
  # Adjust theme and formatting
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 13),  # X-axis label size
    axis.text.y = element_text(size = 13),               # Y-axis label size
    axis.title.x = element_text(size = 13),              # X-axis title size
    axis.title.y = element_text(size = 13),              # Y-axis title size
    panel.background = element_blank(),                  # Remove background
    plot.background = element_blank(),                   # Remove plot background
    panel.border = element_blank(),                      # Remove panel border
    axis.line = element_line(color = "black"))           # Keep axis lines

lean_field_plot_raw


# Step 1: Calculate the residuals from the model
field_fem_comp$residuals_lean <- residuals(field_lean)

# Step 2: Create predicted values for each data point using ggpredict
predicted_field_lean <- ggpredict(field_lean, terms = "Plant")

# Step 3: Calculate the mean of the predicted values
mean_predicted_value <- mean(predicted_field_lean$predicted)

# Step 4: Center the residuals around the mean predicted value
field_fem_comp$centered_residuals_lean <- field_fem_comp$residuals_lean + (mean_predicted_value - mean(field_fem_comp$residuals_lean))

# Step 5: Calculate significance stars based on p-value for Plant
p_field_lean_plant <- anova(field_lean)["Plant", "Pr(>Chisq)"]  # Extract p-value for Plant effect
significance_lean <- ifelse(p_field_lean_plant < 0.001, "***",
                            ifelse(p_field_lean_plant < 0.01, "**",
                                   ifelse(p_field_lean_plant < 0.05, "*", "")))

# Step 6: Create the residuals plot
field_lean_residual_plot <- ggplot(predicted_field_lean, aes(x = x, y = predicted)) +
  geom_point(size = 3) +  # Predicted points
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +  # Error bars
  
  # Add residuals (centered) to the plot, jittered for visibility
  geom_jitter(data = field_fem_comp, aes(x = Plant, y = centered_residuals_lean), 
              width = 0.1, alpha = 0.25, size = 1.5, color = "grey30") +  # Residuals
  
  # Step 7: Add significance annotation (stars)
  annotate("text", x = length(unique(predicted_field_lean$x)) + 0.3, 
           y = max(predicted_field_lean$predicted) + 0.1, 
           label = significance_lean, size = 6) +  # Place significance stars
  
  # Step 8: Customize the plot's appearance
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 13),  # X-axis label size
    axis.text.y = element_text(size = 13),               # Y-axis label size
    axis.title.x = element_text(size = 13),              # X-axis title size
    axis.title.y = element_text(size = 13),              # Y-axis title size
    panel.background = element_blank(),                  # Remove background
    plot.background = element_blank(),                   # Remove plot background
    panel.border = element_blank(),                      # Remove panel border
    axis.line = element_line(color = "black")            # Keep axis lines
  ) +
  # Step 9: Correct the axis labels as needed
  labs(
    y = "Lean mass (g)",  # Correcting y-axis label
    x = NULL  # Remove x-axis label for Plant
  )

# Step 10: Print the final plot
print(field_lean_residual_plot)


########################## 3.4 Water mass ###################################

# Fit the model
field_water <- lmer(WaterMass ~ Plant+TotalSA+RawWeight_day0+(1|EnclID), data = field_fem_comp)
check_model(field_water)

Anova(field_water)

# Predict effects for the model
predicted_field_water <- ggpredict(field_water, terms = "Plant")

# Create the plot for DrywaterMass by ExpLoc
water_field_plot <- ggplot(predicted_field_water, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(y = "Water (g)", x = "", title = "") +
  
  # Add raw data points with jitter
  geom_jitter(data = field_fem_comp, aes(x = Plant, y = WaterMass), 
              width = 0.1, alpha = 0.25, size = 1.5) +
  
  
  # Adjust theme and formatting
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 13),  # X-axis label size
    axis.text.y = element_text(size = 13),               # Y-axis label size
    axis.title.x = element_text(size = 13),              # X-axis title size
    axis.title.y = element_text(size = 13),              # Y-axis title size
    panel.background = element_blank(),                  # Remove background
    plot.background = element_blank(),                   # Remove plot background
    panel.border = element_blank(),                      # Remove panel border
    axis.line = element_line(color = "black"))           # Keep axis lines

water_field_plot

field_multi_plot <- (field_weight_plot + fat_field_plot) / (lean_field_plot + water_field_plot) + 
  plot_annotation(tag_levels = 'a')&
  theme(plot.tag.position = c(0.22, .95))

# Display the combined plot
field_multi_plot

ggsave("plots/comboplot_Field.png", field_multi_plot, width = 27, height = 20, units = "cm", dpi = 300,)


# Step 2: Calculate the residuals from the model
field_fem_comp$residuals_water <- residuals(field_water)

# Step 3: Create predicted values for each data point using ggpredict
predicted_field_water <- ggpredict(field_water, terms = "Plant")

# Step 4: Calculate the mean of the predicted values
mean_predicted_value_water <- mean(predicted_field_water$predicted)

# Step 5: Center the residuals around the mean predicted value
field_fem_comp$centered_residuals_water <- field_fem_comp$residuals_water + (mean_predicted_value_water - mean(field_fem_comp$residuals_water))

# Step 6: Calculate significance stars based on p-value for Plant
p_field_water_plant <- anova(field_water)["Plant", "Pr(>Chisq)"]  # Extract p-value for Plant effect
significance_water <- ifelse(p_field_water_plant < 0.001, "***",
                             ifelse(p_field_water_plant < 0.01, "**",
                                    ifelse(p_field_water_plant < 0.05, "*", "")))

# Step 7: Create the residuals plot
field_water_residual_plot <- ggplot(predicted_field_water, aes(x = x, y = predicted)) +
  geom_point(size = 3) +  # Predicted points
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +  # Error bars
  
  # Add residuals (centered) to the plot, jittered for visibility
  geom_jitter(data = field_fem_comp, aes(x = Plant, y = centered_residuals_water), 
              width = 0.1, alpha = 0.25, size = 1.5, color = "grey30") +  # Residuals
  
  # Step 8: Add significance annotation (stars)
  annotate("text", x = length(unique(predicted_field_water$x)) + 0.3, 
           y = max(predicted_field_water$predicted) + 0.1, 
           label = significance_water, size = 6) +  # Place significance stars
  
  # Step 9: Customize the plot's appearance (same theme)
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 13),  # X-axis label size
    axis.text.y = element_text(size = 13),               # Y-axis label size
    axis.title.x = element_blank(),                      # No X-axis title
    axis.title.y = element_text(size = 13),              # Y-axis title size
    panel.background = element_blank(),                  # Remove background
    plot.background = element_blank(),                   # Remove plot background
    panel.border = element_blank(),                      # Remove panel border
    axis.line = element_line(color = "black")            # Keep axis lines
  ) +
  labs(y = "Water Mass (g)")  # Set the correct y-axis label

# Step 10: Print the final plot
print(field_water_residual_plot)


########################## Make a combination plot for total surface area and each of the responses
# 1. field weight vs total surface area
# 1. For `field_weight`
predicted_field_weight <- ggpredict(field_weight, terms = "TotalSA")
summary_field_weight <- summary(field_weight)

p_field_weight_totalsa <- summary_field_weight$coefficients["TotalSA", "Pr(>|t|)"]
significance_weight <- ifelse(p_field_weight_totalsa < 0.001, "***",
                              ifelse(p_field_weight_totalsa < 0.01, "**",
                                     ifelse(p_field_weight_totalsa < 0.05, "*", "")))

# Plot for `field_weight`
field_weight_sa_plot <- ggplot(predicted_field_weight, aes(x = x, y = predicted)) +
  geom_line(size = 1, aes(group = 1)) +  # Ensures a line is drawn
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) + 
  labs(y = "Weight (g)", x = "Total Surface Area (cm²)", title = "") +
  geom_jitter(data = field_fem_weight, aes(x = TotalSA, y = Weight), 
              width = 0.1, alpha = 0.25, size = 1.5) +
  geom_text(aes(x = max(predicted_field_weight$x), 
                y = max(predicted_field_weight$predicted), 
                label = significance_weight), 
            color = "black", size = 6, hjust = 10, vjust = -11) +
  theme(axis.text.x = element_text(hjust = 0.5, size = 13), 
        axis.text.y = element_text(size = 13), 
        axis.title.x = element_text(size = 13), 
        axis.title.y = element_text(size = 13), 
        panel.background = element_blank(),
        plot.background = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))

print(field_weight_sa_plot)


# 2. For `field_fat`
predicted_field_fat <- ggpredict(field_fat, terms = "TotalSA")
summary_field_fat <- summary(field_fat)
p_field_fat_totalsa <- summary_field_fat$coefficients["TotalSA", "Pr(>|t|)"]
significance_fat <- ifelse(p_field_fat_totalsa < 0.001, "***",
                           ifelse(p_field_fat_totalsa < 0.01, "**",
                                  ifelse(p_field_fat_totalsa < 0.05, "*", "")))

# Plot for `field_fat`
field_fat_sa_plot <- ggplot(predicted_field_fat, aes(x = x, y = predicted)) +
  geom_line(size = 1, aes(group = 1)) +  
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) + 
  labs(y = "Dry Fat Mass (g)", x = "Total Surface Area (cm²)", title = "") +
  geom_jitter(data = field_fem_comp, aes(x = TotalSA, y = DryFatMass), 
              width = 0.1, alpha = 0.25, size = 1.5) +
  geom_text(aes(x = max(predicted_field_fat$x), 
                y = max(predicted_field_fat$predicted), 
                label = significance_fat), 
            color = "black", size = 6, hjust = 10, vjust = -11) +
  theme(axis.text.x = element_text(hjust = 0.5, size = 13), 
        axis.text.y = element_text(size = 13), 
        axis.title.x = element_text(size = 13), 
        axis.title.y = element_text(size = 13), 
        panel.background = element_blank(),
        plot.background = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))


print(field_fat_sa_plot)

# 3. For `field_lean`
predicted_field_lean <- ggpredict(field_lean, terms = "TotalSA")
summary_field_lean <- summary(field_lean)
p_field_lean_totalsa <- summary_field_lean$coefficients["TotalSA", "Pr(>|t|)"]
significance_lean <- ifelse(p_field_lean_totalsa < 0.001, "***",
                            ifelse(p_field_lean_totalsa < 0.01, "**",
                                   ifelse(p_field_lean_totalsa < 0.05, "*", "")))

# Plot for `field_lean`
field_lean_sa_plot <- ggplot(predicted_field_lean, aes(x = x, y = predicted)) +
  geom_line(size = 1, aes(group = 1)) +  
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) + 
  labs(y = "Dry Lean Mass (g)", x = "Total Surface Area (cm²)", title = "") +
  geom_jitter(data = field_fem_comp, aes(x = TotalSA, y = DryLeanMass), 
              width = 0.1, alpha = 0.25, size = 1.5) +
  geom_text(aes(x = max(predicted_field_lean$x), 
                y = max(predicted_field_lean$predicted), 
                label = significance_lean), 
            color = "black", size = 6, hjust = 10, vjust = -11) +
  theme(axis.text.x = element_text(hjust = 0.5, size = 13), 
        axis.text.y = element_text(size = 13), 
        axis.title.x = element_text(size = 13), 
        axis.title.y = element_text(size = 13), 
        panel.background = element_blank(),
        plot.background = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))

print(field_lean_sa_plot)

# 4. Field water vs total SA


# Get predictions from the model for TotalSA effect
predicted_field_water <- ggpredict(field_water, terms = "TotalSA")

# Get the p-value for TotalSA from the model summary
summary_field_water <- summary(field_water)
p_field_totalsa <- summary_field_water$coefficients["TotalSA", "Pr(>|t|)"]

# Determine significance label based on the p-value
significance <- ifelse(p_field_totalsa < 0.001, "***",
                       ifelse(p_field_totalsa < 0.01, "**",
                              ifelse(p_field_totalsa < 0.05, "*", "")))

# Create the plot for TotalSA effect
field_water_sa_plot <- ggplot(predicted_field_water, aes(x = x, y = predicted)) +
  # Line for predicted values (grouping by x)
  geom_line(size = 1, aes(group = 1)) +  # Ensures a line is drawn through the predicted points
  
  # Confidence interval ribbon for the prediction
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) + 
  
  labs(y = "Water (g)", x = "Total Surface Area (cm²)", title = "") +
  
  # Add raw data points with jitter
  geom_jitter(data = field_fem_comp, aes(x = TotalSA, y = WaterMass), 
              width = 0.1, alpha = 0.25, size = 1.5) +
  
  # Add significance label in the top-right corner
  geom_text(aes(x = max(predicted_field_water$x), 
                y = max(predicted_field_water$predicted), 
                label = significance), 
            color = "black", size = 6, hjust = 10, vjust = -11) +
  
  # Adjust theme and formatting
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 13),  # X-axis label size
    axis.text.y = element_text(size = 13),               # Y-axis label size
    axis.title.x = element_text(size = 13),              # X-axis title size
    axis.title.y = element_text(size = 13),              # Y-axis title size
    panel.background = element_blank(),                  # Remove background
    plot.background = element_blank(),                   # Remove plot background
    panel.border = element_blank(),                      # Remove panel border
    axis.line = element_line(color = "black")            # Keep axis lines
  )

# Print the plot
print(field_water_sa_plot)


# Add significance labels and set the color to black for the water plot
field_weight_sa_panel <- field_weight_sa_plot + 
  geom_text(aes(x = max(predicted_field_weight$x), 
                y = max(predicted_field_weight$predicted), 
                label = significance_weight), 
            color = "black", size = 6, hjust = 1.2, vjust = 2)

field_fat_sa_panel <- field_fat_sa_plot + 
  geom_text(aes(x = max(predicted_field_fat$x), 
                y = max(predicted_field_fat$predicted), 
                label = significance_fat), 
            color = "black", size = 6, hjust = 5, vjust = -3.5)

field_lean_sa_panel <- field_lean_sa_plot + 
  geom_text(aes(x = max(predicted_field_lean$x), 
                y = max(predicted_field_lean$predicted), 
                label = significance_lean), 
            color = "black", size = 6, hjust = 5, vjust = -4)

# Add significance labels and set the color to black for the water plot
field_water_sa_panel <- field_water_sa_plot + 
  geom_text(aes(x = max(predicted_field_water$x), 
                y = max(predicted_field_water$predicted), 
                label = significance_water), 
            color = "black", size = 6, hjust = 1.2, vjust = 2)

field_fat_sa_panel <- field_fat_sa_panel + theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank()
)

field_weight_sa_panel <- field_weight_sa_panel + theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank()
)

# Combine the plots with adjustments for the tag position
field_multisa_plot <- (field_weight_sa_panel + field_fat_sa_panel) / 
  (field_lean_sa_panel + field_water_sa_plot) + 
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag.position = c(0.22, 0.95))

# Display the combined plot
field_multisa_plot

ggsave("plots/comboplot_Field_sa.png", field_multisa_plot, width = 27, height = 20, units = "cm", dpi = 300,)




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
p_nolegend + my_plot_theme

mean(GH_data$RelWeightGain_day7)
std.error(GH_data$RelWeightGain_day7)

mean(Field_data$RelWeightGain_day7)
std.error(Field_data$RelWeightGain_day7)

######################## Temperature data ###########################

# Ensure the Date and Time columns are combined into one datetime column
tempdat$DateTime <- as.POSIXct(paste(tempdat$Date, tempdat$Time24hr), format = "%Y-%m-%d %H:%M:%S")

# Add an indicator for day or night
tempdat$DayNight <- ifelse(format(tempdat$DateTime, "%H:%M:%S") >= "08:00:00" & format(tempdat$DateTime, "%H:%M:%S") < "20:00:00", "Day", "Night")

# Calculate daily statistics by location and include JulianDate
daily_summary <- tempdat %>%
  group_by(Date, JulianDate, Location) %>%
  summarize(
    MeanDayTemp = mean(Temperature[DayNight == "Day"], na.rm = TRUE),
    MeanNightTemp = mean(Temperature[DayNight == "Night"], na.rm = TRUE),
    MaxDayTemp = max(Temperature[DayNight == "Day"], na.rm = TRUE),
    MinNightTemp = min(Temperature[DayNight == "Night"], na.rm = TRUE),
    .groups = 'drop'  # Prevents the warning for unused grouping
  )

# Create the plot
mean_day_temp_plot <- ggplot(daily_summary, aes(x = JulianDate, y = MeanDayTemp, color = Location)) +
  geom_line() +  # Add lines for each location
  labs(
    x = "Julian Date",
    y = "Mean Daytime Temperature (°C)",  # Adjust units if necessary
    color = "Location"
  ) +
  theme_minimal() +  # Use a minimal theme for better aesthetics
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(colour = "black", fill = NA),  # Add border
    axis.line = element_line(colour = "black")  # Keep axis lines
  )

# Display the plot
print(mean_day_temp_plot)

# Define the Julian date for September 10th
julian_date_sept10 <- 253  # Replace with the actual Julian date for September 10

# Filter the dataset for plotting between August 31st and September 10th
final_filtered_data <- filtered_daily_summary[filtered_daily_summary$JulianDate >= 243 & 
                                                filtered_daily_summary$JulianDate <= julian_date_sept10, ]

# Create the plot
mean_day_temp_plot <- ggplot(final_filtered_data, aes(x = JulianDate, y = MeanDayTemp, color = Location)) +
  geom_line() +  # Add lines for each location
  labs(
    x = "Julian Date",
    y = "Mean Daytime Temperature (°C)"  # Adjust units if necessary
  ) +
  theme_minimal() +  # Use a minimal theme for better aesthetics
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove box outline
    axis.line = element_line(colour = "black"),  # Keep axis lines
    axis.ticks = element_line(size = 0.5)  # Add major tick marks
  ) +
  scale_y_continuous(limits = c(15, 30)) +  # Set y-axis limits from 15 to 30
  scale_x_continuous(
    limits = c(243, NA),  # Start at 243, end automatically based on data
    breaks = seq(243, max(final_filtered_data$JulianDate), by = 1)  # Only show integer breaks
  )  

# Create the plot for maximum daily temperature
max_temp_plot <- ggplot(final_filtered_data_max, aes(x = JulianDate, y = MaxDayTemp, color = Location)) +
  geom_line() +  # Add lines for each location
  labs(
    x = "Julian Date",
    y = "Maximum Daily Temperature (°C)"  # Adjust units if necessary
  ) +
  theme_minimal() +  # Use a minimal theme for better aesthetics
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove box outline
    axis.line = element_line(colour = "black"),  # Keep axis lines
    axis.ticks = element_line(size = 0.5)  # Add major tick marks
  ) +
  scale_y_continuous(limits = c(15, 30)) +  # Set y-axis limits from 15 to 30
  scale_x_continuous(
    limits = c(243, NA),  # Start at 243, end automatically based on data
    breaks = seq(243, max(final_filtered_data_max$JulianDate), by = 1)  # Only show integer breaks
  )  

# Display the plot
print(max_temp_plot)

# Display the plot
print(mean_day_temp_plot)

lm_daytemp = lm(MeanDayTemp ~ Location, data = daily_summary)
Anova(lm_daytemp)
plot(allEffects(lm_daytemp))
summary(lm_daytemp)

lm_maxtemp = lm(MaxDayTemp ~ Location, data = daily_summary)
Anova(lm_maxtemp)
plot(allEffects(lm_maxtemp))
summary(lm_maxtemp)



##########################RANDOM STUFF I NEED################################

# Count of butterflies that came from couples

counts_GH <- GH_data %>%
  filter(Couple %in% c("Jacob", "MVD", "PCD", "TWA","S47") | is.na(Couple)) %>%
  count(Couple, name = "Count")

counts_field <- Field_data %>%
  filter(Couple %in% c("Jacob", "MVD", "PCD", "TWA","S47") | is.na(Couple)) %>%
  count(Couple, name = "Count")

counts_allLocs <- summarydat %>%
  filter(Couple %in% c("Jacob", "MVD", "PCD", "TWA","S47") | is.na(Couple)) %>%
  count(Couple, name = "Count")

unique_counts <- gh_TreatmentData %>%
  filter(Couple %in% c("Jacob", "MVD", "PCD", "TWA","S47") | is.na(Couple)) %>%
  group_by(Couple) %>%
  summarise(Unique_Count = n_distinct(ID))
unique_counts

# View the result
print(counts_GH)
print(counts_field)
print(counts_allLocs)


