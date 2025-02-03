### Temperature analysis script
### Written by M. Veselovsky
### Last modified 2025-01-31

rm(list=ls())

# load data files
tempdat <-read.csv("processed/tempdat.csv")

# Load necessary libraries
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(performance)
library(ggeffects)
library(splines)  # For natural splines
library(dplyr)
library(tidyr)
library(nlme)



# Preprocess data: Select necessary columns and ensure they're numeric
processed_data <- tempdat %>%
  mutate(
    JulianDate = as.numeric(JulianDate),  # Ensure JulianDate is numeric
    Time24Hr = as.numeric(format(strptime(Time24Hr, format = "%H:%M:%S"), "%H"))  # Time to numeric
  )

# # Fit the lmer model with only fixed effects
# temp_model <- lm(Temperature ~ Location, data = tempdat)
# 
# # Use ggpredict to get model predictions
# predicted_temp <- ggpredict(temp_model, terms = "Location")
# plot(allEffects(temp_model))
# 
# 
# # Update "Fletcher" label to "Field" for the x-axis
# predicted_temp$x <- recode(predicted_temp$x, "Fletcher" = "Field")
# 
# 
# # Plot the predictions with the specified modifications
# ggplot(predicted_temp, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high)) +
#   geom_point(color = "black", size = 3) +  # Mean point as black and bigger
#   geom_errorbar(width = 0.2) +  # Error bars for the 95% CI
#   labs(x = NULL, y = "Temperature (°C)") +  # Y-axis label with degree symbol
#   scale_x_discrete(labels = function(x) gsub("_", " ", x)) +  # Remove location label below x-axis
#   scale_y_continuous(
#     limits = c(17.5, 25),  # Y-axis from 17.5°C to 25°C
#     breaks = seq(18, 24, by = 2)  # Tick marks every 2°C from 18°C to 24°C
#   ) +
#   theme_minimal(base_size = 14) +  # Base size for readability
#   theme(
#     panel.grid = element_blank(),  # Remove gridlines
#     axis.line = element_line(color = "black", size = 0.6),  # Show x and y axis lines
#     axis.ticks = element_line(color = "black"),  # Keep ticks on the axes
#     axis.text.x = element_text(size = 14,color="black"),  # Increase x-axis text size
#     axis.text.y = element_text(size = 12,color="black"),  # Adjust y-axis text size
#     axis.title.x = element_blank()  # Remove x-axis label
#   )
# 
# # Summarize the model to inspect results
# summary(temp_model)
# plot(allEffects(temp_model))

# Adjust Time from GMT to EST (subtract 5 hours)
processed_data <- processed_data %>%
  mutate(
    Time24Hr_EST = Time24Hr - 5  # Convert GMT to EST
  ) %>%
  mutate(
    Time24Hr_EST = ifelse(Time24Hr_EST < 0, Time24Hr_EST + 24, Time24Hr_EST)  # Wrap around if < 0
  )

# Filter for daytime temperatures (8 am EST to 8 pm EST)
daytime_data <- processed_data %>%
  filter(Time24Hr >= 8 & Time24Hr <= 20)  # Keep only daytime hours

# Further restrict data to Julian dates from August 31 to September 10 (243 to 253)
daytime_data <- daytime_data %>%
  filter(JulianDate >= 243 & JulianDate <= 253)

daytime_data <- daytime_data %>%
  mutate(DateTime_combined = paste(Date, Time24Hr, sep = "_"))

daytime_data$Location = as.factor(daytime_data$Location)
# Fit the lmer model with only fixed effects
daytemp_model <- lm(Temperature ~ Location, data = daytime_data)
daytemp_model_2 <- lmer(Temperature ~ Location + (1|DateTime_combined), data = daytime_data)

# Use ggpredict to get model predictions
predicted_temp <- ggpredict(daytemp_model_2, terms = "Location")
plot(allEffects(daytemp_model_2))


# Update "Fletcher" label to "Field" for the x-axis
predicted_temp$x <- recode(predicted_temp$x, "Fletcher" = "Field")


# Plot the predictions with the specified modifications
pred_temp_plot = ggplot(predicted_temp, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high)) +
  geom_point(color = "black", size = 3, position = position_dodge(width = 0.1)) +  # Reduce space between locations
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.1)) +  # Adjust error bars
  labs(x = NULL, y = "Temperature (°C)") +  # Y-axis label with degree symbol
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) +  # Remove location label below x-axis
  scale_y_continuous(
    limits = c(19.5, 26.5),  # Y-axis from 19.5°C to 27°C
    breaks = seq(18, 26, by = 2)  # Tick marks every 2°C from 18°C to 24°C
  ) +
  theme_minimal(base_size = 14) +  # Base size for readability
  theme(
    panel.grid = element_blank(),  # Remove gridlines
    axis.line = element_line(color = "black", size = 0.6),  # Show x and y axis lines
    axis.ticks = element_line(color = "black"),  # Keep ticks on the axes
    axis.text.x = element_text(size = 16, color = "black"),  # Increase x-axis text size
    axis.text.y = element_text(size = 16, color = "black"),  # Adjust y-axis text size
    axis.title.x = element_blank()  # Remove x-axis label
  )

# Display plot
pred_temp_plot

ggsave("plots/pred_temp_plot.jpeg", pred_temp_plot, width = 13, height = 17, units = "cm", dpi = 300,)


# Summarize the model to inspect results
summary(daytemp_model)

check_model(daytemp_model_2)

# Combine the Date and Time24Hr columns to create a datetime column
daytime_data <- daytime_data %>%
  mutate(DateTime_combined = as.POSIXct(paste(Date, Time24Hr), format = "%Y-%m-%d %H:%M"))

# Now reshaping the data to get paired temperatures for both locations at each time point
paired_data <- daytime_data %>%
  filter(Location %in% c("Fletcher", "Greenhouse")) %>%
  drop_na(Temperature) %>%
  pivot_wider(names_from = Location, values_from = Temperature)

# Now calculating the temperature difference between 'Fletcher' and 'Greenhouse'
paired_data <- paired_data %>%
  mutate(Temp_diff = Greenhouse - Fletcher)

t_test_result <- t.test(paired_data$Temp_diff)

# Mixed model to account for random effects due to JulianDate (i.e., individual dates)
paired_model <- lmer(Temp_diff ~ Time24Hr + (1|JulianDate), data = paired_data)

# Plotting the effects
predicted_diff <- ggpredict(paired_model, terms = "Time24Hr")
ggplot(predicted_diff, aes(x = x, y = predicted)) +
  geom_line() +
  geom_point() +
  labs(y = "Temperature Difference (°C)")


# Summarize the average temperature by JulianDate and Location with 95% CI
summary_data <- daytime_data %>%
  group_by(JulianDate, Location) %>%
  summarize(
    mean_temp = mean(Temperature, na.rm = TRUE),
    se_temp = sd(Temperature, na.rm = TRUE) / sqrt(n()),  # Standard Error
    ci_low = mean_temp - qt(0.975, df = n() - 1) * se_temp,  # Lower bound of 95% CI
    ci_high = mean_temp + qt(0.975, df = n() - 1) * se_temp  # Upper bound of 95% CI
  )

# Plotting the data
tempplot = ggplot(summary_data, aes(x = JulianDate, y = mean_temp, color = Location)) +
  geom_line(linewidth = 1.5) +  # Thicker line for visual clarity
  labs(
    x = "Julian Date",
    y = "Temperature (°C)",
    color = "Location"
  ) +
  # Set JulianDate range and tick marks on x axis
  scale_x_continuous(
    limits = c(243, 253),
    breaks = seq(243, 253, by = 2),
    expand = c(0, 0)  # No extra padding
  ) +
  # Set temperature range and tick marks on the y axis
  scale_y_continuous(
    limits = c(13.5, 29),
    breaks = seq(14, 29, by = 5),
    expand = c(0, 0)  # No extra padding
  ) +
  # Customize colors and rename the legend
  scale_color_manual(
    values = c("Greenhouse" = "blue", "Fletcher" = "red"),
    labels = c("Greenhouse", "Field")  # Rename "Fletcher" to "Field"
  ) +
  # Customize the appearance of the plot
  theme_minimal() +
  theme(
    # Set font sizes
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    # Remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Keep just axis lines
    axis.line = element_line(color = "black", size = 1),
    panel.border = element_blank(),
    # Explicitly show external tick marks with color
    axis.ticks.length = unit(0.3, "cm"),
    axis.ticks = element_line(color = "black"),
    # Adjusting plot margins to prevent clipping
    plot.margin = margin(20, 30, 20, 20) 
  ) +
  # Position legend inside graph area at requested spot
  theme(legend.position = c(0.75, 0.9))

# Print the plot
tempplot

# Save the plot with sufficient space to ensure no clipping
ggsave(
  "plots/temp_plot.jpeg",
  tempplot,
  width = 30,  
  height = 22, 
  units = "cm",
  dpi = 300
)



# Compute the range of average daytime temperatures for both locations
temperature_range <- summary_data %>%
  group_by(Location) %>%
  summarize(
    avg_temp = mean(mean_temp, na.rm = TRUE),
    min_temp = min(mean_temp, na.rm = TRUE),
    max_temp = max(mean_temp, na.rm = TRUE),
    range_temp = max_temp - min_temp
  )
print(temperature_range)


# Ensure Time24Hr_EST is numeric (if not already)
daytime_data$Time24Hr_EST <- as.numeric(daytime_data$Time24Hr_EST)

# Create the plot
ggplot(processed_data, aes(x = Time24Hr, y = Temperature, color = Location, group = Location)) +
  geom_point(alpha = 0.4) +  # Scatterplot of temperature readings
  geom_smooth(method = "loess", se = TRUE, span = 0.4) +  # Smoothed trend lines
  scale_x_continuous(breaks = seq(0, 24, by = 2)) +  # Labels every 2 hours
  labs(x = "Time of Day (24hr)", y = "Temperature (°C)", color = "Location") +
  theme_minimal() +
  theme(legend.position = "top")



# Ensure processed_data is ready for modeling
processed_data$Time24Hr <- as.numeric(processed_data$Time24Hr)   # Ensure Time24Hr is numeric
processed_data$Location <- as.factor(processed_data$Location)     # Ensure Location is a factor
processed_data$JulianDate <- as.numeric(processed_data$JulianDate) # JulianDate is numeric


# Fit the spline model with 4 degrees of freedom for Time24Hr, Location as a factor, and JulianDate as a random effect
temp_model <- lmer(Temperature ~ ns(Time24Hr, df = 3) * Location + (1 | JulianDate), data = processed_data)

Anova(temp_model)
summary(temp_model)
diag=check_model(temp_model)
summary(diag)
diag$PP_CHECK

# Use ggpredict to get model predictions
predictions <- ggpredict(temp_model, terms = c("Time24Hr [all]", "Location"))

# Plot the predicted temperature with the CI and location
ggplot(predictions, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(size = 1) +  # Plot the predicted temperature
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +  # Add CI shading
  labs(title = "Predicted Effect of Location on Temperature Over Time of Day",
       x = "Time of Day (EST)",
       y = "Predicted Temperature (°C)",
       color = "Location") +
  theme_minimal() +
  theme(legend.position = "top")

ggplot(predicted_data, aes(x = x, y = predicted, color = factor(group), fill = factor(group))) + 
  geom_line(size = 1) +  # Plot the predicted trend
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +  # Add CI shading
  labs(x = "Time of Day (EST)", y = "Predicted Temperature (°C)", color = "Location") + 
  theme_minimal() +
  theme(legend.position = c(0.1,0.9),  # Position inside the plot area (adjust x and y to move legend)
        legend.title = element_blank(),  # remve legend title
        legend.box.spacing = unit(0, "cm"),  # Tighten the legend spacing
        panel.grid = element_blank(),  # Remove gridlines
        axis.line = element_line(),  # Include x and y axis lines only
        axis.ticks = element_line(),  # Keep ticks at major points
        plot.title = element_blank(),  # Remove plot title
        legend.key = element_rect(fill = "transparent")) +  # Make legend background transparent
  guides(fill = "none")  # Remove the legend for the CI shading (fill)




