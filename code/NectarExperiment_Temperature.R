### Temperature analysis script
### Written by M. Veselovsky
### Last modified 2024-12-10

rm(list=ls())

# load data files
tempdat <-read.csv("processed/tempdat.csv")

# Load necessary libraries
library(lme4)
library(dplyr)
library(ggplot2)

# Preprocess data: Select necessary columns and ensure they're numeric
processed_data <- tempdat %>%
  mutate(
    JulianDate = as.numeric(JulianDate),  # Ensure JulianDate is numeric
    Time24Hr = as.numeric(format(strptime(Time24Hr, format = "%H:%M:%S"), "%H"))  # Time to numeric
  )

# Fit the lmer model with only fixed effects
temp_model <- lm(Temperature ~ Location, data = tempdat)

# Summarize the model to inspect results
summary(temp_model)
plot(allEffects(temp_model))


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
  filter(Time24Hr_EST >= 8 & Time24Hr_EST <= 20)  # Keep only daytime hours

# Further restrict data to Julian dates from August 31 to September 10 (243 to 253)
daytime_data <- daytime_data %>%
  filter(JulianDate >= 243 & JulianDate <= 253)

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

