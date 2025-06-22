# Load required libraries
library(openxlsx)  # For Excel file reading/writing (if needed)
library(dplyr)     # Data manipulation
library(tidyr)     # Data tidying
library(readxl)    # Reading Excel files (if needed)
library(ggplot2)   # Plotting

# Create a dataframe containing Tajima's D values for different populations
tajima_data <- data.frame(
  population = c("Global", "Eur1/Mix", "Eur2/Dairy", "So-HC", "IndNA", "Mosaic", "W29 clade"),
  TajimaD = c(0.375038, 1.18481, -0.0189274, 1.02111, 0.749509, 0.271494, 0.673087)
)

# Define colors for each population
population_colors <- c(
  "Eur1/Mix" = "#F29E6D",
  "Eur2/Dairy" = "#05A6A6",
  "So-HC" = "#03738C",
  "IndNA" = "#F23E2E",
  "W29 clade" = "#152B59",
  "Mosaic" = "grey",
  "Global" = "grey4"
)

# Order populations: place "Global" first, then sort others descending by Tajima's D
rest_of_populations <- tajima_data %>% 
  filter(population != "Global") %>% 
  arrange(desc(TajimaD))

# Convert 'population' to a factor with levels ordered as desired
tajima_data$population <- factor(tajima_data$population, levels = c("Global", rest_of_populations$population))

# Plot ordered barplot of Tajima's D by population with specified colors
ggplot(tajima_data, aes(x = population, y = TajimaD, fill = population)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = population_colors) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    legend.position = "none",
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(color = "black")
  ) +
  labs(x = "Population", y = "Tajima's D")

