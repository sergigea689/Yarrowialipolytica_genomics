#Follow up of pipeline of 01_Data_preparation
#This part is focused on analysis of the data generated.
#Matrices and R obbjects used are created in the previous pipeline

# -----------------------------------------------
# 1. Identifying phenotypes associated to genetic groups in quantitative data
# -----------------------------------------------

# Iterative statistical analysis

# Create an empty dataframe to store results
resultados_kruskal <- data.frame(
  variable = character(),  # Name of the variable (column)
  p_value = numeric(),     # p-value from the test
  stringsAsFactors = FALSE
)

# Combine datasets to be tested
kruskal_cond <- cbind(nut, inhib, lippop_unidades)

# Iterate over each column in the dataframe 'kruskal_cond'
for (columna in colnames(kruskal_cond)) {
  if (columna != "FinePop2") {  # Exclude the grouping variable column if present
    # Perform Kruskal-Wallis test comparing groups defined by 'popgensub$soil'
    kruskal_result <- kruskal.test(kruskal_cond[[columna]] ~ popgensub$soil)
    
    # Store the results in the dataframe
    resultados_kruskal <- rbind(resultados_kruskal, data.frame(
      variable = columna,
      p_value = kruskal_result$p.value
    ))
  }
}

# Filter only the variables with significant p-values (e.g., p < 0.05)
resultados_significativos <- resultados_kruskal[resultados_kruskal$p_value < 0.05, ]

# Display significant results
print(resultados_significativos)

# If you want only the names of significant variables:
variables_significativas <- resultados_significativos$variable
print(variables_significativas)

popconditionsliq <- variables_significativas 


# -----------------------------------------------
# 2. PERMANOVA analysis
# -----------------------------------------------

#Repeat one analysis for each data matrix

# Distance matrix calculation (Euclidean by default)
matriz <- dist(t(nutmat))
matriz <- as.matrix(matriz)

# Prepare dataframes for PERMANOVA with metadata
permainh <- as.data.frame(t(inhmat))
permainh$pop <- popgensub$FinePop2
permainh$morph <- popgensub$morph

permanut <- as.data.frame(t(nutmat))
permanut$pop <- popgensub$FinePop2
permanut$morph <- popgensub$morph

# Perform PERMANOVA for inhibitor dataset
resultados_permanova <- adonis(matriz ~ pop, data = permainh, permutations = 999)

# Perform PERMANOVA for nutrient dataset
resultados_permanova <- adonis(matriz ~ pop, data = permanut, permutations = 999)

# Print PERMANOVA results (ANOVA table)
print(resultados_permanova$aov.tab)

# -----------------------------------------------
# 3. Phylogenetic-phneotypic correlation analysis
# -----------------------------------------------

#Data
paperdata <- kruskal_cond 


# Select columns from 'sol' that start with "growth" and end with "240h"
selected_columns <- grep("^growth.*240h$", colnames(sol), value = TRUE)

# Subset those selected columns from 'sol'
paperdatasol <- sol[, selected_columns]

# Add columns 29 to 34 from 'sol' to 'paperdatasol'
paperdatasol <- cbind(paperdatasol, sol[,29:34])

# Set rownames of 'paperdatasol' equal to those of 'sol'
rownames(paperdatasol) <- rownames(sol)

# Load required libraries for phylogenetic analysis
library(phytools)
library(ape)

# Calculate phenotypic distance matrix using Euclidean distance
phenotypic_distance <- dist(paperdata, method = "euclidean")
print(phenotypic_distance)

# Convert the distance object to a matrix
phenodist <- as.matrix(phenotypic_distance)

# Load phylogenetic tree from file
arbol_filogenetico <- read.tree("biallelicvariants.recode.min4.phy.treefile")

# Calculate phylogenetic distance matrix (cophenetic distances from the tree)
distancia_filogenetica <- cophenetic(arbol_filogenetico)

# Print phylogenetic distances between strains
print(distancia_filogenetica)

# Convert phylogenetic distance object to matrix and print
distancia_filogenetica_matrix <- as.matrix(distancia_filogenetica)
print(distancia_filogenetica_matrix)

# Subset phylogenetic distance matrix to include only strains in 'popgensub$treesamp'
distphylosub <- subset(distancia_filogenetica_matrix, 
                       rownames(distancia_filogenetica_matrix) %in% popgensub$treesamp)

# Reorder the matrix rows and columns to match the order in 'popgensub$treesamp'
distphylosub <- distphylosub[popgensub$treesamp, popgensub$treesamp]

# Replace row and column names with strain names from 'popgensub$Name'
rownames(distphylosub) <- popgensub$Name
colnames(distphylosub) <- popgensub$Name

# Calculate correlation between phenotypic and phylogenetic distance matrices
phenovec <- as.vector(phenodist)
phylovec <- as.vector(distphylosub)
correlation <- cor.test(phenovec, phylovec, method = "pearson")
print(correlation)

# Plot phenotypic vs phylogenetic distances with regression line
plot(phenovec, phylovec, 
     xlab = "Phenotypic distance", ylab = "Phylogenetic distance",
     main = "")
abline(lm(phylovec ~ phenovec), col = "red")

# -----------------------------------------------
# 4. Phylogenetic signal analysis
# -----------------------------------------------

# Calculate phylogenetic signal (lambda) for traits in 'paperdata' using the phylogenetic tree 'arbol_filogenetico'

# Copy the original tree for manipulation
tree <- arbol_filogenetico

# Find common tip labels between the phylogenetic tree and the phenotype sample list
comunes <- intersect(tree$tip.label, popgensub$treesamp)

# Get indices of the common tips in the tree and in the phenotype dataframe
indices_árbol <- match(comunes, tree$tip.label)
indices_dataframe <- match(comunes, popgensub$treesamp)

# Rename the tip labels in the tree with corresponding names from popgensub$Name
tree$tip.label[indices_árbol] <- popgensub$Name[indices_dataframe]

# Match the order of rows in 'paperdata' to the order of tips in the tree
paperdatatree <- paperdata[match(tree$tip.label, rownames(paperdata)), ]

# Initialize vectors to store lambda values and permutation p-values
lambda <- c()
pval <- c()

# Filter the samples that are common between 'paperdatatree' and tree tip labels
common_samples <- intersect(rownames(paperdatatree), tree$tip.label)

# Subset 'paperdatatree' to only keep rows for common samples (avoid NA rows)
paperdatatree <- paperdatatree[common_samples, , drop = FALSE]

# Clean the tree by removing tips with NA labels
tree_clean <- drop.tip(tree, which(is.na(tree$tip.label)))

# Keep only tips in 'tree_clean' that are present in common_samples
tree_common <- drop.tip(tree_clean, setdiff(tree_clean$tip.label, common_samples))

# Loop through each trait (column) in 'paperdatatree' to calculate lambda and test significance
for(i in colnames(paperdatatree)) {
  print(i)  # Print the current trait name
  
  # Extract phenotype vector for the current trait
  fenotipo <- paperdatatree[[i]]
  
  # Remove NA values from phenotype vector to avoid errors
  fenotipo <- fenotipo[!is.na(fenotipo)]
  
  # Calculate observed lambda for the trait using 'phylosig' function
  lambda_value <- phylosig(tree_common, fenotipo, method = 'lambda')$lambda
  lambda <- c(lambda, lambda_value)
  
  # Vector to hold lambda values from permuted datasets
  blambda <- c()
  
  # Perform 200 permutations to generate null distribution of lambda
  for(j in 1:200) {
    # Randomly permute phenotype values without replacement
    t_d <- sample(fenotipo, replace = FALSE)
    
    # Assign names to permuted data to match common_samples order
    names(t_d) <- common_samples
    
    # Calculate lambda for the permuted phenotype vector
    perm_lambda <- phylosig(tree_common, t_d, method = 'lambda')$lambda
    blambda <- c(blambda, perm_lambda)
  }
  
  # Calculate permutation p-value: proportion of permuted lambda greater than observed lambda
  pval_value <- length(blambda[blambda > lambda_value]) / length(blambda)
  pval <- c(pval, pval_value)
}

# Create a dataframe summarizing lambda and permutation p-values per trait
lambda_df <- data.frame(
  trait = colnames(paperdatatree),
  lambda = lambda,
  p_permutation = pval
)

# Save results to CSV file for downstream use
library(data.table)
fwrite(lambda_df, './Phylosig.csv')

# Print results to console
print(lambda_df)

# -----------------------------------------------
# 5. Coefficient of variation analysis
# -----------------------------------------------

#We want to answer if phenotypic diversity is higher in populations than among populations

# PERMANOVA and Phenotypic Coefficient of Variation (CV) analysis

# Load required libraries
library(dplyr)
library(ggplot2)
library(dunn.test)
library(vegan)  # for adonis2 function

# --- PERMANOVA ---

# Calculate distance matrix from phenotypic data
matriz <- dist(paperdata)
matriz <- as.matrix(matriz)

# Prepare metadata dataframe with population information
permapaperdata <- as.data.frame(paperdata)
permapaperdata <- permapaperdata[, 1:34]  # Select first 34 traits
permapaperdata$pop <- popgensub$FinePop2

# Ensure rownames of distance matrix match rows in metadata
rownames(matriz) <- rownames(permapaperdata)

# Perform PERMANOVA to test phenotypic differences among populations
resultados_permanova <- adonis2(matriz ~ pop, data = permapaperdata, permutations = 999)

# Print PERMANOVA results
print(resultados_permanova)

# --- Coefficient of Variation (CV) ---

# Subset paperdata to first 34 phenotypes for CV analysis
paperdata <- paperdata[, 1:34]

# Calculate within-population CV for each phenotype
cv_within <- paperdata %>%
  mutate(Population = popgensub$FinePop2) %>%
  group_by(Population) %>%
  summarise(across(everything(), ~ sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE))) %>%
  pivot_longer(-Population, names_to = "Phenotype", values_to = "CV") %>%
  mutate(Type = "Within Population")

# Calculate between-population CV for each phenotype
cv_between <- paperdata %>%
  summarise(across(everything(), ~ sd(tapply(.x, popgensub$FinePop2, mean, na.rm = TRUE), na.rm = TRUE) / 
                     mean(tapply(.x, popgensub$FinePop2, mean, na.rm = TRUE)))) %>%
  pivot_longer(everything(), names_to = "Phenotype", values_to = "CV") %>%
  mutate(Type = "Between Populations")

# Combine within and between population CV data into one dataframe
cv_data <- bind_rows(cv_within, cv_between)

# Use absolute values for CV (interpretation in absolute terms)
cv_data$CV <- abs(cv_data$CV)

# Replace NA population names with "Interpopulation" for between-population CV
cv_data$Population <- replace(cv_data$Population, is.na(cv_data$Population), "Interpopulation")

# Set factor levels for population order on the x-axis
cv_data$Population <- factor(cv_data$Population, 
                             levels = c("NoAm", "W29 clade", "Hyp", "So-HC", "Eur1", "Eur2", "Interpopulation"))

# --- Plotting ---

p <- ggplot(cv_data, aes(x = Population, y = CV, fill = Type)) +
  geom_boxplot(
    outlier.shape = NA,  # Hide outliers
    alpha = 0.7,         # Transparency for boxplots
    size = 0.25          # Line thickness
  ) +
  geom_jitter(
    width = 0.2, size = 0.8, color = "black", alpha = 0.8  # Jittered points for data
  ) +
  scale_fill_manual(values = c("Within Population" = "skyblue", 
                               "Between Populations" = "salmon")) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    title = "Interpopulation and Intrapopulation phenotypic CV",
    x = "Populations",
    y = "Coefficient of Variation (CV)"
  ) +
  theme_classic(base_size = 8) +  # Smaller base font size
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 9, face = "bold"),
    axis.line = element_line(size = 0.25),  # Thinner axis lines
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5, unit = "mm")
  ) +
  coord_cartesian(clip = "off")  # Prevent clipping of elements outside plot area

print(p)

# Save plot as PDF with publication-quality size and resolution
ggsave("nodiploid_boxplot_CV_pop.pdf", plot = p, 
       width = 24, height = 18, dpi = 300, units = "cm", device = cairo_pdf)

# --- Statistical tests ---

# Shapiro-Wilk test for normality of CV distribution
shapiro_result <- shapiro.test(cv_data$CV)
print(shapiro_result)

# Kruskal-Wallis test to compare CV between 'Type' groups (Within vs Between populations)
kruskal_result <- kruskal.test(CV ~ Type, data = cv_data)
print(kruskal_result)

# Post-hoc Dunn test for pairwise comparisons between populations, with Bonferroni correction
dunn_results <- dunn.test(cv_data$CV, cv_data$Population, method = "bonferroni")
print(dunn_results)

# -----------------------------------------------
# 6. Heatmap representation of significant conditions
# -----------------------------------------------

# Significant conditions:

condiliqfil <- popconditionsliq


# Filter the paperdata dataframe to keep only the selected conditions
sigpaperdata <- paperdata %>%
  select(all_of(condiliqfil))

# Add the Name column from popgensub as an extra column to the filtered data
sigpaperdatastrain <- cbind(sigpaperdata, popgensub$Name)

# Load library to write Excel files
library(readxl)
write.xlsx(sigpaperdatastrain, file  = "sigpaperdataliqpop.xlsx")

# Split the filtered data into two matrices:
# inhmat for environmental stress conditions (inhibitors)
inhmat <- sigpaperdata[, colnames(sigpaperdata) %in% mediosinh]

# nutmat for nutritional conditions
nutmat <- sigpaperdata[, colnames(sigpaperdata) %in% nutrient]

# Heatmap preparation

# Create annotation dataframe from popgensub using Name and FinePop2 columns
anotacion <- as.data.frame(popgensub[, c("Name", "FinePop2")])

# Extract the first column as row names for the annotation dataframe
nombres_fila <- anotacion[, 1]

# Set rownames of annotation dataframe to these names
rownames(anotacion) <- nombres_fila

# Remove the first column (Name) from the annotation dataframe
anotacion <- anotacion[, -1]

# Ensure annotation is a dataframe and convert the annotation column to factor
anotacion <- as.data.frame(anotacion)
anotacion$anotacion <- as.factor(anotacion$anotacion)

# Rename the annotation column
colnames(anotacion) <- c("FinePop2")

# Make sure rownames have no leading/trailing spaces
rownames(anotacion) <- trimws(rownames(anotacion))

# Set rownames of annotation to match collect$Name (assumed master list)
rownames(anotacion) <- collect$Name

# Verify if rownames of nutmat and annotation match exactly
identical(rownames(nutmat), rownames(anotacion))

# Define colors for annotation categories
anno_col <- list(
  FinePop2 = c("Eur1" = "#F29E6D", "Eur2" = "#05A6A6", "So-HC" = "#03738C", 
               "NoAm" = "#F23E2E", "W29 clade" = "#152B59", "Hyp" = "grey") 
)

# Define color palette for heatmap (reverse RdYlBu palette)
mi_paleta <- rev(brewer.pal(n = 10, name = "RdYlBu"))

# Define breaks for color intervals in heatmap
intervalos <- seq(from = -6, to = 4, length.out =10)

# Transpose nutritional matrix for heatmap (rows: conditions, columns: strains)
nutmat <- t(nutmat)

# Generate PDF for nutritional conditions heatmap
pdf("nodiploid_significant_paper_Heatmap_sources_pop_sep.pdf", width = 15, height = 10)
pheatmap(nutmat, cluster_rows = TRUE, cluster_cols = FALSE,
         display_numbers = FALSE, breaks = intervalos, main = "",
         format = "f", color = mi_paleta, annotation_col = anotacion,
         annotation_colors = anno_col,
         fontsize_row = 8,  # Adjust row label font size
         fontsize_col = 7,  # Adjust column label font size
         angle_col = 90,
         annotation_col_margin = 25,
         cellwidth = 8,
         cellheight = 10,
         border_color = "snow3"  # Remove cell border lines
)
dev.off()

# Repeat the same for inhibitors heatmap but with different palette and breaks
mi_paleta <- rev(brewer.pal(n = 7, name = "RdYlBu"))
intervalos <- seq(from = -3, to = 3, length.out =7)

inhmat <- t(inhmat)  # Transpose inhibitors matrix

pdf("nodiploid_significant_paper_Heatmap_inh_pop_nosep_2nd.pdf", width = 15, height = 10)
pheatmap(inhmat, cluster_rows = TRUE, cluster_cols = FALSE,
         display_numbers = FALSE, breaks = intervalos, main = "",
         format = "f", color = mi_paleta, annotation_col = anotacion,
         annotation_colors = anno_col,
         fontsize_row = 8,
         fontsize_col = 7,
         angle_col = 90,
         annotation_col_margin = 25,
         cellwidth = 8,
         cellheight = 10,
         border_color = "snow3"
)
dev.off()

# Lipase activity heatmap preparation

# Select columns starting with "lip" from filtered data
lip <- sigpaperdata %>% 
  select(starts_with("lip"))

# Rename columns for clarity
colnames(lip) <- c("Glucose", "Glycerol", "Urea", "pH 7", "CN High")

# Transpose lipase matrix
tlip <- t(as.matrix(lip))

mi_paleta <- rev(brewer.pal(n = 6, name = "RdYlBu"))
intervalos <- seq(from = -2, to = 4, length.out =6)

# Generate PDF for lipase activity heatmap
pdf("nodiploid_significant_paper_Heatmap_lip_pop_sep.pdf", width = 15, height = 10)
pheatmap(tlip, cluster_rows = TRUE, cluster_cols = FALSE,
         display_numbers = FALSE, breaks = intervalos, main = "Lipase activity on different media",
         format = "f", color = mi_paleta, annotation_col = anotacion,
         annotation_colors = anno_col,
         fontsize_row = 8,
         fontsize_col = 7,
         angle_col = 90,
         annotation_col_margin = 25,
         cellwidth = 6,
         cellheight = 10,
         border_color = "snow3"
)
dev.off()

#Colony morphology and protease activity

sigsoltostay <- c("Proteolitic_SKM") 
sigpaperdatasol <- paperdatasol %>%
  select(all_of(sigsoltostay))

# Creating matrix

protmorph <- as.data.frame(sigpaperdatasol$Proteolitic_SKM)
protmorph$morph <- popgensub$morph
rownames(protmorph) <- rownames(sigpaperdatasol)
colnames(protmorph) <- c("Proteolytic_activity", "Colony_morphology")

#Heatmap representation

anno_col <- list(
  FinePop2 = c("Eur1" = "#F29E6D", "Eur2" = "#05A6A6", "So-HC" = "#03738C", 
               "NoAm" = "#F23E2E", "W29 clade" = "#152B59", "Hyp" = "grey"))

trans<- t(as.matrix(protmorph))
# colnames(trans) <-rownames(solid)
mi_paleta <- rev(brewer.pal(n = 4, name = "RdYlBu"))
intervalos <- c(-1,0,1,2, 3)

pdf("nodiploid_paper_Heatmap_potmorph_pop_sep.pdf", width = 15, height = 10)
pheatmap(trans, cluster_rows = TRUE, cluster_cols = FALSE,
         display_numbers = FALSE, breaks = intervalos, main = "",
         format = "f", color = mi_paleta, annotation_col = anotacion,
         annotation_colors = anno_col,
         fontsize_row = 8,  # Ajusta el tamaño de las etiquetas de las filas según tu preferencia
         fontsize_col = 7,  # Ajusta el tamaño de las etiquetas de las columnas según tu preferencia
         angle_col = 90,
         annotation_col_margin = 25,
         cellwidth = 8,
         cellheight = 10,
         border_color = "snow3"  # Eliminar las líneas de separación entre celdas
)
dev.off()

# -----------------------------------------------
# 7. Wilcoxon tests for trait associations in Eur2/Dairy and So-HC
# -----------------------------------------------

#Do this test for both So-HC strains and Eur2/Dairy strain

# 1. Get all unique non-NA levels of the factor 'Eur2' in popgensub

niveles <- unique(na.omit(popgensub$Eur2))

# Generate all pairwise combinations of these levels
combinaciones <- combn(niveles, 2)

# 2. Apply Wilcoxon test for each pair of levels
resultados <- apply(combinaciones, 2, function(par) {
  # Select values of 'Tributyrin' for each group in the pair
  grupo1 <- paperdata$Tributyrin[popgensub$dairy == par[1]]
  grupo2 <- paperdata$Tributyrin[popgensub$dairy == par[2]]
  
  # Perform two-sided Wilcoxon test between the two groups
  test <- wilcox.test(grupo1, grupo2, alternative = "two.sided")
  
  # Store results in a data frame with medians and p-value
  data.frame(
    Grupo1 = par[1],
    Grupo2 = par[2],
    p_value = test$p.value,
    Mediana1 = median(grupo1),
    Mediana2 = median(grupo2)
  )
})

# 3. Combine all test results into one data frame
resultados_finales <- do.call(rbind, resultados)

# Print the results table
print(resultados_finales)

# Create contingency table between 'dairy' groups and 'Proteolitic_SKM' phenotype
table(popgensub$dairy,paperdatasol$Proteolitic_SKM)

# Perform Chi-square test of independence on this contingency table
chisq.test(table(popgensub$dairy,paperdatasol$Proteolitic_SKM)) 

# Associations with the clade 'So-HC'

# Reorder rows of paperdata to match popgensub$Name order
paperdata <- paperdata[match(popgensub$Name, rownames(paperdata)), ]

# Extract the grouping variable 'So-HC' from popgensub
grupo <- popgensub$`So-HC`

# Check that rownames of paperdata exactly match popgensub$Name
stopifnot(all(rownames(paperdata) == popgensub$Name))

# Apply Wilcoxon test for each column of paperdata comparing the two groups in 'grupo'
wilcoxon_results <- apply(paperdata, 2, function(x) {
  resultado <- wilcox.test(x ~ grupo)
  c(
    p.value = resultado$p.value,
    median_group1 = median(x[grupo == unique(grupo)[1]], na.rm = TRUE),
    median_group2 = median(x[grupo == unique(grupo)[2]], na.rm = TRUE)
  )
})

# Convert results to data.frame and sort by adjusted p-value (FDR correction optional)
wilcoxon_df <- as.data.frame(t(wilcoxon_results))
wilcoxon_df$padj <- p.adjust(wilcoxon_df$p.value, method = "fdr")
wilcoxon_df <- wilcoxon_df[order(wilcoxon_df$padj), ]

# Select condition names with p-value < 0.05
SoHCcond <- rownames(subset(wilcoxon_df, wilcoxon_df$p.value < 0.05)) #for So-HC, can be repeated for Eur2/Dairy


# Combine popgensub with paperdata to create a single dataframe for plotting
popgensub1 <- cbind(popgensub, paperdata)

#Selected conditions significantly associated with So-HC and Eur2/Dairy
#Conditions of interest
cond <- c("Lactate", "Glutamine", "4% NaCl", "Acetate", 
          "Glycerol", "lipase_pH 7","lipase_CN High", "Hexadecane", "Tributyrin", "SDS") #Eur2/Dairy

#Join data and metadata
popgensub1 <- cbind(popgensub, paperdata)

#Subset phenotypic data accordingly for each case:
subpca <- as.matrix(subset(popgensub1, select = colnames(popgensub1) %in% SoHCcond)) #So-HC
subpca <- as.matrix(subset(popgensub1, select = colnames(popgensub1) %in% cond)) #Eur2/Dairy

# Subset and reshape to long format for ggplot
df <- as.data.frame(subpca[, SoHCcond, drop = FALSE]) %>%
  mutate(
    Strain = rownames(subpca),  # Add strain names as a column
    SoHC = popgensub1$`So-HC`   # Group assignment ("So-HC" vs "Other")
  ) %>%
  pivot_longer(
    cols = all_of(SoHCcond), 
    names_to = "Trait", 
    values_to = "Value"
  )

# Reorder levels of 'Trait' by median value within the "So-HC" group
ordered_traits <- df %>%
  filter(SoHC == "So-HC") %>%
  group_by(Trait) %>%
  summarise(median_value = median(Value, na.rm = TRUE)) %>%
  arrange(desc(median_value)) %>%
  pull(Trait)

# Apply ordered levels to the 'Trait' factor
df$Trait <- factor(df$Trait, levels = ordered_traits)

# Define color palettes (unchanged)
coleur2   <- c("Eur2" = "#05A6A6", "Other" = "grey4")
coldairy  <- c("Yes" = "#3EA8FA", "No" = "grey4", "NA" = "white")
col2dairy <- c("Eur2 dairy" = "#05A6A6", "Eur1 dairy" = "#F29E6D", "Other" = "grey")
colsohc   <- c("So-HC" = "#03738C", "Other" = "grey")


#Plot accordingly in each case (this version is prepared for SO-HC:
p <- ggplot(df, aes(x = Trait, y = Value)) +
  geom_jitter(
    aes(fill = SoHC, color = SoHC),
    size = 2.2, alpha = 0.5,
    shape = 21, stroke = 0.3,
    position = position
  ) +
  geom_boxplot(
    aes(group = interaction(Trait, SoHC)), 
    outlier.shape = NA, width = 0.5,
    color = "black", fill = NA, size = 1,
    position = position_dodge(width = 0.7)
  ) +
  scale_fill_manual(values = colsohc) +
  scale_color_manual(values = colsohc) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    y = "Growth",
    fill = "So-HC",
    color = "So-HC"
  ) +
  ggtitle("Phenotypic traits linked So-HC")

print(p)

# -----------------------------------------------
# 8. Pehnotypical analysis of admixed populations
# -----------------------------------------------

# Do admixed groups grow better in more conditions?
# -------------------------------
# Generalism hypothesis analysis
# -------------------------------

# Load required libraries
library(tidyverse)
library(ggpubr)
library(dunn.test)

# Subset growth data (first 29 columns of paperdata)
growth <- paperdata[, 1:29]

# Add strain names as a column
growth$Name <- rownames(growth)

# Reshape to long format and join with population group information
growth_long <- growth %>%
  pivot_longer(-Name, names_to = "Condition", values_to = "Growth") %>%
  left_join(popgensub %>% select(Name, FinePop2), by = "Name")

# Assign broader population groupings
growth_long$FinePop2 <- as.factor(growth_long$FinePop2)
growth_long <- growth_long %>%
  mutate(Groups = case_when(
    FinePop2 %in% c("NoAm", "W29 clade", "Eur2", "So-HC") ~ "Structured",
    FinePop2 %in% c("Hyp", "Eur1") ~ "Admixed",
    TRUE ~ NA_character_  # for unassigned populations
  ))

# Define color palette for population groups
colgroups <- c("Admixed" = "#FADD39", "Structured" = "#6159BA")

# Define jitter position for better visualization
position <- position_jitterdodge(jitter.width = 1, dodge.width = 0.7)

# Plot growth across population groups
p <- ggplot(growth_long, aes(x = Groups, y = Growth)) +
  # Add individual points with jitter
  geom_jitter(
    aes(fill = Groups, color = Groups),
    size = 2.2, alpha = 0.5,
    shape = 21, stroke = 0.3,
    position = position
  ) +
  # Overlay transparent boxplots
  geom_boxplot(
    aes(group = interaction(Groups)), 
    outlier.shape = NA, width = 0.5,
    color = "black", fill = NA, size = 1,
    position = position_dodge(width = 0.7)
  ) +
  # Apply manual colors
  scale_fill_manual(values = colgroups) +
  scale_color_manual(values = colgroups) +
  # Add theme and formatting
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    y = "Growth",
    fill = "Population group",
    color = "Population group"
  ) +
  ggtitle("Growth across conditions by population group") +
  
  # Global Kruskal-Wallis test
  stat_compare_means(
    method = "kruskal.test", 
    label.y = max(growth_long$Growth, na.rm = TRUE) * 1.05
  ) +
  
  # Pairwise Dunn test with Bonferroni correction
  stat_compare_means(
    method = "dunn.test", 
    label = "p.signif", 
    p.adjust.method = "bonferroni",
    comparisons = list(c("Admixed", "Structured")),
    hide.ns = TRUE, 
    label.y = max(growth_long$Growth, na.rm = TRUE) * 1.1
  )

# Print plot
print(p)

#Are the more frequently found among top growers in each condition?

# -----------------------------------------
# Proportion of best growers by population
# -----------------------------------------

# Load required library
library(forcats)

# Calculate the total number of strains per FinePop2 group
group_sizes <- popgensub %>%
  group_by(FinePop2) %>%
  summarise(total_strains = n())

# Define the top performer threshold (top 1%)
top_percent <- 0.01

# Identify top growers per condition (top 1% per condition)
top_growers <- growth_long %>%
  group_by(Condition) %>%
  mutate(rank_growth = rank(-Growth, ties.method = "min")) %>%
  filter(rank_growth <= ceiling(n() * top_percent)) %>%
  ungroup() %>%
  count(FinePop2, name = "top_count")

# Join with total group sizes and normalize frequency
top_growers_norm <- top_growers %>%
  left_join(group_sizes, by = "FinePop2") %>%
  mutate(freq_normalized = top_count / total_strains) %>%
  arrange(desc(freq_normalized))

# Print normalized top grower frequency table
print(top_growers_norm)

# Reorder FinePop2 factor levels to match desired population order
top_growers_norm <- top_growers_norm %>%
  mutate(FinePop2 = factor(FinePop2, levels = orden_pops))

# Define bar plot of normalized top growers per population
p <- ggplot(top_growers_norm, aes(x = FinePop2, y = freq_normalized, fill = FinePop2)) +
  geom_col(color = "black", width = 0.7) +   # Bar height reflects normalized frequency
  scale_fill_manual(values = colpop2) +
  theme_classic(base_size = 15) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  ) +
  labs(
    y = "Proportion of top growers",
    title = "Relative frequency of top growers"
  ) +
  ylim(0, max(top_growers_norm$freq_normalized, na.rm = TRUE) * 1.1)  # Add headroom

# Display the plot
print(p)
