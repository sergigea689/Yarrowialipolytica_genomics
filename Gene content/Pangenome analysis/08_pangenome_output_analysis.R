#Script for the analysis of the pangenome matrix and the accessory gene content
#Load libraries
library(ggplot2)
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tidyr)

#Load data

pangenomecomplete <- read_excel("nonbasalpangenomecomplete.xlsx") #Accesory genome (anything non basal core/not present in every strain
pangenome_matrix_clean <- read_excel("pangenome_matrix_clean_whole.xlsx") #complete genome
popgen <- read_excel("update_manfinetreeposition.xlsx") #metadata
popgensub <- subset(popgen, popgen$treesamp %in% rownames(pangenomecomplete))

# ------------------------------------------------
# 1. Perform Principal Component Analysis (PCA)
# ------------------------------------------------

# Remove strains with extreme outlier values
outliers <- c("ERR5235159_")
filtered_matrix <- pangenomecomplete[!(rownames(pangenomecomplete) %in% outliers), ]

# Perform PCA without scaling variables (set scale. = TRUE to standardize)
pca_result <- prcomp(filtered_matrix, scale. = FALSE)

# 2. Summarize PCA results
summary(pca_result)

# Extract PCA scores for each sample
pca_scores <- pca_result$x

# Calculate proportion of variance explained by each principal component
variance_explained <- (pca_result$sdev)^2 / sum((pca_result$sdev)^2)

# 3. Prepare data for visualization

# Subset metadata to match PCA samples and order accordingly
pca_metadata <- subset(popgensub, treesamp %in% rownames(pca_scores))
pca_metadata <- pca_metadata[order(match(pca_metadata$treesamp, rownames(pca_scores))), ]

# Add population and morphology information to PCA scores
pca_scores <- as.data.frame(pca_scores)
pca_scores$Population <- factor(pca_metadata$FinePop2)
pca_scores$Morphology <- factor(pca_metadata$morph)

# Define custom color palette for populations
population_colors <- c(
  "Eur1" = "#F29E6D", 
  "Eur2" = "#05A6A6", 
  "So-HC" = "#03738C", 
  "NoAm" = "#F23E2E", 
  "W29 clade" = "#152B59", 
  "Hyp" = "grey"
)

# 4. Plot PCA results (PC1 vs PC2)
library(ggplot2)

p <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Population)) +
  geom_point(size = 4, alpha = 0.8) +
  # Uncomment below to add sample labels:
  # geom_text(aes(label = rownames(pca_scores)), vjust = -0.5, hjust = 0.5, size = 2) +
  scale_color_manual(values = population_colors) +
  labs(
    x = sprintf("PC1 (%.1f%%)", variance_explained[1] * 100),
    y = sprintf("PC2 (%.1f%%)", variance_explained[2] * 100),
    color = "Population"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p)

# Save plot as PDF with high resolution for publication
ggsave(
  filename = "paper_PCA_noncorebasalperstudy.pdf",
  plot = p,
  device = "pdf",
  width = 8,
  height = 6,
  dpi = 300,
  units = "in"
)

# 5. Investigate variable loadings on PC2

# Extract PCA loadings (variable contributions)
pca_loadings <- pca_result$rotation

# Select loadings for the second principal component (PC2)
pc2_loadings <- pca_loadings[, 2]

# Sort absolute loadings in descending order to identify top contributors
pc2_loadings_sorted <- sort(abs(pc2_loadings), decreasing = TRUE)

print(pc2_loadings_sorted)

# Visualize variable contributions to PC2 using a barplot
barplot(
  pc2_loadings_sorted,
  main = "Variable Contributions to PC2",
  las = 2,
  col = "skyblue",
  cex.names = 0.7
)

# 6. Notable contributors to PC2

# The top contributors include:
# - YL_90_EKDN230056791.1A_HWT7HDSX7_L4_..Seq68, with 100% homology to mating-type protein A2 in Yarrowia.
# - uniprot|Q6CCQ7 Yarrowia lipolytica YALI0C07436g APN2 DNA-(apurinic or apyrimidinic site) lyase, an endonuclease.
# - uniprot|Q6CCQ6 Yarrowia lipolytica YALI0C07458g MATB2 Mating type B protein 2.
# - uniprot|Q707X9 Yarrowia lipolytica YALI0C07480g MATB1 Mating type B protein 1.
# - ERR5235104_::Seq37 MATA1.
# - uniprot|Q6CCQ4 Yarrowia lipolytica YALI0C07502g SLA2 Sla2p cytoskeleton assembly control protein.
# - Similarities with uniprot|P09230 Yarrowia lipolytica XPR2 Alkaline extracellular protease precursor.

#  ----------------------------------------
# 2.PERMANOVA analysis
#  ----------------------------------------

# Calculate distance matrix
distance_matrix <- dist(pangenomecomplete)

# Filter metadata to match the samples in the distance matrix
pop_metadata <- subset(popgensub, popgensub$treesamp %in% rownames(pangenomecomplete))

# Ensure sample order consistency between distance matrix and metadata
permanova_input <- pangenomecomplete[rownames(pangenomecomplete) %in% pop_metadata$treesamp, ]
permanova_input$Population <- factor(pop_metadata$FinePop2)

# Perform PERMANOVA
library(vegan)
permanova_results <- adonis2(distance_matrix ~ Population, data = permanova_input, permutations = 999)

# Display results
print(permanova_results)
#  ----------------------------------------
# 3.Looking for population patterns
#  ----------------------------------------

# Reorder popgensubperma to match the rows in the pangenome matrix
popgensubperma <- popgensubperma[match(rownames(pangenomecomplete), popgensubperma$treesamp), ]

# Add the population group to the matrix
pangenome_df <- as.data.frame(pangenomecomplete)
pangenome_df$grupo <- popgensubperma$FinePop2

# Sum gene presence (1s) and absence (0s) per group using base R
genes_presentes_por_grupo <- aggregate(. ~ grupo, data = pangenome_df, FUN = sum)
genes_ausentes_por_grupo <- aggregate(. ~ grupo, data = pangenome_df, FUN = function(x) sum(x == 0))

# Preview presence matrix (first 5 columns)
head(genes_presentes_por_grupo[, 1:5]) 

# Count number of strains per group (used to normalize)
cepas_por_grupo <- matriz_con_grupo %>%
  group_by(FinePop2) %>%
  summarise(n = n())
colnames(cepas_por_grupo) <- c("grupo", "n")

# Normalize presence and absence by strain count
genes_presentes_prop <- genes_presentes_por_grupo
genes_presentes_prop[,-1] <- genes_presentes_prop[,-1] / cepas_por_grupo$n

genes_ausentes_prop <- genes_ausentes_por_grupo
genes_ausentes_prop[,-1] <- genes_ausentes_prop[,-1] / cepas_por_grupo$n

# Calculate total number of present genes per group
present_total <- genes_presentes_por_grupo %>%
  rowwise() %>%
  mutate(total_present = sum(c_across(-grupo))) %>%
  select(grupo, total_present)

# Normalize total presence
present_total_norm <- present_total %>%
  left_join(cepas_por_grupo, by = "grupo") %>%
  mutate(norm_present = total_present / n)

# Calculate total number of absent genes per group
ausent_total <- genes_ausentes_por_grupo %>%
  rowwise() %>%
  mutate(total_absent = sum(c_across(-grupo))) %>%
  select(grupo, total_absent)

# Normalize total absence
ausent_total_norm <- ausent_total %>%
  left_join(cepas_por_grupo, by = "grupo") %>%
  mutate(norm_absent = total_absent / n)

# Combine data into a single dataframe for ggplot
df_plot_norm <- left_join(present_total_norm[, c("grupo", "norm_present")],
                          ausent_total_norm[, c("grupo", "norm_absent")],
                          by = "grupo") %>%
  arrange(desc(norm_present)) %>%
  mutate(grupo = factor(grupo, levels = grupo)) %>%
  pivot_longer(cols = starts_with("norm_"),
               names_to = "Status", values_to = "Normalized_count") %>%
  mutate(Status = recode(Status,
                         norm_present = "Genes present",
                         norm_absent = "Genes absent"))

# Plot normalized gene presence and absence per group
ggplot(df_plot_norm, aes(x = grupo, y = Normalized_count, fill = Status)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Genes present" = "#1E8BFA", "Genes absent" = "#FA2B1E")) +
  theme_classic() +
  labs(x = "Group", y = "Normalized gene count (per strain)", fill = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14)
  )


#  ----------------------------------------
# 4.Proportion of types of genes in each genetic group
#  ----------------------------------------

# --- Classification of genes into accessory genome categories ---
# Calculate global frequency of each gene (presence across strains)
# Proportion of each accessory genome class per population

# Classify genes based on their global presence across strains
presencia_gen <- colSums(pangenome_df[,1:384]) / nrow(pangenome_df[,1:384])
clases_gen <- cut(presencia_gen,
                  breaks = c(-Inf, 0.05, 0.95, Inf),
                  labels = c("cloud", "shell", "soft-core"))

# Dataframe mapping each gene to its accessory class
gen_clase_df <- data.frame(
  gen = names(presencia_gen),
  clase = clases_gen,
  presencia_global = presencia_gen
)

# Retrieve gene names for each class
soft_core_genes <- gen_clase_df$gen[gen_clase_df$clase == "soft-core"]
shell_genes     <- gen_clase_df$gen[gen_clase_df$clase == "shell"]
cloud_genes     <- gen_clase_df$gen[gen_clase_df$clase == "cloud"]

# Function to calculate the proportion of genes per class in each strain
proporcion_clase <- function(df, genes) {
  if (length(genes) == 0) return(rep(0, nrow(df)))
  rowSums(df[, genes, drop = FALSE]) / rowSums(df)
}

# Calculate proportions per class for each strain
prop_soft_core <- proporcion_clase(pangenome_df[,1:384], soft_core_genes)
prop_shell     <- proporcion_clase(pangenome_df[,1:384], shell_genes)
prop_cloud     <- proporcion_clase(pangenome_df[,1:384], cloud_genes)

# Create dataframe with class proportions and strain metadata
proporciones_df <- data.frame(
  cepa      = rownames(pangenome_df),
  soft_core = prop_soft_core,
  shell     = prop_shell,
  cloud     = prop_cloud
)

# Reorder rows to match metadata structure
pangenome_df <- pangenome_df[match(popgensubperma$treesamp, rownames(pangenome_df)), ]

# Update class proportions with population group information
proporciones_df <- data.frame(
  cepa      = popgensubperma$treesamp,
  grupo     = popgensubperma$FinePop2,
  soft_core = prop_soft_core,
  shell     = prop_shell,
  cloud     = prop_cloud
)

# Compute mean proportion of each genome class per population
proporciones_grupo <- proporciones_df %>%
  group_by(grupo) %>%
  summarise(
    mean_soft_core = mean(soft_core, na.rm = TRUE),
    mean_shell     = mean(shell, na.rm = TRUE),
    mean_cloud     = mean(cloud, na.rm = TRUE)
  )

# Reshape to long format for plotting
proporciones_long <- proporciones_grupo %>%
  pivot_longer(
    cols = starts_with("mean_"),
    names_to = "tipo_genoma",
    values_to = "proporcion"
  ) %>%
  mutate(
    tipo_genoma = recode(tipo_genoma,
                         mean_soft_core = "Soft-core",
                         mean_shell     = "Shell",
                         mean_cloud     = "Cloud")
  )

# Barplot: mean proportions of accessory genome types per population
orden_grupos <- proporciones_long %>%
  filter(tipo_genoma == "Cloud") %>%
  arrange(desc(proporcion)) %>%
  pull(grupo)

# Order
proporciones_long$grupo <- factor(proporciones_long$grupo, levels = orden_grupos)

ggplot(proporciones_long, aes(x = grupo, y = proporcion, fill = tipo_genoma)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Soft-core" = "#1b9e77", "Shell" = "#7570b3", "Cloud" = "#d95f02")) +
  theme_classic() +
  labs(
    x = "Population group",
    y = "Mean proportion of genes",
    fill = "Accessory genome class",
    title = "Mean proportion of accessory genome classes across populations"
  ) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y     = element_text(size = 12),
    axis.title      = element_text(size = 14, face = "bold"),
    legend.title    = element_text(size = 13),
    legend.text     = element_text(size = 12)
  )

# Boxplot: comparison of cloud genome proportions per population
my_colors_pop <- c(
  "Eur1"      = "#F29E6D",
  "Eur2"      = "#05A6A6",
  "So-HC"     = "#03738C",
  "NoAm"      = "#F23E2E",
  "W29 clade" = "#152B59",
  "Hyp"       = "grey"
)

ggplot(proporciones_df, aes(x = grupo, y = cloud, color = grupo)) +
  geom_boxplot(
    fill = NA,
    color = "black",
    size = 1.2,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15, height = 0,
    size = 2.5,
    alpha = 0.8,
    aes(color = grupo)
  ) +
  # stat_compare_means(...)  # Optional significance testing
  scale_color_manual(values = my_colors_pop) +
  theme_classic(base_size = 16) +
  labs(
    x = "",
    y = "Proportion of cloud genes",
    title = "Cloud gene proportions per genetic group"
  ) +
  ylim(0, 0.11) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
    axis.text.y     = element_text(size = 14, color = "black"),
    axis.title      = element_text(size = 16, face = "bold"),
    legend.position = "none",
    panel.border    = element_rect(color = "black", fill = NA, size = 1)
  )

# Dunn post-hoc test after Kruskal-Wallis to compare cloud gene proportions
kruskal.test(soft_core ~ grupo, data = proporciones_df)

dunn_result <- dunnTest(cloud ~ grupo, data = proporciones_df, method = "bh")
print(dunn_result)


#  ----------------------------------------
# 4.Number of unique exclusive genes per population
#  ----------------------------------------

# Transpose the pangenome matrix to have genes as rows
pangenome_transposed <- as.data.frame(t(pangenome_df[, -ncol(pangenome_df)]))  # exclude FinePop2

# Retrieve unique population groups
pop_list <- unique(popgensubperma$FinePop2)

# Initialize list to store exclusive genes for each population
genes_unicos <- list()

# Loop through each population to identify unique genes
for (pop in pop_list) {
  
  # Identify strains belonging to the current population
  strains_in_pop <- popgensubperma$treesamp[popgensubperma$FinePop2 == pop]
  
  # Identify strains from other populations
  strains_out_pop <- setdiff(rownames(pangenome_df), strains_in_pop)
  
  # Genes present in at least one strain of the current population
  present_in_pop <- rowSums(pangenome_transposed[, strains_in_pop, drop = FALSE]) > 0
  
  # Genes completely absent in all strains of other populations
  absent_out_pop <- rowSums(pangenome_transposed[, strains_out_pop, drop = FALSE]) == 0
  
  # Identify genes exclusively present in the current population
  unique_genes <- names(which(present_in_pop & absent_out_pop))
  
  # Store result
  genes_unicos[[pop]] <- unique_genes
}

# Exclude genes likely to be contaminants (e.g., from ERR5235159)
genes_unicos_filtrados <- lapply(genes_unicos, function(glist) {
  glist[!grepl("ERR5235159_", glist)]
})

# Count the number of unique genes per population
genes_unicos_counts <- sapply(genes_unicos_filtrados, length)

# Format into dataframe for plotting
genes_unicos_df <- data.frame(
  grupo = names(genes_unicos_counts),
  n_genes_unicos = as.numeric(genes_unicos_counts)
)

# Barplot showing the number of unique genes per population
ggplot(genes_unicos_df, aes(x = reorder(grupo, -n_genes_unicos), y = n_genes_unicos, fill = grupo)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", size = 1) +  # bars with black borders
  scale_fill_manual(values = my_colors_pop) +
  theme_classic(base_size = 16) +
  labs(
    x = "",
    y = "Number of genes",
    title = "Unique genes per genetic group"
  ) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
    axis.text.y     = element_text(size = 14, color = "black"),
    axis.title      = element_text(size = 16, face = "bold"),
    legend.position = "none",
    panel.border    = element_rect(color = "black", fill = NA, size = 1)  # border for publication style
  )

# Analysis of the strain ERR5235159_ suspected of contributing spurious genes

# Identify columns associated with the suspect strain
cols_err <- grep("ERR5235159_", colnames(pangenome_df), value = TRUE)

# Filter genes present in this strain
genes_en_ERR <- pangenome_df %>%
  filter(if_any(all_of(cols_err), ~ . == 1)) %>%
  select(-starts_with("gene_id"))  # optional: exclude gene ID column if present

# Count how many times each gene appears as '1' in ERR5235159_
conteo_1_ERR <- rowSums(pangenome_df[, cols_err] == 1)
conteo_1_ERR <- as.data.frame(conteo_1_ERR)

# Display the first rows of the count
head(conteo_1_ERR)


#  ----------------------------------------
# 5.Saturation curve
#  ----------------------------------------

library(vegan)
library(micropan)

# Calculate pangenome saturation curve using rarefaction with 100% confidence interval
saturation_curve <- specaccum(pangenome_matrix_clean, method = "rarefaction", ci = 1)

# Plot saturation curve with confidence interval polygon in grey
plot(saturation_curve, main = "Pangenome Saturation Curve", ci = 4, ci.type = "polygon", col = "grey")

# Model pangenome openness with Heap's law and 1000 permutations
heap_model <- heaps(pangenome_matrix_clean, n.perm = 1000)

# Print fitted parameters of Heap model
print(heap_model)

# Create data frame for pangenome size and standard deviation from saturation curve
pan_genome_df <- data.frame(
  Size = saturation_curve$sites,
  PanGenome = saturation_curve$richness,
  SD_PanGenome = saturation_curve$sd
)

# Plot saturation curve using ggplot2 with confidence ribbons
saturation_plot <- ggplot(pan_genome_df, aes(x = Size, y = PanGenome)) +
  geom_line(color = "red", size = 1.2) +  # Curve line
  geom_ribbon(aes(ymin = PanGenome - SD_PanGenome, ymax = PanGenome + SD_PanGenome),
              fill = "grey", alpha = 0.5) +  # Confidence interval ribbon
  labs(
    title = "",
    x = "Number of samples",
    y = "Number of genes",
    caption = ""
  ) +
  ylim(6300, 6600) +  # Set Y-axis limits
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.8, color = "black")
  )
