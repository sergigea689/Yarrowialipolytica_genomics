#Script for the exploration of haplotypes per population in genes with high mean Fst value

# ------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------

library(readxl)
library(dplyr)
library(tidyr)
library(pheatmap)
library(vcfR)
library(RColorBrewer)

# ------------------------------------------------------------
# 2. Data loading
# ------------------------------------------------------------

# Load top genes identified by Fst analysis
top_genes <- read_excel("top_genes_genome_ann.xlsx")

# Load SNP presence/absence matrix
snp_matrix <- read_xlsx("presenceabscencesnpmatrix.xlsx")

# Load VCF annotations
vcf <- read.vcfR("phenosamples1.ann.vcf")
vcf_data <- vcf@fix

# Extract relevant annotation columns
annots <- vcf_data[, c("CHROM", "POS", "INFO")]
annots_df <- as.data.frame(annots)

# ------------------------------------------------------------
# 3. SNP parsing - annotations
# ------------------------------------------------------------

# Extract ANN field from INFO and split into multiple annotation columns
annots_df <- annots_df %>%
  mutate(ANN = sub(".*ANN=([^;]+);?.*", "\\1", INFO)) %>%
  separate_rows(ANN, sep = ",") %>%
  separate(
    ANN,
    into = c("Allele", "Effect", "Impact", "Gene_Name", "Gene_ID", "Feature_Type", 
             "Feature_ID", "Transcript_Biotype", "Rank", "HGVS.c", "HGVS.p", 
             "cDNA.pos", "CDS.pos", "AA.pos", "Distance", "ERRORS"),
    sep = "\\|",
    fill = "right"
  )

annots_df$POS <- as.numeric(annots_df$POS)

# Clean gene names by removing underscores
annots_df$Gene_Name <- gsub("_", "", annots_df$Gene_Name)

# Filter annotations for genes in the top genes list
gene_pattern <- paste0(top_genes$Name, collapse = "|")
annots_subset <- annots_df[grepl(gene_pattern, annots_df$Gene_Name), ]


# Rename first column to "SNP" if necessary
colnames(snp_matrix)[1] <- "SNP"

# Convert SNP matrix from wide to long format
snp_matrix_long <- snp_matrix %>%
  pivot_longer(-SNP, names_to = "Strain", values_to = "Presence")

# Create unique SNP identifier by combining chromosome and position
annots_subset_1_unique <- annots_subset %>%
  mutate(SNP = paste0(CHROM, "_", POS)) %>%
  select(SNP, Gene_Name) %>%
  distinct(SNP, .keep_all = TRUE)

# Filter SNP matrix to keep only SNPs present in annotation subset
snp_matrix_filtered <- snp_matrix_long %>%
  filter(SNP %in% annots_subset$SNP)

# Join SNP data with gene annotation
snp_annotated_1 <- snp_matrix_filtered %>%
  left_join(annots_subset_1_unique, by = "SNP")

# ------------------------------------------------------------
# 4. Constructing Haplotype matrix
# ------------------------------------------------------------

# Collapse SNP presence/absence by gene and strain into haplotype strings
haplotype_matrix <- snp_annotated_1 %>%
  arrange(Gene_Name, SNP) %>%
  group_by(Gene_Name, Strain) %>%
  summarise(Haplotype = paste0(Presence, collapse = ""), .groups = "drop")

# Convert from long to wide format: genes as rows, strains as columns
haplotype_wide <- haplotype_matrix %>%
  pivot_wider(names_from = Strain, values_from = Haplotype)

# Convert haplotype data frame to have gene names as rownames
haplotype_wide_df <- as.data.frame(haplotype_wide)
rownames(haplotype_wide_df) <- haplotype_wide_df[[1]]
haplotype_wide_df <- haplotype_wide_df[, -1]

# Count unique haplotypes per gene
haplotype_counts <- apply(haplotype_wide_df, 1, function(x) length(unique(x)))

# Create summary dataframe of haplotype diversity per gene
haplotype_summary <- data.frame(
  Gene = rownames(haplotype_wide_df),
  Num_Haplotypes = haplotype_counts
)

# Sort by number of haplotypes in descending order (optional)
haplotype_summary <- haplotype_summary[order(-haplotype_summary$Num_Haplotypes), ]

#Recode Haplotypes as Numeric Factors --------------------------------------

# Convert haplotype strings into numeric factor IDs for each gene and strain
haplotype_ids <- t(apply(haplotype_wide_df, 1, function(x) as.numeric(factor(x))))

haplotype_id_df <- as.data.frame(haplotype_ids)
rownames(haplotype_id_df) <- rownames(haplotype_wide_df)
colnames(haplotype_id_df) <- colnames(haplotype_wide_df)

# ------------------------------------------------------------
# 5. Load metadata
# ------------------------------------------------------------

#Load Strain Metadata and Prepare for Heatmap ------------------------------

strains <- read_excel("Acetate.xlsx")
strains_ordered <- strains %>% arrange(growth)

# Reorder SNP matrix columns to match ordered strains
snp_matrix_subset <- snp_matrix_subset[, strains_ordered$treesamp, drop = FALSE]
snp_matrix_subset[] <- lapply(snp_matrix_subset, as.numeric)

# Set rownames for SNP matrix subset
rownames(snp_matrix_subset) <- rowsub

# Create annotation dataframe for heatmap columns (strains)
annotation_df <- data.frame(group = strains$group)
rownames(annotation_df) <- strains$treesamp

# Define color palette for strain groups
col_acetate <- c("high" = "#F08705", "low" = "#6CD4EF")

# Define heatmap color palette and breakpoints
heatmap_palette <- rev(brewer.pal(n = 2, name = "RdYlBu"))
breaks <- seq(0, 1, length.out = 3)

# Plot Heatmap of Functional Haplotypes -------------------------------------

pheatmap(
  haplotype_id_df,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  display_numbers = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  breaks = breaks,
  main = "Functional Haplotypes by Gene",
  color = heatmap_palette,
  annotation_col = annotation_df,
  annotation_colors = list(group = col_acetate),
  fontsize_row = 6,
  fontsize_col = 7,
  angle_col = 90,
  annotation_col_margin = 25,
  cellwidth = 5
)

# ------------------------------------------------------------
# 6. Identification of common haplotypes
# ------------------------------------------------------------

# Function to assign top N common haplotypes, label others as "rare"
assign_top_haplotypes <- function(haps, top_n = 3) {
  freq <- sort(table(haps), decreasing = TRUE)
  top_haps <- names(freq)[1:min(top_n, length(freq))]
  ifelse(haps %in% top_haps, haps, "rare")
}

# Categorize haplotypes across all genes and strains
haplotype_categorized <- as.data.frame(
  t(apply(haplotype_wide_df, 1, assign_top_haplotypes))
)
rownames(haplotype_categorized) <- rownames(haplotype_wide_df)
colnames(haplotype_categorized) <- colnames(haplotype_wide_df)

# Identify strains labeled "high" in metadata
high_strains <- strains %>%
  filter(group == "high" & !is.na(treesamp)) %>%
  pull(treesamp)

# Recode haplotypes with respect to "high" strains:
# Place the most frequent haplotype in "high" strains first,
# followed by others, and assign "rare" last.
haplotype_ids <- t(apply(haplotype_categorized, 1, function(x) {
  x_high <- x[high_strains]
  x_high_no_rare <- x_high[x_high != "rare"]
  
  if (length(x_high_no_rare) == 0) {
    unique_levels <- unique(x)
    levels_ordered <- c(setdiff(unique_levels, "rare"), "rare")
  } else {
    freq_table <- sort(table(x_high_no_rare), decreasing = TRUE)
    majority_haplo <- names(freq_table)[1]
    others <- setdiff(unique(x), c(majority_haplo, "rare"))
    levels_ordered <- c(majority_haplo, others, "rare")
  }
  
  factorized <- factor(x, levels = levels_ordered)
  nums <- as.numeric(factorized)
  nums[x == "rare"] <- 4  # Explicitly assign numeric 4 to "rare"
  return(nums)
}))

haplotype_id_df <- as.data.frame(haplotype_ids)
rownames(haplotype_id_df) <- rownames(haplotype_categorized)
colnames(haplotype_id_df) <- colnames(haplotype_categorized)

# ------------------------------------------------------------
# 7. PCA Analysis of Functional Haplotype Profiles
# ------------------------------------------------------------

# Load required libraries for data manipulation
library(tidyr)
library(dplyr)
library(tibble)

# Transform the haplotype matrix: genes as rows and strains as columns
haplo_long <- as.data.frame(haplotype_id_df) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Strain", values_to = "Haplotype") %>%
  mutate(Haplotype = paste0("haplo_", Haplotype)) %>%
  unite("Feature", Gene, Haplotype, sep = "_") %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Feature, values_from = value, values_fill = 0)

# Prepare matrix for PCA: strains as rows, features as columns
pca_matrix <- haplo_long %>%
  column_to_rownames("Strain")

# Perform PCA with scaling
pca_res <- prcomp(pca_matrix, scale. = TRUE)

# Extract PCA coordinates for each strain
pca_coords <- as.data.frame(pca_res$x)

# Prepare metadata for merging
rownames(strains_ordenado) <- strains_ordenado$treesamp
pca_coords <- pca_coords %>%
  rownames_to_column(var = "treesamp")

# Join PCA coordinates with strain metadata
pca_coords <- left_join(pca_coords, strains_ordenado, by = "treesamp")

# Load visualization libraries
library(ggplot2)
library(grid)    # For arrow annotation (if needed)
library(tibble)
library(dplyr)

# Extract and scale PCA loadings for the first two principal components
loading_df <- as.data.frame(pca_res$rotation[, 1:2]) %>%
  rownames_to_column(var = "Feature") %>%
  mutate(magnitude = sqrt(PC1^2 + PC2^2)) %>%
  arrange(desc(magnitude)) %>%
  slice(1:10) %>%
  mutate(PC1 = PC1 * 100,   # Scale loadings for visualization purposes
         PC2 = PC2 * 100)

# Calculate variance explained by PC1 and PC2 (percentage)
var_exp <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100
pc1_var <- round(var_exp[1], 1)
pc2_var <- round(var_exp[2], 1)

# Format categorical variables for plotting
pca_coords$pop <- as.factor(pca_coords$pop)
pca_coords$group <- as.character(pca_coords$group)    # Convert to character to handle NA
pca_coords$group[is.na(pca_coords$group)] <- "intermediate"
pca_coords$group <- factor(pca_coords$group)

# Define color palettes and shapes for populations and groups
my_colors_pop <- c(
  "Eur1" = "#F29E6D", "Eur2" = "#05A6A6", "So-HC" = "#03738C",
  "NoAm" = "#F23E2E", "W29 clade" = "#152B59", "Hyp" = "grey"
)
my_shapes_group <- c(
  "high" = 24, "low" = 22, "intermediate" = 21
)

# First PCA scatterplot: colored by population and shaped by group
p <- ggplot(pca_coords, aes(x = PC1, y = PC2, fill = pop, shape = group)) +
  geom_point(size = 3, color = "black", alpha = 0.9) +  # Points with black borders
  scale_fill_manual(values = my_colors_pop) +
  scale_shape_manual(values = my_shapes_group) +
  theme_classic(base_size = 18) +  # Publication-ready font size
  labs(
    title = "PCA of functional haplotype profiles",
    x = paste0("PC1 (", pc1_var, "% variance explained)"),
    y = paste0("PC2 (", pc2_var, "% variance explained)"),
    fill = "Population",
    shape = "Group"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.title = element_text(face = "bold")
  )

print(p)


