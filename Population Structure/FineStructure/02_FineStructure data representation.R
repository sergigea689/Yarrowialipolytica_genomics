##################################################################
## FineSTRUCTURE R Workflow
## Author: Daniel Lawson (dan.lawson@bristol.ac.uk)
## Website: www.paintmychromosomes.com (see "R Library" section)
## Date: 14/02/2012
## 
## Description:
## Helper script to process FineSTRUCTURE output.
## NOTE: This is not a robust or generalized R package. It is provided 
## as-is for use with specific FineSTRUCTURE outputs. Use with caution.
##
## License: GPL v3
##################################################################

### Load Required Packages
library(XML)
library(ape)
library(pheatmap)
library(dendextend)
library(ggtree)
library(readxl)
library(openxlsx)

### Load Custom FineSTRUCTURE Library
source("FinestructureLibrary.R")

### Set Working Directory
setwd("C:/Users/sergi/Desktop/finestructure_whole")

### Define Input File Names
chunkfile <- "output.chunkcounts.out"        # ChromoPainter chunkcounts file
mcmcfile  <- "out.YALI.mcmc.xml"             # FineSTRUCTURE MCMC XML output
treefile  <- "structure_tree.out"            # FineSTRUCTURE Newick tree
coincidence_file <- "structure_meancoincidence.csv"  # Mean coincidence matrix

### File Existence Check
stopifnot(file.exists(chunkfile),
          file.exists(treefile),
          file.exists(mcmcfile),
          file.exists(coincidence_file))

### Load Chunkcount Matrix (pairwise coancestry)
dataraw <- as.matrix(read.table(chunkfile, row.names = 1, header = TRUE, skip = 1))

### Load MCMC Output
mcmcxml  <- xmlTreeParse(mcmcfile)
mcmcdata <- as.data.frame.myres(mcmcxml)

### Load and Visualize Mean Coincidence Matrix
mat <- as.matrix(read.csv(coincidence_file, row.names = 1))
pheatmap(mat, 
         main = "Mean Coancestry Matrix", 
         clustering_distance_rows = "euclidean", 
         show_colnames = FALSE)

### Load FineSTRUCTURE Tree
newick_string <- readLines(treefile)
newick_string <- newick_string[grep(";", newick_string)][1]
tree <- read.tree(text = newick_string)

### Basic Tree Plot
plot(tree, main = "FineSTRUCTURE Tree", cex = 0.5, no.margin = TRUE)

### Load Metadata and Assign Tip Colors
metadata <- read_excel("update_manfinetreeposition.xlsx")
subset <- subset(metadata, metadata$treesamp %in% tree$tip.label)

my_colors_pop <- c("Eur1" = "#F29E6D", "Eur2" = "#05A6A6", "So-HC" = "#03738C",
                   "NoAm" = "#F23E2E", "W29 clade" = "#152B59", "Hyp" = "grey")

# Match metadata order with tree tips
tip_labels <- tree$tip.label
tip_colors <- my_colors_pop[subset$FinePop[match(tip_labels, subset$treesamp)]]
subset_for_plot <- subset[match(tip_labels, subset$treesamp), ]
stopifnot(all(tree$tip.label == subset_for_plot$treesamp))

### Visualize Tree with Metadata (ggtree)
p <- ggtree(tree) %<+% subset_for_plot +
  geom_tippoint(aes(color = FinePop2), size = 1.5) +
  geom_tiplab(aes(label = Name), size = 2, hjust = -0.1) +
  theme_tree2() +
  scale_color_manual(values = my_colors_pop) +
  theme(legend.position = "right")
print(p)


### Prepare and Visualize Coancestry Matrix
replace_dot_with_dash <- function(col_names) gsub("\\.", "-", col_names)
colnames(dataraw) <- replace_dot_with_dash(colnames(dataraw))

# Load tree order
fullorder <- scan("paperlabels.txt", what = character(), sep = "\t")
datamatrix <- dataraw[fullorder, fullorder]

# Cap matrix values to improve contrast in heatmap
tmatmax <- 500
tmpmat <- datamatrix
tmpmat[tmpmat > tmatmax] <- tmatmax

# Create annotation for rows/columns
subset_for_plot_ordered <- subset_for_plot[match(fullorder, subset_for_plot$treesamp), ]
stopifnot(all(subset_for_plot_ordered$treesamp == fullorder))
annotation <- data.frame(FinePop2 = subset_for_plot_ordered$FinePop2)
rownames(annotation) <- subset_for_plot_ordered$treesamp

# Plot coancestry heatmap
pheatmap(tmpmat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         breaks = seq(200, 300, length.out = 20),
         annotation_row = annotation,
         annotation_col = annotation,
         annotation_colors = list(FinePop2 = my_colors_pop),
         show_rownames = FALSE,
         show_colnames = FALSE,
         color = colorRampPalette(c("yellow", "red"))(19))

### Compute and Visualize Average Coancestry Between Populations
pop_labels <- annotation$FinePop2
names(pop_labels) <- rownames(annotation)
unique_pops <- unique(pop_labels)

# Initialize population-level coancestry matrix
popmat <- matrix(NA, nrow = length(unique_pops), ncol = length(unique_pops),
                 dimnames = list(unique_pops, unique_pops))

# Calculate average pairwise coancestry between populations
for (i in unique_pops) {
  for (j in unique_pops) {
    rows_i <- names(pop_labels)[pop_labels == i]
    cols_j <- names(pop_labels)[pop_labels == j]
    submatrix <- tmpmat[rows_i, cols_j, drop = FALSE]
    popmat[i, j] <- mean(submatrix)
  }
}

# Visualize raw and symmetrized matrices
pheatmap(popmat,
         color = colorRampPalette(c("yellow", "red"))(5),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.1f",
         main = "Average Coancestry Between Populations")

mean_matrix_sym <- (popmat + t(popmat)) / 2  # Symmetrize the matrix

pheatmap(mean_matrix_sym,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("yellow", "red"))(7),
         display_numbers = TRUE,
         number_format = "%.1f",
         fontsize_number = 10,
         border_color = NA,
         angle_col = 45,
         main = "Symmetrized Average Coancestry Between Populations")
