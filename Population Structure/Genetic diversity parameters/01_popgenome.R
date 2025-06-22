# PopGenome Analysis Script
# Note: PopGenome must be installed from a local .tar.gz file as it is no longer maintained on CRAN.

# Set working directory to where the PopGenome tar.gz file is located
setwd("C:/Users/USUARIO/Downloads")

# Specify path to the PopGenome archive (adjust to your local path)
ruta_archivo_tar <- "C:/Users/USUARIO/Downloads/PopGenome_2.7.5.tar.gz"

# Install PopGenome from source
install.packages(ruta_archivo_tar, repos = NULL, type = "source")

# Install required dependencies
install.packages("vcfR")            # For working with VCF files
install.packages("ff", type = "binary")  # Dependency for PopGenome

# Load the libraries
library(PopGenome)
library(vcfR)

# Confirm working directory
getwd()

# Set working directory to location of the VCF file and input data
setwd("C:/Users/USUARIO/Desktop/tree/nomutants/admixture")

# Define genomic regions (start and end positions per chromosome/contig)
posfrom <- c(1, 1, 1, 1, 1, 1)
posto <- c(2303260, 3066374, 3272609, 3661369, 4224102, 4003361)
chr <- c("YALI0A", "YALI0B", "YALI0C", "YALI0D", "YALI0E", "YALI0F")

# Initialize list to store per-contig VCF data
vcf_data_list <- list()

# Loop through each contig and extract the corresponding region from the compressed and indexed VCF
# Note: Ensure the VCF file has been bgzipped and indexed using tabix prior to this step
for (i in 1:length(chr)) {
  vcf_data <- readVCF(filename = "C:/Users/USUARIO/Desktop/tree/nomutants/admixture/biallelicvariantspure.recode.vcf.gz",
                      numcols = 1000, tid = chr[i], frompos = posfrom[i], topos = posto[i], approx = TRUE)
  vcf_data_list[[i]] <- vcf_data
}

# If needed, combine the VCF objects (optional)
# all_vcf_data <- rbind.vcf(vcf_data_list)

## Genetic Analysis Section

# (Re-)Install required packages for analysis and plotting
install.packages("../../../../../../../Downloads/PopGenome", repos = NULL, type = "source")
install.packages("tidyverse")
install.packages("dplyr")

# Load analysis libraries
library(PopGenome)
library(tidyverse)
library(dplyr)
# Load custom plotting functions (if needed)
source("./plotting_funcs.R")


# Concatenate VCF objects into a single genome object
Whole_genome <- concatenate.classes(list(vcf_data_list[[1]], vcf_data_list[[2]], vcf_data_list[[3]],
                                         vcf_data_list[[4]], vcf_data_list[[5]], vcf_data_list[[6]]))

# Display summary statistics
get.sum.data(Whole_genome)
Whole_genome@populations  # Check if population data is present

# Load population assignments from a text file
# The file should have sample names in the first column and population labels in the second
pop <- read.table("C:/Users/USUARIO/Desktop/tree/nomutants/admixture/poppure.txt")

# Create population groupings
Eur1 <- subset(pop, pop$V2 == "Eur1")$V1
Eur2 <- subset(pop, pop$V2 == "Eur2")$V1
SoHC <- subset(pop, pop$V2 == "So-HC")$V1
NoAm <- subset(pop, pop$V2 == "NoAm")$V1

# Assign populations to the genome object
Whole_genome <- set.populations(Whole_genome, list(Eur1, Eur2, SoHC, NoAm))
# Optional: set outgroup if applicable
# Whole_genome <- set.outgroup(Whole_genome, outgroup)

# Concatenate all regions for genome-wide analysis
Whole_genome <- concatenate.regions(Whole_genome)

# Display object structure and SNP count
Whole_genome
Whole_genome@region.names
Whole_genome@n.biallelic.sites

# Extract individual sample names
sample_names <- get.individuals(Whole_genome)

# Compute FST between populations using nucleotide diversity
Whole_genome <- F_ST.stats(Whole_genome, mode = "nucleotide")
pairwise.FST <- t(Whole_genome@nuc.F_ST.pairwise)
head(pairwise.FST)

# Save FST matrix to file
write.table(pairwise.FST, file = "Fst_pure", sep = "\t", row.names = F, col.names = T)

# Compute nucleotide diversity (π)
Whole_genome <- diversity.stats(Whole_genome)
nucdiv <- Whole_genome@nuc.diversity.within / 20000000
head(nucdiv)

# Save π values to file
write.table(nucdiv, file = "Pi_pure", sep = "\t", row.names = F, col.names = T)

#### Visualization Section: Pi and FST

library(tidyverse)
library(ggplot2)

# Example data for visualization of pairwise FST
datos_fst <- data.frame(
  poblacion_origen = c("Eur1", "Eur1", "Eur1", "Eur2", "Eur2", "So-HC"),
  poblacion_destino = c("Eur2", "So-HC", "NoAm", "So-HC", "NoAm", "NoAm"),
  fst = c(0.76, 0.57, 0.61, 0.76, 0.83, 0.64)
)

# Reshape data for ggplot
library(reshape2)
datos_fst_melt <- melt(datos_fst, id.vars = c("poblacion_origen", "poblacion_destino"))

# Define factor levels for consistent plotting order
datos_fst_melt$poblacion_origen <- factor(datos_fst_melt$poblacion_origen, levels = unique(datos_fst_melt$poblacion_origen))
datos_fst_melt$poblacion_destino <- factor(datos_fst_melt$poblacion_destino, levels = unique(datos_fst_melt$poblacion_destino))

# Create heatmap of pairwise FST values
plot <- ggplot(datos_fst_melt, aes(x = poblacion_destino, y = poblacion_origen, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "orange", high = "red") +
  geom_text(aes(label = round(value, 2)), color = "black") +
  theme_minimal() +
  labs(x = "", y = "", title = "Fixation Index") +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 13, face = "bold")
  )
print(plot)

# Export FST heatmap to PDF
ggsave(
  filename = "paper_fixation_index_plot.pdf",
  plot = plot,
  device = "pdf",
  width = 8, height = 6,
  dpi = 300
)

# Create a bubble chart of nucleotide diversity (π)
datos_pi <- data.frame(
  Pop = c("Eur1", "Eur2", "So-HC", "NoAm"),
  NucDiv = c("0.001", "0.00004", "0.0008", "0.0006")
)
datos_pi$NucDiv <- as.numeric(datos_pi$NucDiv)
datos_pi$log10 <- log10(datos_pi$NucDiv)

# Define custom colors for populations
colfinepop <- c("So-HC" = "#03738C", "Eur2" = "#05A6A6",
                "Eur1" = "#F29E6D", "NoAm" = "#F23E2E")


# Reorder populations for bar chart
datos_pi$Pop <- factor(datos_pi$Pop, levels = datos_pi$Pop)
order <- c("NoAm", "So-HC", "Eur1", "Eur2")
datos_pi <- datos_pi[order(match(datos_pi$Pop, order)),]
datos_pi$Pop <- factor(datos_pi$Pop, levels = unique(datos_pi$Pop))

# Generate bar chart of nucleotide diversity
bar_chart <- ggplot(datos_pi, aes(x = Pop, y = NucDiv, fill = Pop)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = colfinepop) +
  labs(title = "Nucleotide diversity", x = "", y = "") +
  scale_y_continuous(limits = c(0.0000, 0.0015), breaks = seq(0, 0.002, by = 0.0005)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 0.1, colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold"))

