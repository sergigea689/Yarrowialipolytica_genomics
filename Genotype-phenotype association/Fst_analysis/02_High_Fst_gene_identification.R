# Script for the analysis of global Fst values

# Set working directory
setwd("C:/Users/USUARIO/Desktop/Acetate/fst-eur2")

# Load required libraries
library(readxl)     # For reading Excel files
library(dplyr)      # For data manipulation
library(vcfR)       # For handling VCF data
library(openxlsx)   # For exporting to Excel
library(readr)      # For reading delimited files
library(stringr)    # For string operations
library(tidyr)      # For reshaping data

# ------------------------------------------------------------
# 1. Annotation data
# ------------------------------------------------------------


#We want to annotate the genes as well, we use anno_yali0.tsv.
anno <- read_tsv("anno_yali0.tsv")

tsv_data <- read.table("anno_yali0.tsv", header = TRUE, sep = "\t", quote = "", fill = TRUE, comment.char = "")
length(unique(tsv_data$Name))
#Fix tsv in R

# Create columns if they don't exist
if (!"Contig" %in% colnames(tsv_data)) tsv_data$Contig <- NA
if (!"Start" %in% colnames(tsv_data)) tsv_data$Start <- NA
if (!"End" %in% colnames(tsv_data)) tsv_data$End <- NA
if (!"Name" %in% colnames(tsv_data)) tsv_data$Name <- NA
if (!"Note" %in% colnames(tsv_data)) tsv_data$Note <- NA


# Modificy rows where start is an empty chain
rows_to_modify <- tsv_data$Start == "" | is.na(tsv_data$Start)

#Divide V1 column by spaces for selected columns, keeping only first 6 splits
split_data <- strsplit(tsv_data$Nº.gen.chr[rows_to_modify], " +", perl = TRUE)

# Temporal lists to store separated columns
Contig <- sapply(split_data, function(x) if (length(x) >= 2) x[2] else NA)
Start <- sapply(split_data, function(x) if (length(x) >= 3) as.integer(x[3]) else NA)
End <- sapply(split_data, function(x) if (length(x) >= 4) as.integer(x[4]) else NA)
Kind <- sapply(split_data, function(x) if (length(x) >= 5) x[5] else NA)
Name <- sapply(split_data, function(x) if (length(x) >= 6) x[6] else NA)

#Join remaining elements in a string for Note column
Note <- sapply(split_data, function(x) if (length(x) > 6) paste(x[7:length(x)], collapse = " ") else NA)

#Replace first row number by the original row number
OriginalRowNumber <- sapply(split_data, function(x) if (length(x) >= 1) as.integer(x[1]) else NA)

# Change selected rows in the original dataframe
tsv_data$Contig[rows_to_modify] <- Contig
tsv_data$Start[rows_to_modify] <- Start
tsv_data$End[rows_to_modify] <- End
tsv_data$Kind[rows_to_modify] <- Kind
tsv_data$Name[rows_to_modify] <- Name
tsv_data$Note[rows_to_modify] <- Note

tsv_data$Nº.gen.chr[rows_to_modify] <- OriginalRowNumber

#Just take CDS
tsv_data <- subset(tsv_data, tsv_data$Kind == "CDS")


tsv_data$Start <- as.numeric(tsv_data$Start)
tsv_data$End <- as.numeric(tsv_data$End)
any_na_in_start <- any(is.na(tsv_data$Start))
any_na_in_end <- any(is.na(tsv_data$End))

#Replace missing information in END
# Replace NAs in End column by value Start + 1000
tsv_data$End <- ifelse(is.na(tsv_data$End), tsv_data$Start + 1000, tsv_data$End)
tsv_data$Start <- as.numeric(tsv_data$Start)
tsv_data$End <- as.numeric(tsv_data$End)

anno <- tsv_data

# Filter genes on the list
anno_filtered <- anno

#Create columns for positions including 2000 bases upstream and downstream
anno_filtered$reg_Start <- anno_filtered$Start - 2000
anno_filtered$reg_End <- anno_filtered$End + 2000

# Select columns
result <- anno_filtered[, c("Contig", "Start", "End", "reg_Start", "reg_End", "Name", "Note")]

result$reg_Start <- pmax(0, result$reg_Start)
# Crear un dataframe BED con las columnas requeridas (sin encabezado)
bed_df <- result[, c("Contig", "reg_Start", "reg_End")]

# ------------------------------------------------------------
# 2. Load genome-wide Fst values and calculate percentiles
# ------------------------------------------------------------

# Read genome-wide Fst values (SNP-level)
fst_global <- read.table("fst_results.weir.fst", header = TRUE)

# Remove SNPs with missing Fst values
fst_global <- fst_global[!is.na(fst_global$WEIR_AND_COCKERHAM_FST), ]

#Load data
# Cargar FST de los SNPs en tus genes de interés

## Add a column for the mean of Fst
result$FST_mean <- NA
result$Num_SNPs <- 0
#Iterate for each gene
for (i in 1:nrow(result)) {
  chrom <- result$Contig[i]
  start <- result$reg_Start[i]
  end <- result$reg_End[i]
  
  # Subset SNPs dentro del rango del gen
  snps_in_range <- subset(fst_global, CHROM == chrom & POS >= start & POS <= end)
  
  # Calcular mediana y contar SNPs
  if (nrow(snps_in_range) > 0) {
    result$FST_mean[i] <- mean(snps_in_range$WEIR_AND_COCKERHAM_FST, na.rm = TRUE)
    result$Num_SNPs[i] <- nrow(snps_in_range)
  }
}


# Define a function to calculate the percentile of a given Fst value
get_percentile <- function(fst_value, global_values) {
  ecdf(global_values)(fst_value)
}

# Add a column with percentile values to the result table
# NOTE: make sure 'result' and its column FST_mean are defined before this step
result$FST_percentile <- sapply(result$FST_mean, get_percentile, global_values = fst_global$WEIR_AND_COCKERHAM_FST)

# ------------------------------------------------------------
# 3. Compute mean Fst value per gene
# ------------------------------------------------------------

# Read gene-level Fst values (each SNP assigned to a gene region)
result_genome <- result
fst_genes_global <- read.table("fst_results.weir.fst", header = TRUE)
head(fst_genes_global)

# Initialize columns for Fst mean and number of SNPs per gene
result_genome$FST_mean <- NA
result_genome$Num_SNPs <- 0

# Loop over each gene to calculate mean Fst and SNP count
for (i in 1:nrow(result_genome)) {
  chrom <- result_genome$Contig[i]
  start <- result_genome$reg_Start[i]
  end <- result_genome$reg_End[i]
  
  # Subset SNPs within the gene region
  snps_in_range <- subset(fst_genes_global, CHROM == chrom & POS >= start & POS <= end)
  
  # Compute mean Fst and SNP count if SNPs are present
  if (nrow(snps_in_range) > 0) {
    result_genome$FST_mean[i] <- mean(snps_in_range$WEIR_AND_COCKERHAM_FST, na.rm = TRUE)
    result_genome$Num_SNPs[i] <- nrow(snps_in_range)
  }
}

# ------------------------------------------------------------
# 4. Identify top 2% of genes based on Fst values
# ------------------------------------------------------------

# Define the 98th percentile threshold for gene-level Fst mean
threshold_genome <- quantile(result_genome$FST_mean, 0.98, na.rm = TRUE)

# Get the names of genes with Fst above the threshold
top_genes_genome <- result_genome$Name[result_genome$FST_mean >= threshold_genome]

# Save the list of top genes to a text file
write.table(top_genes_genome, file = "top2percentFstgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Extract annotations for the top genes
top_genes_genome_ann <- result_genome %>%
  filter(FST_mean >= threshold_genome) %>%
  select(Name, Note)

# Export annotations to Excel
write.xlsx(top_genes_genome_ann, "top_genes_genome_ann.xlsx")

# ------------------------------------------------------------
# 5. Extract SNPs located within the top 2% genes
# ------------------------------------------------------------

# Subset genes that are in the top 2%
subset_genome <- result_genome[result_genome$Name %in% top_genes_genome, ]

# Initialize empty dataframe to store SNPs
snps_top_genome <- data.frame()

# Loop over each top gene and collect SNPs in the corresponding region
for (i in 1:nrow(subset_genome)) {
  chr <- subset_genome$Contig[i]
  start <- subset_genome$reg_Start[i]
  end <- subset_genome$reg_End[i]
  gene_name <- subset_genome$Name[i]
  
  # Subset SNPs within the gene coordinates
  snps_in_gene <- subset(fst_genes_global, CHROM == chr & POS >= start & POS <= end)
  
  # Add gene name to each SNP
  snps_in_gene$Gene <- gene_name
  
  # Append SNPs to the output dataframe
  snps_top_genome <- rbind(snps_top_genome, snps_in_gene)
}

# Export SNP list with gene information
write.table(snps_top_genome, "snps_top2pct_genes_with_gene.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Also export the top gene list again (redundant with previous but included here)
write.table(top_genes_genome, "list_top2pct_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------
# 6. Annotate the SNPs located in the top genes using VCF annotation
# ------------------------------------------------------------

# Load the VCF file with annotations
vcf <- read.vcfR("phenosamples1.ann.vcf")

# Extract fixed fields (including CHROM, POS, INFO)
vcf_data <- vcf@fix

# Select relevant columns for annotation: chromosome, position, and INFO field
annots <- vcf_data[, c("CHROM", "POS", "INFO")]

# Convert to a data frame for easier manipulation
annots_df <- as.data.frame(annots)

# Extract annotation details from the INFO column
annots_df <- annots_df %>%
  # Extract the ANN field from INFO (everything between 'ANN=' and the next semicolon)
  mutate(ANN = sub(".*ANN=([^;]+);?.*", "\\1", INFO)) %>%
  # Separate multiple annotations per SNP (comma-separated)
  tidyr::separate_rows(ANN, sep = ",") %>%
  # Split the ANN field into its subfields based on the '|' delimiter
  tidyr::separate(ANN, into = c("Allele", "Effect", "Impact", "Gene_Name", "Gene_ID", 
                                "Feature_Type", "Feature_ID", "Transcript_Biotype", 
                                "Rank", "HGVS.c", "HGVS.p", "cDNA.pos", "CDS.pos", 
                                "AA.pos", "Distance", "ERRORS"),
                  sep = "\\|", fill = "right")

# Convert POS column to numeric type for proper joining
annots_df$POS <- as.numeric(annots_df$POS)

# Join SNP data with annotation information by chromosome and position
snps_top_genome_annot <- snps_top_genome %>%
  left_join(annots_df[, c("CHROM", "POS", "Effect", "Impact", "Allele", "Gene_Name")],
            by = c("CHROM", "POS"))

# Clean up gene names by removing underscores
annots_df$Gene_Name <- gsub("_", "", annots_df$Gene_Name)


# -----------------------------------------------
# 7. Obtain Uniprot codes for the genes
# -----------------------------------------------

# Load the list of top 2% Fst genes
top20pcfst <- read.table("So-HC_unique_gains.txt")

# Load necessary libraries for API requests and JSON parsing
library(httr)
library(jsonlite)

# Extract gene identifiers from the loaded file
gene_ids <- top20pcfst$V1  # Gene identifiers

# Initialize empty lists to store results
kegg_ids <- list()
uniprot_ids <- list()

# Query UniProt for each gene to get UniProt IDs
for (gene_id in gene_ids) {
  # Skip genes with names 'gene' or 'CHREND'
  if (gene_id %in% c("gene", "CHREND")) {
    next  # Skip to next iteration
  }
  
  # Construct UniProt API URL to search by gene_id
  url <- paste0("https://rest.uniprot.org/uniprotkb/stream?query=", gene_id, "&format=json")
  
  # Perform GET request
  response <- GET(url)
  
  # Parse JSON content from response
  data <- fromJSON(content(response, "text"))
  
  # Check if any UniProt entries were returned
  if (length(data$results) > 0) {
    uniprot_ids[[gene_id]] <- data$results$primaryAccession  # Store UniProt ID
  } else {
    uniprot_ids[[gene_id]] <- NA  # Assign NA if no UniProt ID found
  }
}

# Print the UniProt IDs found
print(uniprot_ids)

# Save UniProt IDs as a backup file
save(uniprot_ids, file = "SoHC_unique_gains_uniprots_ids_backup.RData")

# Extract just the UniProt IDs (values) from the list
uniprots_ids_only <- unlist(uniprot_ids)

# Filter out NA values and convert to a data frame
uniprot_no_na <- uniprots_ids_only[!is.na(uniprots_ids_only)]
uniprot_df <- data.frame(Uniprot_IDs = uniprot_no_na)

# Save the UniProt IDs to a tab-separated file for clusterProfiler input
write_tsv(uniprot_df, "uniprot_ids_for_clusterProfiler_SoHC_unique_gains.txt", col_names = FALSE, quote_escape = "none")

