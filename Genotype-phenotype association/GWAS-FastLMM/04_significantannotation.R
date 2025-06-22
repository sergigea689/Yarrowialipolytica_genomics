# Script for extracting and annotating genes associated with FastLMM results

# Load required libraries
library(dplyr)
library(vcfR)
library(openxlsx)
library(readr)
library(stringr)
library(tidyr)
library(httr)
library(jsonlite)
library(readxl)

# Set working directory
setwd("C:/Users/USUARIO/Desktop/Acetate/Fast-LMM")

# ---------------------
# Load FastLMM results
# ---------------------
results <- read.csv("results_df.csv")  # Use read.csv2() if the file uses semicolons as separators

# Count total number of variants
total_variants <- nrow(results)

# Filter significant variants based on p-value threshold
significant_results <- results %>%
  filter(PValue < 0.005)

# ---------------------
# Load VCF and extract annotation fields
# ---------------------
vcf <- read.vcfR("phenosamples1.ann.vcf")
vcf_data <- vcf@fix
annots <- vcf_data[, c("CHROM", "POS", "INFO")]

# Convert to data frame and extract annotations from INFO column
annots_df <- as.data.frame(annots) %>%
  mutate(ANN = sub(".*ANN=([^;]+);?.*", "\\1", INFO)) %>%
  separate_rows(ANN, sep = ",") %>%
  separate(
    ANN,
    into = c("Allele", "Effect", "Impact", "Gene_Name", "Gene_ID", "Feature_Type", 
             "Feature_ID", "Transcript_Biotype", "Rank", "HGVS.c", "HGVS.p", 
             "cDNA.pos", "CDS.pos", "AA.pos", "Distance", "ERRORS"),
    sep = "\\|", fill = "right"
  )

# Rename chromosomes to match the format in the VCF
significant_results <- significant_results %>%
  mutate(Chr = recode(Chr, 
                      `1` = "YALI0A", 
                      `2` = "YALI0B", 
                      `3` = "YALI0C", 
                      `4` = "YALI0D", 
                      `5` = "YALI0E", 
                      `6` = "YALI0F"))

# Join significant results with VCF annotations
significant_annots <- annots_df %>%
  mutate(POS = as.numeric(POS)) %>%
  semi_join(significant_results, by = c("CHROM" = "Chr", "POS" = "ChrPos"))

# ---------------------
# Gene extraction
# ---------------------
genes_significativos <- significant_annots %>%
  select(CHROM, POS, Gene_Name, Effect, Impact) %>%
  distinct()

# Filter cases where genes are affected by both coding and upstream/intergenic variants
# Step 1: Identify genes with at least one missense variant
genes_con_missense <- genes_significativos %>%
  filter(Effect == "missense_variant") %>%
  pull(Gene_Name) %>%
  unique()

# Step 2: Keep only missense variants for those genes
#         For genes without missense variants, keep all entries
genes_filtrados <- genes_significativos %>%
  filter((Gene_Name %in% genes_con_missense & Effect == "missense_variant") |
           !(Gene_Name %in% genes_con_missense))

# ---------------------
# Annotate final table with extra variant information
# ---------------------
genes_anotados <- genes_filtrados %>%
  left_join(
    significant_results %>%
      select(Chr, ChrPos, 6:12),
    by = c("CHROM" = "Chr", "POS" = "ChrPos")
  )

# Save final annotated variants
write.xlsx(genes_anotados, file = "acetate_genes_anotados.xlsx")

# ---------------------
# Extract unique gene names
# ---------------------
genes_lista <- genes_anotados %>%
  select(Gene_Name) %>%
  separate_rows(Gene_Name, sep = "-") %>%                 # For intergenic SNPs
  mutate(Gene_Name = str_replace_all(Gene_Name, "_", "")) %>%  # Remove underscores
  distinct()

# Save as plain text file for downstream use
write_tsv(genes_lista, "acetate_gene_names.txt", col_names = FALSE, quote_escape = "none")

# ---------------------
# Query UniProt API for gene annotations
# ---------------------
gene_ids <- genes_lista$Gene_Name  # Gene identifiers
uniprot_ids <- list()              # Dictionary to store results

for (gene_id in gene_ids) {
  # Skip invalid or ambiguous gene names
  if (gene_id %in% c("gene", "CHREND")) {
    next
  }
  
  # Query UniProt
  url <- paste0("https://rest.uniprot.org/uniprotkb/stream?query=", gene_id, "&format=json")
  response <- GET(url)
  data <- fromJSON(content(response, "text"))
  
  # Extract UniProt ID if available
  if (length(data$results) > 0) {
    uniprot_ids[[gene_id]] <- data$results$primaryAccession
  } else {
    uniprot_ids[[gene_id]] <- NA
  }
}

# Print the retrieved UniProt IDs
print(uniprot_ids)

# Save backup
save(uniprot_ids, file = "uniprots_ids_backup.RData")

# Flatten list and filter non-missing values
uniprot_no_na <- unlist(uniprot_ids)[!is.na(unlist(uniprot_ids))]
uniprot_df <- data.frame(Uniprot_IDs = uniprot_no_na)

# Save as TSV file for enrichment analysis (e.g., clusterProfiler)
write_tsv(uniprot_df, "uniprot_ids_for_clusterProfiler.txt", col_names = FALSE, quote_escape = "none")

# ---------------------
# Count gene name repetitions
# ---------------------
rep <- read_excel("acetate_genes_anotados.xlsx")
repetitions <- as.data.frame(table(rep$Gene_Name))

# Save frequency of genes across variants
write.xlsx(repetitions, file = "genes_repetitions.xlsx")