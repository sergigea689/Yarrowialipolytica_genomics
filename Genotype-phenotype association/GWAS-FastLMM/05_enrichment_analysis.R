# Script for the enrichment analysis with clusterProfiler (KEGG) and TopGO (GO)

# ---------------------
#1. Directories, installation and libraries
# ---------------------

getwd()
setwd("C:/Users/USUARIO/Desktop/Acetate/Fast-LMM")

# Install packages
install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install(c("enrichplot", "DOSE", "pathview", "ggplot2"))
install.packages("biomaRt")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
install.packages("httr")
install.packages("jsonlite")

# Libraries
library(httr)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(pathview)
library(ggplot2)
library(biomaRt)
library(jsonlite)
library(stringr)

# ---------------------
#2. Transform uniprot codes into KEGG codes
# ---------------------



# Load data
uniprot_file <- "uniprot_ids_for_clusterProfiler.txt"
uniprot_ids <- readLines(uniprot_file)
kegg_ids <- list()

# Iterate over each uniprot_id to get KEGG ids from UniProt API
for (uniprot_id in uniprot_ids) {
  url <- paste0("https://www.uniprot.org/uniprot/", uniprot_id, ".json")
  response <- GET(url)
  json_text <- content(response, "text")
  kegg_id_match <- str_extract(json_text, '"KEGG","id":"(yli:[^"]+)"')
  if (!is.na(kegg_id_match)) {
    kegg_id <- str_remove(kegg_id_match, '"database": "KEGG", "id": "')
    kegg_id <- str_remove(kegg_id, '"')
  } else {
    kegg_id <- NA
  }
  kegg_ids[[uniprot_id]] <- kegg_id
}

print(kegg_ids)

# Clean KEGG ids
kegg_ids_clean <- na.omit(sapply(kegg_ids, function(x) {
  result <- str_extract(x, "(yli:[^\",]+)")
  if (!is.na(result)) {
    return(result)
  } else {
    return(NA)
  }
}))

kegg <- as.data.frame(kegg_ids_clean)
kegg_ids_clean_numeric <- gsub("yli:", "", kegg$kegg_ids_clean)

# ---------------------
#3. # KEGG Enrichment analysis with clusterProfiler
# ---------------------

enrichment_results <- enrichKEGG(gene = kegg_ids_clean_numeric, organism = "yli", pvalueCutoff = 0.05)
results <- enrichment_results@result
summary(enrichment_results)

# ---------------------
#4. # GO Enrichment analysis preparatoion
# ---------------------

# Load GO terms

# Install tidyverse if needed
install.packages("tidyverse")

library(clusterProfiler)
library(tidyverse)
library(tidyr)
library(readxl)
library(dplyr)

# Read GO database and GO to test
go_mapping <- read_excel("go_mapping.xlsx")
go_final <- read_excel("go_final.xlsx")

# Separate Go terms - one per row
term2gene_all <- go_mapping %>%
  separate_rows(GOs, sep = ",") %>%
  mutate(GOs = trimws(GOs)) %>%
  select(GOs, YALI_ID)

# Extract genes of interest
genes_interes <- go_final$YALI_ID

# Alternative: download full gene ontology file and load
goa <- read.delim("20011.Y_lipolytica.goa", comment.char = "!", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Extract important columns
term2gene <- goa[, c(5, 2)]  # UniProt ID - GO Term
colnames(term2gene) <- c("GO", "UniProt ID")
term2name <- unique(goa[, c(5, 10)])  # GO Term - GO Name
colnames(term2name) <- c("GO", "TermName")

# TopGO strategy
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("topGO")
library(topGO)

# Data prep
all_genes <- unique(term2gene$`UniProt ID`)
geneList <- rep(0.5, length(all_genes))
names(geneList) <- all_genes
geneList[names(geneList) %in% uniprot_ids] <- 0.01

topDiffGenes <- function(geneList) {
  return(geneList < 0.05)
}

gene2GO <- split(term2gene$GO, term2gene$`UniProt ID`)


# Table with top 40 terms according to Fisher test
allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "classicFisher", topNodes = 40)
print(allRes)

# Check genes with a specific GO term
genes_con_termino <- names(Filter(function(x) "GO:0042908" %in% x, gene2GO))
print(genes_con_termino)

#Other ways of checking
# genesGO <- genesInTerm(sampleGOdata, whichGO = "GO:0035825")[[1]]
# print(genesGO)
# 
# coincidencias <- intersect(genesGO, uniprot_ids)
# print(coincidencias)
# 
# unidf <- as.data.frame(uniprot_ids)
# g <- subset(unidf, unidf$uniprot_ids %in% genes_con_termino)

# ---------------------
#5. # GO enrichment analysis with TopGo
# ---------------------

run_topGO_weight01 <- function(ont, geneList, gene2GO) {
  sampleGOdata <- new("topGOdata",
                      description = paste("GO analysis -", ont),
                      ontology = ont,
                      allGenes = geneList,
                      geneSel = topDiffGenes,
                      nodeSize = 10,
                      annot = annFUN.gene2GO,
                      gene2GO = gene2GO)
  
  resultWeight01 <- runTest(sampleGOdata, algorithm = "weight01", statistic = "fisher")
  
  resTable <- GenTable(sampleGOdata,
                       weight01Fisher = resultWeight01,
                       orderBy = "weight01Fisher",
                       topNodes = 50)
  
  resTable$ontology <- ont
  resTable$pvalue <- as.numeric(sub("< ", "", resTable$weight01Fisher))
  return(resTable)
}

# Run the function for all ontologies
go_bp <- run_topGO_weight01("BP", geneList, gene2GO)
go_cc <- run_topGO_weight01("CC", geneList, gene2GO)
go_mf <- run_topGO_weight01("MF", geneList, gene2GO)

go_all <- dplyr::bind_rows(go_bp, go_cc, go_mf)

go_all$Term <- as.character(go_all$Term)
go_all$Term <- as.factor(go_all$Term)

str(go_all$Term)

go_all_df <- as.data.frame(go_all)

library(dplyr)

# ---------------------
#6. # Representation of results
# ---------------------

top_terms <- read_excel("Eur2unique_enrichment_cat_unadjusted.xlsx")
top_terms <- top_terms %>%
  arrange(ontology, pvalue) %>%
  mutate(Term = factor(Term, levels = rev(unique(Term))))

ontology_palette <- c(
  "BP" = "#FF8C00",  # naranja intenso (biological process)
  "MF" = "#00BFA0",  # verde-azulado vibrante (molecular function)
  "CC" = "#8A2BE2"   # púrpura brillante (cellular component)
)
library(ggplotw)
ggplot(top_terms, aes(x = -log10(pvalue), y = Term, 
                      size = Significant, color = ontology)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = ontology_palette) +
  scale_size_continuous(range = c(3, 10)) +
  labs(
    x = expression(-log[10](pvalue)),
    y = "GO Term",
    size = "Number of genes",
    color = "Ontology",
    title = "Top enriched GO terms (p < 0.05)"
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", color = "black")
  )

# Export results
library(openxlsx)
write.xlsx(top_terms, file = "enrichment_cat_unadjusted.xlsx")

# Examples for working with GO offspring terms
offspring <- GOBPOFFSPRING[["GO:0035825"]]
head(offspring)

any(offspring %in% gene2GO)
head(gene2GO)

genes_matching <- names(gene2GO)[sapply(gene2GO, function(go_list) any(go_list %in% offspring))]
print(genes_matching)
