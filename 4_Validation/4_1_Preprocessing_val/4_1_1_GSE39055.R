
#############################################################################
#> Script to perform preprocessing steps on GSE39055 
#> 
#> Inputs: gene_signature
#> 
#> Outputs: metadata_gse39055, gene_expression_matrix
#> 
#############################################################################

library(GEOquery)
library(tidyverse)
library(Biobase)
library(illuminaHumanWGDASLv4.db)


# 1.--------------- Load processed expression data -----------------

#NOTE: No normalization step was needed because GSE39055 had already been
# normalized using VST and quantile normalization

gse39055 <- getGEO("GSE39055", GSEMatrix = TRUE, getGPL = FALSE)

length(gse39055)

gse39055_data <- gse39055[[1]]

expr_matrix <- exprs(gse39055_data) #Extract expression matrix

dim(expr_matrix)


# Map Illumina probe ID's to Gene Symbols 
probe_gene <- AnnotationDbi::select(
  illuminaHumanWGDASLv4.db,
  keys = rownames(expr_matrix),
  columns = "SYMBOL",
  keytype = "PROBEID"
)

# Clean probe annotation
# Remove probes without a Gene Symbol and remove identical probe-to-gene mappings
probe_gene_clean <- probe_gene %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  distinct(PROBEID, SYMBOL)

# Identify ambiguous probes that map to more than one Gene Symbol
ambiguous_probes <- probe_gene_clean %>%
  group_by(PROBEID) %>%
  summarise(n_symbols = n_distinct(SYMBOL), .groups = "drop") %>%
  filter(n_symbols > 1)

# Remove ambiguous probes and keep only probes that map to one Gene Symbol
probe_gene_clean <- probe_gene_clean %>%
  group_by(PROBEID) %>%
  filter(n_distinct(SYMBOL) == 1) %>%
  dplyr::slice(1) %>%
  ungroup()


# Probe ID to Symbol 
counts_data_gse39055 <- expr_matrix %>%
  as.data.frame() %>%
  rownames_to_column("PROBEID") %>% # Move probe ID to column
  inner_join(probe_gene_clean, by = "PROBEID") %>% # Add Gene Symbolss
  mutate(variance = apply(dplyr::select(., -PROBEID, -SYMBOL), 1, var)) %>%  # Calculate variance of each probe across all sample
  group_by(SYMBOL) %>% slice_max(order_by = variance,
                                 n = 1,
                                 with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(-PROBEID, -variance) %>% # Remove unwanted columns
  column_to_rownames("SYMBOL") %>% # SYMBOL as rowname
  as.matrix()

dim(counts_data_gse39055)
head(rownames(counts_data_gse39055))
anyDuplicated(rownames(counts_data_gse39055))
sum(is.na(counts_data_gse39055))




# 2.- Load Metadata GSE39055 ----------------------------------------------

pre_metadata <- pData(gse39055_data)

metadata_gse39055_pre <- pre_metadata

dim(metadata_gse39055_pre)
colnames(metadata_gse39055_pre)


# Verify that metadata is identical to the expression matrix
identical(colnames(expr_matrix), rownames(metadata_gse39055_pre))


# 3.- Preprocessing metadata GSE39055 ------------------------------------

metadata_gse39055 <- metadata_gse39055_pre %>%
  mutate(
    REC_STAT = case_when(`recurrence:ch1` == "Y" ~ 1, `recurrence:ch1` == "N" ~ 0),
    # 1: recurrence, 0: no recurrence/last follow-up
    REC_MON = as.numeric(`time until first recurrence or latest follow-up (months):ch1`),
    SURV_STAT = case_when(`death:ch1` == "Y" ~ 1, `death:ch1` == "N" ~ 0),
  )

# Add GEO sample accession as a simplified sample identifier
metadata_gse39055 <- metadata_gse39055 %>%
  mutate(id = geo_accession)

# Verify that expression data and metadata contain the same samples in the same order
identical(colnames(counts_data_gse39055), metadata_gse39055$id)


# Match common genes from TARGET-OS signature with GSE39055

gene_signature <- scan("output_data/gene_signature.csv", sep = ",", what = character())

# Confirm all genes in signature are found in new gse

common_genes_TARGET_GSE39055 <- intersect(gene_signature, rownames(counts_data_gse39055))

length(common_genes_TARGET_GSE39055)

saveRDS(metadata_gse39055, "output_data/metadata_gse39055.RDS")
saveRDS(counts_data_gse39055, "output_data/counts_data_gse39055.RDS")

rm(list = ls())
