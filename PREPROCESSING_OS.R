library(TCGAbiolinks)
library(data.table)
library(dplyr)
library(biomaRt)
library(UCSCXenaTools)


os_directory <- "~/Documents/OSTEOSARCOMA/R.project/Hueso"


#1 - DOWNLOAD COUNTS DATA

query_os <- GDCquery(
  project = "TARGET-OS",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts")

# GDCdownload(query_os, directory = os_directory)


os_data <- GDCprepare(
  query = query_os,
  directory = os_directory,
  summarizedExperiment = FALSE)

# Raw counts
counts_col <- grep(
  pattern = "^fpkm_uq_unstranded_",
  x = names(os_data),
  value = TRUE)


length(counts_col)
head(counts_col)
tail(counts_col)

counts_raw <- os_data %>% 
  dplyr::select(all_of(counts_col),
         gene_name)


counts_raw$variance <- apply(counts_raw %>% dplyr::select(-gene_name),1 , var, na.rm = TRUE)

counts_raw <- counts_raw %>% 
  group_by(gene_name) %>%
  slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% 
  ungroup %>% 
  tibble::column_to_rownames("gene_name") %>% 
  dplyr::select(-variance)
  


class(counts_raw)
dim(counts_raw)


# Use ENSEML ID as identifier of each raw

colnames(counts_raw) <- sub(
  pattern = "^fpkm_uq_unstranded_", # remove "unstranded"
  replacement = "",
  x = colnames(counts_raw))

head(colnames(counts_raw))

# Check for duplicates

# Ensembl IDs

sum(duplicated(rownames(counts_raw)))

# Duplicate samples

sum(duplicated(colnames(counts_raw)))

# Symbol

sum(duplicated(gene_annotation$gene_name))

# Search for Na in counts

anyNA(counts_raw)

# Negative values

any(counts_raw < 0)


mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl")


genes <- biomaRt::getBM(
  attributes = c(
    "hgnc_symbol",
    "transcript_biotype"
  ),
  filters = c("transcript_biotype"),
  values = list("protein_coding"),
  mart = mart
)


counts_data <- counts_raw[rownames(counts_raw) %in% genes$hgnc_symbol ,]

colnames(counts_data) <- sub(
  pattern = "^fpkm_uq_unstranded_",
  replacement = "",
  x = colnames(counts_data))

colnames(counts_data) <- sub(
  pattern = "-01R", 
  replacement = "",
  x = colnames(counts_data))

head(colnames(counts_data))


# ---------- 2 - METADATA PREPREOCCESSING ----------------

metadata_query <- XenaGenerate(subset = XenaDatasets == "TARGET-OS.clinical.tsv") %>% 
  XenaQuery()

# 2.2 Download 

xe_download <- XenaDownload(metadata_query, destdir = os_directory)

metadata_raw <- XenaPrepare(xe_download)

metadata_os <- metadata_raw

metadata_os <- metadata_os %>% 
  filter(sample %in% colnames(counts_data))

metadata_os <- metadata_os %>% 
  mutate(
    metastasis_at_diagnosis = ifelse(
      metastasis_at_diagnosis.diagnoses == "Metastasis, NOS",
      yes = 1,
      no = 0 
    ),
    survival_status = ifelse(
      vital_status.demographic == "Dead",
      yes = 1,
      no = 0
    )
  )




