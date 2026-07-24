# This script preprocesses counts data and metadata

#> Main output objects:
#> counts_data: Unnormalized object for future differential expression
#> vst_counts: VST normalized object for most downstream analysis
#> metadata_os: metadata


library(TCGAbiolinks)
library(data.table)
library(dplyr)
library(biomaRt)
library(UCSCXenaTools)
library(DESeq2)


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
  pattern = "^unstranded_",
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


# Check for duplicates

# Ensembl IDs

sum(duplicated(rownames(counts_raw)))

# Duplicate samples

sum(duplicated(colnames(counts_raw)))

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
  pattern = "^unstranded_",
  replacement = "",
  x = colnames(counts_data))

colnames(counts_data) <- sub(
  pattern = "-01R", 
  replacement = "",
  x = colnames(counts_data))

colnames(counts_data) <- sub(
  pattern = "-01A", 
  replacement = "",
  x = colnames(counts_data))

head(colnames(counts_data))

vst_counts <- vst(as.matrix(counts_data), blind = TRUE)[-1,]

# ---------- 2 - METADATA PREPREOCCESSING ----------------

# install.packages('BiocManager')
# BiocManager::install("seandavi/TargetOsteoAnalysis")

library(TargetOsteoAnalysis)


metadata_raw <- TargetOsteoAnalysis::target_load_clinical()

names(metadata_raw)[names(metadata_raw) == "TARGET USI"] <- "sample"

# Check for sample duplicates and keep only one

sum(duplicated(metadata_raw$sample))

metadata_raw$sample[duplicated(metadata_raw$sample)]


metadata_raw <- metadata_raw %>%
  distinct(sample, .keep_all = TRUE) # Keep only 1 sample per patient 


# Keep only the patients with available counts
metadata_os <- metadata_raw %>% 
  filter(sample %in% colnames(vst_counts))


# Add modified columns needed for further analysis (metastasis and survival)

# Complete metadata
metadata_os <- metadata_os %>%
  mutate(
    metastasis_at_diagnosis = ifelse(
      `Disease at diagnosis` == "Metastatic (confirmed)" |
        `Disease at diagnosis` == "Metastatic",
      yes = 1,
      no = 0
    ),
    survival_stat = ifelse(`Vital Status` == "Dead", yes = 1, no = 0),
    survival_time = `Overall Survival Time in Days`,
    relapse_stat = ifelse(`First Event` == "Relapse", yes = 1, no = 0), # 1 = relapse, 0 = death, censored, SMN or no EVENT
    time_to_first_event = `Time to First Event in Days`,
  )




# # Survival analysis 
# metadata_os_surv <- metadata_os %>% 
#   as.data.frame()
#   mutate(EVENT_STAT = as.numeric(survival_status),
#          EVENT_DAYS = as.numeric (survival_time))



