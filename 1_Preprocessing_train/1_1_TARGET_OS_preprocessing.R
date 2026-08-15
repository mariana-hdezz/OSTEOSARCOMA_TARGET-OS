# This script preprocesses counts data and metadata

#> Main output objects:
#> counts_data: Unnormalized object for future differential expression
#> vst_counts: VST normalized object for most downstream analysis
#> metadata_os: complete metadata
#> metadata_os_surv: metadata useful for survival analysis
#> metadata_os_rec: metadata useful for recurrence analysis
#> metadata_os_met: metadata useful for metastasis analysis


library(TCGAbiolinks)
library(data.table)
library(dplyr)
library(biomaRt)
library(UCSCXenaTools)
library(DESeq2)


os_directory <- "~/Documents/OSTEOSARCOMA/R.project/Hueso" # Directory for download and retrieval of data


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


counts_raw <- os_data %>% 
  dplyr::select(all_of(counts_col),
                gene_name,
                gene_type)

gene_dist <- os_data %>% 
  dplyr::select(gene_name,
                gene_type,
                gene_id)


gene_dist$ensembl <- gsub("\\..", "", gene_dist$gene_id)

counts_raw$variance <- apply(counts_raw %>% dplyr::select(-gene_name),1 , var, na.rm = TRUE)

counts_raw <- counts_raw %>% 
  group_by(gene_name) %>%
  slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% 
  ungroup %>% 
  tibble::column_to_rownames("gene_name") %>% 
  dplyr::select(-variance)



class(counts_raw)
dim(counts_raw)


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
  filter(sample %in% colnames(counts_data))


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
    survival_stat = case_when(`Vital Status` == "Dead" ~ 1,
                              `Vital Status` == "Alive" ~ 0,
                              .default = NA ),
    survival_time = `Overall Survival Time in Days`,
    relapse_stat = ifelse(`First Event` == "Relapse", yes = 1, no = 0), # 1 = relapse, 0 = death, censored, SMN or no EVENT
    time_to_first_event = `Time to First Event in Days`,
    hist_res = gsub(".*\\b(\\d+/\\d+)\\b.*", "\\1", `Histologic response`),
    necrosis_at_surg = as.numeric(gsub("< ", "", `Percent necrosis`)),
    huvos_bin = case_when(
      hist_res == "1/2" | hist_res == "0-90" | necrosis_at_surg < 90 ~ "Bad response",
      hist_res == "3/4" | hist_res == "91-100" | necrosis_at_surg > 90 ~ "Good response"
    ),
    age_yr = `Age at Diagnosis in Days` / 365,
    ethnicity = Ethnicity,
    race = Race,
    tumor_site = `Primary tumor site`,
    met_site = `Metastasis site`,
    gender = Gender,
    tumor_side = `Specific tumor side`,
    tumor_site_spec = `Specific tumor site`,
    tumor_reg = `Specific tumor region`,
    rel_type = `Relapse Type`,
    time_to_f_rel = `Time to first relapse in days`,
    time_to_smn = `Time to first SMN in days`,
    surgery = `Definitive Surgery`,
    site_progression = `Primary site progression`,
    treatment = Therapy,
    first_event = "First Event"
  ) %>% 
  dplyr::select(- c(
    "Gender",
    "Race",
    "Ethnicity",
    "Age at Diagnosis in Days",
    "First Event",
    "Time to First Event in Days",
    "Vital Status",
    "Overall Survival Time in Days",
    "Year of Diagnosis",
    "Year of Last Follow Up",
    "Protocol",
    "Disease at diagnosis",
    "Metastasis site",
    "Primary tumor site",
    "Specific tumor site",
    "Specific tumor side",
    "Specific tumor region",
    "Definitive Surgery",
    "Primary site progression",
    "Site of initial relapse",
    "Time to first relapse in days",
    "Time to first enrollment on relapse protocol in days",
    "Time to first SMN in days",
    "Time to death in days",
    "Histologic response",
    "Percent necrosis",
    "Relapse Type",
    "Therapy",
    "Comment",
    "cohort",
    "Percent necrosis at Definitive Surgery"
  ))


# Survival analysis
metadata_os_surv <- metadata_os %>% 
  as.data.frame() %>% 
  mutate(EVENT_STAT = as.numeric(survival_stat),
         EVENT_DAYS = as.numeric(survival_time)) %>% 
  dplyr::select(-c(survival_stat, survival_time, metastasis_at_diagnosis, relapse_stat))


# Relapse analysis
metadata_os_rec <- metadata_os %>%
  as.data.frame() %>%
  mutate(EVENT_STAT = as.numeric(relapse_stat),
         EVENT_DAYS = as.numeric(time_to_first_event)) %>% 
  dplyr::select(-c(survival_stat, survival_time, relapse_stat, metastasis_at_diagnosis))


# Metastasis analysis
metadata_os_met <- metadata_os %>% 
  as.data.frame() %>% 
  mutate(EVENT_STAT = metastasis_at_diagnosis) %>% 
  dplyr::select(-c(survival_stat, survival_time, relapse_stat, metastasis_at_diagnosis))



# Counts data final adjust ------------------------------------------------

counts_data <- counts_data[, colnames(counts_data) %in% metadata_os$sample]


# VST and FPKM ------------------------------------------------------------


vst_counts <- vst(as.matrix(counts_data), blind = TRUE)


# FPKM Preprocessing ------------------------------------------------------


counts_col_fpkm <- grep(
  pattern = "^fpkm_uq_unstranded_",
  x = names(os_data),
  value = TRUE)


counts_raw_fpkm <- os_data %>% 
  dplyr::select(all_of(counts_col_fpkm),
                gene_name,
                gene_type)

counts_raw_fpkm$variance <- apply(counts_raw_fpkm %>% dplyr::select(-gene_name),1 , var, na.rm = TRUE)

counts_raw_fpkm <- counts_raw_fpkm %>% 
  group_by(gene_name) %>%
  slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% 
  ungroup %>% 
  tibble::column_to_rownames("gene_name") %>% 
  dplyr::select(-variance)



class(counts_raw_fpkm)
dim(counts_raw_fpkm)


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


fpkm_data <- counts_raw_fpkm[rownames(counts_raw_fpkm) %in% genes$hgnc_symbol ,]



colnames(fpkm_data) <- sub(
  pattern = "^fpkm_uq_unstranded_",
  replacement = "",
  x = colnames(fpkm_data))

colnames(fpkm_data) <- sub(
  pattern = "-01R", 
  replacement = "",
  x = colnames(fpkm_data))

colnames(fpkm_data) <- sub(
  pattern = "-01A", 
  replacement = "",
  x = colnames(fpkm_data))


fpkm_data <- fpkm_data[, colnames(fpkm_data) %in% metadata_os$sample]


fpkm_data_log <- log(fpkm_data + 1)


saveRDS(vst_counts, "./output_data/vst_counts.RDS")

saveRDS(metadata_os, "./output_data/metadata_os.RDS")

saveRDS(counts_data, "./output_data/counts_data.RDS")

saveRDS(fpkm_data, "./output_data/fpkm_data.RDS")

saveRDS(gene_dist, "./output_data/gene_dist.RDS")

saveRDS(metadata_os_surv,"./output_data/metadata_os_surv.RDS" )

saveRDS(metadata_os_rec,"./output_data/metadata_os_rec.RDS" )

saveRDS(metadata_os_met,"./output_data/metadata_os_met.RDS" )



rm(list = ls())
