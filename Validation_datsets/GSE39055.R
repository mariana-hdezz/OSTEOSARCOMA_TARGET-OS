# -------- PREPROCESSING OF GSE39055 DATABASE FOR EXTERNAL VALIDATION -------

library(GEOquery)
library(tidyverse)
library(Biobase)
library(limma)
library(illuminaHumanWGDASLv4.db)

# Set directory
gse39055_directory <- "~/Documents/OSTEOSARCOMA/R.project/Hueso/GSE39055"


# 1.--------------- Load processed expression data -----------------

#NOTE: No normalization step was needed because GSE39055 had already been 
# normalized using VST and quantile normalization

getGEOSuppFiles("GSE39055", baseDir = gse39055_directory)

gse39055 <- getGEO("GSE39055", GSEMatrix = TRUE, getGPL = FALSE)

length(gse39055)

gse39055_data <- gse39055[[1]]

expr_matrix <- exprs(gse39055_data ) #Extract expression matrix

dim(expr_matrix)


# Annotation

probe_gene <- AnnotationDbi::select(
  illuminaHumanWGDASLv4.db,
  keys = rownames(expr_matrix),
  columns = "SYMBOL",
  keytype = "PROBEID"
)

dim(probe_gene)
sum(is.na(probe_gene$SYMBOL))
length(unique(probe_gene$SYMBOL[!is.na(probe_gene$SYMBOL)]))




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
    REC_STAT = case_when(`recurrence:ch1` == "Y" ~ 1, `recurrence:ch1` == "N" ~ 0), # 1: recurrence, 0: no recurrence/last follow-up
    REC_MON = as.numeric(`time until first recurrence or latest follow-up (months):ch1`),
    SURV_STAT = case_when(`death:ch1` == "Y" ~ 1, `death:ch1` == "N" ~ 0),
  )



