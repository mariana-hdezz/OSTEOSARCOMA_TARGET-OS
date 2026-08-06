library(oligo)
library(GEOquery)
library(dplyr)
library(tibble)
library(limma)
library(hgu133a.db)

# getGEOSuppFiles(GEO = "GSE21257", makeDirectory = TRUE, baseDir = "d:/")



# getGEOSuppFiles(GEO = "GSE33382", makeDirectory = TRUE, baseDir = "d:/")



gpl <- getGEO("GPL10295", AnnotGPL = TRUE)

annot <- Table(gpl) %>% 
  dplyr::select(ID,
                Symbol,
                Illumina_ProbeID,
                Illumina_Gene,
                GB_ACC)


# GSE21257 ----------------------------------------------------------------

# Counts

gset_gse21257 <- getGEO("GSE21257", GSEMatrix = TRUE)

pheno_gse21257 <- pData(gset_gse21257[[1]])

unique(pheno_gse21257$data_processing)

counts_data_norm_gse21257 <- exprs(gset_gse21257[[1]]) %>% 
  as.data.frame()

counts_data_norm_gse21257 <- counts_data_norm_gse21257[rownames(counts_data_norm_gse21257) %in% annot$ID, ]

annot_gse21257 <- annot[annot$ID %in% rownames(counts_data_norm_gse21257), ]

counts_data_gse21257 <- counts_data_norm_gse21257 %>% 
  rownames_to_column("ID") %>% 
  left_join(annot_gse21257, by = "ID") %>% 
  select(- c(ID,
             Illumina_ProbeID,
             Illumina_Gene,
             GB_ACC
             )) 

counts_data_gse21257$var <- apply(counts_data_gse21257[, !(colnames(counts_data_gse21257) %in% "Symbol")], 1, function(i){var(i, na.rm = TRUE)})

counts_data_gse21257 <- counts_data_gse21257[!is.na(counts_data_gse21257$Symbol) & counts_data_gse21257$Symbol != "", ]

counts_data_gse21257 <- 
  counts_data_gse21257 %>% 
  group_by(Symbol) %>% 
  slice_max(order_by = var, n = 1, with_ties = FALSE) %>% 
  ungroup() %>% 
  select(- var) %>% 
  column_to_rownames("Symbol")

# Metadata


# GSE33382 ----------------------------------------------------------------


gset_gse33382 <- getGEO("GSE33382", GSEMatrix = TRUE)

pheno_gse33382 <- pData(gset_gse33382[[1]])

unique(pheno_gse33382$data_processing)

counts_data_norm_gse33382 <- exprs(gset_gse33382[[1]]) %>% 
  as.data.frame()


counts_data_norm_gse33382 <- counts_data_norm_gse33382[rownames(counts_data_norm_gse33382) %in% annot$ID, ]

annot_gse33382 <- annot[annot$ID %in% rownames(counts_data_norm_gse33382), ]

counts_data_gse33382 <- counts_data_norm_gse33382 %>% 
  rownames_to_column("ID") %>% 
  left_join(annot_gse33382, by = "ID") %>% 
  select(- c(ID,
             Illumina_ProbeID,
             Illumina_Gene,
             GB_ACC
  )) 

counts_data_gse33382$var <- apply(counts_data_gse33382[, !(colnames(counts_data_gse33382) %in% "Symbol")], 1, function(i){var(i, na.rm = TRUE)})

counts_data_gse33382 <- counts_data_gse33382[!is.na(counts_data_gse33382$Symbol) & counts_data_gse33382$Symbol != "", ]

counts_data_gse33382 <- 
  counts_data_gse33382 %>% 
  group_by(Symbol) %>% 
  slice_max(order_by = var, n = 1, with_ties = FALSE) %>% 
  ungroup() %>% 
  select(- var) %>% 
  column_to_rownames("Symbol")


# Metadata