library(GEOquery)
library(dplyr)
library(tibble)
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
  dplyr::select(- c(ID,
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
  dplyr::select(- var) %>% 
  column_to_rownames("Symbol")

# Metadata

metadata_gse21257 <- pheno_gse21257 %>%  
  mutate(age_yr = round(as.numeric(as.character(gsub("^age: (.*) months$", "\\1", pheno_gse21257$characteristics_ch1))) / 12),
         gender = `gender:ch1`,
         hist_sub = `histological subtype:ch1`,
         tumor_loc = `tumor location:ch1`,
         huvos = `huvos grade:ch1`,
         survival_stat = ifelse(substr(`status:ch1`, 1,8) == "Deceased",
                                1,
                                0),
         survival_time = ifelse(substr(`status:ch1`, 1,8) == "Deceased",
                                gsub("Deceased (at|after) (.*) (=?months).*", "\\2", `status:ch1`),
                                gsub("Alive at (.*) (=?months).*", "\\1", `status:ch1`)),
         relapse_stat = as.numeric(as.character(ifelse(`group:ch1` == "No metastases",
                               0,
                               1))),
         relapse_time = case_when(
           relapse_stat == 0 ~ as.numeric(as.character(survival_time)),
           (relapse_stat == 1) & (substr(`group:ch1`, 1, 13) == "Metastases at") ~ as.numeric(gsub("^Metastases at (.*) (=?months).*", "\\1" , `group:ch1`)),
           (relapse_stat == 1) & (substr(`group:ch1`, 1, 13) == "Metastases pr") ~ 0
         )
         ) %>% 
  
  dplyr::select(-c(status,
                   submission_date,
                   last_update_date,
                   channel_count,
                   organism_ch1,
                   supplementary_file,
                   `contact_zip/postal_code`,
                   contact_country,
                   contact_laboratory,
                   contact_address,
                   contact_city,
                   contact_institute,
                   platform_id, data_processing,
                   description,
                   scan_protocol,
                   hyb_protocol,
                   taxid_ch1,
                   label_ch1,
                   label_protocol_ch1,
                   extract_protocol_ch1,
                   molecule_ch1,
                   contact_name,
                   supplementary_file.1,
                   data_row_count,
                   source_name_ch1,
                   characteristics_ch1.4,
                   `histological subtype:ch1`,
                   characteristics_ch1.3,
                   `huvos grade:ch1`,
                   characteristics_ch1.2,
                   `gender:ch1`,
                   characteristics_ch1,
                   characteristics_ch1.6,
                   characteristics_ch1.1,
                   `age:ch1`,
                   `tumor location:ch1`,
                   characteristics_ch1.5,
                   `status:ch1`,
                   `group:ch1`
                   
                   )
                )

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
  dplyr::select(- c(ID,
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
  dplyr::select(- var) %>% 
  column_to_rownames("Symbol")


# Metadata

pheno_gse33382

