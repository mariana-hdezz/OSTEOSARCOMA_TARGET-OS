
#############################################################################
#> Script to perform preprocessing steps on GSE21257 and GSE33382 
#> 
#> Inputs: None
#> 
#> Outputs: metadata_gse21257, metadata_33382, counts_data_gse21257, counts_data_gse33382, annot,
#> counts_merged
#> 
#############################################################################

library(GEOquery)
library(dplyr)
library(tibble)
library(readr)
library(stringr)
library(hgu133a.db)


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
  mutate(age_yr = round(parse_number(as.character(characteristics_ch1)) / 12),
         gender = `gender:ch1`,
         hist_sub = tolower(`histological subtype:ch1`),
         hist_sub = gsub("^possibly ", "", hist_sub),
         hist_sub = gsub("fibroma like", "fibroma-like", hist_sub),
         tumor_loc = `tumor location:ch1`,
         huvos = `huvos grade:ch1`,
         survival_stat = ifelse(str_starts(`status:ch1`, "Deceased"), 
                                yes = 1, 
                                no = 0),
         survival_time = as.numeric(parse_number(`status:ch1`)),
         relapse_stat = as.numeric(as.character(ifelse(`group:ch1` == "No metastases",
                               0,
                               1))),
         relapse_time = case_when(
           relapse_stat == 0 ~ survival_time,
           str_starts(`group:ch1`, "Metastases pr") ~ 0,
           str_starts(`group:ch1`, "Metastases at") ~ parse_number(
             ifelse(str_starts(`group:ch1`, "Metastases at"), `group:ch1`, NA_character_)
           ),
           TRUE ~ NA_real_
         ),
         cohort = "GSE21257"
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

metadata_33382 <- pheno_gse33382 %>% 
  mutate(
    age = gsub("^(.*) months" , "\\1",  `age:ch1`),
    gender = factor(`gender:ch1`),
    hist_sub_com = factor(gsub(" osteosarcoma$", "" , `histological subtype:ch1`)),
    huvos = factor(`huvos grade:ch1`),
    metastasis_5y = factor(`metastasis within 5yrs:ch1`),
    tumor_loc = `tumor location:ch1`,
    tissue = `type:ch1`,
    hist_sub = ifelse(hist_sub_com == "fibroblastic MFH-like" | hist_sub_com == "fibroblastic giant cell rich" | hist_sub_com == "gibroblastic",
                          "fibroblastic",
                          as.character(hist_sub_com)),
    cohort = "GSE33382"
  ) %>% 
  dplyr::select(- c(
    extract_protocol_ch1,
    `matching copy number data:ch1`,
    label_ch1,
    label_protocol_ch1,
    hyb_protocol,
    scan_protocol,
    description,
    data_processing,
    contact_name,
    contact_laboratory,
    contact_institute,
    contact_address,
    contact_city,
    `contact_zip/postal_code`,
    contact_country,
    supplementary_file,
    supplementary_file.1,
    data_row_count,
    `growth protocol:ch1`,
    title,
    status,
    submission_date,
    last_update_date,
    source_name_ch1,
    channel_count,
    organism_ch1,
    characteristics_ch1.1,
    characteristics_ch1.2,
    characteristics_ch1.3,
    characteristics_ch1.4,
    characteristics_ch1.5,
    characteristics_ch1.6,
    characteristics_ch1.7,
    taxid_ch1,
    `age:ch1`,
    `gender:ch1`,
    `histological subtype:ch1`,
    `huvos grade:ch1`,
    `matching copy number data:ch1`,
    `metastasis within 5yrs:ch1`,
    `tumor location:ch1`,
    type,
    characteristics_ch1
  )) %>% 
  filter(tissue == "biopsy")

counts_data_gse33382 <- counts_data_gse33382[, colnames(counts_data_gse33382) %in% metadata_33382$geo_accession]


# Batch correction


 
metadata_cohort <- bind_rows(metadata_33382 %>% 
              dplyr::select(geo_accession,
                     cohort),
            metadata_gse21257 %>% 
              dplyr::select(geo_accession,
                            cohort) 
            )


counts_merged <- merge(counts_data_gse33382, counts_data_gse21257, by = 0 , all = TRUE) %>% 
  tibble::column_to_rownames("Row.names")

pca_data <- prcomp(t(counts_merged), center = TRUE, scale. = FALSE)


pca_df <- as_tibble(pca_data$x) %>%
  mutate(geo_accession = rownames(pca_data$x)) %>%
  left_join(metadata_cohort, by = "geo_accession")


ggplot(pca_df) +
  aes(x = PC1, y = PC2, color = cohort) +
  geom_point(size = 4, alpha = 0.8) +
  labs(
    title = "PCA",
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal()

counts_merged <- counts_merged[, metadata_cohort$geo_accession]

counts_batch <- limma::removeBatchEffect(counts_merged, batch = metadata_cohort$cohort)


pca_data <- prcomp(t(counts_batch), center = TRUE, scale. = FALSE)


pca_df <- as_tibble(pca_data$x) %>%
  mutate(geo_accession = rownames(pca_data$x)) %>%
  left_join(metadata_cohort, by = "geo_accession")


ggplot(pca_df) +
  aes(x = PC1, y = PC2, color = cohort) +
  geom_point(size = 4, alpha = 0.8) +
  labs(
    title = "PCA",
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal()



saveRDS(metadata_gse21257, "output_data/metadata_gse21257.RDS")
saveRDS(metadata_33382, "output_data/metadata_33382.RDS")
saveRDS(counts_data_gse21257, "output_data/counts_data_gse21257.RDS")
saveRDS(counts_data_gse33382, "output_data/counts_data_gse33382.RDS")
saveRDS(annot, "output_data/annot.RDS")
saveRDS(counts_merged, "output_data/counts_merged.RDS")


rm(list = ls())