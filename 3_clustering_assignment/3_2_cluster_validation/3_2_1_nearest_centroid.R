
#############################################################################
#> Nearest centroid calculation and assigning clusters to the test method
#> 
#> Inputs: counts_data_gse21257, counts_data_gse33382, counts_data_gse39055,
#> metadata_gse21257, metadata_33382, metadata_gse39055, vst_counts, metadata_os,
#> gene_signature, gene_signature_gse
#> 
#> Outputs: Overwrite metadata_gse21257, metadata_33382, metadata_gse39055
#> to contain clusters. Creates metadata_gse33382_for_merge and metadata_gse21257_for_merge
#> where it contains the objects for limma diff ex metadata
#> 
#############################################################################

library(dplyr)
library(tibble)

dictionary_genes <- 
  list(
    "GRAMD1B" = "GRAMD1B",
    "KIF25" = "KIF25",
    "GRRP1" = "FAM110D",
    "ITGA10" = "ITGA10",
    "GBP1" = "GBP1",
    "MXI1" = "MXI1",
    "TRIM68" = "TRIM68",
    "C9orf100" = "ARHGEF39",
    "TPD52" = "TPD52",
    "CGREF1" = "CGREF1",
    "KERA" = "KERA",
    "UBE2D4" = "UBE2D4",
    "BNIP3" = "BNIP3",
    "MKL2" = "MRTFB",
    "COL22A1" = "COL22A1",
    "CPE" = "CPE",
    "LSG1" = "LSG1",
    "CTNNBIP1" = "CTNNBIP1",
    "LGR6" = "LGR6",
    "RHBDL2" = "RHBDL2",
    "TMEM49" = "VMP1",
    "HSD11B2" = "HSD11B2",
    "EID-3" = "EID2B",
    "UHRF2" = "UHRF2",
    "SLC12A4" = "SLC12A4",
    "MEF2A" = "MEF2A",
    "ACTA2" = "ACTA2",
    "MSC" = "MSC",
    "PIP5K1C" = "PIP5K1C",
    "F13A1" = "F13A1",
    "ARHGAP28" = "ARHGAP28",
    "FKBP11" = "FKBP11",
    "PSMC4" = "PSMC4",
    "RGS9" = "RGS9",
    "SH3PXD2A" = "SH3PXD2A",
    "CD163" = "CD163",
    "SF3B3" = "SF3B3"
  )


counts_data_gse21257 <- readRDS("output_data/counts_data_gse21257.RDS")
counts_data_gse33382 <- readRDS("output_data/counts_data_gse33382.RDS")
counts_data_gse39055 <- readRDS("output_data/counts_data_gse39055.RDS")
metadata_gse21257 <- readRDS("output_data/metadata_gse21257.RDS")
metadata_33382 <- readRDS("output_data/metadata_33382.RDS")
metadata_gse39055 <- readRDS("output_data/metadata_gse39055.RDS")

vst_counts <- readRDS("output_data/vst_counts.RDS")
metadata_os <- readRDS("output_data/metadata_os.RDS")

gene_signature <- scan("output_data/gene_signature.csv", sep = ",", what = character())
gene_signature_gse <- scan("output_data/gene_signature_gse.csv", sep = ",", what = character()) 

counts_data_list <- c("counts_data_gse21257", "counts_data_gse33382", "counts_data_gse39055")

metadata_list <- c("metadata_gse21257", "metadata_33382", "metadata_gse39055")

for (t in 1:3) {


    
  # Load needed objects (counts_data_test_centr, metadata_centroids, gene_signature, gene_signature_gse)
  
  counts_data_test_centr <- get(counts_data_list[t])
  
  metadata_centroids <- get(metadata_list[t]) 
  
  
  # Scale vst and convert to df
  scaled_counts <- scale(t(vst_counts)[metadata_os$sample, ])
  
  scaled_counts_df <- as.data.frame(scaled_counts)
  
  scaled_counts_gse <- scale(t(counts_data_test_centr)[metadata_centroids$geo_accession, ])
  
  scaled_counts_gse_df <- as.data.frame(scaled_counts_gse)
  
  # keep the genes in signature
  
  scaled_counts_df <- scaled_counts_df[, gene_signature] 
  
  if(all(gene_signature %in% colnames(scaled_counts_gse_df)) ){
    
    scale_counts_gse_test <- scaled_counts_gse_df[, gene_signature]
    
  }else{
    scale_counts_gse_df_41 <- scaled_counts_gse_df[, gene_signature_gse] 
    
    all(colnames(scale_counts_gse_df_41) %in% names(dictionary_genes))
    
    # Position of the colnames on the dictionary
    
    idx <- match(colnames(scale_counts_gse_df_41), names(dictionary_genes))
    
    # Asign the other term of the dictionary to the previously determined position
    
    colnames(scale_counts_gse_df_41)[!is.na(idx)] <- unlist(dictionary_genes)[idx[!is.na(idx)]]
    
    scale_counts_gse_test <-  scale_counts_gse_df_41
    
    scaled_counts_df <- scaled_counts_df[, colnames(scale_counts_gse_test)]

    
  }
  
  # Centroids from training
  
  all(rownames(scaled_counts_df) == metadata_os$sample) # Same order
  
  scaled_counts_df$clusters <- metadata_os$clusters # Assign clusters
  
  train_centroids <- aggregate(. ~ clusters, data = scaled_counts_df, FUN = mean) %>% 
    column_to_rownames("clusters")
  
  
  # Convert both to matrices
  
  train_centroids_mat <- as.matrix(train_centroids[, colnames(train_centroids)])
  
  test_mat <- as.matrix(scale_counts_gse_test[, colnames(train_centroids)])
  
  all(colnames(train_centroids_mat) == colnames(test_mat))
  
  # Eucliean distances
  
  cluster_list <- list() # output list
  
  for (i in 1:nrow(test_mat)) { # Repeat to all rows on test matrix
    x <- as.numeric(test_mat[i, ]) # Convert that row into numeric vector
    sample_name <- rownames(test_mat)[i] # maintain the patient name
    
    distances <- numeric(nrow(train_centroids_mat)) # to register the distance with respect to patient
    names(distances) <- rownames(train_centroids_mat) #assign row names as names
    
    for (e in 1:nrow(train_centroids_mat)) { # Repeat the distance calculation for each cluster
      centroid <- as.numeric(train_centroids_mat[e, ]) # convert the train cluster row to numbers 
      distances[e] <- sqrt(sum((x - centroid)^2)) # calculate euclidean distance
    }
    
    
    cluster_list[[sample_name]] <- names(distances)[which.min(distances)]
  }
  
  
  cluster_val <- data.frame(
    geo_accession = names(cluster_list),
    clusters = unlist(cluster_list),
    stringsAsFactors = FALSE
  )
  
  assign(
    x = metadata_list[t], 
    value = metadata_centroids %>% left_join(cluster_val, by = "geo_accession")
  )

  print(metadata_list[t])
}

metadata_gse33382_for_merge <- 
  metadata_33382 %>% 
  dplyr::select(geo_accession,
                cohort, 
                hist_sub,
                clusters)

metadata_gse21257_for_merge <- 
  metadata_gse21257 %>% 
  dplyr::select(geo_accession,
                cohort, 
                hist_sub,
                clusters)

saveRDS(metadata_gse21257, "output_data/metadata_gse21257.RDS")
saveRDS(metadata_33382, "output_data/metadata_33382.RDS")
saveRDS(metadata_gse39055, "output_data/metadata_gse39055.RDS")
saveRDS(metadata_gse33382_for_merge, "output_data/metadata_gse33382_for_merge.RDS")
saveRDS(metadata_gse21257_for_merge, "output_data/metadata_gse21257_for_merge.RDS")

rm(list = ls())
