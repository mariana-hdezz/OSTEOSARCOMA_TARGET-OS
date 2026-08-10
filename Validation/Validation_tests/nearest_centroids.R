library(dplyr)
library(tibble)

dictionary_genes <- 
  list(
    "APEX2" = "APEX2",
    "ARHGAP1" = "ARHGAP1",
    "C9orf100" = "ARHGEF39",
    "MGC20255" = "CCDC97",
    "CGREF1" = "CGREF1",
    "CLUAP1" = "CLUAP1",
    "COL22A1" = "COL22A1",
    "CPE" = "CPE",
    "CTNNBIP1" = "CTNNBIP1",
    "CYFIP1" = "CYFIP1",
    "MGC4172" = "DHRS11",
    "DLX1" = "DLX1",
    "ERCC4" = "ERCC4",
    "F13A1" = "F13A1",
    "GRRP1" = "FAM110D",
    "FAT" = "FAT1",
    "FKBP11" = "FKBP11",
    "GALNT14" = "GALNT14",
    "GBP1" = "GBP1",
    "GMIP" = "GMIP",
    "GRAMD1B" = "GRAMD1B",
    "HSD11B2" = "HSD11B2",
    "INPP4A" = "INPP4A",
    "KERA" = "KERA",
    "KIF25" = "KIF25",
    "LGR6" = "LGR6",
    "C9orf150" = "LURAP1L",
    "MEF2A" = "MEF2A",
    "MKL2" = "MRTFB",
    "MXI1" = "MXI1",
    "NUBP1" = "NUBP1",
    "SF3B3" = "SF3B3",
    "SLC12A4" = "SLC12A4",
    "SLC45A4" = "SLC45A4",
    "SLC8A3" = "SLC8A3",
    "STAT5B" = "STAT5B",
    "TIMM50" = "TIMM50",
    "TPD52" = "TPD52",
    "TRIM68" = "TRIM68",
    "ZNF537" = "TSHZ3",
    "UBE2D4" = "UBE2D4",
    "UNC5B" = "UNC5B",
    "TMEM49" = "VMP1"
  )

# Load needed objects (counts_data_test_centr, metadata_centroids, pucky, pucky_gse)

counts_data_test_centr <- counts_data_gse21257

metadata_centroids <- metadata_gse21257 


# Scale vst and convert to df
scaled_counts <- scale(t(vst_counts)[metadata_os$sample, ])

scaled_counts_df <- as.data.frame(scaled_counts)

scaled_counts_gse <- scale(t(counts_data_test_centr)[metadata_centroids$geo_accession, ])

scaled_counts_gse_df <- as.data.frame(scaled_counts_gse)

# keep the genes in signature

scaled_counts_df <- scaled_counts_df[, pucky] 

if(all(pucky %in% colnames(scaled_counts_gse_df)) ){
  
  scaled_counts_gse_test <- scaled_counts_gse_df[, pucky]
  
}else{
  scale_counts_gse_df_41 <- scaled_counts_gse_df[, pucky_gse] 
  
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
  x <- as.numeric(test_mat[i, ]) # Convert that row intonumeric vector
  sample_name <- rownames(test_mat)[i] # maintain the patient name
  
  distances <- numeric(nrow(train_centroids_mat)) # to register the distance with respect to patient
  names(distances) <- rownames(train_centroids_mat) #assign rownames as names
  
  for (e in 1:nrow(train_centroids_mat)) { # Repeat the distance calculation for each cluster
    centroid <- as.numeric(train_centroids_mat[e, ]) # convert the train cluster row to numbers 
    distances[e] <- sqrt(sum((x - centroid)^2)) # calculate euqlidean distance
  }
  
  
  cluster_list[[sample_name]] <- names(distances)[which.min(distances)]
}


cluster_val <- data.frame(
  geo_accession = names(cluster_list),
  clusters = unlist(cluster_list),
  stringsAsFactors = FALSE
)


metadata_centroids <- metadata_centroids %>% 
  left_join(cluster_val, by = "geo_accession")
