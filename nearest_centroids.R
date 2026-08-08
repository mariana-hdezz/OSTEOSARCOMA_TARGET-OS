
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
    "VMP1" = "VMP1"
  )

all(colnames(counts_gse_centr) %in% names(dictionary_genes))

# Scale vst and convert to df

scaled_counts <- scale(vst_counts_t_43)[metadata_os$sample, ]

scaled_counts_df <- as.data.frame(scaled_counts)

scaled_counts_df$clusters <- metadata_os$clusters

# Obtain centroids

train_centroids <- aggregate(. ~ clusters, data = scaled_counts_df, FUN = mean) %>% 
  column_to_rownames("clusters")

# Kepp only the counts of the gene list (of the previously mapped genes on common_siugnature.R)

counts_gse_centr <- t(counts_data_gse21257[rownames(counts_data_gse21257) %in% pucky_gse, ])

# Position of the colnames on the dictionary

idx <- match(colnames(counts_gse_centr), names(dictionary_genes))

#

colnames(counts_gse_centr)[!is.na(idx)] <- unlist(dictionary_genes)[idx[!is.na(idx)]]

counts_gse_centr_scaled <- scale(counts_gse_centr)

common_genes <- intersect(colnames(train_centroids), colnames(counts_gse_centr_scaled))

train_centroids_mat <- as.matrix(train_centroids[, common_genes])
test_mat <- as.matrix(counts_gse_centr_scaled[, common_genes])

all(colnames(train_centroids_mat) == colnames(test_mat))

cluster_list <- list()

for (i in 1:nrow(test_mat)) {
  x <- as.numeric(test_mat[i, ])
  sample_name <- rownames(test_mat)[i]
  
  distances <- numeric(nrow(train_centroids_mat))
  
  for (e in 1:nrow(train_centroids_mat)) {
    centroid <- as.numeric(train_centroids_mat[e, ])
    distances[e] <- sqrt(sum((x - centroid)^2))
  }
  
  names(distances) <- rownames(train_centroids_mat)
  
  cluster_list[[sample_name]] <- names(distances)[which.min(distances)]
}


cluster_val <- data.frame(
  geo_accession = names(cluster_list),
  clusters = unlist(cluster_list),
  stringsAsFactors = FALSE
)
