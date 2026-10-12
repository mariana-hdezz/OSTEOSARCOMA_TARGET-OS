library(cluster)
library(Boruta)
library(factoextra)

#############################################################################
#> Hierarchical clustering on train set (GSE39055" patients), utilizing
#> the gene list obtained from Boruta.
#> 
#> 
#> Inputs: gene_signature,counts_data_gse39055, metadata_os
#> 
#> Outputs: Overwrite metadata_os to inclide clusters column for future analysis
#> 
#############################################################################

# List obtained from Boruta

counts_data_gse39055 <- readRDS("./output_data/counts_data_gse39055.RDS")

metadata_gse39055 <- readRDS("./output_data/metadata_gse39055.RDS")

gene_signature <- readRDS("results/boruta/boruta_conf_gse39055.RDS")
gene_signature <- unlist(gene_signature)


# Keep only genes for clustering

vst_counts_hc <- counts_data_gse39055[rownames(counts_data_gse39055) %in% gene_signature, ]

# Transpose so pt in rows and genes in columns

vst_counts_hc <-
  vst_counts_hc %>%
  t()

# Clustering ----------------------------------------------------------

fviz_nbclust(vst_counts_hc, FUN = hcut, method = "silhouette")

# Distance matrix between samples

dist_counts <- get_dist(vst_counts_hc, method = "manhattan")

# Clustering

hc_counts <- hclust(dist_counts, method = "ward.D2")


# Clustering characteristics ------------------------------------------

# Dendrogram

plot(hc_counts,
     labels = FALSE,
     hang = -1,
     main = "CLUSTERING JERÁRQUICO GSEO39055")

rect.hclust(hc_counts, k = 3, border = "purple") # bottom up approach

# Tree of clusters

clusters <- cutree(hc_counts, k = 3)

sil <- silhouette(clusters, dist_counts)

mean(sil[, "sil_width"])

# 6. Plot the silhouette profile
plot(sil, col = 2:(length(unique(as.integer(clusters))) + 1), main = "Hierarchical Silhouette Plot")

# Observe how many patients in each cluster

table(clusters)

# Create object as data frame

cluster_df <- data.frame(clusters = clusters)

# Assign a column named sample with the row names such that we can then merge based on that column

cluster_df$sample <- rownames(vst_counts_hc)

# Merge object so that in the metadata there is a column corresponding to that patients cluster

metadata_gse39055 <-
  metadata_gse39055 %>%
  left_join(cluster_df, by = "sample")

saveRDS(metadata_gse39055, "./output_data/metadata_gse39055.RDS")

rm(list = ls())
gc()
