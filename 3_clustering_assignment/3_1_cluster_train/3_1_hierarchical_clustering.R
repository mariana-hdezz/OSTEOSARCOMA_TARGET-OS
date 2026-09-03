library(cluster)
library(Boruta)

#############################################################################
#> Hierarchical clustering on train set (TARGET_OS patients), utilizing
#> the gene list obtained from Boruta.
#> 
#> It first creates the signature with the boruta output
#> 
#> Inputs: survival_signature, vst_counts, metadata_os, boruta_signature
#> 
#> Outputs: Overwrite metadata_os to inclide clusters column for future analysis,
#> gene_signature
#> 
#############################################################################

# List obtained from Boruta

vst_counts <- readRDS("./output_data/vst_counts.RDS")

metadata_os <- readRDS("./output_data/metadata_os.RDS")

gene_signature <- scan("output_data/gene_signature.csv", sep = ",", what = character())


# Keep only genes for clustering

vst_counts_hc <- vst_counts[rownames(vst_counts) %in% gene_signature, ]

# Transpose so pt in rows and geenes in columns

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
     main = "CLUSTERING JERÁRQUICO TARGET-OS")

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

metadata_os <-
  metadata_os %>%
  left_join(cluster_df, by = "sample")

saveRDS(metadata_os, "./output_data/metadata_os.RDS")

rm(list = ls())
gc()
