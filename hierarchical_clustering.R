library(factoextra)

# List obtained from Boruta

gene_list <- c("APEX2", "ARHGAP1", "ARHGEF39", "CCDC97", "CGREF1", "CLUAP1", "COL22A1", "CPE", "CTNNBIP1", "CYFIP1", "DHRS11", "DLX1", "ERCC4", "F13A1", "FAM110D", "FAT1", "FKBP11", "GALNT14", "GBP1", "GMIP", "GRAMD1B", "HSD11B2", "INPP4A", "KERA", "KIF25", "LGR6", "LURAP1L", "MEF2A", "MRTFB", "MXI1", "NUBP1", "SF3B3", "SLC12A4", "SLC45A4", "SLC8A3", "STAT5B", "TIMM50", "TPD52", "TRIM68", "TSHZ3", "UBE2D4", "UNC5B", "VMP1")

# Keep only genes for clustering

vst_counts_hc <- vst_counts[rownames(vst_counts) %in% gene_list, ]

# Transpose so pt in rows and geenes in columns

vst_counts_hc <-
  vst_counts_hc %>%
  t()

# Clustering ----------------------------------------------------------

fviz_nbclust(vst_counts_hc, FUN = hcut, method = "silhouette")

# Distance matrix between samples

dist_counts <- get_dist(vst_counts_hc, method = "pearson")

# Clustering 

hc_counts <- hclust(dist_counts, method = "ward.D2")

# Clustering characteristics ------------------------------------------

# Dendrogram

plot(
  hc.out_microarray.brca,
  labels = FALSE,
  hang = -1,
  main = "CLUSTERING JERÁRQUICO METABRIC"
)

rect.hclust(hc.out_microarray.brca,
            k = 3,
            border = "purple") # bottom up approach

# Tree of clusters

clusters <- cutree(hc.out_microarray.brca, k = 3)

# Observe how many patients in each cluster

table(clusters)

# Create object as data frame

cluster_df <- data.frame(clusters = clusters)

# Asign a column named sample with the rownames such that we can then merge based on that column

cluster_df$sample <- rownames(vst_counts_hc)

# Merge object so that in the metadata there is a column corresponding to that patients cluster

metadata_os <- 
  metadata_os %>% 
  left_join(cluster_df, by = "sample")


 rm(list = setdiff(ls(), c("vst_counts", "counts_data", "metadata_os")))
