library(dplyr)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
library(msigdbr)
library(aplot)
library(ggtree)

#############################################################################
#> Script to perform differential expression analysis on TARGET-OS patients 
#> utilizing the means of each gene in each cluster.
#> Steps: 
##> Calculate mean of each gene of the patients in each cluster (similar to
##> nearest centroids)
##> Rank genes based on mean
##> GSEA on ranked list
#> 
#>  
#> Inputs: metadatao_os, vst_counts
#> 
#> Outputs: No outputs used in further analsis
#> 
#> Results: 
##> Objects with GSEA results:
###> c1_cent_GO
###> c2_cent_hallmark
###> c3_cent_GO
###> c1_cent_hallmark
###> c2_cent_GO
###> c3_cent_hallmark
#
#############################################################################

# Load data

metadata_os <- readRDS("./output_data/metadata_os.RDS")
vst_counts <- readRDS("./output_data/vst_counts.RDS")


# Scale counts

scaled_counts <- scale(t(vst_counts)[metadata_os$sample, ])

scaled_counts_df <- as.data.frame(scaled_counts)

all(rownames(scaled_counts_df) == metadata_os$sample) # Same order

# Asign clusters

scaled_counts_df$clusters <- metadata_os$clusters # Assign clusters

# Obtain means

cluster_gene_mean <- scaled_counts_df %>%
  group_by(clusters) %>%
  summarise(
    across(everything(), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  tibble::column_to_rownames("clusters")

# Rank genes

for (i in 1:nrow(cluster_gene_mean)) {

  rownames(cluster_gene_mean[i,])
  
  res <- unlist((cluster_gene_mean[i,]))
  
  names(res) <- colnames(cluster_gene_mean)
  
  gene_list <- na.omit(sort(res, decreasing = TRUE))
  
  # GSEA with GO terms
  
  gse <- gseGO(
    geneList = gene_list, 
    ont = "ALL", 
    keyType = "SYMBOL", 
    pvalueCutoff = 0.01, 
    OrgDb = "org.Hs.eg.db",
    minGSSize = 30, 
    maxGSSize = 100 
  )
  
  
  gse_df <- as.data.frame(gse)
  
  # GSEA with msigdbr
  
  msigdbr_collections()
  
  m_df <- msigdbr(species = "Homo sapiens", category = "H")
  
  msig_t2g <- m_df %>% dplyr::select(gs_name, gene_symbol)
  
  set.seed(123) # For reproducibility
  
  gsea_res <- GSEA(
    geneList     = gene_list,
    TERM2GENE    = msig_t2g,
    pvalueCutoff = 0.01,
    pAdjustMethod = "BH",
    verbose      = FALSE
  )
  
  gsea_df <- as.data.frame(gsea_res)
  
  if (i == "1") {
    
    # Results for C1
    
    c1_cent_GO <- data.frame(row.names = gse$Description, c1  = gse$NES)
    
    
    c1_cent_hallmark <- data.frame(row.names = gsea_res$Description, c1  = gsea_res$NES)
    
  } else if (i == "2") {
    
    # C2
    
    c2_cent_GO <- data.frame(row.names = gse$Description, c2  = gse$NES)
    
    
    c2_cent_hallmark <- data.frame(row.names = gsea_res$Description, c2  = gsea_res$NES)
    
    
    
  } else if (i == "3") {
    
    # C3
    
    c3_cent_GO <- data.frame(row.names = gse$Description, c3  = gse$NES)
    
    
    c3_cent_hallmark <- data.frame(row.names = gsea_res$Description, c3  = gsea_res$NES)
  }
}

# Keep top and lower 10% of C1, C2 and C3

c1_cent_GO_10 <- c1_cent_GO %>% 
  filter(c1 > quantile(c1, 0.9) | c1 < quantile(c1, 0.1))

c2_cent_GO_10 <- c2_cent_GO %>% 
  filter(c2 > quantile(c2, 0.9) | c2 < quantile(c2, 0.1))

c3_cent_GO_10 <- c3_cent_GO %>% 
  filter(c3 > quantile(c3, 0.9) | c3 < quantile(c3, 0.1))

# Join and convert to matrix

gsea_GO_mat <- merge(c3_cent_GO_10, (merge(c1_cent_GO_10, c2_cent_GO_10, by = 0, all = TRUE) %>% tibble::column_to_rownames("Row.names")), by = 0, all = TRUE)

gsea_GO_mat <- gsea_GO_mat %>% 
  tibble::column_to_rownames("Row.names") %>% 
  as.matrix()

gsea_GO_mat[is.na(gsea_GO_mat)] <- 0

# Distances for dendofram

row_hc <- hclust(dist(gsea_GO_mat))
col_hc <- hclust(dist(t(gsea_GO_mat)))


# Heatmap

heatmap_gsea <- gsea_GO_mat %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("path") %>% 
  pivot_longer(cols = c("c1", "c2", "c3"), names_to = "clusters") %>% 
  ggplot(aes(x = clusters, y = path, fill = value)) +
  geom_tile() +
  scale_fill_distiller(palette = "Spectral", direction = -1) +
  scale_x_discrete(expand = c(0, 0), labels = c("c1" = "C1", "c3" = "C3", "c2" = "C2")) +  
  scale_y_discrete(expand = c(0, 0)) +  
  theme_classic() +
  labs(x = "Clusters", y = "Pathway", fill = "NES", title = "GSEA from clusters means") +
  theme(
  axis.text.y = element_text(size = 5),
  axis.text.x = element_text(size = 9),
  plot.title = element_text(size = 11)
)


# Dendograms

tree_right <- ggtree(row_hc) + 
  scale_x_reverse() + 
  scale_y_continuous(expand = c(0, 0))

tree_top <- ggtree(col_hc, hang = -1) + 
  layout_dendrogram() + 
  scale_y_reverse(expand = c(0, 0))

# Full heatmap

heatmap_gsea %>% 
  insert_right(tree_right, width = 0.1) %>% 
  insert_top(tree_top, height = 0.1)

################################################################################

c1_cent_hallmark_10 <- c1_cent_hallmark %>% 
  filter(c1 > quantile(c1, 0.75) | c1 < quantile(c1, 0.75))

c2_cent_hallmark_10 <- c2_cent_hallmark %>% 
  filter(c2 > quantile(c2, 0.75) | c2 < quantile(c2, 0.75))

c3_cent_hallmark_10 <- c3_cent_hallmark %>% 
  filter(c3 > quantile(c3, 0.75) | c3 < quantile(c3, 0.75))

# Join and convert to matrix

gsea_hallmark_mat <- merge(c3_cent_hallmark_10, (merge(c1_cent_hallmark_10, c2_cent_hallmark_10, by = 0, all = TRUE) %>% tibble::column_to_rownames("Row.names")), by = 0, all = TRUE)

gsea_hallmark_mat <- gsea_hallmark_mat %>% 
  tibble::column_to_rownames("Row.names") %>% 
  as.matrix()

gsea_hallmark_mat[is.na(gsea_hallmark_mat)] <- 0


row_hc <- hclust(dist(gsea_hallmark_mat))
col_hc <- hclust(dist(t(gsea_hallmark_mat)))


heatmap_gsea_hm <- gsea_hallmark_mat %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("path") %>% 
  pivot_longer(cols = c("c1", "c2", "c3"), names_to = "clusters") %>% 
  ggplot(aes(x = clusters, y = path, fill = value)) +
  geom_tile() +
  scale_fill_distiller(palette = "Spectral", direction = -1) +
  scale_x_discrete(expand = c(0, 0), labels = c("c1" = "C1", "c3" = "C3", "c2" = "C2")) +  
  scale_y_discrete(expand = c(0, 0)) +  
  theme_classic() +
  labs(x = "Clusters", y = "Pathway", fill = "NES", title = "GSEA from clusters means")


tree_right_hm <- ggtree(row_hc) + 
  scale_x_reverse() + 
  scale_y_continuous(expand = c(0, 0))

tree_top_hm <- ggtree(col_hc, hang = -1) + 
  layout_dendrogram() + 
  scale_y_reverse(expand = c(0, 0))

heatmap_gsea_hm %>% 
  insert_right(tree_right_hm, width = 0.1) %>% 
  insert_top(tree_top_hm, height = 0.1)



# Save results for GSEA GO
write.csv(c1_cent_GO, "./results/diffex_gsea_target/c1_cent_GO.csv")
write.csv(c2_cent_GO, "./results/diffex_gsea_target/c2_cent_GO.csv")
write.csv(c3_cent_GO, "./results/diffex_gsea_target/c3_cent_GO.csv")

# Save results for GSEA HM
write.csv(c1_cent_hallmark, "./results/diffex_gsea_target/c1_cent_hallmark.csv")
write.csv(c2_cent_hallmark, "./results/diffex_gsea_target/c2_cent_hallmark.csv")
write.csv(c3_cent_hallmark, "./results/diffex_gsea_target/c3_cent_hallmark.csv")


