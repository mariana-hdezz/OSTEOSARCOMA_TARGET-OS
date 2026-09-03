library(dplyr)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
library(msigdbr)
library(aplot)
library(ggtree)

#############################################################################
#> Script to perform differential expression analysis on batch corrected
#> GSE patients utilizing the means of each gene in each cluster. Same as
#> done on 3_3_GSEA_means.R
#> 
#> Steps: 
##> Calculate mean of each gene of the patients in each cluster (similar to
##> nearest centroids)
##> Rank genes based on mean
##> GSEA on ranked list
#> 
#>  
#> Inputs: metadata_gse33382_for_merge, metadata_gse21257_for_merge, counts_merged
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
###> Heatmaps
#
#############################################################################

# Load data


metadata_gse33382_for_merge <- readRDS("output_data/metadata_gse33382_for_merge.RDS")
metadata_gse21257_for_merge <- readRDS("output_data/metadata_gse21257_for_merge.RDS")

counts_merged <- readRDS("output_data/counts_batch.RDS")

gene_signature_gse <- scan("output_data/gene_signature_gse.csv", sep = ",", what = character())

metadata_merged <- bind_rows(metadata_gse33382_for_merge,
                             metadata_gse21257_for_merge)


# Scale counts

scaled_counts <- scale(t(counts_merged)[metadata_merged$geo_accession, ])

scaled_counts_df <- as.data.frame(scaled_counts)

all(rownames(scaled_counts_df) == metadata_merged$geo_accession) # Same order

# Asign clusters

scaled_counts_df$clusters <- metadata_merged$clusters # Assign clusters

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



write.csv(c1_cent_GO      , "./results/diffex_gsea_gse/c1_cent_GO.csv"      )
write.csv(c1_cent_hallmark, "./results/diffex_gsea_gse/c1_cent_hallmark.csv")
write.csv(c2_cent_GO      , "./results/diffex_gsea_gse/c2_cent_GO.csv"      )
write.csv(c2_cent_hallmark, "./results/diffex_gsea_gse/c2_cent_hallmark.csv")
write.csv(c3_cent_GO      , "./results/diffex_gsea_gse/c3_cent_GO.csv"      )
write.csv(c3_cent_hallmark, "./results/diffex_gsea_gse/c3_cent_hallmark.csv")



rm(list = ls())
gc()