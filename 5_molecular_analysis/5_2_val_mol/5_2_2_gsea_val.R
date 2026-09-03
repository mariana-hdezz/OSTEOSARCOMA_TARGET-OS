library(stats)
library(fgsea)


res$entrez <- mapIds(org.Hs.eg.db,
                     keys=rownames(res),
                     column="ENTREZID", 
                     keytype="SYMBOL", 
                     multiVals="first") 

res$name <-   mapIds(org.Hs.eg.db,
                     keys=rownames(res), 
                     column="GENENAME",
                     keytype="SYMBOL",
                     multiVals="first")



gene_symbols <- rownames(count_data)


res_lfc <- res$logFC

# Asign the gene names to log fold change object

names(res_lfc) <- rownames(res)

# Delete NA

gene_list <- na.omit(res_lfc)


# Sort descending

gene_list <- sort(gene_list, decreasing = TRUE)


# GSEA object GO

gse <- gseGO(
  geneList = gene_list, 
  ont = "ALL", 
  keyType = "SYMBOL", 
  pvalueCutoff = 0.01, 
  OrgDb = "org.Hs.eg.db",
  minGSSize = 30, 
  maxGSSize = 100,
  eps = 0
)
  
gsea_go_df <- as.data.frame(gse)


# gsea hallmarks

msigdbr_collections()

m_df <- msigdbr(species = "Homo sapiens", category = "H")

msig_t2g <- m_df %>% dplyr::select(gs_name, gene_symbol)

set.seed(123) # For reproducibility

gsea_res <- GSEA(
  geneList     = gene_list,
  TERM2GENE    = msig_t2g,
  pvalueCutoff = 0.01,
  pAdjustMethod = "BH",
  verbose      = FALSE,
  eps = 0
)

gsea_go_hm <- as.data.frame(gsea_res)

# Create objects for heatmaps ---------------------------------------------

# if else conditional logic where if in the contrast.matrix object generated on 5_2_1diffex_val.R the contrast is between cluster i and j 
# then that section of the script will run
# For example, since the contrast matrix looks like this in c3 vs c1 
#> contrast.matrix
#> cluster 1          -1
#> cluster 2           0
#> cluster 3           1
#> Then the sum of (cluster 1)^2 and (clsuter 3)^2 = 2 meanwhile the sum of (cluster 2)^2 with any of the other 2 will yield 1 so that means 
#> cluster 2 is not taken into account in this round

if((contrast.matrix[2,])^2 + (contrast.matrix[3,])^2 == 2){

  
  out_path_3v2_go <- paste0("results/diffex_gsea_gse/", "gsea_c3_vs_c2_GO_", t, ".csv")
  out_path_3v2_hm <- paste0("results/diffex_gsea_gse/", "gsea_c3_vs_c2_hm_", t, ".csv")
  
  write.csv(gsea_go_df, out_path_3v2_go)
  write.csv(gsea_go_hm, out_path_3v2_hm)
  
}else if((contrast.matrix[2,])^2 + (contrast.matrix[1,])^2 == 2){

  
  out_path_1v2_go <- paste0("results/diffex_gsea_gse/", "gsea_c1_vs_c2_GO_", t, ".csv")
  out_path_1v2_hm <- paste0("results/diffex_gsea_gse/", "gsea_c1_vs_c2_hm_", t, ".csv")
  
  write.csv(gsea_go_df, out_path_1v2_go)
  write.csv(gsea_go_hm, out_path_1v2_hm)
  
}else if((contrast.matrix[3,])^2 + (contrast.matrix[1,])^2 == 2){


  
  out_path_3v1_go <- paste0("results/diffex_gsea_gse/", "gsea_c3_vs_c1_GO_", t, ".csv")
  out_path_3v1_hm <- paste0("results/diffex_gsea_gse/", "gsea_c3_vs_c1_hm_", t, ".csv")
  
  write.csv(gsea_go_df, out_path_3v1_go)
  write.csv(gsea_go_hm, out_path_3v1_hm)

}
