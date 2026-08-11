


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

# 9.1 Asignar al objeto de Log Fold los nombres de ENSEMBL

names(res_lfc) <- rownames(res)

# 9.2 Eliminar los NA

gene_list <- na.omit(res_lfc)


# 9.3 Ordena en orden descendente

gene_list <- sort(gene_list, decreasing = TRUE)


# 9.4 Objeto GSEA de Gene ontology

gse <- gseGO(
  geneList = gene_list, 
  ont = "ALL", 
  keyType = "SYMBOL", 
  pvalueCutoff = 0.01, 
  OrgDb = "org.Hs.eg.db",
  minGSSize = 30, 
  maxGSSize = 100 
)
  
# 9.4.1 Observar como data frame

gse_df <- as.data.frame(gse)



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


# Create objects for heatmaps ---------------------------------------------

# if else conditional logic where if in the contrast.matrix object generated on diff_expr_val_sets.R the contrast is between cluster i and j 
# then that section of the script will run
# For example, since the contrast matrix looks like this in c3 vs c1 
#> contrast.matrix
#> cluster 1          -1
#> cluster 2           0
#> cluster 3           1
#> Then the sum of (cluster 1)^2 and (clsuter 3)^2 = 2 meanwhile the sum of (cluster 2)^2 with any of the other 2 will yield 1 so that means 
#> cluster 2 is not taken into account in this round

if((contrast.matrix[2,])^2 + (contrast.matrix[3,])^2 == 2){
  c3_vs_c2_GO <- data.frame(row.names = gse$Description, # GO cluster 3 vs cluster 2
                            c3_vs_c2  = gse$NES)
  
  
  c3_vs_c2_hallmark <- data.frame(row.names = gsea_res$Description, # hallmarks cluster 3 vs cluster 2
                                  c3_vs_c2  = gsea_res$NES)
  
}else if((contrast.matrix[2,])^2 + (contrast.matrix[1,])^2 == 2){
  c1_vs_c2_GO <- data.frame(row.names = gse$Description,
                            c1_vs_c2  = gse$NES)
  
  c1_vs_c2_hallmark <- data.frame(row.names = gsea_res$Description,
                                  c1_vs_c2  = gsea_res$NES)
  
}else if((contrast.matrix[3,])^2 + (contrast.matrix[1,])^2 == 2){

c3_vs_c1_GO <- data.frame(row.names = gse$Description,
                          c3_vs_c1  = gse$NES)


c3_vs_c1_hallmark <- data.frame(row.names = gsea_res$Description,
                                c3_vs_c1  = gsea_res$NES)
}
