library(EnhancedVolcano)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggridges)
library(enrichplot)
library(msigdbr)
library(ComplexHeatmap)
library(circilize)


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


fit$genes <- data.frame(Gene.symbol = gene_symbols)

toptable <- topTable(fit, n = Inf)

EnhancedVolcano(toptable,
                lab = toptable$Gene.symbol,
                x = 'logFC',
                y = 'P.Value')




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

View(gse_df)




# Convert Ensembl to Entrez (KEGG format)

entrez_df <- bitr(names(gene_list), 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db)

# Keep available genes

gene_list_kegg <- gene_list[names(gene_list) %in% entrez_df$SYMBOL]

# Keep unique entrez

entrez_df <- entrez_df[entrez_df$SYMBOL %in% names(gene_list_kegg), ] %>% 
  group_by(SYMBOL) %>% 
  slice_head()

# Join entrez with value

x <- as.data.frame(gene_list_kegg) %>% 
  rownames_to_column("SYMBOL") %>% 
  filter(SYMBOL %in% entrez_df$SYMBOL) %>% 
  left_join(entrez_df, by = "SYMBOL")

genes_values <- x$gene_list_kegg

names(genes_values) <- make.names(x$ENTREZID)

names(genes_values) <- gsub("^X", "", names(genes_values))

gse_kegg <- gseKEGG(
  geneList = genes_values,
  keyType = "kegg", 
  pvalueCutoff = 0.01,
  organism = "hsa",  
  minGSSize = 30, 
  maxGSSize = 500,
  BPPARAM = BiocParallel::SerialParam()
)


View(as.data.frame(gse_kegg))

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
View(gsea_df)


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
  
  
  c3_vs_c2_KEGG <- data.frame(row.names = gse_kegg$Description, # KEGG cluster 3 vs cluster 2
                              c3_vs_c2  = gse_kegg$NES)
  
  
  c3_vs_c2_hallmark <- data.frame(row.names = gsea_res$Description, # hallmarks cluster 3 vs cluster 2
                                  c3_vs_c2  = gsea_res$NES)
  
}else if((contrast.matrix[2,])^2 + (contrast.matrix[1,])^2 == 2){
  c1_vs_c2_GO <- data.frame(row.names = gse$Description,
                            c1_vs_c2  = gse$NES)
  
  
  c1_vs_c2_KEGG <- data.frame(row.names = gse_kegg$Description,
                              c1_vs_c2  = gse_kegg$NES)
  
  
  c1_vs_c2_hallmark <- data.frame(row.names = gsea_res$Description,
                                  c1_vs_c2  = gsea_res$NES)
  
}else if((contrast.matrix[3,])^2 + (contrast.matrix[1,])^2 == 2){

c3_vs_c1_GO <- data.frame(row.names = gse$Description,
                          c3_vs_c1  = gse$NES)


c3_vs_c1_KEGG <- data.frame(row.names = gse_kegg$Description,
                            c3_vs_c1  = gse_kegg$NES)


c3_vs_c1_hallmark <- data.frame(row.names = gsea_res$Description,
                                c3_vs_c1  = gsea_res$NES)
}

# Keep only the top


c3_vs_c2_GO_10 <- c3_vs_c2_GO %>% 
  filter(c3_vs_c2 > quantile(c3_vs_c2, 0.95) | c3_vs_c2 < quantile(c3_vs_c2, 0.05))

c1_vs_c2_GO_10 <- c1_vs_c2_GO %>% 
  filter(c1_vs_c2_GO > quantile(c1_vs_c2, 0.95) | c1_vs_c2_GO < quantile(c1_vs_c2, 0.05))

c3_vs_c1_GO_10 <- c3_vs_c1_GO %>% 
  filter(c3_vs_c1_GO > quantile(c3_vs_c1, 0.95) | c3_vs_c1_GO < quantile(c3_vs_c1, 0.05))


# Join and convert to matrix

gsea_GO_mat <- merge(c3_vs_c2_GO_10, (merge(c1_vs_c2_GO_10, c3_vs_c1_GO_10, by = 0, all = TRUE) %>% column_to_rownames("Row.names")), by = 0, all = TRUE)

gsea_GO_mat <- gsea_GO_mat %>% 
  column_to_rownames("Row.names") %>% 
  as.matrix()
  
gsea_GO_mat[is.na(gsea_GO_mat)] <- 0

# Plot

col_fun = colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red"))

Heatmap(
  gsea_GO_mat,
  name = "NES",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE, 
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_title = "GSEA Normalized Enrichment Scores Across Contrasts"
)



###############################################################################

c3_vs_c2_KEGG_10 <- c3_vs_c2_KEGG %>% 
  filter(c3_vs_c2 > quantile(c3_vs_c2, 0.95) | c3_vs_c2 < quantile(c3_vs_c2, 0.05))

c1_vs_c2_KEGG_10 <- c1_vs_c2_KEGG %>% 
  filter(c1_vs_c2_KEGG > quantile(c1_vs_c2, 0.95) | c1_vs_c2_KEGG < quantile(c1_vs_c2, 0.05))

c3_vs_c1_KEGG_10 <- c3_vs_c1_KEGG %>% 
  filter(c3_vs_c1_KEGG > quantile(c3_vs_c1, 0.95) | c3_vs_c1_KEGG < quantile(c3_vs_c1, 0.05))


gsea_KEGG_mat <- merge(c3_vs_c2_KEGG_10, (merge(c1_vs_c2_KEGG_10, c3_vs_c1_KEGG_10, by = 0, all = TRUE) %>% column_to_rownames("Row.names")), by = 0, all = TRUE)

gsea_KEGG_mat <- gsea_KEGG_mat %>% 
  column_to_rownames("Row.names") %>% 
  as.matrix()

gsea_KEGG_mat[is.na(gsea_KEGG_mat)] <- 0



col_fun = colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red"))

Heatmap(
  gsea_KEGG_mat,
  name = "NES",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE, 
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_title = "GSEA Normalized Enrichment Scores Across Contrasts"
)


##############################################################################

c3_vs_c2_hallmark_10 <- c3_vs_c2_hallmark %>% 
  filter(c3_vs_c2 > quantile(c3_vs_c2, 0.95) | c3_vs_c2 < quantile(c3_vs_c2, 0.05))

c1_vs_c2_hallmark_10 <- c1_vs_c2_hallmark %>% 
  filter(c1_vs_c2_hallmark > quantile(c1_vs_c2, 0.95) | c1_vs_c2_hallmark < quantile(c1_vs_c2, 0.05))

c3_vs_c1_hallmark_10 <- c3_vs_c1_hallmark %>% 
  filter(c3_vs_c1_hallmark > quantile(c3_vs_c1, 0.95) | c3_vs_c1_hallmark < quantile(c3_vs_c1, 0.05))


gsea_hallmark_mat <- merge(c3_vs_c2_hallmark_10, (merge(c1_vs_c2_hallmark_10, c3_vs_c1_hallmark_10, by = 0, all = TRUE) %>% column_to_rownames("Row.names")), by = 0, all = TRUE)

gsea_hallmark_mat <- gsea_hallmark_mat %>% 
  column_to_rownames("Row.names") %>% 
  as.matrix()

gsea_hallmark_mat[is.na(gsea_hallmark_mat)] <- 0


col_fun = colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red"))

Heatmap(
  gsea_hallmark_mat,
  name = "NES",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE, 
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_title = "GSEA Normalized Enrichment Scores Across Contrasts"
)

