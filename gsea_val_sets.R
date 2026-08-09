library(EnhancedVolcano)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggridges)
library(enrichplot)
library(msigdbr)


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
  pvalueCutoff = 0.05, 
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
  pvalueCutoff = 0.05,
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
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose      = FALSE
)

gsea_df <- as.data.frame(gsea_res)
View(gsea_df)
