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

gsea_go_df
gsea_go_hm
