library(DESeq2)
library(apeglm)
library(tidyverse)
library(ggiraph)
library(clusterProfiler)
library(enrichplot)
library(igraph)
library(ggplot2)
library(tidyr)

# Add clusters column to each sample of the metadata
metadata_os_clusters <- metadata_os %>% 
  dplyr::select(sample, clusters) %>% 
  mutate(clusters = factor(clusters))


keep_expression <- rowSums(counts_data >= 10) >= ceiling(0.25 * ncol(counts_data)) # genes with more than 1 count on more than 10% of patients


counts_data_dseq <- counts_data[rownames(counts_data) %in% names(keep_expression)[keep_expression == TRUE], colnames(counts_data) %in% metadata_os_clusters$sample]

all(colnames(counts_data_dseq) == metadata_os$sample) # check if samples in counts and metadata are in the same order

counts_data_dseq <- counts_data_dseq[ , metadata_os$sample]


#Pt in columns
#-------------------- DIFFERENTIAL EXPRESSION DESeq2 --------------------------

dds <- DESeqDataSetFromMatrix(countData = counts_data_dseq,
                              colData = metadata_os_clusters,
                              design= ~ clusters)

dds$clusters <- relevel(dds$clusters, ref = "3")

dds <- DESeq(dds)

resultsNames(dds) # lists the coefficients

res_1_vs_2 <- results(dds, name = "clusters_1_vs_2", lfcThreshold = 0.58)

res_sig_1v2 <- subset(res_1_vs_2, res_1_vs_2$padj < 0.05) %>% 
  as.data.frame()

res_3_vs_2 <- results(dds, name = "clusters_3_vs_2", lfcThreshold = 0.58)

res_sig_3v2 <- subset(res_3_vs_2, res_3_vs_2$padj < 0.05) %>% 
  as.data.frame()

res_1_vs_3 <- results(dds, name = "clusters_1_vs_3", lfcThreshold = 0.58)

res_sig_1v3 <- subset(res_1_vs_3, res_1_vs_3$padj < 0.05) %>% 
  as.data.frame()

# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef = "clusters_1_vs_3", type = "apeglm")



#------- PREPARE INPUT DATA ---------

# Read results from DESeq2 analysis

df <- res

# Extract log2foldchange (column we are interested in)
original_gene_list <- df$log2FoldChange
names(original_gene_list) <- rownames(df)

sum(is.na(original_gene_list))

# Omit the missing values if present
gene_list <- na.omit(original_gene_list)


# Sorting in a descending order (required for clusterProfiler)
gene_list <- sort(gene_list, decreasing = TRUE)

# Set the organism 
library(org.Hs.eg.db)


#---------- GENE SET ENRICHMENT ------------

gse <- gseGO(
  geneList = gene_list,
  ont = "ALL", # One of "BP", "MF, and "CC subontologies or "ALL"
  OrgDb = "org.Hs.eg.db",
  keyType = "SYMBOL",
  minGSSize = 30, # Minimum number of genes in set (gene sets with lower than this many genes in your dataset will be ignored).
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,  # Print message or not
  seed = TRUE,
  nPermSimple = 10000,
  eps = 0
  )

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
  geneList = genes_values, # Objeto creado anteriormente que ordena de manera descendente
  keyType = "kegg", # El codigo de identificación de los genes
  pvalueCutoff = 0.05, 
  organism = "hsa", # Organismo en este caso Homo sapiens (a diferencia de org.Mm.eg.db que es de raton o org.Sc.eg.db que es de levadura) 
  minGSSize = 30, # Tamaño minimo de ada set para analisis
  maxGSSize = 500, # Tamaño maximo de genes anotados para analizar
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
