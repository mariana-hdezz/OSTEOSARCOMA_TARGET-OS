#############################################################################
#> Script to perform differential expression analysis on TARGET-OS patients 
#> comparing the clusters established on 3_clustering_assignment/3_1_cluster_train/3_1_hierarchical_clustering.R
#> 
#> Inputs: metadatao_os, counts_data
#> 
#> Outputs: No outputs used in further analysis
#> 
#> Results: 
##> Objects with GSEA results:
###> gsea_df_GO3.vs.1
###> gsea_df_GO1.vs.2
###> gsea_df_GO3.vs.2
###> gsea_df_GO3.vs.1
###> gsea_df_HM3.vs.1
###> gsea_df_HM3.vs.2
###> gsea_df_HM1.vs.2
#
##> Objects with results of diffexp
###> res_1_vs_2
###> res_3_vs_2
###> res_3_vs_1
#
##> Diffexp objects filtered by p value
###> res_sig_1v2
###> res_sig_3v2
###> res_sig_3v1
#############################################################################


library(DESeq2)
library(apeglm)
library(tidyverse)
library(ggiraph)
library(clusterProfiler)
library(enrichplot)
library(igraph)
library(ggplot2)
library(tidyr)
library(msigdbr)
library(ggtree)
library(aplot)
library(AnnotationDbi)


# Load data

metadata_os <- readRDS("./output_data/metadata_os.RDS")
counts_data <- readRDS("./output_data/counts_data.RDS")
counts_data <- counts_data[-1, ]

# Create output directory


if(dir.exists("./results/diffex_gsea_target/")){
  "differential expr4ession and gsea directory already exists"
}else{
  dir.create("./results/diffex_gsea_target/")
} 


# Add clusters column to each sample of the metadata

metadata_os_clusters <- metadata_os %>% 
  dplyr::select(sample, clusters) %>% 
  mutate(clusters = factor(clusters))


keep_expression <- rowSums(counts_data >= 10) >= ceiling(0.25 * ncol(counts_data)) # genes with more than 1 count on more than 10% of patients


counts_data_dseq <- counts_data[rownames(counts_data) %in% names(keep_expression)[keep_expression == TRUE], colnames(counts_data) %in% metadata_os_clusters$sample]

all(colnames(counts_data_dseq) == metadata_os$sample) # check if samples in counts and metadata are in the same order

counts_data_dseq <- counts_data_dseq[ , metadata_os$sample]


#Pt in columns

dds <- DESeqDataSetFromMatrix(countData = counts_data_dseq,
                              colData = metadata_os_clusters,
                              design= ~ clusters)

rm(list = setdiff(ls(), "dds"))

results_list_go <- list()
results_list_hm <- list()

for (i in 1:(length(unique(dds$clusters)) - 1)) {
  
  dds$clusters <- factor(dds$clusters)
  
  dds$clusters <- relevel(dds$clusters, ref = as.character(i))
  
  dds <- DESeq(dds)
  
  resultsNames(dds) # lists the coefficients
  
  for (e in 1:i) {
    
    
    if(i == 1){
      res_3_vs_1 <- results(dds, name = "clusters_3_vs_1", lfcThreshold = 0.58)
      
      res_sig_3v1 <- subset(res_3_vs_1, res_3_vs_1$padj < 0.05) %>% 
        as.data.frame()
      
      # or to shrink log fold changes association with condition:
      
      res <- lfcShrink(dds, coef = "clusters_3_vs_1", type = "apeglm")
      
    }else if(i == 2){
      
      if(e == 1){
        
        res_1_vs_2 <- results(dds, name = "clusters_1_vs_2", lfcThreshold = 0.58)
        
        res_sig_1v2 <- subset(res_1_vs_2, res_1_vs_2$padj < 0.05) %>% 
          as.data.frame()
        
        res <- lfcShrink(dds, coef = "clusters_1_vs_2", type = "apeglm")
        
      }else if(e == 2){
        
        
        res_3_vs_2 <- results(dds, name = "clusters_3_vs_2", lfcThreshold = 0.58)
        
        res_sig_3v2 <- subset(res_3_vs_2, res_3_vs_2$padj < 0.05) %>% 
          as.data.frame()
        
        res <- lfcShrink(dds, coef = "clusters_3_vs_2", type = "apeglm")
        
      }
      
    }
    
    
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
      pvalueCutoff = 0.01,
      pAdjustMethod = "BH",
      verbose = TRUE,  # Print message or not
      seed = TRUE,
      nPermSimple = 10000,
      eps = 0
    )
    
    gsea_go_df <- as.data.frame(gse)
    
    name_go <- paste0("gsea_df_GO", gsub("^X", "", make.names(gsub("^Wald test p-value: clusters (.*)", "\\1", mcols(res)$description[4]))))
    
    assign(name_go, gsea_go_df)
    
    results_list_go[[paste(i, e)]] <- name_go
    
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
    
    gsea_hm_df <- as.data.frame(gsea_res)
    
    name_hm <- paste0("gsea_df_HM", gsub("^X", "", make.names(gsub("^Wald test p-value: clusters (.*)", "\\1", mcols(res)$description[4]))))
    
    assign(name_hm, gsea_hm_df)
    
    results_list_hm[[paste(i, e)]] <- name_hm
    
  }
}

# Save results


write.csv(gsea_df_GO3.vs.1, "results/diffex_gsea_target/gsea_df_GO3_vs_1.csv")
write.csv(gsea_df_GO3.vs.2, "results/diffex_gsea_target/gsea_df_GO3_vs_2.csv")
write.csv(gsea_df_GO1.vs.2, "results/diffex_gsea_target/gsea_df_GO1_vs_2.csv")
write.csv(gsea_df_GO3.vs.1, "results/diffex_gsea_target/gsea_df_GO3_vs_1.csv")
write.csv(gsea_df_HM3.vs.1, "results/diffex_gsea_target/gsea_df_HM3_vs_1.csv")
write.csv(gsea_df_HM3.vs.2, "results/diffex_gsea_target/gsea_df_HM3_vs_2.csv")
write.csv(gsea_df_HM1.vs.2, "results/diffex_gsea_target/gsea_df_HM1_vs_2.csv")

write.csv(res_sig_1v2, "results/diffex_gsea_target/res_sig_1v2.csv")
write.csv(res_sig_3v2, "results/diffex_gsea_target/res_sig_3v2.csv")
write.csv(res_sig_3v1, "results/diffex_gsea_target/res_sig_3v1.csv")

write.csv(res_1_vs_2 , "results/diffex_gsea_target/res_1_vs_2.csv")
write.csv(res_3_vs_2 , "results/diffex_gsea_target/res_3_vs_2.csv")
write.csv(res_3_vs_1 , "results/diffex_gsea_target/res_3_vs_1.csv")



rm(list = ls())
gc()
