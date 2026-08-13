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

#############################################################################
#> Script to perform differential expression analysis on TARGET-OS patients 
#> comparing the clusters established on 3_Clustering_and_enrichment/3_1_hierarchical_clustering.R
#> 
#> Inputs: metadatao_os, counts_data
#> 
#> Outputs: No outputs used in further analsis
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

# Load data

metadata_os <- readRDS("./output_data/metadata_os.RDS")
counts_data <- readRDS("./output_data/counts_data.RDS")


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

rm(list = setdiff(ls(), "dds"))

results_list_go <- list()
results_list_hm <- list()

for (i in 1:(length(unique(dds$clusters)) - 1)) {
  dds$clusters <- relevel(dds$clusters, ref = i)
  
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

write.csv(res_1_vs_2, "results/diffex_gsea_target/res_1_vs_2.csv")
write.csv(res_3_vs_2, "results/diffex_gsea_target/res_3_vs_2.csv")
write.csv(res_3_vs_1, "results/diffex_gsea_target/res_3_vs_1.csv")



rm(list = ls())


# Load objects for heatmap

gsea_df_GO3_vs_1 <- read_csv("results/diffex_gsea_target/gsea_df_GO3_vs_1.csv")
gsea_df_GO3_vs_2 <- read_csv("results/diffex_gsea_target/gsea_df_GO3_vs_2.csv")
gsea_df_GO1_vs_2 <- read_csv("results/diffex_gsea_target/gsea_df_GO1_vs_2.csv")
gsea_df_GO3_vs_1 <- read_csv("results/diffex_gsea_target/gsea_df_GO3_vs_1.csv")
gsea_df_HM3_vs_1 <- read_csv("results/diffex_gsea_target/gsea_df_HM3_vs_1.csv")
gsea_df_HM3_vs_2 <- read_csv("results/diffex_gsea_target/gsea_df_HM3_vs_2.csv")
gsea_df_HM1_vs_2 <- read_csv("results/diffex_gsea_target/gsea_df_HM1_vs_2.csv")



# Objectts for heatmaps

c1_v_c2_GO <- data.frame(row.names = gsea_df_GO1_vs_2$Description, # Rownames contains tha pathway
                         C1_vs_C2 = gsea_df_GO1_vs_2$NES) # Column named after the comparison contains NES

c3_v_c1_GO <- data.frame(row.names = gsea_df_GO3_vs_1$Description,
                         C3_vs_C1 = gsea_df_GO3_vs_1$NES)

c3_v_c2_GO <- data.frame(row.names = gsea_df_GO3_vs_2$Description,
                         C3_vs_C2 = gsea_df_GO3_vs_2$NES)


c1_v_c2_HM <- data.frame(row.names = gsea_df_HM1_vs_2$Description,
                         C1_vs_C2 = gsea_df_HM1_vs_2$NES)

c3_v_c1_HM <- data.frame(row.names = gsea_df_HM3_vs_1$Description,
                         C3_vs_C1 = gsea_df_HM3_vs_1$NES)

c3_v_c2_HM <- data.frame(row.names = gsea_df_HM3_vs_2$Description,
                         C3_vs_C2 = gsea_df_HM3_vs_2$NES)

# limit how many paths t plot

c1_v_c2_GO <- c1_v_c2_GO %>% 
  filter(C1_vs_C2 > quantile(C1_vs_C2, 0.7) | C1_vs_C2 < quantile(C1_vs_C2, 0.25))

c3_v_c1_GO <- c3_v_c1_GO %>% 
  filter(C3_vs_C1 > quantile(C3_vs_C1, 0.85) | C3_vs_C1 < quantile(C3_vs_C1, 0.15))

c3_v_c2_GO <- c3_v_c2_GO %>% 
  filter(C3_vs_C2 > quantile(C3_vs_C2, 0.75) | C3_vs_C2 < quantile(C3_vs_C2, 0.75))


# Merge the isolated columns to be able tp plot

gsea_GO_heatmap_obj <-  merge(c1_v_c2_GO, c3_v_c1_GO , by = 0, all = TRUE) %>%
  column_to_rownames("Row.names") %>%
  merge(c3_v_c2_GO, by = 0, all = TRUE) %>% 
  column_to_rownames("Row.names") 

#Na to 0

gsea_GO_heatmap_obj[is.na(gsea_GO_heatmap_obj )] <- 0



gsea_HM_heatmap_obj <-  merge(c1_v_c2_HM, c3_v_c1_HM , by = 0, all = TRUE) %>%
  column_to_rownames("Row.names") %>%
  merge(c3_v_c2_HM, by = 0, all = TRUE) %>%
  column_to_rownames("Row.names") 

gsea_HM_heatmap_obj[is.na(gsea_HM_heatmap_obj )] <- 0

# Calculate distances for dendogram

row_go <- hclust(dist(gsea_GO_heatmap_obj))
col_go <- hclust(dist(t(gsea_GO_heatmap_obj)))

# Dndograms

tree_right_go <- ggtree(row_go) + 
  scale_x_reverse() + 
  scale_y_continuous(expand = c(0, 0))

tree_top_go <- ggtree(col_go, hang = -1) + 
  layout_dendrogram() + 
  scale_y_reverse(expand = c(0, 0))

# Pivot longer to be able to plot

obj_for_GO <- gsea_GO_heatmap_obj %>% 
  rownames_to_column("path") %>% 
  pivot_longer(cols = c(C1_vs_C2, C3_vs_C1, C3_vs_C2),
               names_to = "clusters")

# Prepare object

heatmap_go <- obj_for_GO %>%
  ggplot(aes(x = clusters, y = path, fill = value)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Spectral") +
  theme_classic() + 
  scale_x_discrete(labels = c("C3_vs_C1" = "C3 vs C1", "C1_vs_C2" = "C1 vs C2", "C3_vs_C2" = "C3 vs C2")) +
  labs(title = "GSEA between clusters Gene Ontology. TARGET-OS", 
       x = "Clusters")

# Plot

heatmap_go %>% 
  insert_right(tree_right_go, width = 0.1) %>% 
  insert_top(tree_top_go, height = 0.1)

# Repeat for Hallmarks

row_hm <- hclust(dist(gsea_HM_heatmap_obj))
col_hm <- hclust(dist(t(gsea_HM_heatmap_obj)))

tree_right_hm <- ggtree(row_hm) + 
  scale_x_reverse() + 
  scale_y_continuous(expand = c(0, 0))

tree_top_hm <- ggtree(col_hm, hang = -1) + 
  layout_dendrogram() + 
  scale_y_reverse(expand = c(0, 0))


obj_for_hm <- gsea_HM_heatmap_obj %>% 
  rownames_to_column("path") %>% 
  pivot_longer(cols = c(C1_vs_C2, C3_vs_C1, C3_vs_C2),
               names_to = "clusters")

heatmap_hm <- obj_for_hm %>%
  ggplot(aes(x = clusters, y = path, fill = value)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Spectral")+
  theme_classic() +
  scale_x_discrete(labels = c("C3_vs_C1" = "C3 vs C1", "C1_vs_C2" = "C1 vs C2", "C3_vs_C2" = "C3 vs C2")) +
  labs(title = "GSEA between clusters Hallmarks of cancer TARGET-OS", 
       x = "Clusters")

heatmap_hm %>% 
  insert_right(tree_right_hm, width = 0.1) %>% 
  insert_top(tree_top_hm, height = 0.1)


rm(list = ls())
