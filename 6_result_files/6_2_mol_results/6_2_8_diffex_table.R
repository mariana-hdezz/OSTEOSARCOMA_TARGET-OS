
#############################################################################
#> Script to obtain the heatmaps for the gsea of GSE analysis 
#> 
#> Inputs: res_sig_1v2, res_sig_3v1, res_sig_3v2, res_sig_c3_vs_c2_gse33382, res_sig_c3_vs_c1_gse33382
#> res_sig_c1_vs_c2_gse33382, res_sig_c3_vs_c2_gse_21257, res_sig_c3_vs_c1_gse_21257
#> res_sig_c1_vs_c2_gse_21257
#> 
#> 
#> Results: 
##> Table with differentially expressed genes in each set
##> 
##>  
#############################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(flextable)

res_sig_1v2 <- read.csv("results/diffex_gsea_target/res_sig_1v2.csv")
res_sig_3v1 <- read.csv("results/diffex_gsea_target/res_sig_3v1.csv")
res_sig_3v2 <- read.csv("results/diffex_gsea_target/res_sig_3v2.csv")

res_sig_c3_vs_c2_gse33382 <- read.csv("results/diffex_gsea_gse/gse33382_res_sig_3v2_.csv") %>% 
  column_to_rownames("X")
res_sig_c3_vs_c1_gse33382 <- read.csv("results/diffex_gsea_gse/gse33382_res_sig_3v1_.csv") %>% 
  column_to_rownames("X")
res_sig_c1_vs_c2_gse33382 <- read.csv("results/diffex_gsea_gse/gse33382_res_sig_1v2_.csv") %>% 
  column_to_rownames("X")

res_sig_c3_vs_c2_gse_21257 <- read.csv("results/diffex_gsea_gse/gse21257_res_sig_3v2_.csv")
res_sig_c3_vs_c1_gse_21257 <- read.csv("results/diffex_gsea_gse/gse21257_res_sig_3v1_.csv")
res_sig_c1_vs_c2_gse_21257 <- read.csv("results/diffex_gsea_gse/gse21257_res_sig_1v2_.csv")

common_3v2 <- rownames(res_sig_c3_vs_c2_gse33382)[rownames(res_sig_c3_vs_c2_gse33382) %in% res_sig_c3_vs_c2_gse_21257$X]
common_3v1 <- rownames(res_sig_c3_vs_c1_gse33382)[rownames(res_sig_c3_vs_c1_gse33382) %in% res_sig_c3_vs_c1_gse_21257$X]
common_1v2 <- rownames(res_sig_c1_vs_c2_gse33382)[rownames(res_sig_c1_vs_c2_gse33382) %in% res_sig_c1_vs_c2_gse_21257$X]


# For each result keep only top 5 and botom 5 genes, if log fold chang eis positive oit is upregulated.


t1 <- (res_sig_1v2 %>%
    arrange(desc(log2FoldChange))) %>%
  slice(1:5, (n() - 4):n()) %>%
  mutate(
    Cluster = ifelse(log2FoldChange > 0, yes = "Upregulated", no = "Downregulated"),
    Comparison = "Cluster 1 vs Cluster 2",
    database = "TARGET-OS",
    Reference_cluster = "Cluster 1"
         ) %>% 
  dplyr::select(Cluster, X, database, Reference_cluster,
                Comparison)


t2 <- res_sig_3v2  %>% 
  arrange(desc(log2FoldChange)) %>% 
  slice(1:5, (n() - 4):n()) %>%
  mutate(
    Cluster = ifelse(log2FoldChange > 0, yes = "Upregulated", no = "Downregulated"),
    Comparison = "Cluster 3 vs Cluster 2" ,
    database = "TARGET-OS"    ,
    Reference_cluster = "Cluster 3"
  ) %>% 
  dplyr::select(Cluster, X, database, Reference_cluster,
                Comparison)

t3 <- res_sig_3v1  %>% 
  arrange(desc(log2FoldChange)) %>% 
  slice(1:5, (n() - 4):n()) %>%
  mutate(
    Cluster = ifelse(log2FoldChange > 0, yes = "Upregulated", no = "Downregulated"),
    Comparison = "Cluster 3 vs Cluster 1" ,
    database = "TARGET-OS"    ,
    Reference_cluster = "Cluster 3"
  ) %>% 
  dplyr::select(Cluster, X, database, Reference_cluster,
                Comparison)


t4 <- res_sig_c1_vs_c2_gse33382[common_1v2, ] %>% 
  arrange(desc(logFC) )%>%
  mutate(
    Cluster = ifelse(logFC > 0, yes = "Upregulated", no = "Downregulated"),
    Comparison = "Cluster 1 vs Cluster 2"   ,
    database = "Validation set"  ,
    Reference_cluster = "Cluster 1"
  ) %>% 
  rownames_to_column("X") %>% 
  dplyr::select(Cluster, X, database, Reference_cluster,
                Comparison)

t5 <- res_sig_c3_vs_c2_gse33382[common_3v2, ] %>% 
  arrange(desc(logFC))  %>% 
  slice(1:5) %>%
  mutate(
    Cluster = ifelse(logFC > 0, yes = "Upregulated", no = "Downregulated"),
    Comparison = "Cluster 3 vs Cluster 2"    ,
    database = "Validation set"    ,
    Reference_cluster = "Cluster 3"
  ) %>% 
  rownames_to_column("X") %>% 
  dplyr::select(Cluster, X, database, Reference_cluster,
                Comparison)

t6 <- res_sig_c3_vs_c1_gse33382[common_3v1, ] %>% 
  arrange(desc(logFC)) %>% 
  slice(1:5, (n() - 4):n()) %>%
  mutate(
    Cluster = ifelse(logFC > 0, yes = "Upregulated", no = "Downregulated"),
    Comparison = "Cluster 3 vs Cluster 1"     ,
    database = "Validation set"   ,
    Reference_cluster = "Cluster 3"
  ) %>% 
  rownames_to_column("X") %>% 
  dplyr::select(Cluster, X, database, Reference_cluster,
                Comparison)

binded <- bind_rows(t1, t2, t3, t4, t5, t6)
binded <- binded[, c("Comparison", "Reference_cluster", "Cluster", "database", "X")] %>% 
  dplyr::rename("Regulation" = "Cluster",
                "Reference cluster" = "Reference_cluster"
                )


binded %>% 
  tidyr::pivot_wider(
    names_from = "database",
    values_from = "X",
    values_fn = ~ paste(unique(.x), collapse = ", ")
  ) %>% 
  flextable() %>% 
  autofit() 

