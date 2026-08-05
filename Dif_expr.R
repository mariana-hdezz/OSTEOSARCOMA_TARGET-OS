library(dplyr)
library(DESeq2)
library(apeglm)
library(tidyverse)
library(ggiraph)
library(clusterProfiler)
library(enrichplot)

# Add clusters column to each sample of the metadata
metadata_os_clusters <- metadata_os %>% 
  dplyr::select(sample, clusters) %>% 
  mutate(clusters = factor(clusters))


keep_expression <- rowSums(counts_data > 10) >= ceiling(0.10 * ncol(counts_data)) # genes with more than 1 count on more than 10% of patients


counts_data_dseq <- counts_data[rownames(counts_data) %in% names(keep_expression)[keep_expression == TRUE], colnames(counts_data) %in% metadata_os_clusters$sample]


#Pt in columns
#-------------------- DIFFERENTIAL EXPRESSION DESeq2 --------------------------

dds <- DESeqDataSetFromMatrix(countData = counts_data_dseq[-1, ],
                              colData = metadata_os_clusters,
                              design= ~ clusters)

dds$clusters <- relevel(dds$clusters, ref = "2")

dds <- DESeq(dds)

resultsNames(dds) # lists the coefficients

# res_1_vs_2 <- results(dds, name = "clusters_1_vs_2")
# 
# res_3_vs_2 <- results(dds, name = "clusters_3_vs_2")


# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef = "clusters_1_vs_2", type = "apeglm")
