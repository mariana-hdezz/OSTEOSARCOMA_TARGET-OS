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



#---------- GENE SET ENRICHMENT ------------

gse <- gseGO(
  geneList = gene_list,
  ont = "All", # One of "BP", "MF, and "CC subontologies or "ALL"
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  minGSSize = 50, # Minimum number of genes in set (gene sets with lower than this many genes in your dataset will be ignored).
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,  # Print message or not
  seed = TRUE,
  nPermSimple   = 10000,
  eps = 0
  )


# nPerm = 1000, # The higher the no. of permutations, the more accurate the result, but the longer the analysis will take

# ------------- CREATE DOTPLOT -----------
require(DOSE)


dotplot(gse, showCategory = 15,  split = ".sign") + 
  facet_grid(.~.sign) + # apply "facet_grid" to the plot by specified variable (e.g., "ONTOLOGY", "category" and "intersect")
  theme(
    axis.text.x = element_text(angle = 10, hjust = 1),
    axis.text.y = element_text(size = 7)
  )

# --------- CREATE CNETPLOT -----------
  # depicts the linkages of genes and biological concepts (e.g. GO terms or KEGG pathways) as a network (helpful 
  # to see which genes are involved in enriched pathways and genes that may belong to multiple annotation categories).

cnetplot(gse, foldChange = gene_list, showCategory = 3)



# --------- CREATE HEATMAP -----------

p1 <- heatplot(gse, showCategory = 3)
p2 <- heatplot(gse, foldChange = gene_list, showCategory = 5)

cowplot::plot_grid(p1, p2, ncol = 1, labels = LETTERS[1:2])



# -------- RIDGEPLOT ------------
# Helpful to interpret up/down-regulated pathways.

library(ggridges)

ridgeplot(gse) +
  labs(x = "Enrichment Distribution") +
  theme(axis.text.y = element_text(size = 8))


#--------- CREATE GSEAPLOT -----------

# For a gene set as the analysis walks down the ranked gene list, 
# including the location of the maximum enrichment score (the red line).

# The black lines in the Running Enrichment Score show where the members 
# of the gene set appear in the ranked list of genes, indicating the leading edge subset.

# The Ranked list metric shows the value of the ranking metric (log2 fold change) as you
# move down the list of ranked genes. The ranking metric measures a gene’s correlation with a phenotype.


gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)


gseaplot2(gse, title = gse$Description[1], geneSetID = 1)


gseaplot2(gse, geneSetID = 1:3,
          color = c("lightpink", "lightgreen", "lightblue"))



#--------- CREATE PMCPLOT -----

library(europepmc)

terms <- gse$Description[1:3]
pmcplot(terms, 2010:2024, proportion = FALSE)

