setwd("~/Documents/OSTEOSARCOMA/R.project/Hueso/Networks/WGCNA_OS")

library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(igraph)

options(stringsAsFactors = FALSE)

allowWGCNAThreads() # Allow multi-threading 


# WGCNA requires:
# rows = samples
# columns = genes

#----- Transposition needed for WGCNA -------
datExpr <- t(vst_counts)

dim(datExpr)

#------- Quality check -------
good_sample_genes <- goodSamplesGenes(datExpr, verbose = 3) 
good_sample_genes$allOK
## Note: it founded 148 genes without useful info (missing samples or zero variance)

# Eliminate unuseful samples
datExpr <- datExpr[good_sample_genes$goodSamples, good_sample_genes$goodGenes]
dim(datExpr) # look how many samples we got left


#---- Choose soft thresholding power ---
powers <- c(1:10, seq(from = 12, to = 30, by = 2)) # Create a list of values to test 
# It tests powers from 1 to 10 and then pairs from 12 to 30 to find the best one

soft_threshold <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  verbose = 5
)

soft_threshold$fitIndices

soft_power <- soft_threshold$powerEstimate


# ---- Construct network and identify modules ----

set.seed(123)

os_network <- blockwiseModules(
  datExpr,
  power = soft_power,
  networkType = "signed",
  maxBlockSize = 5000,
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  deepSplit = 2,
  numericLabels = FALSE,
  pamRespectsDendro = FALSE,
  verbose = 3
)




# Extract colors from the modules
module_colors <- os_network$colors

table(module_colors)

length(unique(module_colors)) # See how many modules were found

# Total analyzed genes: 18,834
# Genes assigned to modules: 13,759 
# Not assigned genes (grey module): 5,075
# Total module/colors: 44 


#------ Plot dendrograms and module colors ------

# Show the 4 blocks together
par(mfrow = c(2, 2), mar = c(5, 4, 4, 2))

for (block_number in seq_along(os_network$dendrograms)) {
  genes_in_block <- os_network$blockGenes[[block_number]]
  WGCNA::plotDendroAndColors(
    os_network$dendrograms[[block_number]],
    cbind(os_network$unmergedColors[genes_in_block], os_network$colors[genes_in_block]),
    c("Unmerged", "Merged"),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = paste("TARGET-OS gene dendrogram - Block", block_number)
  )
}


# Select module eigengene

module_eigengenes <- orderMEs(os_network$MEs)
dim(module_eigengenes)
colnames(module_eigengenes)


# Remove grey module 
module_eigengenes <- module_eigengenes %>%
  dplyr::select(-MEgrey)

dim(module_eigengenes)

identical(rownames(module_eigengenes), rownames(datExpr)) # Confirm that order is preserved



# Samples present in both objects
common_samples <- intersect(rownames(module_eigengenes), metadata_os$sample)

length(common_samples) # There's a patient with missing metadata

missing_metadata <- setdiff(rownames(module_eigengenes), metadata_os$sample)

# Keep the 87 samples with clinical metadata
module_eigengenes_WGCNA <- module_eigengenes[common_samples, , drop = FALSE]

# Align metadata in exactly the same order
metadata_WGCNA <- metadata_os %>%
  filter(sample %in% common_samples) %>%
  arrange(match(sample, common_samples)) %>%
  as.data.frame()

rownames(metadata_WGCNA) <- metadata_WGCNA$sample

dim(module_eigengenes_WGCNA)

dim(metadata_WGCNA)

identical(rownames(module_eigengenes_WGCNA), rownames(metadata_WGCNA))

# Clinical variables for WGCNA

metadata_analysis <- metadata_WGCNA %>%
  dplyr::select(
    metastasis_at_diagnosis,
    survival_stat,
    relapse_stat,
    survival_time,
    time_to_first_event
  ) %>%
  as.data.frame()

rownames(metadata_analysis) <- rownames(metadata_WGCNA)

dim(module_eigengenes_WGCNA)
dim(metadata_analysis)

identical(rownames(module_eigengenes_WGCNA),
          rownames(metadata_analysis))

summary(metadata_analysis)

#------ Clinical correlation --------

traits_correlation <- metadata_analysis %>%
  dplyr::select(metastasis_at_diagnosis, survival_stat, relapse_stat) # leave out time


#----- Module-trait heatmap -----

heatmap_data <- merge(as.data.frame(module_eigengenes_WGCNA),
                      as.data.frame(traits_correlation),
                      by = "row.names")

heatmap_data <- heatmap_data %>%
  tibble::column_to_rownames("Row.names")

CorLevelPlot(
  heatmap_data,
  x = colnames(traits_correlation),
  y = colnames(module_eigengenes_WGCNA),
  col = c("blue1", "skyblue", "white", "pink", "red")
)



module_trait_test <- WGCNA::corAndPvalue(
  x = as.matrix(module_eigengenes_WGCNA),
  y = as.matrix(traits_correlation),
  use = "pairwise.complete.obs"
)

module_trait_cor <- module_trait_test$cor
module_trait_p <- module_trait_test$p
module_trait_n <- module_trait_test$nObs


# Correct p-values for multiple comparisons
module_trait_FDR <- apply(module_trait_p,
                          MARGIN = 2,
                          FUN = p.adjust,
                          method = "BH")

# Keep the same row and column names
dimnames(module_trait_FDR) <- dimnames(module_trait_p)


# Create data frame
module_trait_results <- data.frame(
  module = rep(rownames(module_trait_cor), times = ncol(module_trait_cor)),
  
  trait = rep(colnames(module_trait_cor), each = nrow(module_trait_cor)),
  
  correlation = as.vector(module_trait_cor),
  p_value = as.vector(module_trait_p),
  FDR = as.vector(module_trait_FDR),
  n = as.vector(module_trait_n)
) %>%
  arrange(FDR)

head(module_trait_results, 20)

module_trait_results %>%
  filter(FDR < 0.05)

# 5. Extract genes from the significant red module

red_module_genes <- names(module_colors)[module_colors == "red"]

length(red_module_genes)

all(red_module_genes %in% colnames(datExpr))

# Expression of the 642 genes from the red module
red_expression <- datExpr[, red_module_genes, drop = FALSE]

# Red module eigengene, preserving sample order
red_eigengene <- module_eigengenes[rownames(red_expression), "MEred", drop = FALSE]

# Correlation of each gene with the red eigengene
red_kME_test <- WGCNA::corAndPvalue(
  x = as.matrix(red_expression),
  y = as.matrix(red_eigengene),
  use = "pairwise.complete.obs"
)

red_kME_results <- data.frame(
  gene = rownames(red_kME_test$cor),
  kME_red = as.numeric(red_kME_test$cor[, 1]),
  p_kME_red = as.numeric(red_kME_test$p[, 1]),
  stringsAsFactors = FALSE
) %>%
  dplyr::arrange(desc(kME_red))

head(red_kME_results, 20) # Genes that represent the general behavior of the red module


#------ Construct the red-module co-expression network ------

# Calculate adjacency matrix
red_gene_adjacency <- WGCNA::adjacency(red_expression, power = soft_power, type = "signed")

dim(red_gene_adjacency)

# Calculate Topological Overlap Matrix 
## Note: it measures genes relationship considering their direct relationship and the "neighbors" they share
red_gene_TOM <- WGCNA::TOMsimilarity(red_gene_adjacency, TOMType = "signed")

dimnames(red_gene_TOM) <- list(colnames(red_expression), colnames(red_expression)) # Add gene names

dim(red_gene_TOM)
red_gene_TOM[1:5, 1:5]


#------ Calculate red-module gene connectivity ------

red_gene_connectivity <- rowSums(red_gene_TOM) - 1

red_connectivity_results <- data.frame(
  gene = names(red_gene_connectivity),
  TOM_connectivity = as.numeric(red_gene_connectivity),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(TOM_connectivity))

head(red_connectivity_results, 20)

#------ Combine kME and TOM connectivity ------

red_hub_results <- red_kME_results %>%
  left_join(red_connectivity_results, by = "gene") %>%
  arrange(desc(kME_red), desc(TOM_connectivity))

head(red_hub_results, 20)

# TOM_connectivity: intramodular connectivity
# kME_red: how much does a gene represents the red module


#------ Transform TOM into an igraph network ------

red_gene_graph_all <- igraph::graph_from_adjacency_matrix(red_gene_TOM,
                                                          mode = "undirected",
                                                          weighted = TRUE,
                                                          diag = FALSE)

igraph::vcount(red_gene_graph_all)
igraph::ecount(red_gene_graph_all)


#------ Distribution of TOM edge weights ------

tom_weight_quantiles <- quantile(igraph::E(red_gene_graph_all)$weight,
                                 probs = c(0, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1))

tom_weight_quantiles

# Note:
# Median;value below which 50% of the edges fall.
# 90th percentile retains approximately the strongest 10% of edges
# 95th percentile retains approximately the strongest 5% of edges
# 99th percentile retains approximately the strongest 1% of edges
# Edge filtering is used only to obtain an interpretable visualization
# kME and TOM connectivity were calculated using the complete network


#------ Filter the network for visualization ------

# Use the 99th percentile as the visualization threshold
minimum_edge_weight <- unname(tom_weight_quantiles["99%"])

minimum_edge_weight


# Make copy of the complete graph
red_gene_graph <- red_gene_graph_all


# Remove edges below selected threshold
red_gene_graph <- igraph::delete_edges(red_gene_graph, igraph::E(red_gene_graph)[igraph::E(red_gene_graph)$weight < minimum_edge_weight])


# Remove genes without remaining strong connections
red_gene_graph <- igraph::delete_vertices(red_gene_graph, igraph::V(red_gene_graph)[igraph::degree(red_gene_graph) == 0])

# Check filtered network
igraph::vcount(red_gene_graph)
# [1] 354
igraph::ecount(red_gene_graph)
# [1] 2058


#------ Visualize  filtered red-module network ------

set.seed(123)

red_network_layout <- igraph::layout_with_fr(red_gene_graph)

plot(
  red_gene_graph,
  layout = red_network_layout,
  vertex.label = igraph::V(red_gene_graph)$name,
  vertex.size = 5,
  vertex.label.cex = 0.45,
  edge.width = igraph::E(red_gene_graph)$weight * 5,
  main = "Red-module co-expression network"
)

#------ Identify highly connected genes in the filtered graph ------

# Sum of strong connection's weight of each gene
red_gene_strength <- igraph::strength(red_gene_graph,
                                      mode = "all",
                                      weights = igraph::E(red_gene_graph)$weight)
 
red_graph_hub_results <- data.frame(
  gene = names(red_gene_strength),
  graph_strength = as.numeric(red_gene_strength),
  stringsAsFactors = FALSE
) %>%
  dplyr::arrange(desc(graph_strength))

head(red_graph_hub_results, 20)

#------ Combine all hub-gene measures ------

red_hub_final <- red_hub_results %>%
  dplyr::left_join(red_graph_hub_results, by = "gene") %>%
  dplyr::arrange(desc(kME_red), desc(TOM_connectivity))

head(red_hub_final, 20)









