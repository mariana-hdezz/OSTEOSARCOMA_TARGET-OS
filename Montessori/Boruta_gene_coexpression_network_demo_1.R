# Intento de red con 43 genes de primer boruta
## Falta comentar, tal vez hago enrichment, pero más bien voy a intentar otra cosa con todos los genes de TARGET
## Demo para aprender, pero interesante
#### Hola


library(igraph)
library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot) 
library(gridExtra)

# List of 43 Boruta genes
gene_list <- c(
  "APEX2", "ARHGAP1", "ARHGEF39", "CCDC97", "CGREF1", "CLUAP1",
  "COL22A1", "CPE", "CTNNBIP1", "CYFIP1", "DHRS11", "DLX1",
  "ERCC4", "F13A1", "FAM110D", "FAT1", "FKBP11", "GALNT14",
  "GBP1", "GMIP", "GRAMD1B", "HSD11B2", "INPP4A", "KERA",
  "KIF25", "LGR6", "LURAP1L", "MEF2A", "MRTFB", "MXI1",
  "NUBP1", "SF3B3", "SLC12A4", "SLC45A4", "SLC8A3", "STAT5B",
  "TIMM50", "TPD52", "TRIM68", "TSHZ3", "UBE2D4", "UNC5B", "VMP1"
)


genes_boruta <- intersect(
  gene_list,
  rownames(counts_data)
)


# Check results
length(genes_boruta)

# Extract expression values for the 43 genes
boruta_gene_counts <- counts_data[genes_boruta, , drop = FALSE]

dim(boruta_gene_counts)

boruta_gene_counts[1:5, 1:5]

dim(boruta_gene_counts)

# WGCNA requires:
# rows = samples
# columns = genes
# Transposition needed for the network


datExpr <- boruta_gene_counts %>% 
  as.data.frame() %>% 
  t()


dim(datExpr)

datExpr <- log2(datExpr + 1)

length(genes_boruta)
dim(boruta_gene_counts)

# Quality control
# Check for genes or samples with excessive missing values.
gsg <- goodSamplesGenes(
  datExpr,
  verbose = 3
)

gsg$allOK # We need for it to be TRUE


# 3. Choose soft-thresholding power 

# Test a range of powers to evaluate which one best approximates a scale free co-expression network

powers <- c(1:10, seq(from = 12, to = 30, by = 2))

sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  verbose = 5
)

sft.data <- sft$fitIndices

# Plot scale-free topology fit

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.05) +
  geom_hline(yintercept = 0.80) +
  labs(x = "Power", y = "Scale-free topology model fit, signed R²") +
  theme_classic()


# Plot mean connectivity

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.05) +
  labs(x = "Power", y = "Mean connectivity") +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)

sft$powerEstimate


# Based on the scale-free topology analysis a soft-thresholding power of 10 was selected
## Highest signed R² while maintaining an acceptable mean connectivity


# --------------------- Construct co-expression network ------------------------

# Soft-thresholding power previously chosen
soft_power <- 10


# Calculate adjacency matrix
## Note: Higher values indicate stronger co-expression between genes
gene.adjacency <- adjacency(datExpr, power = soft_power, type = "signed")


# Calculate topological overlap matrix (TOM)
## TOM considers both direct connections and shared neighbors
gene.TOM <- TOMsimilarity(gene.adjacency, TOMType = "signed")


# Add gene names to rows and columns
dimnames(gene.TOM) <- list(colnames(datExpr), colnames(datExpr))


dim(gene.TOM)
gene.TOM[1:5, 1:5]


# Calculate gene connectivity 

# Sum the TOM connections for each gene and take 1 to remove the connection of each gene with itself

gene.connectivity <- rowSums(gene.TOM) - 1


# Create results table
connectivity.results <- data.frame(gene = names(gene.connectivity),
                                   connectivity = as.numeric(gene.connectivity)) %>%
  arrange(desc(connectivity))

head(connectivity.results, 15) # Genes with higher connectivity are potential hub genes

# Visualize network

# Transform TOM matrix into a graph
gene.graph <- graph_from_adjacency_matrix(
  gene.TOM,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)


# Remove weak connections
minimum.edge.weight <- 0.20

gene.graph <- delete_edges(gene.graph, E(gene.graph)[E(gene.graph)$weight < minimum.edge.weight])


# Remove genes without remaining connections
gene.graph <- delete_vertices(gene.graph, degree(gene.graph) == 0)


# Check number of genes and connections remaining
vcount(gene.graph)

ecount(gene.graph)

# Plot
plot(
  gene.graph,
  vertex.label = V(gene.graph)$name,
  vertex.size = 8,
  vertex.label.cex = 0.7,
  edge.width = E(gene.graph)$weight * 3,
  layout = layout_with_fr(gene.graph),
  main = "Boruta gene co-expression network"
)


# 6. Visualize all 43 genes 


# Convert complete TOM matrix into a graph
gene.graph.all <- graph_from_adjacency_matrix(gene.TOM,
                                              mode = "undirected",
                                              weighted = TRUE,
                                              diag = FALSE)

# Use a lower threshold to retain more connections
minimum.edge.weight <- 0.02

gene.graph.all <- delete_edges(gene.graph.all, E(gene.graph.all)[E(gene.graph.all)$weight < minimum.edge.weight])

# Don't delete isolated vertices because we want to display all 43 genes

vcount(gene.graph.all)
ecount(gene.graph.all)

set.seed(123)

plot(
  gene.graph.all,
  vertex.label = V(gene.graph.all)$name,
  vertex.size = 8,
  vertex.label.cex = 0.6,
  edge.width = E(gene.graph.all)$weight * 4,
  layout = layout_with_fr(gene.graph.all),
  main = "Co-expression network of all 43 Boruta genes"
)


# ANALYZE 

# Connected components of the graph
components.results <- components(gene.graph.all)

# Size of each component
table(components.results$membership)

# Identify the largest connected component
largest.component <- which.max(components.results$csize)

# Genes in the largest component
main.network.genes <- names(components.results$membership[components.results$membership == largest.component])

length(main.network.genes)


# IDENTIFY POSSIBLE HUB GENES
# Weighted connectivity of genes in the graph
gene.strength <- strength(gene.graph.all,
                          mode = "all",
                          weights = E(gene.graph.all)$weight)

hub.results <- data.frame(gene = names(gene.strength),
                          connectivity = as.numeric(gene.strength)) %>%
  arrange(desc(connectivity))

head(hub.results, 20)
