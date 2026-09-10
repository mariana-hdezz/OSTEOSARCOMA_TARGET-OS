#############################################################################
#> Script to observe attempt to replicate results from https://doi.org/10.1371/journal.pone.0295364
#> 
#> Inputs: fpkm_data and metadata_os
#> 
#> Outputs: No outputs utilized in other analysis
#> 
#> Results
#> 
############################################################################


# Libraries

library(dplyr)
library(tibble)
library(tidyr)
library(ConsensusClusterPlus)

# Data load

fpkm <- readRDS("./output_data/fpkm_data.RDS")

fpkm_log <- log(fpkm + 1)

metadata_os <- readRDS("./output_data/metadata_os.RDS")

gene_signature <- c("TOR1A", "PDK2", "PML", "MAPKAPK1", "G6PD", "ARNTL", "FBXO7", "TLDC2", "GFER", "CAT", "PXN", "TXN2", "ERCC2", "MAP3K5", "GSKIP", "TPM1", "SIRPA", "MAPK1", "HGF", "PRNP", "PAWR", "TREM2", "EGFR", "HMOX1", "MGST1", "AXL", "APOE", "CA3", "SMPD3", "GPR37", "CD36", "CYGB", "STC2", "PSIP1", "NOL3", "NUDT1", "BNIP3", "PDK1", "SIGMAR1", "ERMP1", "ANKZF1", "RPS3", "ATF4", "GPX7", "FANCC")


fpkm_log_hc <- fpkm_log[rownames(fpkm_log) %in% gene_signature, ]

# Prepare for clustering

fpkm_matrix <- as.matrix(fpkm_log_hc)


fpkm_scaled <- t(scale(t(fpkm_matrix)))

# 3.3 Clustering

results <- ConsensusClusterPlus(
  fpkm_scaled,
  maxK = 10, 
  reps = 500, 
  clusterAlg = "hc", 
  distance = "pearson",
  seed = 500,
  plot = NULL 
)

# Metadata for patients new cluster

metadata_filtro <- metadata_os[metadata_os$sample %in% names(results[[2]]$consensusClass),]

# Delete our cluster

metadata_filtro$clusters <- NULL
 
#  Clusters and names of patients

x <- as.data.frame(results[[2]]$consensusClass) %>% 
  dplyr::rename(cluster = "results[[2]]$consensusClass") %>% 
  rownames_to_column("sample")
  
table(x$cluster)

# Join with metadata

metadata_os_f <- 
  metadata_filtro %>% 
  left_join(x, by = "sample")

# Table where we join our clusters and the simulation of theur clusters for the same patient
# We also binarize our clsuters into bad repsonse (cluster 1)and good response (c2 and c3)
# so as to compare theri good prognoisis clsuter (c2) with our biunarizewd version

comp <- tibble(
  lab_comp = ifelse(metadata_os$clusters == 1, 1, 2),
  lab = metadata_os$clusters,
  ros = metadata_os_f$cluster,
  y = seq_along(metadata_os_f$cluster)
) %>% 
  mutate(change = case_when(
    (lab == 1 & ros == 1) | (lab == 2 & ros == 2) | (lab == 3 & ros == 2) ~ "no_change",
    (lab == 1 & ros == 2) | (lab == 1 & ros == 3) ~"change_to_2",
    (lab == 2 & ros == 1) | (lab == 3 & ros == 1) ~ "change_to_1"
  ))

# Object for plot

comp_long <- 
  comp  %>% 
  pivot_longer(cols = c(lab, ros),
               names_to = "vers",
               values_to = "cluster") %>% 
  mutate(cluster = factor(cluster))

# Plot

comp_long %>% 
  ggplot(aes(x = vers, y = y, colour = cluster)) + 
  geom_point() + 
  theme_classic(base_size = 25) + 
  geom_line(aes(group = y, colour = change), alpha = 0.2, arrow = grid::arrow(ends = "last"))

# Proportion of patients that changed cluster (in the binarized form)

sum(comp$lab_comp == comp$ros) / length(comp$lab_comp)

# To where did they change

table(comp$change, comp$lab)

