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
  dplyr::arrange(FDR)

head(module_trait_results, 20)

module_trait_results %>%
  dplyr::filter(FDR < 0.05)
