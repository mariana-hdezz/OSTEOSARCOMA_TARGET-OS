#############################################################################
#> Script to perform differential expression analysis on the clusters
#> assigned to GSE patients 
#> 
#> Note. it also runs GSE calling the script "4_Validation/4_2_Validation_tests/4_2_3_gsea_val.R"
#> inside the diff expression loop
#>  
#> Inputs: respective GSE metadata and counts
#> 
#> Results consist on res objects which are the full differential expression results, res_sig
#> which are thje genes filtered by a p value of 0.05 and a LFC of 1.5 and GSEA 
#> results
#> 
##> gse21257_res_1v2_.csv
##> gse21257_res_3v1_.csv
##> gse21257_res_3v2_.csv
##> gse33382_res_1v2_.csv
##> gse33382_res_3v1_.csv
##> gse33382_res_3v2_.csv
##> 
##> 
##> gse21257_res_sig_1v2_.csv
##> gse21257_res_sig_3v1_.csv
##> gse21257_res_sig_3v2_.csv
##> gse33382_res_sig_1v2_.csv
##> gse33382_res_sig_3v1_.csv
##> gse33382_res_sig_3v2_.csv
##> 
##> gsea_c1_vs_c2_GO_gse21257.csv
##> gsea_c1_vs_c2_GO_gse33382.csv
##> gsea_c1_vs_c2_hm_gse21257.csv
##> gsea_c1_vs_c2_hm_gse33382.csv
##> gsea_c3_vs_c1_GO_gse21257.csv
##> gsea_c3_vs_c1_GO_gse33382.csv
##> gsea_c3_vs_c1_hm_gse21257.csv
##> gsea_c3_vs_c1_hm_gse33382.csv
##> gsea_c3_vs_c2_GO_gse21257.csv
##> gsea_c3_vs_c2_GO_gse33382.csv
##> gsea_c3_vs_c2_hm_gse21257.csv
##> gsea_c3_vs_c2_hm_gse33382.csv
##> 
############################################################################

library(limma)
library(dplyr)
library(tibble)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggridges)
library(enrichplot)
library(msigdbr)
library(ComplexHeatmap)
library(circlize)
library(ggtree)
library(ggplot2)
library(aplot)
library(tidyr)

# Load data

metadata_gse21257 <- readRDS("output_data/metadata_gse21257.RDS")
metadata_gse33382 <- readRDS("output_data/metadata_33382.RDS")

counts_data_gse21257 <- readRDS("output_data/counts_data_gse21257.RDS")
counts_data_gse33382 <- readRDS("output_data/counts_data_gse33382.RDS")

gene_signature_gse <- scan("output_data/gene_signature_gse.csv", sep = ",", what = character()) 

# Create output directory


if(dir.exists("./results/diffex_gsea_gse/")){
  "differential expr4ession and gsea directory already exists"
}else{
  dir.create("./results/diffex_gsea_gse/")
} 

# Loop to apply analysis to both datasets


for(t in c("gse21257", "gse33382")) {
  
  if(t == "gse21257"){
  
    metadata_difex <- metadata_gse21257
    
    counts_data_difex <- counts_data_gse21257
    
  }else if(t == "gse33382"){
    
    metadata_difex <- metadata_gse33382
    
    counts_data_difex <- counts_data_gse33382
    
  }
  
  # 1.1 Generate column corresponding to clusters nodes, those that have 0 in one group and those with more than 0 in another
  
  col_data <- metadata_difex %>% 
    dplyr::select(geo_accession, clusters) %>% # Create only the object to use for Limma
    column_to_rownames("geo_accession")
  
  
  
  
  # 2.- Differential expression -----------------------------------------------
  
  
  # 2.2 Data counts of the patients that had clusters node information in the metadata
  
  count_data <- counts_data_difex[colnames(counts_data_difex) %in% rownames(col_data)]
  
  # 2.2.2 Making shure they are in the same order
  
  count_data <- count_data[match(rownames(col_data), colnames(count_data))]
  
  all(colnames(count_data) == rownames(col_data))
  
  
  # 2.4 Design based on object separating on clusters nodes
  
  design <- model.matrix(~ 0 + clusters, data = col_data)
  
  # 2.4.2 Asign make.names objects as colnames
  
  colnames(design) <- make.names(colnames(design)) 
  
  # 2.5 Fit
  
  fit <- lmFit(count_data, design)
  
  
  for (i in 1:ncol(design)) {
    fit <- lmFit(count_data, design)
    
    if(i == 1){
      
      # 2.5.2 Contrast matrix comparing clusters 0 to > 0 clusters
      
      contrast.matrix <- makeContrasts(clusters1 - clusters2,
                                       levels = design)
      
    }else if(i == 2){
      
      contrast.matrix <- makeContrasts(clusters3 - clusters2,
                                       levels = design)
    }else if(i == 3){
      
      contrast.matrix <- makeContrasts(clusters3 - clusters1,
                                       levels = design)
    }
  
  
  
  # 2.5.3 Fit based on contrasts
  
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit)
  
  # 2.6 Results
  
  res <- topTable(fit, coef = 1, number = Inf)
  
  # 2.6.2 Results that correspond to a signfiicant p value and log fold change
  
  res_sig <- res %>%
    filter(adj.P.Val < 0.05 & abs(logFC) > 1.5) # 0.1
  
  if(i == 1){
    
    write.csv(res, paste0("results/diffex_gsea_gse/", t, "_res_1v2_.csv"))
    
    write.csv(res_sig, paste0("results/diffex_gsea_gse/", t, "_res_sig_1v2_.csv"))
  
    }else if(i == 2){
    
    write.csv(res, paste0("results/diffex_gsea_gse/", t, "_res_3v2_.csv"))
    
    write.csv(res_sig, paste0("results/diffex_gsea_gse/", t, "_res_sig_3v2_.csv"))
  
    }else if(i == 3){
    
    write.csv(res, paste0("results/diffex_gsea_gse/", t, "_res_3v1_.csv"))
    
    write.csv(res_sig, paste0("results/diffex_gsea_gse/", t, "_res_sig_3v1_.csv"))
  
    }
  
  source("5_molecular_analysis/5_2_val_mol/5_2_2_gsea_val.R")
  
  
  }
  
  
}

rm(list = ls())
gc()