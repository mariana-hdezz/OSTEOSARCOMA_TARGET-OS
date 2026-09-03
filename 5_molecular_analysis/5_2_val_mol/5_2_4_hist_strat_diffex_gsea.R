library(limma)
library(dplyr)
library(ggplot2)
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
library(aplot)


metadata_gse33382_for_merge <- readRDS("output_data/metadata_gse33382_for_merge.RDS")
metadata_gse21257_for_merge <- readRDS("output_data/metadata_gse21257_for_merge.RDS")
counts_merged <- readRDS("output_data/counts_merged.RDS")
gene_signature_gse <- scan("output_data/gene_signature_gse.csv", sep = ",", what = character())

metadata_merged <- bind_rows(metadata_gse33382_for_merge,
                             metadata_gse21257_for_merge)

table_cluster <- (table(metadata_merged$clusters, metadata_merged$hist_sub))

table_cluster <- table_cluster[, colSums(table_cluster >= 3) >= 2]


for(t in colnames(table_cluster)) {
  
  metadata_merged_hist_sub <- metadata_merged %>% 
    filter(
      hist_sub == as.character(t)
    )
  
  
  col_data <- metadata_merged_hist_sub %>% 
    column_to_rownames("geo_accession")
  
  counts_data_difex <- counts_merged
  
  
  # 2.- Differential expression -----------------------------------------------
  
  
  # 2.2 Data counts of the patients that had clusters node information in the metadata
  
  count_data <- counts_data_difex[colnames(counts_data_difex) %in% rownames(col_data)]
  
  # 2.2.2 Making shure they are in the same order
  
  count_data <- count_data[match(rownames(col_data), colnames(count_data))]
  
  all(colnames(count_data) == rownames(col_data))
  
  
  # 2.4 Design based on object separating on clusters nodes
  
  design <- model.matrix(~ 0  + clusters + cohort, data = col_data)
  
  # 2.4.2 Asign make.names objects as colnames
  
  colnames(design) <- make.names(colnames(design)) 
  
  # 2.5 Fit
  
  fit <- lmFit(count_data, design)
  
  
  for (i in (1:(ncol(design)-1))) {
    
    if(i == 1 & t != "telangiectatic"){
      
      print(t)
      print("1-2")
      # 2.5.2 Contrast matrix comparing clusters 0 to > 0 clusters
      
      contrast.matrix <- makeContrasts(clusters1 - clusters2,
                                       levels = design)
      
    }else if(i == 2 & t != "telangiectatic"){
      print(t)
      print("3-2")
      contrast.matrix <- makeContrasts(clusters3 - clusters2,
                                       levels = design)
    }else if(i == 3 | (i == 1 & t == "telangiectatic")){
      print(t)
      print("3-1")
      contrast.matrix <- makeContrasts(clusters3 - clusters1,
                                       levels = design)
    }else if(i == 2 & t == "telangiectatic"){
      break
    }
    
    
    
    # 2.5.3 Fit based on contrasts
    
    fit_2 <- contrasts.fit(fit, contrast.matrix)
    fit_2 <- eBayes(fit_2)
    
    # 2.6 Results
    
    res <- topTable(fit_2, coef = 1, number = Inf)
    
    # 2.6.2 Results that correspond to a signfiicant p value and log fold change
    
    res_sig <- res %>%
      filter(adj.P.Val < 0.05 & abs(logFC) > 1.5) # 0.1
    
    if(i == 1 & t != "telangiectatic"){
      
      write.csv(res, paste0("results/diffex_gsea_gse/", t, "_hist_strat_res_1v2_.csv"))
      
      write.csv(res_sig, paste0("results/diffex_gsea_gse/", t, "_hist_strat_res_sig_1v2_.csv"))
      
    }else if(i == 2 & t != "telangiectatic"){
      
      write.csv(res, paste0("results/diffex_gsea_gse/", t, "_hist_strat_res_3v2_.csv"))
      
      write.csv(res_sig, paste0("results/diffex_gsea_gse/", t, "_hist_strat_res_sig_3v2_.csv"))
      
    }else if(i == 3 | (i == 1 & t == "telangiectatic")){
      
      write.csv(res, paste0("results/diffex_gsea_gse/", t, "_hist_strat_res_3v1_.csv"))
      
      write.csv(res_sig, paste0("results/diffex_gsea_gse/", t, "_hist_strat_res_sig_3v1_.csv"))
      
    }
    
    source("5_molecular_analysis/5_2_val_mol/5_2_2_gsea_val.R")
    
    
  }
  
  
}

rm(list = ls())

gc()
