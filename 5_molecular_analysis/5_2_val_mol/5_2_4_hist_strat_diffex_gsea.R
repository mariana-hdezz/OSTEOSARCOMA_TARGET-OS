#############################################################################
#> Script that perforrms differential expression stratified by histologic subtypes
#> 
#> Steps: 
##> 
#>  
#> Inputs: metadata_gse33382_for_merge, metadata_gse21257_for_merge, counts_merged
#> 
#> Outputs: No outputs used in further analsis
#> 
#> Results: 
##> Objects with GSEA results:
###> c1_cent_GO
###> c2_cent_hallmark
###> c3_cent_GO
###> c1_cent_hallmark
###> c2_cent_GO
###> c3_cent_hallmark
###> Heatmaps
#
#############################################################################

# Libraries

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

# Load data

metadata_gse33382_for_merge <- readRDS("output_data/metadata_gse33382_for_merge.RDS")
metadata_gse21257_for_merge <- readRDS("output_data/metadata_gse21257_for_merge.RDS")
counts_merged <- readRDS("output_data/counts_merged.RDS")
gene_signature_gse <- scan("output_data/gene_signature_gse.csv", sep = ",", what = character())

# Merge metadata

metadata_merged <- bind_rows(metadata_gse33382_for_merge,
                             metadata_gse21257_for_merge)

# Create table

table_cluster <- (table(metadata_merged$clusters, metadata_merged$hist_sub))

# Filter histologic subtypes that hace at least 3 patietns in at least 2 clusters

table_cluster <- table_cluster[, colSums(table_cluster >= 3) >= 2]

# For loop

for(t in colnames(table_cluster)) { # Filtered histologic subtypes as input for loop
  # Notice that t will be a character vector corresponding to a histologic subtype
  # This is important for when the GSEA script is sourced since it names objects based on t 
  # so the reuslt will have the hist subtype name
  
  metadata_merged_hist_sub <- metadata_merged %>%  # Keep histologic subtype selected for that loop
    filter(
      hist_sub == as.character(t)
    )
  

  col_data <- metadata_merged_hist_sub %>% 
    column_to_rownames("geo_accession")
  
  counts_data_difex <- counts_merged
  
  
  # Differential expression -----------------------------------------------
  
  # Data counts 
  
  count_data <- counts_data_difex[colnames(counts_data_difex) %in% rownames(col_data)]
  
  # Making shure they are in the same order
  
  count_data <- count_data[match(rownames(col_data), colnames(count_data))]
  
  all(colnames(count_data) == rownames(col_data))
  
  
  # Design based on clusters adjusting by cohort batch effect
  
  design <- model.matrix(~ 0  + clusters + cohort, data = col_data)
  
  # Asign make.names objects as colnames
  
  colnames(design) <- make.names(colnames(design)) 
  
  # Fit
  
  fit <- lmFit(count_data, design)
  
  
  for (i in (1:(ncol(design)-1))) { # For loop that creates the contrast matrix based on hist subtype
    
    if(i == 1 & t != "telangiectatic"){ # SInce telangiectatic only compares c3 vs c1 then we dont want for the 1st loop of telangiectatic to go along, for the other subtypes we do
      
      print(t)
      print("1-2")

      contrast.matrix <- makeContrasts(clusters1 - clusters2,
                                       levels = design)
      
    }else if(i == 2 & t != "telangiectatic"){ # SInce telangiectatic only compares c3 vs c1 then we dont want for the 2nd loop of telangiectatic to go along
      print(t)
      print("3-2")
      contrast.matrix <- makeContrasts(clusters3 - clusters2,
                                       levels = design)
    }else if(i == 3 | (i == 1 & t == "telangiectatic")){ # SInce telangiectatic only compares c3 vs c1 then we do want for the 1st loop of telangiectatic to go along here
      print(t)
      print("3-1")
      contrast.matrix <- makeContrasts(clusters3 - clusters1,
                                       levels = design)
    }else if(i == 2 & t == "telangiectatic"){
      break
    }
    
    
    # Fit based on contrasts
    
    fit_2 <- contrasts.fit(fit, contrast.matrix)
    fit_2 <- eBayes(fit_2)
    
    # Results
    
    res <- topTable(fit_2, coef = 1, number = Inf)
    
    #  Results that correspond to a signfiicant p value and log fold change
    
    res_sig <- res %>%
      filter(adj.P.Val < 0.05 & abs(logFC) > 1.5) # 0.1
    
    if(i == 1 & t != "telangiectatic"){ # Once again takes into account telangiectatic to not save that under the comparisons of c1 vc c2 and c2 and c3
                                        # For the other subtypes the if runs nirmally based on the previous loop
      
      write.csv(res, paste0("results/diffex_gsea_gse/", t, "_hist_strat_res_1v2_.csv"))
      
      write.csv(res_sig, paste0("results/diffex_gsea_gse/", t, "_hist_strat_res_sig_1v2_.csv"))
      
    }else if(i == 2 & t != "telangiectatic"){
      
      write.csv(res, paste0("results/diffex_gsea_gse/", t, "_hist_strat_res_3v2_.csv"))
      
      write.csv(res_sig, paste0("results/diffex_gsea_gse/", t, "_hist_strat_res_sig_3v2_.csv"))
      
    }else if(i == 3 | (i == 1 & t == "telangiectatic")){ # Save telangiectatic comparison when due
      
      write.csv(res, paste0("results/diffex_gsea_gse/", t, "_hist_strat_res_3v1_.csv"))
      
      write.csv(res_sig, paste0("results/diffex_gsea_gse/", t, "_hist_strat_res_sig_3v1_.csv"))
      
    }
    
    source("7_isolated_functions/5_2_2_gsea_val.R")    
    
  }
  
  
}

rm(list = ls())

gc()
