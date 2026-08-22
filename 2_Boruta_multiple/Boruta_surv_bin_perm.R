
print("all_instaled")

library(Boruta)
library(survival)
library(ranger)
library(dplyr)
library(tibble)
library(tidyr)



fpkm_data <- readRDS("output_data/fpkm_data.RDS")
metadata_os <- readRDS("output_data/metadata_os.RDS")

fpkm_data <- fpkm_data[-1,]

metadata_os_bor <- metadata_os %>% 
  dplyr::select(sample, 
                survival_stat) 

boruta_df_small <- 
  fpkm_data %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  left_join(metadata_os_bor, by = "sample") %>% 
  dplyr::rename(survival_status = survival_stat) %>% 
  column_to_rownames("sample") %>% 
  drop_na(survival_status) 


# 5.3 Fix any non-standard gene names 

colnames(boruta_df_small) <- make.names(colnames(boruta_df_small))

boruta_df_small <- boruta_df_small %>% 
  mutate(survival_status = factor(survival_status))

dim(boruta_df_small)

# 6. Parallel Function  ---------------------------------

impRangerSurv <- function(x, y, ...) {
  temp_df <- as.data.frame(x) # 1. Combine X and Y into one dataframe for ranger, x corresponds to the gene names and y to the surv obj
  temp_df$target_surv <- y
  
  # 2. Use 'dependent.variable.name' instead of 'formula'
  # This is the "fast" interface for high-dimensional data
  
  res <- ranger::ranger(
    dependent.variable.name = "target_surv", 
    data = temp_df, # Data frame created before
    importance = "permutation", 
    num.trees = 500,            
    num.threads = 6,            # Ensure threads are passed here
    ...
  )
  return(res$variable.importance)
}

boruta_list <- list()
boruta_tent <- list()
boruta_res <- list()

for (i in 1:100) {
  
  set.seed(110 + i)
  
  if(i < 25){
    leave_5_out <- sample(c(1:85), ceiling(85 * 0.05), replace = FALSE)
    
    boruta_df_small_pert <-  boruta_df_small[- leave_5_out, ]
    
    print(dim(boruta_df_small_pert))
    
  }else if(i >= 25 & i < 50){
    leave_10_out <- sample(c(1:85), ceiling(85 * 0.1), replace = FALSE)
    
    boruta_df_small_pert <-  boruta_df_small[- leave_10_out, ]
    
    print(dim(boruta_df_small_pert))
    
  }else if(i >= 50 & i < 75){
    
    leave_15_out <- sample(c(1:85), ceiling(85 * 0.15), replace = FALSE)
    
    boruta_df_small_pert <-  boruta_df_small[- leave_15_out, ]
    
    print(dim(boruta_df_small_pert))
    
  }else if(i >= 75){
    leave_20_out <- sample(c(1:85), ceiling(85 * 0.2), replace = FALSE)
    
    boruta_df_small_pert <-  boruta_df_small[- leave_20_out, ]
    
    print(dim(boruta_df_small_pert))
  }
  
  
  # Prepare clean X and Y
  
  x_data <- boruta_df_small_pert[, setdiff(colnames(boruta_df_small_pert), "survival_status")]
  
  y_data <- boruta_df_small_pert$survival_status
  
  # log2 transformation (FPKM + 1)
  x_log2 <- log2(x_data + 1)
  
  
  # Genes expresados en al menos 10% de los pacientes
  keep_expression <- colSums(x_data > 1) >= ceiling(0.10 * nrow(x_data))
  
  x_log2_filtered <- x_log2[, keep_expression, drop = FALSE]
  
  # Eliminate genes with zero variance
  gene_variance <- apply(x_log2_filtered, 2, var, na.rm = TRUE)
  
  x_log2_filtered <- x_log2_filtered[
    ,
    is.finite(gene_variance) & gene_variance > 0,
    drop = FALSE
  ]
  
  # Round for reproducibility since when testing on different computers, some decimals where different on positions >15 after log transformation
  
  x_log2_filtered <- round(x_log2_filtered, 10)
  
  print(i)
  
  # Run boruta
  
  set.seed(110 + i)
  
  boruta.signature <- Boruta(
    x = x_log2_filtered,
    y = y_data,
    getImp = impRangerSurv, 
    doTrace = 1,
    maxRuns = 500
  )
  
  
  boruta_list[[i]] <- names(boruta.signature$finalDecision)[boruta.signature$finalDecision == "Confirmed"]
  boruta_tent[[i]] <- names(TentativeRoughFix(boruta.signature)$finalDecision)[TentativeRoughFix(boruta.signature)$finalDecision == "Confirmed" &  boruta.signature$finalDecision == "Tentative"]
  
  boruta_res[[i]] <- boruta.signature
  
}

# 2. Save the final model object

saveRDS(boruta_res, "results/boruta/boruta_signature_pert.RDS")

saveRDS(boruta_list, "results/boruta/boruta_conf_pert.RDS")

saveRDS(boruta_tent, "results/boruta/boruta_tent_pert.RDS")


