library(Boruta)
library(survival)
library(parallel)
library(ranger)
library(tidyverse)
library(dplyr)



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

# 7. Run Parallelized Boruta ---------------------------

# Prepare clean X and Y

x_data <- boruta_df_small[, setdiff(colnames(boruta_df_small), "survival_status")]

y_data <- boruta_df_small$survival_status

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

dim(x_log2_filtered)

# Run boruta

set.seed(111)

boruta.signature <- Boruta(
  x = x_log2_filtered,
  y = y_data,
  getImp = impRangerSurv, 
  doTrace = 3,
  maxRuns = 500
)

#  8. View Results ------------------

print(boruta.signature)

# 2. Save the final model object

saveRDS(boruta.signature, "results/boruta/boruta_signature.RDS")


