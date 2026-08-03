library(Boruta)
library(survival)
library(parallel)
library(ranger)
library(tidyverse)
library(dplyr)

boruta_df_small <- readRDS("/datos/home/marh/OSTEOSARCOMA/boruta_df_1.rds")

boruta_df_small <- boruta_df_small %>% 
  dplyr::select(-survival_status)


# 5.3 Fix any non-standard gene names 

colnames(boruta_df_small) <- make.names(colnames(boruta_df_small))

boruta_df_small <- boruta_df_small %>% 
  mutate(metastasis_at_diagnosis = factor(metastasis_at_diagnosis))

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
    num.threads = 16,            # Ensure threads are passed here
    ...
  )
  return(res$variable.importance)
}

# 7. Run Parallelized Boruta ---------------------------

# Prepare clean X and Y

x_data <- boruta_df_small[, setdiff(colnames(boruta_df_small), "metatastasis_at_diagnosis")]

x_data <- as.data.frame(apply(x_data, 2, as.numeric))

y_data <- boruta_df_small$metastasis_at_diagnosis

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

boruta.signature <- boruta.signature

# Only the 11 confirmed genes
confirmed_only <- getSelectedAttributes(boruta.signature, withTentative = FALSE)

# The 22 tentative  genes
# (We find these by taking the full list and removing the confirmed ones)

all_selected <- getSelectedAttributes(boruta.signature, withTentative = TRUE)
tentative_only <- setdiff(all_selected, confirmed_only)

cat(confirmed_only, sep = ", ")
cat(tentative_only, sep = ", ")


# Define path
out_path <- "/datos/home/marh/OSTEOSARCOMA/Final_signatures/Metastasis_bin"

# 1. Force the decision on those 22 tentative genes
final_boruta_decided <- TentativeRoughFix(boruta.signature)

final_boruta_decided <- final_boruta_decided$finalDecision[final_boruta_decided$finalDecision == "Confirmed"]

saveRDS(boruta.signature, paste0(out_path, "/boruta_met_bin_18.rds"))

# 2. Save the final model object
saveRDS(final_boruta_decided, paste0(out_path, "/final_boruta_decided_met_bin.rds"))

stats <- attStats(boruta.signature)


# 3. Get the names of all selected genes (Confirmed + Fixed Tentatives)
final_gene_names <- names(final_boruta_decided)


# 4. Save the gene list as a CSV
write.csv(final_gene_names, paste0(out_path, "/final_gene_met_bin_18.csv"), row.names = FALSE)

# 5. Export the importance values (the numbers used in your plot)
stats <- attStats(boruta.signature)
write.csv(stats, paste0(out_path, "/gene_importance_full_stats_met_bin_18.csv"))


# 9.- Load data ---------------------------

# Load the final decided Boruta object
final_boruta <- readRDS(paste0(out_path, "final_boruta_decided.rds"))

# Verify it loaded correctly
print(final_boruta)

final_boruta_decided$finalDecision

cat(final_gene_names, sep = ", ")

cat(paste0('"', names(final_boruta_decided) , '"'), sep = ", ")
