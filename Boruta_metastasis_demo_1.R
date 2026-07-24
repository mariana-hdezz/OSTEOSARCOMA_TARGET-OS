library(Boruta)
library(survival)
library(parallel)
library(ranger)
library(tidyverse)
library(dplyr)

boruta_df_small <- readRDS("/datos/home/marh/OSTEOSARCOMA/boruta_df_1.rds")

# Load and transpose data

vst_counts <- t(vst_counts)

counts_data <- t(counts_data)

# Filter counts

keep_expression <- colSums(counts_data > 1) >= ceiling(0.10 * nrow(x_data)) # genes with more than 1 count on more than 10% of patients

boruta_counts <- vst_counts[, keep_expression, drop = FALSE]


# Join with metadata

# create x data and y data

x_data <-  boruta_df[, setdiff(colnames(boruta_df), "surv_obj")]

y_data

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

set.seed(111)

# Run boruta

boruta.signature <- Boruta(
  x = x_data,
  y = y_data,
  getImp = impRangerSurv, 
  doTrace = 3,
  maxRuns = 500
)

#  8. View Results ------------------

print(boruta.signature)
plot(boruta.signature, las = 2, cex.axis = 0.7)

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
out_path <- "/datos/home/marh/OSTEOSARCOMA/sign_demos"

# 1. Force the decision on those 22 tentative genes
final_boruta_decided <- TentativeRoughFix(boruta.signature)

final_boruta_decided$finalDecision[final_boruta_decided$finalDecision == "Confirmed"]

saveRDS(boruta.signature, paste0(out_path, "/boruta_met_sign-1.rds"))

# 2. Save the final model object
saveRDS(final_boruta_decided, paste0(out_path, "/final_boruta_met_decided-1.rds"))

stats <- attStats(boruta.signature)


# 3. Get the names of all selected genes (Confirmed + Fixed Tentatives)
final_gene_names <- getSelectedAttributes(final_boruta_decided)

# 4. Save the gene list as a CSV
write.csv(final_gene_names, paste0(out_path, "/final_gene_met_sign-1.csv"), row.names = FALSE)

# 5. Export the importance values (the numbers used in your plot)
stats <- attStats(final_boruta_decided)
write.csv(stats, paste0(out_path, "/gene_met_importance_full_stats-1.csv"))


# 9.- Load data ---------------------------

# Load the final decided Boruta object
final_boruta <- readRDS(paste0(out_path, "final_boruta_decided.rds"))

# Verify it loaded correctly
print(final_boruta)

final_boruta_decided$finalDecision

cat(final_gene_names, sep = ", ")
