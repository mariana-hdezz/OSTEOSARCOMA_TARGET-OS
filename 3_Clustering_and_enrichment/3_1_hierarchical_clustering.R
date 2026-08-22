library(factoextra)
library(dplyr)
library(tidyr)
library(survival)
library(coxphf)
library(cluster)
library(Boruta)

#############################################################################
#> Hierarchical clustering on train set (TARGET_OS patients), utilizing
#> the gene list obtained from Boruta.
#> 
#> It first creates the signature with the boruta output
#> 
#> Inputs: survival_signature, vst_counts, metadata_os, boruta_signature
#> 
#> Outputs: Overwrite metadata_os to inclide clusters column for future analysis,
#> gene_signature
#> 
#> Results:------------------------------------------------------------------ 
#> 
##> Plots:------------------------------------------------------------------ %>% 
###> surv_plot: Kaplan Meier curve with survival as outcome
###> surv_plot_rec: Kaplan Meier curve with recurrence as outcome
#
##> Statistics:-------------------------------------------------------------
###> summary_cox: Cox analysis of survival as outcome
###> summary_cox_rec: Cox analysis of recurrence as outcome
#
##> Accompanying plots and results (not as objects):------------------------
###> Silhouette plots
###> Silhouette mean
###> Dendogram
#############################################################################

# List obtained from Boruta

vst_counts <- readRDS("./output_data/vst_counts.RDS")

metadata_os <- readRDS("./output_data/metadata_os.RDS")

gene_signature <- scan("output_data/gene_signature.csv", sep = ",", what = character())

# Keep only genes for clustering

vst_counts_hc <- vst_counts[rownames(vst_counts) %in% gene_signature, ]

# Transpose so pt in rows and geenes in columns

vst_counts_hc <-
  vst_counts_hc %>%
  t()

# Clustering ----------------------------------------------------------

fviz_nbclust(vst_counts_hc, FUN = hcut, method = "silhouette")

# Distance matrix between samples

dist_counts <- get_dist(vst_counts_hc, method = "manhattan")

# Clustering

hc_counts <- hclust(dist_counts, method = "ward.D2")


# Clustering characteristics ------------------------------------------

# Dendrogram

plot(hc_counts,
     labels = FALSE,
     hang = -1,
     main = "CLUSTERING JERÁRQUICO TARGET-OS")

rect.hclust(hc_counts, k = 3, border = "purple") # bottom up approach

# Tree of clusters

clusters <- cutree(hc_counts, k = 3)

sil <- silhouette(clusters, dist_counts)

mean(sil[, "sil_width"])

# 6. Plot the silhouette profile
plot(sil, col = 2:(length(unique(as.integer(clusters))) + 1), main = "Hierarchical Silhouette Plot")

# Observe how many patients in each cluster

table(clusters)

# Create object as data frame

cluster_df <- data.frame(clusters = clusters)

# Asign a column named sample with the rownames such that we can then merge based on that column

cluster_df$sample <- rownames(vst_counts_hc)

# Merge object so that in the metadata there is a column corresponding to that patients cluster

metadata_os <-
  metadata_os %>%
  left_join(cluster_df, by = "sample")


metadata_os %>%
  group_by(clusters) %>%
  dplyr::count(survival_stat)


metadata_os %>%
  mutate(
    survival_stat = factor(survival_stat),
    clusters = factor(clusters),
    metastasis_at_diagnosis = factor(metastasis_at_diagnosis)
  ) %>%
  tidyr::drop_na(survival_stat) %>%
  ggplot(aes(x = survival_stat, fill = factor(relapse_stat))) +
  geom_histogram(stat = "count") +
  facet_wrap( ~ clusters) + 
  scale_fill_manual(values = c("#8d79dd", "#55c3fcfb"), labels = c("FALSE" = "No relapse", "TRUE" = "Relapse")) +
  scale_x_discrete(labels = c("0" = "Alive", "1" = "Deceased")) + 
  labs(x = "Survival", y = "N. of patients", fill = "Relapse") +
  theme_classic()

metadata_os$clusters <- relevel(factor(metadata_os$clusters), 2)

cox <- survival::coxph(
  Surv(metadata_os$survival_time, metadata_os$survival_stat) ~ clusters,
  metadata_os %>% mutate(clusters = factor(clusters))
)

summary_cox <- summary(cox)

pha <- survival::cox.zph(cox)
fit_km <- survival::survfit(Surv(metadata_os$survival_time, metadata_os$survival_stat) ~ clusters,
                            metadata_os)

surv_plot <- survminer::ggsurvplot(
  fit_km,
  data = metadata_os,
  pval = TRUE,
  risk.table = TRUE,
  
  xlim = c(0, 6000),
  break.time.by = 500,
  ggtheme = theme_minimal(),
  
  
  
  linewidth = 3,
  # Line size
  palette = c("#c380d3" , "#ff89d4", "#33ccff"),
)

surv_plot$plot <- surv_plot$plot + 
  annotate(
    geom = "text", 
    x = 500,          # X-axis position
    y = 0.10,         # Y-axis position
    label = paste0("PH assumption ", round(pha$table[1,3], 2)), 
    color = "black", 
    size = 5, 
    fontface = "bold"
  )

surv_plot

# ------ Recurrence -----------


metadata_os %>%
  group_by(clusters) %>%
  dplyr::count(relapse_stat)

metadata_os$clusters <- relevel(factor(metadata_os$clusters), 2)

cox_rec <- survival::coxph(
  Surv(metadata_os$time_to_first_event, metadata_os$relapse_stat) ~ clusters,
  metadata_os %>% mutate(clusters = factor(clusters))
)

summary_cox_rec <- summary(cox_rec)

pha_rec <- survival::cox.zph(cox_rec)
fit_km_rec <- survival::survfit(Surv(metadata_os$time_to_first_event, metadata_os$relapse_stat) ~ clusters,
                                metadata_os)

surv_plot_rec <- survminer::ggsurvplot(
  fit_km_rec,
  data = metadata_os,
  pval = TRUE,
  risk.table = TRUE,
  
  xlim = c(0, 6000),
  break.time.by = 500,
  ggtheme = theme_minimal(),
  
  
  
  linewidth = 3,
  # Line size
  palette = c("#c380d3" , "#ff89d4", "#33ccff"),
)

surv_plot_rec$plot <- surv_plot_rec$plot + 
  annotate(
    geom = "text", 
    x = 500,          # X-axis position
    y = 0.10,         # Y-axis position
    label = paste0("PH assumption ", round(pha_rec$table[1,3], 2)), 
    color = "black", 
    size = 5, 
    fontface = "bold"
  )

surv_plot_rec


# Huvos chi squared -------------------------------------------------------

huvos_meta <- metadata_os %>% 
  tidyr::drop_na(huvos_bin)

huvos_chi <- chisq.test(table(huvos_meta$huvos_bin, huvos_meta$clusters), simulate.p.value = TRUE, B = 10000)


saveRDS(metadata_os, "./output_data/metadata_os.RDS")



cat("--------------Cox results--------------\n")
cat("#####################\nCluster", levels(metadata_os$clusters)[1], "as reference:\n\n", "Cox survival analysis:\n#####################\n\n\n"); print(summary_cox)
cat("#####################\nCluster", levels(metadata_os$clusters)[1], "as reference:\n\n", "Cox recurrence analysis:\n#####################\n\n\n"); print(summary_cox_rec)

cat("Comparing response to neoadjuvant treatment between clusters")
cat("Chi squared results\n"); print(huvos_chi)
cat("Residuals\n"); print(huvos_chi$residuals)


rm(list = ls())
