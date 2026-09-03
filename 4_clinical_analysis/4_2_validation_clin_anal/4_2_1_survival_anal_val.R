#############################################################################
#> Script for clinical analysis in GSE clusters
#> 
#> Results saved:------------------------------------------------------------------ 
#> 
##> Plots:------------------------------------------------------------------ 
###> surv_plot_gse21257: Kaplan Meier curve with survival as outcome for gse21257
###> surv_plot_gse21257_rec: Kaplan Meier curve with recurrence as outcome
###> heat_hist: heatmap for histologic subtypes across clusters
###> heat_huvos: huvos grade across clusters
#
##> Statistics:-------------------------------------------------------------
###> cox_gse21257:  cox analysis of survival as outcome
###> cox_gse21257_rec: Cox analysis of recurrence as outcome
###> hist_chi_sqr: chi squared test comparing hist subtypes across clusters
###> huvos_val_chsq: chi squared test comparing binary huvos response across clusters
#
#############################################################################

# Libraries

library(patchwork)
library(ggplot2)
library(tidyr)
library(ggtree)
library(aplot)

# Load data

metadata_gse21257 <- readRDS("output_data/metadata_gse21257.RDS")
metadata_33382 <- readRDS("output_data/metadata_33382.RDS")
counts_data_gse21257 <- readRDS("output_data/counts_data_gse21257.RDS")
counts_data_gse33382 <- readRDS("output_data/counts_data_gse33382.RDS")

# Assure clusters are a factor and relevel to 3 as reference

metadata_gse21257$clusters <- as.factor(metadata_gse21257$clusters)

metadata_gse21257$clusters <- relevel(metadata_gse21257$clusters, 3)

# Cox model for relapse

cox_gse21257_rec <- survival::coxph(Surv(as.numeric(metadata_gse21257$relapse_time), metadata_gse21257$relapse_stat) ~ clusters, metadata_gse21257)

# PH asumption

pha <- survival::cox.zph(cox_gse21257_rec)

# Fit for KM plots

fit_km <- survival::survfit(Surv(as.numeric(metadata_gse21257$relapse_time), metadata_gse21257$relapse_stat) ~ clusters, metadata_gse21257)

# Prepare relapse plot

surv_plot_gse21257_rec <- survminer::ggsurvplot(
  fit_km,
  data = metadata_gse21257,
  pval = TRUE,
  risk.table = TRUE,
  
  xlim = c(0, 300),
  break.time.by = 50,
  ggtheme = theme_minimal(),
  
  linewidth = 3,
  
  palette = c("#c380d3" , "#ff89d4", "#33ccff"),
)

# Add ph asumption to plot

surv_plot_gse21257_rec$plot <- surv_plot_gse21257_rec$plot + 
  annotate(
    geom = "text", 
    x = 30,      
    y = 0,           
    label = paste0("PH assumption ", round(pha$table[1,3], 2)), 
    color = "black", 
    size = 5, 
    fontface = "bold"
  )

# Cox for survival

cox_gse21257 <- survival::coxph(Surv(as.numeric(metadata_gse21257$survival_time), metadata_gse21257$survival_stat) ~ clusters, metadata_gse21257)

# PH asumption

pha <- survival::cox.zph(cox_gse21257)

fit_km <- survival::survfit(Surv(as.numeric(metadata_gse21257$survival_time), metadata_gse21257$survival_stat) ~ clusters, metadata_gse21257)

# Prepare plot for survival

surv_plot_gse21257 <- survminer::ggsurvplot(
  fit_km,
  data = metadata_gse21257,
  pval = TRUE,
  risk.table = TRUE,
  
  xlim = c(0, 300),
  break.time.by = 50,
  ggtheme = theme_minimal(),
  
  
  
  linewidth = 3,

  palette = c("#c380d3" , "#ff89d4", "#33ccff"),
)

surv_plot_gse21257$plot <- surv_plot_gse21257$plot + 
  annotate(
    geom = "text", 
    x = 30,         
    y = 0,         
    label = paste0("PH assumption ", round(pha$table[1,3], 2)), 
    color = "black", 
    size = 5, 
    fontface = "bold"
  )


# Hist sub and huvos analysis ---------------------------------------------

# Select the columns to fuse metadata objects


metad_heatmap_33382 <- metadata_33382 %>% 
  dplyr::select(clusters, hist_sub, huvos, response, metastasis_5y, cohort) %>% 
  dplyr::rename(hist_sub = "hist_sub" ,
                relapse_stat = "metastasis_5y")


metad_heatmap_21257 <- metadata_gse21257 %>% 
  dplyr::select(clusters, hist_sub, huvos, response, relapse_stat, cohort) %>% 
  mutate(hist_sub = tolower(hist_sub),
         huvos = tolower(huvos))

# Bind the metadata objects

metad_heatmap <- bind_rows(metad_heatmap_21257, metad_heatmap_33382) %>% 
  mutate(hist_sub = paste0(toupper(substr(hist_sub, 1, 1)), substr(hist_sub, 2, nchar(hist_sub))))

# Distances for dendograms in heatmap

row_hc_hist <- hclust(dist(table(metad_heatmap$hist_sub, metad_heatmap$clusters)))
col_hc_hist <- hclust(dist(t(table(metad_heatmap$hist_sub, metad_heatmap$clusters))))

# Obtain counts and proportions of hist subtype by cluster

metad_table_hist <- 
  metad_heatmap %>% 
  group_by(clusters) %>% 
  dplyr::count(hist_sub) %>% 
  ungroup() %>% 
  group_by(hist_sub) %>% 
  mutate(prop = n / sum(n))

# Initiall heatmap

heat_hist <- metad_table_hist %>% 
  ggplot(aes(x = clusters, y = hist_sub, fill = prop)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) + 
  scale_fill_distiller(palette = "Spectral", direction = -1) + 
  theme_classic(base_size = 30) +
  labs(x = "Clusters", y = "Histologic subtype", fill = "Freq", title = "Proportion of cluster for each histologic subtype")

# Prepare dendograms

tree_right_hist <- ggtree(row_hc_hist, hang = -1) + 
  scale_x_reverse() + 
  scale_y_continuous(expand = c(0, 0))

tree_top_hist <- ggtree(col_hc_hist, hang = -1) + 
  layout_dendrogram() + 
  scale_y_reverse(expand = c(0, 0))

# Add dendograms

heat_hist <- 
  heat_hist %>% 
  insert_right(tree_right_hist, width = 0.1) %>% 
  insert_top(tree_top_hist, height = 0.1)

# Chi squared test comparing hist subtypes across clusters

hist_chi_sqr <- chisq.test(table(metad_heatmap$clusters, metad_heatmap$hist_sub), simulate.p.value = TRUE, B = 10000)

# Huvos heatmap ignoring unknown

huvos_clean <- metad_heatmap %>%
  filter(huvos != "unknown") %>%
  mutate(
    huvos = factor(huvos, levels = c("1", "2", "3", "4")),
    clusters = factor(clusters)
  )

# Convert to binary

huvos_binary <- huvos_clean %>%
  mutate(response = ifelse(huvos %in% c("1", "2"), "Poor (<90%)", "Good (>=90%)"))

tab_binary <- table(huvos_binary$clusters, huvos_binary$response)

# Chi squared test comparing huvos between clusters

huvos_val_chsq <- chisq.test(tab_binary)

# Prepare for hatmap

metad_table_huvos <- 
  metad_heatmap %>% 
  filter(!(huvos == "unknown")) %>% 
  group_by(clusters) %>% 
  dplyr::count(huvos) %>% 
  ungroup() %>% 
  group_by(huvos) %>% 
  mutate(prop = n / sum(n))

# Create heatmap

heat_huvos <- metad_table_huvos %>% 
  ggplot(aes(x = clusters, y = huvos, fill = prop)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) + 
  scale_fill_distiller(palette = "Spectral", direction = -1) + 
  theme_classic(base_size = 30) +
  labs(x = "Clusters", y = "Histologic subtype", fill = "Freq", title = "Proportion of cluster for each histologic subtype")


saveRDS(surv_plot_gse21257     , "./results/clinical_res/surv_plot_gse21257.RDS")
saveRDS(surv_plot_gse21257_rec , "./results/clinical_res/surv_plot_gse21257_rec.RDS")
saveRDS(cox_gse21257           , "./results/clinical_res/cox_gse21257.RDS")
saveRDS(cox_gse21257_rec       , "./results/clinical_res/cox_gse21257_rec.RDS")
saveRDS(heat_hist              , "./results/clinical_res/heat_hist.RDS")
saveRDS(heat_huvos             , "./results/clinical_res/heat_huvos.RDS")
saveRDS(hist_chi_sqr           , "./results/clinical_res/hist_chi_sqr.RDS")
saveRDS(huvos_val_chsq         , "./results/clinical_res/huvos_val_chsq.RDS")

rm(list = ls())

gc()