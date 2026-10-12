#############################################################################
#> Script for clinical analysis in GSE38055 clusters
#> 
#> Results saved:------------------------------------------------------------------ 
#> 
##> Plots:------------------------------------------------------------------ 
###> surv_plot: Kaplan Meier curve with survival as outcome
###> surv_plot_rec: Kaplan Meier curve with recurrence as outcome
#
##> Statistics:-------------------------------------------------------------
###> summary_cox: Cox analysis of survival as outcome
###> summary_cox_rec: Cox analysis of recurrence as outcome
###> huvos_chi: Chi squared test comparing huvos grade between clusters
#
##> Accompanying plots and results (not as objects):------------------------
###> Silhouette plots
###> Silhouette mean
###> Dendogram
#>
#############################################################################

library(factoextra)
library(dplyr)
library(tidyr)
library(survival)
library(coxphf)

# Load data

metadata_gse39055 <- readRDS("./output_data/metadata_gse39055.RDS")


# Survival ----------------------------------------------------------------

metadata_gse39055 %>%
  group_by(clusters) %>%
  dplyr::count(SURV_STAT)

table(metadata_gse39055$clusters, metadata_gse39055$SURV_STAT)

x <- chisq.test(table(metadata_gse39055$clusters, metadata_gse39055$SURV_STAT), simulate.p.value = TRUE, B = 2000)



# Recurrence --------------------------------------------------------------

table(metadata_gse39055$clusters, metadata_gse39055$REC_STAT)

y <- chisq.test(table(metadata_gse39055$clusters, metadata_gse39055$REC_STAT), simulate.p.value = TRUE, B = 2000)

cox_rec <- survival::coxph(
  Surv(metadata_gse39055$REC_MON, metadata_gse39055$REC_STAT) ~ clusters,
  metadata_gse39055 %>% mutate(clusters = factor(clusters))
)

summary_cox_rec <- summary(cox_rec)

pha_rec <- survival::cox.zph(cox_rec)
fit_km_rec <- survival::survfit(Surv(metadata_gse39055$REC_MON, metadata_gse39055$REC_STAT) ~ clusters,
                                metadata_gse39055)

surv_plot_rec <- survminer::ggsurvplot(
  fit_km_rec,
  data = metadata_gse39055,
  pval = TRUE,
  risk.table = TRUE,
  
  xlim = c(0, 200),
  break.time.by = 50,
  ggtheme = theme_classic(base_size = 15),
  
  
  
  linewidth = 3,
  # Line size
  palette = c("#c380d3" , "#ff89d4", "#33ccff"),
)

surv_plot_rec$plot <- surv_plot_rec$plot + 
  annotate(
    geom = "text", 
    x = 150,          # X-axis position
    y = 0.10,         # Y-axis position
    label = paste0("PH assumption ", round(pha_rec$table[1,3], 2)), 
    color = "black", 
    size = 5, 
    fontface = "bold"
  ) +
  labs(y = "Recurrence")


# Huvos chi squared -------------------------------------------------------

huvos_meta <- metadata_gse39055 %>% 
  tidyr::drop_na(huvos_bin)

huvos_chi <- chisq.test(table(huvos_meta$huvos_bin, huvos_meta$clusters), simulate.p.value = TRUE, B = 10000)

huvos_meta_plot <- 
  huvos_meta %>% 
  group_by(clusters) %>% 
  count(huvos_bin) %>% 
  ggplot(aes(x = clusters, y = huvos_bin, fill = n))+
  geom_tile(color = "white", lwd = 0.5, linetype = 1) + 
  scale_fill_distiller(palette = "Spectral", direction = -1) + 
  theme_classic(base_size = 30) +
  labs(x = "Clusters", y = "Huvos grade", fill = "Freq", title = "Proportion of patients with each Huvos grade for each cluster")

if(dir.exists("./results/clinical_res/")){
  "Directory already exists"
}else{
  dir.create("./results/clinical_res/")
}



saveRDS(cox , "./results/clinical_res/cox.RDS")
saveRDS(cox_rec , "./results/clinical_res/cox_rec.RDS")
saveRDS(surv_plot , "./results/clinical_res/surv_plot.RDS")
saveRDS(surv_plot_rec , "./results/clinical_res/surv_plot_rec.RDS")
saveRDS(huvos_chi, "./results/clinical_res/huvos_chisqr_target.RDS")
saveRDS(huvos_meta_plot, "./results/clinical_res/huvos_meta_plot.RDS")

rm(list = ls())
gc()
