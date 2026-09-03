#############################################################################
#> Script for clinical analysis in TARGET-OS clusters
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
#
##> Accompanying plots and results (not as objects):------------------------
###> Silhouette plots
###> Silhouette mean
###> Dendogram
#############################################################################

library(factoextra)
library(dplyr)
library(tidyr)
library(survival)
library(coxphf)

# Load data

metadata_os <- readRDS("./output_data/metadata_os.RDS")


# Survival ----------------------------------------------------------------


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

metadata_os$clusters <- relevel(factor(metadata_os$clusters), 3)

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



# Recurrence --------------------------------------------------------------


metadata_os %>%
  group_by(clusters) %>%
  dplyr::count(relapse_stat)

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


# Huvos chi squared -------------------------------------------------------

huvos_meta <- metadata_os %>% 
  tidyr::drop_na(huvos_bin)

huvos_chi <- chisq.test(table(huvos_meta$huvos_bin, huvos_meta$clusters), simulate.p.value = TRUE, B = 10000)


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