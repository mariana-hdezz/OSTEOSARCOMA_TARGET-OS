#############################################################################
#> Script for clinical results for both TARGET and validation
#> 
#> Results:------------------------------------------------------------------ 
#> 
##> Plots:------------------------------------------------------------------ 
###> surv_plot: Kaplan Meier curve with survival as outcome for TARGET
###> surv_plot_rec: Kaplan Meier curve with recurrence as outcome for TARGET
###> 
###> 
#
##> Statistics:-------------------------------------------------------------
###> summary_cox: Cox analysis of survival as outcome
###> summary_cox_rec: Cox analysis of recurrence as outcome
###> huvos_chi and its residuals
#>
#> At the end of the script prints results, it also prints a table and the KM plots
#>
#############################################################

# Libraries

library(flextable)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)

# Load data

cox <- readRDS("./results/clinical_res/cox.RDS")
cox_rec <- readRDS("./results/clinical_res/cox_rec.RDS")
surv_plot <- readRDS("./results/clinical_res/surv_plot.RDS")
surv_plot_rec <- readRDS("./results/clinical_res/surv_plot_rec.RDS")
huvos_chi <- readRDS("./results/clinical_res/huvos_chisqr_target.RDS")

metadata_33382 <- readRDS("./output_data/metadata_33382.RDS")

surv_plot_gse21257     <- readRDS("./results/clinical_res/surv_plot_gse21257.RDS")
surv_plot_gse21257_rec <- readRDS("./results/clinical_res/surv_plot_gse21257_rec.RDS")
cox_gse21257           <- readRDS("./results/clinical_res/cox_gse21257.RDS")
cox_gse21257_rec       <- readRDS("./results/clinical_res/cox_gse21257_rec.RDS")
heat_hist              <- readRDS("./results/clinical_res/heat_hist.RDS")
heat_huvos             <- readRDS("./results/clinical_res/heat_huvos.RDS")
hist_chi_sqr           <- readRDS("./results/clinical_res/hist_chi_sqr.RDS")
huvos_val_chsq         <- readRDS("./results/clinical_res/huvos_val_chsq.RDS")

# Summary of the survival cox

cox_sum <- summary(cox)

# Similar for gse results of survival

cox_sum_gse <- summary(cox_gse21257_rec)

# Similar but for recurrence analysis for both target and validation sets

cox_rec_sum <- summary(cox_rec)

cox_rec_sum_gse21257 <- summary(cox_gse21257_rec)

# Prepare table for survival target

surv_tab <- 
  data.frame(
    Cohort = "TARGET-OS",
    Comparison = c("Cluster 3 vs Cluster 1", "Cluster 3 vs Cluster 2"),
    Odds_ratio = c(round(cox_sum$coefficients[1, 2], 3), round(cox_sum$coefficients[2, 2], 3)),
    Conf_int_low = c(round(cox_sum$conf.int[1, 3], 3), round(cox_sum$conf.int[2, 3], 3)),
    Conf_int_high = c(round(cox_sum$conf.int[1, 4], 3), round(cox_sum$conf.int[2, 4], 3)),
    p_value = c(round(cox_sum$coefficients[1, 5], 4), round(cox_sum$coefficients[2, 5], 4))
  )

# Prepare table for survival validation

surv_tab_gse <- 
  data.frame(
    Cohort = "GSE21257",
    Comparison = c("Cluster 3 vs Cluster 1", "Cluster 3 vs Cluster 2"),
    Odds_ratio = c(round(cox_sum_gse$coefficients[1, 2], 3), round(cox_sum_gse$coefficients[2, 2], 3)),
    Conf_int_low = c(round(cox_sum_gse$conf.int[1, 3], 3), round(cox_sum_gse$conf.int[2, 3], 3)),
    Conf_int_high = c(round(cox_sum_gse$conf.int[1, 4], 3), round(cox_sum_gse$conf.int[2, 4], 3)),
    p_value = c(round(cox_sum_gse$coefficients[1, 5], 4), round(cox_sum_gse$coefficients[2, 5], 4))
  )
  
# Prepare surv table

table_cox_surv <- bind_rows(surv_tab, surv_tab_gse)

table_cox_surv <- table_cox_surv %>% 
  janitor::clean_names(case = "sentence")

# Same as previous lines but for recurrence

rec_tab <- 
  data.frame(
    Cohort = "TARGET-OS",
    Comparison = c("Cluster 3 vs Cluster 1", "Cluster 3 vs Cluster 2"),
    Odds_ratio = c(round(cox_rec_sum$coefficients[1, 2], 3), round(cox_rec_sum$coefficients[2, 2], 3)),
    Conf_int_low = c(round(cox_rec_sum$conf.int[1, 3], 3), round(cox_rec_sum$conf.int[2, 3], 3)),
    Conf_int_high = c(round(cox_rec_sum$conf.int[1, 4], 3), round(cox_rec_sum$conf.int[2, 4], 3)),
    p_value = c(round(cox_rec_sum$coefficients[1, 5], 5), round(cox_rec_sum$coefficients[2, 5], 4))
  )


rec_tab_gse <- 
  data.frame(
    Cohort = "GSE21257",
    Comparison = c("Cluster 3 vs Cluster 1", "Cluster 3 vs Cluster 2"),
    Odds_ratio = c(round(cox_rec_sum_gse21257$coefficients[1, 2], 3), round(cox_rec_sum_gse21257$coefficients[2, 2], 3)),
    Conf_int_low = c(round(cox_rec_sum_gse21257$conf.int[1, 3], 3), round(cox_rec_sum_gse21257$conf.int[2, 3], 3)),
    Conf_int_high = c(round(cox_rec_sum_gse21257$conf.int[1, 4], 3), round(cox_rec_sum_gse21257$conf.int[2, 4], 3)),
    p_value = c(round(cox_rec_sum_gse21257$coefficients[1, 5], 4), round(cox_rec_sum_gse21257$coefficients[2, 5], 4))
  )


table_cox_rec <- bind_rows(rec_tab, rec_tab_gse)

table_cox_rec <- table_cox_rec %>%  
  janitor::clean_names(case = "sentence")

# Headers for dividing between outcome in table

upper_header <-  data.frame(Cohort = "Survival")
middle_header <- data.frame(Cohort = "Recurrence")

# Create table

survival_analysis_table <- 
  bind_rows(
  upper_header,
  table_cox_surv,
  middle_header,
  table_cox_rec
) %>% 
  flextable() %>% 
  merge_at(i = 1, j = 1:6, part = "body") %>% 
  merge_at(i = 6, j = 1:6, part = "body") %>% 
  hline(i = c(1, 5, 6), part = "body") %>%  
  autofit()

cat("--------------Cox results--------------\n")

survival_analysis_table

cat("Comparing response to neoadjuvant treatment between clusters in TARGET-OS")
cat("Chi squared results\n"); print(huvos_chi)
cat("Residuals\n"); print(huvos_chi$residuals)

cat("\nComparing response to neoadjuvant treatment between clusters in merged data of GSE21257 and GSE33382")
cat("Chi squared results\n"); print(huvos_val_chsq)
cat("Residuals\n"); print(huvos_val_chsq$residuals)

cat("\nComparing histologic subtype between clusters in merged data of GSE21257 and GSE33382")
cat("Chi squared results\n"); print(hist_chi_sqr)
cat("Residuals\n"); print(hist_chi_sqr$residuals)



((surv_plot$plot / surv_plot_rec$plot) | (surv_plot_gse21257$plot / surv_plot_gse21257_rec$plot)) +
  patchwork::plot_annotation(tag_levels = "A")

heat_hist <- ggplotify::as.ggplot(heat_hist)

((heat_hist) / (heat_huvos)) +
  plot_layout(heights = c(5, 2), tag_level = "new") +
  plot_annotation(tag_levels = "A")  &
  theme(plot.tag = element_text(size = 20))


# Drop NA in metastasis for an analysisi in GSE33382

metadata_33382_no_na <- metadata_33382 %>% 
  drop_na(metastasis_5y) 


chisq.test(table(metadata_33382_no_na$clusters, metadata_33382_no_na$metastasis_5y), simulate.p.value = 2000)


rm(list = ls())
gc()


