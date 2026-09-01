
library(flextable)

cox <- readRDS("results/res_3_1_hc/cox.RDS")
cox_rec <- readRDS("results/res_3_1_hc/cox_rec.RDS")
surv_plot <- readRDS("results/res_3_1_hc/surv_plot.RDS")
surv_plot_rec <- readRDS("results/res_3_1_hc/surv_plot_rec.RDS")

cox
cox_gse21257

cox_rec
cox_gse21257_rec

cox_sum <- summary(cox)
cox_sum_gse <- summary(cox_gse21257_rec)

cox_sum$coefficients
cox_sum$conf.int

cox_rec_sum <- summary(cox_rec)

cox_rec_sum_gse21257 <- summary(cox_gse21257_rec)


surv_tab <- 
  data.frame(
    Cohort = "TARGET-OS",
    Comparison = c("Cluster 3 vs Cluster 1", "Cluster 3 vs Cluster 2"),
    Odds_ratio = c(round(cox_sum$coefficients[1, 2], 3), round(cox_sum$coefficients[2, 2], 3)),
    Conf_int_low = c(round(cox_sum$conf.int[1, 3], 3), round(cox_sum$conf.int[2, 3], 3)),
    Conf_int_high = c(round(cox_sum$conf.int[1, 4], 3), round(cox_sum$conf.int[2, 4], 3)),
    p_value = c(round(cox_sum$coefficients[1, 5], 4), round(cox_sum$coefficients[2, 5], 4))
  )


surv_tab_gse <- 
  data.frame(
    Cohort = "GSE21257",
    Comparison = c("Cluster 3 vs Cluster 1", "Cluster 3 vs Cluster 2"),
    Odds_ratio = c(round(cox_sum_gse$coefficients[1, 2], 3), round(cox_sum_gse$coefficients[2, 2], 3)),
    Conf_int_low = c(round(cox_sum_gse$conf.int[1, 3], 3), round(cox_sum_gse$conf.int[2, 3], 3)),
    Conf_int_high = c(round(cox_sum_gse$conf.int[1, 4], 3), round(cox_sum_gse$conf.int[2, 4], 3)),
    p_value = c(round(cox_sum_gse$coefficients[1, 5], 4), round(cox_sum_gse$coefficients[2, 5], 4))
  )
  
table_cox_surv <- bind_rows(surv_tab, surv_tab_gse)

table_cox_surv <- table_cox_surv %>% 
  janitor::clean_names(case = "sentence")


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



upper_header <-  data.frame(Cohort = "Survival")
middle_header <- data.frame(Cohort = "Recurrence")

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
