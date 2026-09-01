library(patchwork)
library(ggplot2)
library(tidyr)

metadata_gse21257 <- readRDS("output_data/metadata_gse21257.RDS")
metadata_33382 <- readRDS("output_data/metadata_33382.RDS")

metadata_gse21257$clusters <- as.factor(metadata_gse21257$clusters)

metadata_gse21257$clusters <- relevel(metadata_gse21257$clusters, 3)

cox_gse21257_rec <- survival::coxph(Surv(as.numeric(metadata_gse21257$relapse_time), metadata_gse21257$relapse_stat) ~ clusters, metadata_gse21257)


pha <- survival::cox.zph(cox_gse21257_rec)

fit_km <- survival::survfit(Surv(as.numeric(metadata_gse21257$relapse_time), metadata_gse21257$relapse_stat) ~ clusters, metadata_gse21257)

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

surv_plot_gse21257_rec

summary(cox_gse21257_rec)



cox_gse21257 <- survival::coxph(Surv(as.numeric(metadata_gse21257$survival_time), metadata_gse21257$survival_stat) ~ clusters, metadata_gse21257)


pha <- survival::cox.zph(cox_gse21257)

fit_km <- survival::survfit(Surv(as.numeric(metadata_gse21257$survival_time), metadata_gse21257$survival_stat) ~ clusters, metadata_gse21257)

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

surv_plot_gse21257

summary(cox_gse21257)


metadata_gse21257 %>% dplyr::group_by(clusters) %>% dplyr::count(survival_stat)

metadata_gse21257 %>% group_by(clusters) %>% count(relapse_stat)


metadata_33382_no_na <- metadata_33382%>% 
  drop_na(metastasis_5y) 


chisq.test(table(metadata_33382_no_na$clusters, metadata_33382_no_na$metastasis_5y), simulate.p.value = 2000)


x <- chisq.test(table(metad_heatmap$clusters, metad_heatmap$hist_sub), simulate.p.value = TRUE, B = 10000)
x$residuals

