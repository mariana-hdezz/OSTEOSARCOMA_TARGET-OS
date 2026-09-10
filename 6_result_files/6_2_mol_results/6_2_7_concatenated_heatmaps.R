
library(patchwork)
library(ggplotify)

heatmap_hm_target           <- readRDS("./results/mul_hm/heatmap_go_target.RDS")
heatmap_go_target           <- readRDS("./results/mul_hm/heatmap_go_target.RDS")
mean_heatmap_gsea_target    <- readRDS("./results/mul_hm/mean_heatmap_gsea_target.RDS")
mean_heatmap_gsea_hm_target <- readRDS("./results/mul_hm/mean_heatmap_gsea_hm_target.RDS")
mean_heatmap_gsea_hm_gse    <- readRDS("./results/mul_hm/mean_heatmap_gsea_hm_gse.RDS")
mean_heatmap_gsea_gse       <- readRDS("./results/mul_hm/mean_heatmap_gsea_gse.RDS")
gse21257_plot_gse_go        <- readRDS("./results/mul_hm/gse21257_plot_gse_go.RDS")
gse21257_plot_gse_hm        <- readRDS("./results/mul_hm/gse21257_plot_gse_hm.RDS")
gse33382_plot_gse_go        <- readRDS("./results/mul_hm/gse33382_plot_gse_go.RDS")
gse33382_plot_gse_hm        <- readRDS("./results/mul_hm/gse33382_plot_gse_hm.RDS")
heatmap_hist_hm             <- readRDS("./results/mul_hm/heatmap_hist_hm.RDS")


heatmap_hm_target <- as.ggplot(heatmap_hm_target)
heatmap_go_target <- as.ggplot(heatmap_go_target)
gse21257_plot_gse_go <- as.ggplot(gse21257_plot_gse_go)
gse21257_plot_gse_hm <- as.ggplot(gse21257_plot_gse_hm)
gse33382_plot_gse_go <- as.ggplot(gse33382_plot_gse_go)
gse33382_plot_gse_hm <- as.ggplot(gse33382_plot_gse_hm)

(((heatmap_go_target / gse33382_plot_gse_go) + patchwork::plot_layout(heights = c(4,3))) | (heatmap_hm_target / gse33382_plot_gse_hm) ) +
  plot_annotation(tag_levels = "A")



mean_heatmap_gsea_hm_gse    <- as.ggplot(mean_heatmap_gsea_hm_gse   )
mean_heatmap_gsea_gse       <- as.ggplot(mean_heatmap_gsea_gse      )
mean_heatmap_gsea_target    <- as.ggplot(mean_heatmap_gsea_target   )
mean_heatmap_gsea_hm_target <- as.ggplot(mean_heatmap_gsea_hm_target)
heatmap_hist_hm             <- as.ggplot(heatmap_hist_hm            )

(((mean_heatmap_gsea_target / mean_heatmap_gsea_gse) + plot_layout(heights = c(5, 3))) | 
    (mean_heatmap_gsea_hm_target / mean_heatmap_gsea_hm_gse) | 
    heatmap_hist_hm) + 
  plot_layout(widths = c(3, 3, 4)) +
  plot_annotation(tag_levels = "A")


