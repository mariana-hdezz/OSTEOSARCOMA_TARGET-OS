library(factoextra)
library(tidyverse)

# List obtained from Boruta


gene_list_bor_sur_22 <- c("ACTA2", "ACTG2", "ATG4A", "CGREF1", "CORO6", "DENND2C",
                          "EDIL3", "EFTUD2", "GBP1", "GRAMD1B", "MINDY2", "MMP27", 
                          "MRC2", "MXI1", "NSUN6", "PODXL2", "RHBDL2", "SMPD3", 
                          "TAC4", "TTC9B", "TUBA1A", "ZNF587B")

gene_list_bor_rec_21 <- c("ACTA2", "ACTG2", "ADAM10", "BBOX1", "BMP8B", "CLK3", 
                          "COL13A1", "COL22A1", "ERICH1", "GBP1", "GCNT4", "PCDHB6", 
                          "PIP5K1C", "PPIL2", "PUM3", "RHBDL2", "SLC12A4", "SNAPC3", 
                          "SPICE1", "TAF5L", "TAS2R10")

gene_list_bor_survbin_43 <- c("APEX2", "ARHGAP1", "ARHGEF39", "CCDC97", "CGREF1", 
                              "CLUAP1", "COL22A1", "CPE", "CTNNBIP1", "CYFIP1", 
                              "DHRS11", "DLX1", "ERCC4", "F13A1", "FAM110D", 
                              "FAT1", "FKBP11", "GALNT14", "GBP1", "GMIP", 
                              "GRAMD1B", "HSD11B2", "INPP4A", "KERA", "KIF25",
                              "LGR6", "LURAP1L", "MEF2A", "MRTFB", "MXI1",
                              "NUBP1", "SF3B3", "SLC12A4", "SLC45A4", "SLC8A3", 
                              "STAT5B", "TIMM50", "TPD52", "TRIM68", "TSHZ3", 
                              "UBE2D4", "UNC5B", "VMP1")

gene_list_bor_metbin_19 <- c("AFMID", "ARHGEF2", "ATP6V0D1", "ATRX", "CKMT2", 
                             "DOCK8", "FBXW8", "GADD45GIP1", "GPC1", "GSTCD", 
                             "GUF1", "IL17RA", "ITPR3", "MAP7D1", "MYO19", 
                             "RAP1B", "TCAP", "TFDP1", "TNK2")

# Keep only genes for clustering

vst_counts_hc <- vst_counts[rownames(vst_counts) %in% gene_list, ]

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

rect.hclust(hc_counts, k = 2, border = "purple") # bottom up approach

# Tree of clusters

clusters <- cutree(hc_counts, k = 4)

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
  mutate(survival_stat = factor(survival_stat),
         clusters = factor(clusters)) %>%
  tidyr::drop_na(survival_stat) %>%
  ggplot(aes(x = survival_stat, fill = clusters)) +
  geom_histogram(stat = "count")

metadata_os %>%
  mutate(
    survival_stat = factor(survival_stat),
    clusters = factor(clusters),
    metastasis_at_diagnosis = factor(metastasis_at_diagnosis)
  ) %>%
  tidyr::drop_na(survival_stat) %>%
  ggplot(aes(x = survival_stat, fill = `First Event` == "Relapse")) +
  geom_histogram(stat = "count") +
  facet_wrap( ~ clusters)


metadata_os$clusters <- relevel(factor(metadata_os$clusters), 2)

cox <- survival::coxph(
  Surv(metadata_os$survival_time, metadata_os$survival_stat) ~ clusters,
  metadata_os %>% mutate(clusters = factor(clusters))
)

firth_fit <- coxphf(
  formula = Surv(survival_time, survival_stat) ~ factor(clusters), 
  data = metadata_os %>% drop_na(survival_stat)
)

summary(firth_fit)

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

metadata_os$clusters <- NULL

#rm(list = setdiff(ls(), c("vst_counts", "counts_data", "metadata_os")))
