library(ggplot2)
library(ggtree)


metadata_gse21257 <- readRDS("output_data/metadata_gse21257.RDS")
metadata_33382 <- readRDS("output_data/metadata_33382.RDS")

counts_data_gse21257 <- readRDS("output_data/counts_data_gse21257.RDS")
counts_data_gse33382 <- readRDS("output_data/counts_data_gse33382.RDS")


metad_heatmap_33382 <- metadata_33382 %>% 
  mutate(response = ifelse(huvos %in% c("1", "2"), "Poor (<90%)", "Good (>=90%)"),
         metastasis_5y = case_when(
           metastasis_5y == "no"~ 0,
           metastasis_5y == "yes" ~ 1,
           is.na(metastasis_5y) ~ NA
         ),
         cohort = "GSE33382") %>% 
  dplyr::select(clusters, hist_sub, huvos, response, metastasis_5y, cohort) %>% 
  rename(hist_sub = "hist_sub" ,
         relapse_stat = "metastasis_5y")


metad_heatmap_21257 <- metadata_gse21257 %>%
  mutate(response = ifelse(huvos %in% c("1", "2"), "Poor (<90%)", "Good (>=90%)"),
         cohort = "GSE21257") %>% 
  dplyr::select(clusters, hist_sub, huvos, response, relapse_stat, cohort) %>% 
  mutate(hist_sub = tolower(hist_sub),
         huvos = tolower(huvos))


metad_heatmap <- bind_rows(metad_heatmap_21257, metad_heatmap_33382) %>% 
  mutate(hist_sub = paste0(toupper(substr(hist_sub, 1, 1)), substr(hist_sub, 2, nchar(hist_sub))))

row_hc_hist <- hclust(dist(table(metad_heatmap$hist_sub, metad_heatmap$clusters)))
col_hc_hist <- hclust(dist(t(table(metad_heatmap$hist_sub, metad_heatmap$clusters))))

metad_table_hist <- 
  metad_heatmap %>% 
  group_by(clusters) %>% 
  dplyr::count(hist_sub) %>% 
  ungroup() %>% 
  group_by(hist_sub) %>% 
  mutate(prop = n / sum(n))

heat_hist <- metad_table_hist %>% 
  ggplot(aes(x = clusters, y = hist_sub, fill = prop)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) + 
  scale_fill_distiller(palette = "Spectral", direction = -1) + 
  theme_classic(base_size = 30) +
  labs(x = "Clusters", y = "Histologic subtype", fill = "Freq", title = "Proportion of cluster for each histologic subtype")

tree_right_hist <- ggtree(row_hc_hist, hang = -1) + 
  scale_x_reverse() + 
  scale_y_continuous(expand = c(0, 0))

tree_top_hist <- ggtree(col_hc_hist, hang = -1) + 
  layout_dendrogram() + 
  scale_y_reverse(expand = c(0, 0))

heat_hist %>% 
  insert_right(tree_right_hist, width = 0.1) %>% 
  insert_top(tree_top_hist, height = 0.1)

x <- chisq.test(table(metad_heatmap$clusters, metad_heatmap$hist_sub), simulate.p.value = TRUE, B = 10000)
x$residuals


huvos_clean <- metad_heatmap %>%
  filter(huvos != "unknown") %>%
  mutate(
    huvos = factor(huvos, levels = c("1", "2", "3", "4")),
    clusters = factor(clusters)
  )

tab_clean <- table(huvos_clean$clusters, huvos_clean$huvos)

print(tab_clean)


huvos_binary <- huvos_clean %>%
  mutate(response = ifelse(huvos %in% c("1", "2"), "Poor (<90%)", "Good (>=90%)"))

tab_binary <- table(huvos_binary$clusters, huvos_binary$response)

y <- chisq.test(tab_binary)
y$residuals



metad_table_huvos <- 
  metad_heatmap %>% 
  filter(!(huvos == "unknown")) %>% 
  group_by(clusters) %>% 
  dplyr::count(huvos) %>% 
  ungroup() %>% 
  group_by(huvos) %>% 
  mutate(prop = n / sum(n))

heat_huvos <- metad_table_huvos %>% 
  ggplot(aes(x = clusters, y = huvos, fill = prop)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) + 
  scale_fill_distiller(palette = "Spectral", direction = -1) + 
  theme_classic(base_size = 30) +
  labs(x = "Clusters", y = "Histologic subtype", fill = "Freq", title = "Proportion of cluster for each histologic subtype")


