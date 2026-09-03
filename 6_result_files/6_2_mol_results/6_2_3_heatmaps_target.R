#############################################################################
#> Script to obtain the heatmaps for the gsea of TARGET-OS analysis 
#> 
#> Inputs: All of the GSEA results from 5_molecular_analysis/5_1_target_mol/5_1_1_diffEx_gsea_target.R
#> and 5_molecular_analysis/5_1_target_mol/5_1_2_mean_ranked_gsea.R
#> 
#> Outputs: No outputs used in further analysis
#> 
#> Results: 
##> Heaatmap for the GSEA comparing between all 3 clusters for both GO
##> and Hallmark terms
##> 
##> Heatmap for the mean ranked GSEA with both GO and Hallmark terms
##>  
#############################################################################



# Diff_expr and gsea ------------------------------------------------------

# Load objects for heatmap

gsea_df_GO3_vs_1 <- read_csv("results/diffex_gsea_target/gsea_df_GO3_vs_1.csv")
gsea_df_GO3_vs_2 <- read_csv("results/diffex_gsea_target/gsea_df_GO3_vs_2.csv")
gsea_df_GO1_vs_2 <- read_csv("results/diffex_gsea_target/gsea_df_GO1_vs_2.csv")
gsea_df_GO3_vs_1 <- read_csv("results/diffex_gsea_target/gsea_df_GO3_vs_1.csv")
gsea_df_HM3_vs_1 <- read_csv("results/diffex_gsea_target/gsea_df_HM3_vs_1.csv")
gsea_df_HM3_vs_2 <- read_csv("results/diffex_gsea_target/gsea_df_HM3_vs_2.csv")
gsea_df_HM1_vs_2 <- read_csv("results/diffex_gsea_target/gsea_df_HM1_vs_2.csv")



# Objects for heatmaps

c1_v_c2_GO <- data.frame(row.names = gsea_df_GO1_vs_2$Description, # Rownames contains tha pathway
                         C1_vs_C2 = gsea_df_GO1_vs_2$NES) # Column named after the comparison contains NES

c3_v_c1_GO <- data.frame(row.names = gsea_df_GO3_vs_1$Description,
                         C3_vs_C1 = gsea_df_GO3_vs_1$NES)

c3_v_c2_GO <- data.frame(row.names = gsea_df_GO3_vs_2$Description,
                         C3_vs_C2 = gsea_df_GO3_vs_2$NES)


c1_v_c2_HM <- data.frame(row.names = gsea_df_HM1_vs_2$Description,
                         C1_vs_C2 = gsea_df_HM1_vs_2$NES)

c3_v_c1_HM <- data.frame(row.names = gsea_df_HM3_vs_1$Description,
                         C3_vs_C1 = gsea_df_HM3_vs_1$NES)

c3_v_c2_HM <- data.frame(row.names = gsea_df_HM3_vs_2$Description,
                         C3_vs_C2 = gsea_df_HM3_vs_2$NES)

# limit how many paths t plot

c1_v_c2_GO <- c1_v_c2_GO %>% 
  filter(C1_vs_C2 > quantile(C1_vs_C2, 0.9) | C1_vs_C2 < quantile(C1_vs_C2, 0.05))

c3_v_c1_GO <- c3_v_c1_GO %>% 
  
  filter(C3_vs_C1 > quantile(C3_vs_C1, 0.9) | C3_vs_C1 < quantile(C3_vs_C1, 0.05))

c3_v_c2_GO <- c3_v_c2_GO %>% 
  filter(C3_vs_C2 > quantile(C3_vs_C2, 0.9) | C3_vs_C2 < quantile(C3_vs_C2, 0.05))


# Merge the isolated columns to be able tp plot

gsea_GO_heatmap_obj <-  merge(c1_v_c2_GO, c3_v_c1_GO , by = 0, all = TRUE) %>%
  column_to_rownames("Row.names") %>%
  merge(c3_v_c2_GO, by = 0, all = TRUE) %>% 
  column_to_rownames("Row.names") 

#Na to 0

gsea_GO_heatmap_obj[is.na(gsea_GO_heatmap_obj )] <- 0



gsea_HM_heatmap_obj <-  merge(c1_v_c2_HM, c3_v_c1_HM , by = 0, all = TRUE) %>%
  column_to_rownames("Row.names") %>%
  merge(c3_v_c2_HM, by = 0, all = TRUE) %>%
  column_to_rownames("Row.names") 

gsea_HM_heatmap_obj[is.na(gsea_HM_heatmap_obj )] <- 0

# Calculate distances for dendogram

row_go <- hclust(dist(gsea_GO_heatmap_obj))
col_go <- hclust(dist(t(gsea_GO_heatmap_obj)))

# Dendograms

tree_right_go <- ggtree(row_go) + 
  scale_x_reverse() + 
  scale_y_continuous(expand = c(0, 0))

tree_top_go <- ggtree(col_go, hang = -1) + 
  layout_dendrogram() + 
  scale_y_reverse(expand = c(0, 0))

# Pivot longer to be able to plot

obj_for_GO <- gsea_GO_heatmap_obj %>% 
  rownames_to_column("path") %>% 
  pivot_longer(cols = c(C1_vs_C2, C3_vs_C1, C3_vs_C2),
               names_to = "clusters")

# Prepare object

heatmap_go <- obj_for_GO %>% 
  ggplot(aes(x = clusters, y = path, fill = value)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Spectral") + 
  theme_classic() + 
  scale_x_discrete(labels = c(
    "C3_vs_C1" = "C3 vs C1",
    "C1_vs_C2" = "C1 vs C2",
    "C3_vs_C2" = "C3 vs C2")) + 
  labs(title = "GSEA between clusters Gene Ontology. TARGET-OS", x = "Clusters") + 
  theme(
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 9),
    plot.title = element_text(size = 11)
  )



# Plot

heatmap_go %>% 
  insert_right(tree_right_go, width = 0.1) %>% 
  insert_top(tree_top_go, height = 0.1)

# Repeat for Hallmarks

row_hm <- hclust(dist(gsea_HM_heatmap_obj))
col_hm <- hclust(dist(t(gsea_HM_heatmap_obj)))

tree_right_hm <- ggtree(row_hm) + 
  scale_x_reverse() + 
  scale_y_continuous(expand = c(0, 0))

tree_top_hm <- ggtree(col_hm, hang = -1) + 
  layout_dendrogram() + 
  scale_y_reverse(expand = c(0, 0))


obj_for_hm <- gsea_HM_heatmap_obj %>% 
  rownames_to_column("path") %>% 
  pivot_longer(cols = c(C1_vs_C2, C3_vs_C1, C3_vs_C2),
               names_to = "clusters")

heatmap_hm <- obj_for_hm %>%
  ggplot(aes(x = clusters, y = path, fill = value)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Spectral")+
  theme_classic() +
  scale_x_discrete(labels = c("C3_vs_C1" = "C3 vs C1", "C1_vs_C2" = "C1 vs C2", "C3_vs_C2" = "C3 vs C2")) +
  labs(title = "GSEA between clusters Hallmarks of cancer TARGET-OS", 
       x = "Clusters")

heatmap_hm %>% 
  insert_right(tree_right_hm, width = 0.1) %>% 
  insert_top(tree_top_hm, height = 0.1)


rm(list = ls())



# Mean_ranked GSEA --------------------------------------------------------

# Load objects for heatmap
# Save results for GSEA GO

c1_cent_GO       <- read.csv("./results/diffex_gsea_target/c1_cent_GO.csv") %>% 
  column_to_rownames("X")
c2_cent_GO       <- read.csv("./results/diffex_gsea_target/c2_cent_GO.csv") %>% 
  column_to_rownames("X")
c3_cent_GO       <- read.csv("./results/diffex_gsea_target/c3_cent_GO.csv") %>% 
  column_to_rownames("X")
c1_cent_hallmark <- read.csv("./results/diffex_gsea_target/c1_cent_hallmark.csv") %>% 
  column_to_rownames("X")
c2_cent_hallmark <- read.csv("./results/diffex_gsea_target/c2_cent_hallmark.csv") %>% 
  column_to_rownames("X")
c3_cent_hallmark <- read.csv("./results/diffex_gsea_target/c3_cent_hallmark.csv") %>% 
  column_to_rownames("X")

# Keep top and lower 10% of C1, C2 and C3

c1_cent_GO_10 <- c1_cent_GO %>% 
  filter(c1 > quantile(c1, 0.9) | c1 < quantile(c1, 0.1))

c2_cent_GO_10 <- c2_cent_GO %>% 
  filter(c2 > quantile(c2, 0.9) | c2 < quantile(c2, 0.1))

c3_cent_GO_10 <- c3_cent_GO %>% 
  filter(c3 > quantile(c3, 0.9) | c3 < quantile(c3, 0.1))

# Join and convert to matrix

gsea_GO_mat <- merge(c3_cent_GO_10, (merge(c1_cent_GO_10, c2_cent_GO_10, by = 0, all = TRUE) %>% tibble::column_to_rownames("Row.names")), by = 0, all = TRUE)

gsea_GO_mat <- gsea_GO_mat %>% 
  tibble::column_to_rownames("Row.names") %>% 
  as.matrix()

gsea_GO_mat[is.na(gsea_GO_mat)] <- 0

# Distances for dendofram

row_hc <- hclust(dist(gsea_GO_mat))
col_hc <- hclust(dist(t(gsea_GO_mat)))


# Heatmap

heatmap_gsea <- gsea_GO_mat %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("path") %>% 
  pivot_longer(cols = c("c1", "c2", "c3"), names_to = "clusters") %>% 
  ggplot(aes(x = clusters, y = path, fill = value)) +
  geom_tile() +
  scale_fill_distiller(palette = "Spectral", direction = -1) +
  scale_x_discrete(expand = c(0, 0), labels = c("c1" = "C1", "c3" = "C3", "c2" = "C2")) +  
  scale_y_discrete(expand = c(0, 0)) +  
  theme_classic() +
  labs(x = "Clusters", y = "Pathway", fill = "NES", title = "GSEA from clusters means") +
  theme(
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 9),
    plot.title = element_text(size = 11)
  )


# Dendograms

tree_right <- ggtree(row_hc) + 
  scale_x_reverse() + 
  scale_y_continuous(expand = c(0, 0))

tree_top <- ggtree(col_hc, hang = -1) + 
  layout_dendrogram() + 
  scale_y_reverse(expand = c(0, 0))

# Full heatmap

heatmap_gsea %>% 
  insert_right(tree_right, width = 0.1) %>% 
  insert_top(tree_top, height = 0.1)

# Heatmap for hallamrks

c1_cent_hallmark_10 <- c1_cent_hallmark %>% 
  filter(c1 > quantile(c1, 0.75) | c1 < quantile(c1, 0.75))

c2_cent_hallmark_10 <- c2_cent_hallmark %>% 
  filter(c2 > quantile(c2, 0.75) | c2 < quantile(c2, 0.75))

c3_cent_hallmark_10 <- c3_cent_hallmark %>% 
  filter(c3 > quantile(c3, 0.75) | c3 < quantile(c3, 0.75))

# Join and convert to matrix

gsea_hallmark_mat <- merge(c3_cent_hallmark_10, (merge(c1_cent_hallmark_10, c2_cent_hallmark_10, by = 0, all = TRUE) %>% tibble::column_to_rownames("Row.names")), by = 0, all = TRUE)

gsea_hallmark_mat <- gsea_hallmark_mat %>% 
  tibble::column_to_rownames("Row.names") %>% 
  as.matrix()

gsea_hallmark_mat[is.na(gsea_hallmark_mat)] <- 0


row_hc <- hclust(dist(gsea_hallmark_mat))
col_hc <- hclust(dist(t(gsea_hallmark_mat)))


heatmap_gsea_hm <- gsea_hallmark_mat %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("path") %>% 
  pivot_longer(cols = c("c1", "c2", "c3"), names_to = "clusters") %>% 
  ggplot(aes(x = clusters, y = path, fill = value)) +
  geom_tile() +
  scale_fill_distiller(palette = "Spectral", direction = -1) +
  scale_x_discrete(expand = c(0, 0), labels = c("c1" = "C1", "c3" = "C3", "c2" = "C2")) +  
  scale_y_discrete(expand = c(0, 0)) +  
  theme_classic() +
  labs(x = "Clusters", y = "Pathway", fill = "NES", title = "GSEA from clusters means")


tree_right_hm <- ggtree(row_hc) + 
  scale_x_reverse() + 
  scale_y_continuous(expand = c(0, 0))

tree_top_hm <- ggtree(col_hc, hang = -1) + 
  layout_dendrogram() + 
  scale_y_reverse(expand = c(0, 0))

heatmap_gsea_hm %>% 
  insert_right(tree_right_hm, width = 0.1) %>% 
  insert_top(tree_top_hm, height = 0.1)

rm(list = ls())
gc()

