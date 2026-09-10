#############################################################################
#> Script to obtain the heatmaps for the gsea of GSE analysis 
#> 
#> Inputs: All of the GSEA results from 5_molecular_analysis/5_2_val_mol/5_2_1_diffEx_val.R and
#> 5_molecular_analysis/5_2_val_mol/5_2_3_mean_ranked_gsea_val.R
#> 
#> 
#> Results: 
##> Heaatmap for the GSEA comparing between all 3 clusters for both GO
##> and Hallmark terms
##> 
##> Heatmap for the mean ranked GSEA with both GO and Hallmark terms
##>  
#############################################################################

# Load data

## Diff expr gsea

gsea_c1_vs_c2_GO_gse21257 <- read.csv("results/diffex_gsea_gse/gsea_c1_vs_c2_GO_gse21257.csv", row.names = 1)
gsea_c1_vs_c2_GO_gse33382 <- read.csv("results/diffex_gsea_gse/gsea_c1_vs_c2_GO_gse33382.csv", row.names = 1)
gsea_c1_vs_c2_hm_gse21257 <- read.csv("results/diffex_gsea_gse/gsea_c1_vs_c2_hm_gse21257.csv", row.names = 1)
gsea_c1_vs_c2_hm_gse33382 <- read.csv("results/diffex_gsea_gse/gsea_c1_vs_c2_hm_gse33382.csv", row.names = 1)
gsea_c3_vs_c2_GO_gse21257 <- read.csv("results/diffex_gsea_gse/gsea_c3_vs_c2_GO_gse21257.csv", row.names = 1)
gsea_c3_vs_c2_GO_gse33382 <- read.csv("results/diffex_gsea_gse/gsea_c3_vs_c2_GO_gse33382.csv", row.names = 1)
gsea_c3_vs_c2_hm_gse21257 <- read.csv("results/diffex_gsea_gse/gsea_c3_vs_c2_hm_gse21257.csv", row.names = 1)
gsea_c3_vs_c2_hm_gse33382 <- read.csv("results/diffex_gsea_gse/gsea_c3_vs_c2_hm_gse33382.csv", row.names = 1)
gsea_c3_vs_c1_GO_gse33382 <- read.csv("results/diffex_gsea_gse/gsea_c3_vs_c1_GO_gse33382.csv", row.names = 1)
gsea_c3_vs_c1_hm_gse33382 <- read.csv("results/diffex_gsea_gse/gsea_c3_vs_c1_hm_gse33382.csv", row.names = 1)
gsea_c3_vs_c1_GO_gse21257 <- read.csv("results/diffex_gsea_gse/gsea_c3_vs_c1_GO_gse21257.csv", row.names = 1)
gsea_c3_vs_c1_hm_gse21257 <- read.csv("results/diffex_gsea_gse/gsea_c3_vs_c1_hm_gse21257.csv", row.names = 1)


## Mean ranked gsea

c1_cent_GO_gse       <- read.csv("./results/diffex_gsea_gse/c1_cent_GO.csv"      ) %>% 
  column_to_rownames("X")
c1_cent_hallmark_gse <- read.csv("./results/diffex_gsea_gse/c1_cent_hallmark.csv") %>% 
  column_to_rownames("X")
c2_cent_GO_gse       <- read.csv("./results/diffex_gsea_gse/c2_cent_GO.csv"      ) %>% 
  column_to_rownames("X")
c2_cent_hallmark_gse <- read.csv("./results/diffex_gsea_gse/c2_cent_hallmark.csv") %>% 
  column_to_rownames("X")
c3_cent_GO_gse       <- read.csv("./results/diffex_gsea_gse/c3_cent_GO.csv"      ) %>% 
  column_to_rownames("X")
c3_cent_hallmark_gse <- read.csv("./results/diffex_gsea_gse/c3_cent_hallmark.csv") %>% 
  column_to_rownames("X")


# Heatmaps for differential expression ranekd GSEA -----------------------------


heatmap_names <- list()

#> How the loop works:
#> 
#> for e in (names of validation dataset): thus the first run will correspond to the objects of 
#> gse21257 and the second to objects corresponding to gse33382
#> 
#> If statement: here it evaluates which of the gse sets is being analyzed
#> (first gse21257 and then gse33382) and assigns the set sepcfic object to
#> a generall object that will be used in the rest of the script
#> 
#> So for example in the first run c3_vs_c2_GO will correspond to gsea_c3_vs_c2_GO_gse21257
#> and in the second run c3_vs_c2_GO will correpsond to gsea_c3_vs_c2_GO_gse33382 but
#> 
#> The next section is justthe preparing of the heatmaps for both GO and
#> hallmark terms
#> 
#> Finally the name of the object is created pasting the name of the cohort 
#> and the gsea term utilized (go or hm), the name is added to the list of names 
#> for future visualization and the heatmap object is assigned to the name

for (e in c("gse21257", "gse33382")) {
  
  
  if(e == "gse21257"){
    
    c3_vs_c2_GO <- gsea_c3_vs_c2_GO_gse21257
    c1_vs_c2_GO <- gsea_c1_vs_c2_GO_gse21257
    c3_vs_c1_GO <- gsea_c3_vs_c1_GO_gse21257
    c3_vs_c2_hallmark <- gsea_c3_vs_c2_hm_gse21257
    c1_vs_c2_hallmark <- gsea_c1_vs_c2_hm_gse21257
    c3_vs_c1_hallmark <- gsea_c3_vs_c1_hm_gse21257
    
  }else{
    
    c3_vs_c2_GO <- gsea_c3_vs_c2_GO_gse33382
    c1_vs_c2_GO <- gsea_c1_vs_c2_GO_gse33382
    c3_vs_c1_GO <- gsea_c3_vs_c1_GO_gse33382
    c3_vs_c2_hallmark <- gsea_c3_vs_c2_hm_gse33382
    c1_vs_c2_hallmark <- gsea_c1_vs_c2_hm_gse33382
    c3_vs_c1_hallmark <- gsea_c3_vs_c1_hm_gse33382
    
  }
  
  
  c3_vs_c2_GO <- data.frame(row.names = c3_vs_c2_GO$Description, # GO cluster 3 vs cluster 2
                            c3_vs_c2  = c3_vs_c2_GO$NES)
  
  
  c3_vs_c2_hallmark <- data.frame(row.names = c3_vs_c2_hallmark$Description, # hallmarks cluster 3 vs cluster 2
                                  c3_vs_c2  = c3_vs_c2_hallmark$NES)
  
  c1_vs_c2_GO <- data.frame(row.names = c1_vs_c2_GO$Description,
                            c1_vs_c2  = c1_vs_c2_GO$NES)
  
  c1_vs_c2_hallmark <- data.frame(row.names = c1_vs_c2_hallmark$Description,
                                  c1_vs_c2  = c1_vs_c2_hallmark$NES)
  
  c3_vs_c1_GO <- data.frame(row.names = c3_vs_c1_GO$Description,
                            c3_vs_c1  = c3_vs_c1_GO$NES)
  
  c3_vs_c1_hallmark <- data.frame(row.names = c3_vs_c1_hallmark$Description,
                                  c3_vs_c1  = c3_vs_c1_hallmark$NES)
  
  # Keep only the top
  
  
  c3_vs_c2_GO_10 <- c3_vs_c2_GO %>% 
    filter(c3_vs_c2 > quantile(c3_vs_c2, 0.95) | c3_vs_c2 < quantile(c3_vs_c2, 0.05))
  
  c1_vs_c2_GO_10 <- c1_vs_c2_GO %>% 
    filter(c1_vs_c2 > quantile(c1_vs_c2, 0.95) | c1_vs_c2 < quantile(c1_vs_c2, 0.05))
  
  c3_vs_c1_GO_10 <- c3_vs_c1_GO %>% 
    filter(c3_vs_c1 > quantile(c3_vs_c1, 0.95) | c3_vs_c1 < quantile(c3_vs_c1, 0.05))
  
  
  # Join and convert to matrix
  
  gsea_GO_mat <- merge(c3_vs_c2_GO_10, (merge(c1_vs_c2_GO_10, c3_vs_c1_GO_10, by = 0, all = TRUE) %>% column_to_rownames("Row.names")), by = 0, all = TRUE)
  
  gsea_GO_mat <- gsea_GO_mat %>% 
    column_to_rownames("Row.names") %>% 
    as.matrix()
  
  gsea_GO_mat[is.na(gsea_GO_mat)] <- 0
  
  row_hc_gse_GO <- hclust(dist(gsea_GO_mat))
  col_hc_gse_GO <- hclust(dist(t(gsea_GO_mat)))
  
  
  # This one changes the position of C3 V C2 and C3 V C1 so that it is the same order in TARGET 
  # for visualization we leave it commented for replication porpuses
  # 
  # col_hc_gse_GO$order
  # col_hc_gse_GO$labels <- col_hc_gse_GO$labels[c(3, 2, 1)]
  # 
  gsea_GO_mat <-
    gsea_GO_mat %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("path") %>% 
    pivot_longer(cols = c("c3_vs_c2", "c1_vs_c2", "c3_vs_c1"), names_to = "clusters") 
  
  
  heat_gse_GO <- gsea_GO_mat %>% 
    mutate(clusters = factor(clusters)) %>% 
    ggplot(aes(x = clusters, y = path, fill = value)) +
    geom_tile(color = "white", lwd = 0.5, linetype = 1) + 
    scale_fill_distiller(palette = "Spectral", direction = -1) + 
    theme_classic() +
    labs(x = "Clusters", y = "Pathways", fill = "NES", title = paste("GSEA between clusters Gene Ontology:", e)) +
    scale_x_discrete(labels = c(
      "c3_vs_c1" = "C3 vs C1",
      "c1_vs_c2" = "C1 vs C2",
      "c3_vs_c2" = "C3 vs C2")) +
    theme(
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 9),
      plot.title = element_text(size = 11),
      legend.position = "none"
    )
  
  tree_right_gse_GO <- ggtree(row_hc_gse_GO, hang = -1) + 
    scale_x_reverse() + 
    scale_y_continuous(expand = c(0, 0))
  
  tree_top_gse_GO <- ggtree(col_hc_gse_GO, hang = -1) + 
    layout_dendrogram() + 
    scale_y_reverse(expand = c(0, 0))
  
  heat_gse_GO[1]$data <- heat_gse_GO[1]$data %>% 
    mutate(clusters = factor(clusters, levels = c("c1_vs_c2", "c3_vs_c1", "c3_vs_c2")))
  
  plot_gse_go <- heat_gse_GO %>% 
    insert_right(tree_right_gse_GO, width = 0.1) %>% 
    insert_top(tree_top_gse_GO, height = 0.1)
  
  
  
  ##############################################################################
  
  c3_vs_c2_hallmark_10 <- c3_vs_c2_hallmark %>% 
    filter(c3_vs_c2 > quantile(c3_vs_c2, 0.85) | c3_vs_c2 < quantile(c3_vs_c2, 0.85))
  
  c1_vs_c2_hallmark_10 <- c1_vs_c2_hallmark %>% 
    filter(c1_vs_c2 > quantile(c1_vs_c2, 0.85) | c1_vs_c2 < quantile(c1_vs_c2, 0.85))
  
  c3_vs_c1_hallmark_10 <- c3_vs_c1_hallmark %>% 
    filter(c3_vs_c1 > quantile(c3_vs_c1, 0.85) | c3_vs_c1 < quantile(c3_vs_c1, 0.85))
  
  
  gsea_hallmark_mat <- merge(c3_vs_c2_hallmark_10, (merge(c1_vs_c2_hallmark_10, c3_vs_c1_hallmark_10, by = 0, all = TRUE) %>% column_to_rownames("Row.names")), by = 0, all = TRUE)
  
  gsea_hallmark_mat <- gsea_hallmark_mat %>% 
    column_to_rownames("Row.names") %>% 
    as.matrix()
  
  gsea_hallmark_mat[is.na(gsea_hallmark_mat)] <- 0
  
  row_hc_gse_hallmark <- hclust(dist(gsea_hallmark_mat))
  col_hc_gse_hallmark <- hclust(dist(t(gsea_hallmark_mat)))
  
  # This one changes the position of C3 V C2 and C3 V C1 so that it is the same order in TARGET 
  # for visualization we leave it commented for replication porpuses
  
  # col_hc_gse_hallmark$order
  # col_hc_gse_hallmark$labels <- col_hc_gse_hallmark$labels[c(3, 2, 1)]

  
    
  gsea_hallmark_mat <-
    gsea_hallmark_mat %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("path") %>% 
    pivot_longer(cols = c("c3_vs_c2", "c1_vs_c2", "c3_vs_c1"), names_to = "clusters") 
  
  
  
  heat_gse_hallmark <- 
    gsea_hallmark_mat %>% 
    mutate(clusters = factor(clusters)) %>% 
    ggplot(aes(x = clusters, y = path, fill = value)) +
    geom_tile(color = "white", lwd = 0.5, linetype = 1) + 
    scale_fill_distiller(palette = "Spectral", direction = -1) + 
    theme_classic() +
    labs(x = "Clusters", y = "Pathway", fill = "NES", title = paste("GSEA between clusters Hallmarks:", e)) +
    scale_x_discrete(labels = c(
      "c3_vs_c1" = "C3 vs C1",
      "c1_vs_c2" = "C1 vs C2",
      "c3_vs_c2" = "C3 vs C2"))
  
  tree_right_gse_hallmark <- ggtree(row_hc_gse_hallmark, hang = -1) + 
    scale_x_reverse() + 
    scale_y_continuous(expand = c(0, 0))
  
  tree_top_gse_hallmark <- ggtree(col_hc_gse_hallmark, hang = -1) + 
    layout_dendrogram() + 
    scale_y_reverse(expand = c(0, 0))
  
  heat_gse_hallmark[1]$data <- heat_gse_hallmark[1]$data %>% 
    mutate(clusters = factor(clusters, levels = c("c1_vs_c2", "c3_vs_c1", "c3_vs_c2")))
  
  plot_gse_hm <- heat_gse_hallmark %>% 
    insert_right(tree_right_gse_hallmark, width = 0.1) %>% 
    insert_top(tree_top_gse_hallmark, height = 0.1)
  
  
  # Create names: example in the first run e = gse21257 so name_plot_go will be
  # gse21257_plot_gse_go
  
  name_plot_go <- paste0(e, "_plot_gse_go")
  
  name_plot_hm <- paste0(e, "_plot_gse_hm")
  
  # Then add that to the list of names
  
  heatmap_names[[paste(e, "go")]] <- name_plot_go
  heatmap_names[[paste(e, "hm")]] <- name_plot_hm
  
  # Finally assign the heatmap to the objecy name
  
  assign(name_plot_go, plot_gse_go)
  
  assign(name_plot_hm, plot_gse_hm)
  
}


heatmap_names

gse21257_plot_gse_go
gse21257_plot_gse_hm
gse33382_plot_gse_go
gse33382_plot_gse_hm



# Mean ranked gsea --------------------------------------------------------


# Keep top and lower 25% of C1 and all of c2 and c3

c1_cent_GO_gse_10 <- c1_cent_GO_gse %>% 
  filter(c1 > quantile(c1, 0.98) | c1 < quantile(c1, 0.02))

c2_cent_GO_gse_10 <- c2_cent_GO_gse %>% 
  filter(c2 > quantile(c2, 0.90) | c2 < quantile(c2, 0.13))

c3_cent_GO_gse_10 <- c3_cent_GO_gse %>% 
  filter(c3 > quantile(c3, 0.98) | c3 < quantile(c3, 0.02))

# Join and convert to matrix

gsea_GO_mat_mean <- merge(c3_cent_GO_gse_10, (merge(c1_cent_GO_gse_10, c2_cent_GO_gse_10, by = 0, all = TRUE) %>% tibble::column_to_rownames("Row.names")), by = 0, all = TRUE)

gsea_GO_mat_mean <- gsea_GO_mat_mean %>% 
  tibble::column_to_rownames("Row.names") %>% 
  as.matrix()

gsea_GO_mat_mean[is.na(gsea_GO_mat_mean)] <- 0

# Distances for dendofram

mean_row_hc <- hclust(dist(gsea_GO_mat_mean))
mean_col_hc <- hclust(dist(t(gsea_GO_mat_mean)))


# Heatmap

mean_heatmap_gsea <- gsea_GO_mat_mean %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("path") %>% 
  pivot_longer(cols = c("c1", "c2", "c3"), names_to = "clusters") %>% 
  ggplot(aes(x = clusters, y = path, fill = value)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) + 
  scale_fill_distiller(palette = "Spectral", direction = -1) +
  scale_x_discrete(expand = c(0, 0), labels = c("c1" = "C1", "c3" = "C3", "c2" = "C2")) +  
  scale_y_discrete(expand = c(0, 0)) +  
  theme_classic() +
  labs(x = "Clusters", y = "Pathway", fill = "NES", title = "GSEA GO means merged GSE21257 and GSE33382") +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 9),
    plot.title = element_text(size = 9),
    legend.position = "none"
  )


# Dendograms

mean_tree_right <- ggtree(mean_row_hc) + 
  scale_x_reverse() + 
  scale_y_continuous(expand = c(0, 0))

mean_tree_top <- ggtree(mean_col_hc, hang = -1) + 
  layout_dendrogram() + 
  scale_y_reverse(expand = c(0, 0))

# Full heatmap

mean_heatmap_gsea_gse <- 
  mean_heatmap_gsea %>% 
  insert_right(mean_tree_right, width = 0.1) %>% 
  insert_top(mean_tree_top, height = 0.1)

################################################################################

c1_cent_hallmark_gse_10 <- c1_cent_hallmark_gse %>% 
  filter(c1 > quantile(c1, 0.75) | c1 < quantile(c1, 0.75))

c2_cent_hallmark_gse_10 <- c2_cent_hallmark_gse %>% 
  filter(c2 > quantile(c2, 0.75) | c2 < quantile(c2, 0.75))

c3_cent_hallmark_gse_10 <- c3_cent_hallmark_gse %>% 
  filter(c3 > quantile(c3, 0.75) | c3 < quantile(c3, 0.75))

# Join and convert to matrix

mean_gsea_hallmark_mat <- merge(c3_cent_hallmark_gse_10, (merge(c1_cent_hallmark_gse_10, c2_cent_hallmark_gse_10, by = 0, all = TRUE) %>% tibble::column_to_rownames("Row.names")), by = 0, all = TRUE)

mean_gsea_hallmark_mat <- mean_gsea_hallmark_mat %>% 
  tibble::column_to_rownames("Row.names") %>% 
  as.matrix()

mean_gsea_hallmark_mat[is.na(mean_gsea_hallmark_mat)] <- 0


mean_row_hc_hm <- hclust(dist(mean_gsea_hallmark_mat))
mean_col_hc_hm <- hclust(dist(t(mean_gsea_hallmark_mat)))


mean_heatmap_gsea_hm <- mean_gsea_hallmark_mat %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("path") %>% 
  pivot_longer(cols = c("c1", "c2", "c3"), names_to = "clusters") %>% 
  ggplot(aes(x = clusters, y = path, fill = value)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) + 
  scale_fill_distiller(palette = "Spectral", direction = -1) +
  scale_x_discrete(expand = c(0, 0), labels = c("c1" = "C1", "c3" = "C3", "c2" = "C2")) +  
  scale_y_discrete(expand = c(0, 0)) +  
  theme_classic() +
  labs(x = "Clusters", y = "Pathway", fill = "NES", title = "GSEA Hallmark means merged GSE21257 and GSE33382") + 
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 7)
  )


mean_tree_right_hm <- ggtree(mean_row_hc_hm) + 
  scale_x_reverse() + 
  scale_y_continuous(expand = c(0, 0))

mean_tree_top_hm <- ggtree(mean_col_hc_hm, hang = -1) + 
  layout_dendrogram() + 
  scale_y_reverse(expand = c(0, 0))

mean_heatmap_gsea_hm_gse <- 
  mean_heatmap_gsea_hm %>% 
  insert_right(mean_tree_right_hm, width = 0.1) %>% 
  insert_top(mean_tree_top_hm, height = 0.1)

saveRDS(mean_heatmap_gsea_hm_gse   , "results/mul_hm/mean_heatmap_gsea_hm_gse.RDS")
saveRDS(mean_heatmap_gsea_gse      , "results/mul_hm/mean_heatmap_gsea_gse.RDS")
saveRDS(gse21257_plot_gse_go       , "results/mul_hm/gse21257_plot_gse_go.RDS")
saveRDS(gse21257_plot_gse_hm       , "results/mul_hm/gse21257_plot_gse_hm.RDS")
saveRDS(gse33382_plot_gse_go       , "results/mul_hm/gse33382_plot_gse_go.RDS")
saveRDS(gse33382_plot_gse_hm       , "results/mul_hm/gse33382_plot_gse_hm.RDS")


rm(list = ls())
gc()
