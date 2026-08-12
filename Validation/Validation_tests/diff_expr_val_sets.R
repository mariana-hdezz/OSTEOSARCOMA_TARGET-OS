library(limma)
library(dplyr)
library(tibble)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggridges)
library(enrichplot)
library(msigdbr)
library(ComplexHeatmap)
library(circlize)


metadata_difex <- metadata_gse21257

counts_data_difex <- counts_data_gse21257

# 1.1 Generate column corresponding to clusters nodes, those that have 0 in one group and those with more than 0 in another

col_data <- metadata_difex %>% 
  dplyr::select(geo_accession, clusters) %>% # Create only the object to use for Limma
  column_to_rownames("geo_accession")




# 2.- Differential expression -----------------------------------------------


# 2.2 Data counts of the patients that had clusters node information in the metadata

count_data <- counts_data_difex[colnames(counts_data_difex) %in% rownames(col_data)]

# 2.2.2 Making shure they are in the same order

count_data <- count_data[match(rownames(col_data), colnames(count_data))]

all(colnames(count_data) == rownames(col_data))


# 2.4 Design based on object separating on clusters nodes

design <- model.matrix(~ 0 + clusters, data = col_data)

# 2.4.2 Asign make.names objects as colnames

colnames(design) <- make.names(colnames(design)) 

# 2.5 Fit

fit <- lmFit(count_data, design)


for (i in 1:ncol(design)) {
  fit <- lmFit(count_data, design)
  if(i == 1){
    # 2.5.2 Contrast matrix comparing clusters 0 to > 0 clusters
    
    contrast.matrix <- makeContrasts(clusters1 - clusters2,
                                     levels = design)
    
  }else if(i == 2){
    
    contrast.matrix <- makeContrasts(clusters3 - clusters2,
                                     levels = design)
  }else if(i == 3){
    
    contrast.matrix <- makeContrasts(clusters3 - clusters1,
                                     levels = design)
  }




# 2.5.3 Fit based on contrasts

fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

topTable(fit)

# 2.6 Results

res <- topTable(fit, coef = 1, number = Inf)

# 2.6.2 Results that correspond to a signfiicant p value and log fold change

res_sig <- res %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1.5) # 0.1
print(res_sig)

print(intersect(pucky_gse, rownames(res_sig)))


source("Validation/Validation_tests/gsea_val_sets.R")


}



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

row_hc_hist_GO <- hclust(dist(gsea_GO_mat))
col_hc_hist_GO <- hclust(dist(t(gsea_GO_mat)))

gsea_GO_mat <-
  gsea_GO_mat %>% 
  as.data.frame() %>% 
    tibble::rownames_to_column("path") %>% 
    pivot_longer(cols = c("c3_vs_c2", "c1_vs_c2", "c3_vs_c1"), names_to = "clusters") 


heat_hist_GO <- gsea_GO_mat %>% 
  ggplot(aes(x = clusters, y = path, fill = value)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) + 
  scale_fill_distiller(palette = "Spectral", direction = -1) + 
  theme_classic() +
  labs(x = "Clusters", y = "Pathways", fill = "NES", title = "GSEA between clusters Gene Ontology")

tree_right_hist_GO <- ggtree(row_hc_hist_GO, hang = -1) + 
  scale_x_reverse() + 
  scale_y_continuous(expand = c(0, 0))

tree_top_hist_GO <- ggtree(col_hc_hist_GO, hang = -1) + 
  layout_dendrogram() + 
  scale_y_reverse(expand = c(0, 0))

heat_hist_GO %>% 
  insert_right(tree_right_hist_GO, width = 0.1) %>% 
  insert_top(tree_top_hist_GO, height = 0.1)



##############################################################################

c3_vs_c2_hallmark_10 <- c3_vs_c2_hallmark %>% 
  filter(c3_vs_c2 > quantile(c3_vs_c2, 0.85) | c3_vs_c2 < quantile(c3_vs_c2, 0.85))

c1_vs_c2_hallmark_10 <- c1_vs_c2_hallmark %>% 
  filter(c1_vs_c2_hallmark > quantile(c1_vs_c2, 0.85) | c1_vs_c2_hallmark < quantile(c1_vs_c2, 0.85))

c3_vs_c1_hallmark_10 <- c3_vs_c1_hallmark %>% 
  filter(c3_vs_c1 > quantile(c3_vs_c1, 0.85) | c3_vs_c1 < quantile(c3_vs_c1, 0.85))


gsea_hallmark_mat <- merge(c3_vs_c2_hallmark_10, (merge(c1_vs_c2_hallmark_10, c3_vs_c1_hallmark_10, by = 0, all = TRUE) %>% column_to_rownames("Row.names")), by = 0, all = TRUE)

gsea_hallmark_mat <- gsea_hallmark_mat %>% 
  column_to_rownames("Row.names") %>% 
  as.matrix()

gsea_hallmark_mat[is.na(gsea_hallmark_mat)] <- 0

row_hc_hist_hallmark <- hclust(dist(gsea_hallmark_mat))
col_hc_hist_hallmark <- hclust(dist(t(gsea_hallmark_mat)))


gsea_hallmark_mat <-
  gsea_hallmark_mat %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("path") %>% 
  pivot_longer(cols = c("c3_vs_c2", "c1_vs_c2", "c3_vs_c1"), names_to = "clusters") 




heat_hist_hallmark <- gsea_hallmark_mat %>% 
  ggplot(aes(x = clusters, y = path, fill = value)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) + 
  scale_fill_distiller(palette = "Spectral", direction = -1) + 
  theme_classic() +
  labs(x = "Clusters", y = "Pathway", fill = "NES", title = "GSEA between clusters Gene Ontology")

tree_right_hist_hallmark <- ggtree(row_hc_hist_hallmark, hang = -1) + 
  scale_x_reverse() + 
  scale_y_continuous(expand = c(0, 0))

tree_top_hist_hallmark <- ggtree(col_hc_hist_hallmark, hang = -1) + 
  layout_dendrogram() + 
  scale_y_reverse(expand = c(0, 0))

heat_hist_hallmark %>% 
  insert_right(tree_right_hist_hallmark, width = 0.1) %>% 
  insert_top(tree_top_hist_hallmark, height = 0.1)
