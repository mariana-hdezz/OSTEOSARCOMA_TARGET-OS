library(ggplot2)
library(dplyr)
library(tibble)
library(factoextra)
library(gridExtra)
library(sigPCA)
library(RColorBrewer)


vst_counts <- readRDS("output_data/vst_counts.RDS")

metadata <- readRDS("output_data/metadata_os.RDS")

gene_signature <- scan("output_data/gene_signature.csv", sep = ",", what = character()) 

vst_counts <- vst_counts[rownames(vst_counts) %in% gene_signature, ]

vst_counts <- t(vst_counts)

pca <- prcomp(x = vst_counts)

pca_full <- pca$x

pca_plot <- merge(pca_full, metadata %>% dplyr::select(clusters, sample) %>% tibble::column_to_rownames("sample"), by = 0) %>% 
  tibble::column_to_rownames("Row.names") %>% 
  mutate(clusters = factor(clusters))

eigens <- get_eigenvalue(pca)

sig_pc <- sigPCA::sigPCA(vst_counts)

pca_plot <- 
  pca_plot %>% 
  ggplot(aes(x = PC1, y = PC2, colour = clusters)) + 
  geom_point(size = 5) + 
  labs(
    x = paste0("PC1 ",  round(eigens[1, 1], 2), "%"),
    y =   paste0("PC2 ", round(eigens[2, 1], 2), "%"),
    title = "Principal component analysis of the gene signature in TARGET-OS"
  ) +
  theme_classic(base_size = 18) + 
  scale_color_discrete(palette = brewer.pal(n = 12, name = "Set3")[c(8, 5, 10)])


a <- fviz_contrib(pca, choice = "var", axes = 1, top = 30, fill = "#ff89d4", color = "#ff89d4") 


b <- fviz_contrib(pca, choice = "var", axes = 2, top = 30, fill = "#c380d3", color = "#c380d3") 


components <- 
  grid.arrange(a, b, ncol = 2,
               top = paste0("Top Gene Contributors to PC1 and PC2 from gene signature in TARGET-OS"))


batch_counts <- readRDS("output_data/counts_batch.RDS")
metad_gse33382 <- readRDS("output_data/metadata_gse33382_for_merge.RDS")
metad_gse21257 <- readRDS("output_data/metadata_gse21257_for_merge.RDS")
gene_signature_gse <- scan("output_data/gene_signature_gse.csv", sep = ",", what = character()) 


metad_com <- bind_rows(metad_gse33382,
                       metad_gse21257)

batch_counts <- batch_counts[rownames(batch_counts) %in% gene_signature_gse, ]

batch_counts <- t(batch_counts)

pca_gse <- prcomp(x = batch_counts)

pca_gse_full <- pca_gse$x

pca_gse_plot <- merge(pca_gse_full, metad_com %>% dplyr::select(clusters, geo_accession) %>% tibble::column_to_rownames("geo_accession"), by = 0) %>% 
  tibble::column_to_rownames("Row.names") %>% 
  mutate(clusters = factor(clusters))

eigens <- get_eigenvalue(pca_gse)

sig_pc <- sigPCA::sigPCA(batch_counts)

pca_gse_plot <- 
  pca_gse_plot %>% 
  ggplot(aes(x = PC1, y = PC2, colour = clusters)) + 
  geom_point(size = 5) + 
  labs(
    x = paste0("PC1 ",  round(eigens[1, 1], 2), "%"),
    y =   paste0("PC2 ", round(eigens[2, 1], 2), "%"),
    title = "Principal component analysis of the gene signature in validation set"
  ) +
  theme_classic(base_size = 18) + 
  scale_color_discrete(palette = brewer.pal(n = 12, name = "Set3")[c(8, 5, 10)])



a_gse <- fviz_contrib(pca_gse, choice = "var", axes = 1, top = 30, fill = "#ff89d4", color = "#ff89d4") 


b_gse <- fviz_contrib(pca_gse, choice = "var", axes = 2, top = 30, fill = "#c380d3", color = "#c380d3") 


components_gse <- 
  grid.arrange(a_gse, b_gse, ncol = 2,
               top = paste0("Top Gene Contributors to PC1 and PC2 from signature in validation set"))



((pca_plot / components) | (pca_gse_plot / components_gse)) + 
  patchwork::plot_annotation(tag_levels = "A")
