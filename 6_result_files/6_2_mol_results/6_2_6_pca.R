library(ggplot2)
library(dplyr)
library(factoextra)
library(gridExtra)
library(sigPCA)

counts <- readRDS("output_data/vst_counts.RDS")

metadata <- readRDS("output_data/metadata_os.RDS")

gene_signature <- scan("output_data/gene_signature.csv", sep = ",", what = character()) 

counts <- counts[rownames(counts) %in% gene_signature, ]

counts <- t(counts)

pca <- prcomp(x = counts)

pca_full <- pca$x

pca_plot <- merge(pca_full, metadata %>% dplyr::select(clusters, sample) %>% tibble::column_to_rownames("sample"), by = 0) %>% 
  tibble::column_to_rownames("Row.names") %>% 
  mutate(clusters = factor(clusters))

eigens <- get_eigenvalue(pca)

sig_pc <- sigPCA::sigPCA(counts)

pca_plot <- 
  pca_plot %>% 
  ggplot(aes(x = PC1, y = PC2, colour = clusters)) + 
  geom_point(size = 5) + 
  labs(
    x = paste0("PC1 ",  round(eigens[1, 1], 2), "%"),
    y =   paste0("PC2 ", round(eigens[2, 1], 2), "%"),
    title = "Principal component analysis of the gene signature"
    ) +
  theme_classic(base_size = 15)


a <- fviz_contrib(pca, choice = "var", axes = 1, top = 30, fill = "#ff89d4", color = "#ff89d4") 


b <- fviz_contrib(pca, choice = "var", axes = 2, top = 30, fill = "#c380d3", color = "#c380d3") 


components <- 
  grid.arrange(a, b, ncol = 2,
             top = paste0("Top Gene Contributors to PC1 and PC2 from gene _sign"))




pca_plot / components
