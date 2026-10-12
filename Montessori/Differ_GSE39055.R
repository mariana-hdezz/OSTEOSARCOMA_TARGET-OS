#############################################################################
#> Script to perform differential expression analysis on the clusters
#> assigned to GSE patients
#>
#> Note. it also runs GSE calling the script "4_Validation/4_2_Validation_tests/4_2_3_gsea_val.R"
#> inside the diff expression loop
#>

############################################################################

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
library(ggtree)
library(ggplot2)
library(aplot)
library(tidyr)


# Load data

metadata_gse39055 <- readRDS("output_data/metadata_gse39055.RDS")

counts_data_gse39055 <- readRDS("output_data/counts_data_gse39055.RDS")

gene_signature <- readRDS("results/boruta/boruta_conf_gse39055.RDS")
gene_signature <- unlist(gene_signature)


metadata_difex <- metadata_gse39055

counts_data_difex <- counts_data_gse39055


# 1.1 Generate column corresponding to clusters nodes, those that have 0 in one group and those with more than 0 in another

col_data <- metadata_difex %>%
  dplyr::select(geo_accession, clusters) %>% # Create only the object to use for Limma
  column_to_rownames("geo_accession")



# 2.- Differential expression -----------------------------------------------

# 2.2 Data counts of the patients that had clusters node information in the metadata

count_data <- counts_data_difex[, colnames(counts_data_difex) %in% rownames(col_data)]

# 2.2.2 Making shure they are in the same order

count_data <- count_data[,match(rownames(col_data), colnames(count_data))]

all(colnames(count_data) == rownames(col_data))


# 2.4 Design based on object separating on clusters nodes

col_data$clusters <- as.factor(col_data$clusters)

design <- model.matrix( ~ 0 + clusters, data = col_data)

# 2.4.2 Asign make.names objects as colnames

colnames(design) <- make.names(colnames(design))

# 2.5 Fit

fit <- lmFit(count_data, design)


# 2.5.2 Contrast matrix comparing clusters 0 to > 0 clusters

contrast.matrix <- makeContrasts(clusters1 - clusters2, levels = design)


contrast.matrix <- makeContrasts(clusters3 - clusters2, levels = design)


contrast.matrix <- makeContrasts(clusters3 - clusters1, levels = design)




# 2.5.3 Fit based on contrasts

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit)

# 2.6 Results

res <- topTable(fit2, coef = 1, number = Inf)

# 2.6.2 Results that correspond to a signfiicant p value and log fold change

res_sig <- res %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1.5) # 0.1


