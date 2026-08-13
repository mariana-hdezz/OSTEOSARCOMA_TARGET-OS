metadata_33382_osteoblastic <- metadata_33382 %>% 
  filter(hist_sub_sim == "osteoblastic",
         !(clusters == 3))

table(metadata_33382$hist_sub_sim, metadata_33382$clusters)



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
library(aplot)


metadata_difex <- metadata_33382_osteoblastic

counts_data_difex <- counts_data_gse33382

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


contrast.matrix <- makeContrasts(clusters3 - clusters2,
                                 levels = design)


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

print(intersect(gene_signature_gse, rownames(res_sig)))


res$entrez <- mapIds(org.Hs.eg.db,
                     keys=rownames(res),
                     column="ENTREZID", 
                     keytype="SYMBOL", 
                     multiVals="first") 

res$name <-   mapIds(org.Hs.eg.db,
                     keys=rownames(res), 
                     column="GENENAME",
                     keytype="SYMBOL",
                     multiVals="first")

gene_symbols <- rownames(count_data)

res_lfc <- res$logFC

# 9.1 Assign ENSEMBL names to Log Fold object

names(res_lfc) <- rownames(res)

# 9.2 Eliminate NA

gene_list <- na.omit(res_lfc)


# 9.3 Descendant order

gene_list <- sort(gene_list, decreasing = TRUE)


# 9.4 Objeto GSEA de Gene ontology

gse <- gseGO(
  geneList = gene_list, 
  ont = "ALL", 
  keyType = "SYMBOL", 
  pvalueCutoff = 0.05, 
  OrgDb = "org.Hs.eg.db",
  minGSSize = 30, 
  maxGSSize = 100 
)

# 9.4.1 Observar como data frame

gse_df <- as.data.frame(gse)


msigdbr_collections()

m_df <- msigdbr(species = "Homo sapiens", category = "H")

msig_t2g <- m_df %>% dplyr::select(gs_name, gene_symbol)

gsea_res <- GSEA(
  geneList     = gene_list,
  TERM2GENE    = msig_t2g,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose      = FALSE
)

gsea_df <- as.data.frame(gsea_res)

# Create objects for heatmaps
c3_vs_c2_GO <- data.frame(row.names = gse$Description,
                          c3_vs_c2  = gse$NES)


c3_vs_c2_hallmark <- data.frame(row.names = gsea_res$Description, 
                                c3_vs_c2  = gsea_res$NES)


c3_vs_c2_GO[is.na(c3_vs_c2_GO)] <- 0

c3_vs_c2_hallmark[is.na(c3_vs_c2_hallmark)] <- 0


ggplot(c3_vs_c2_GO, aes(x = colnames(c3_vs_c2_GO), y = rownames(c3_vs_c2_GO), fill = c3_vs_c2)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Spectral") +
  theme_classic() + 
  scale_x_discrete(labels = "C3 vs C2") +
  labs(title = "GSEA between clusters Gene Ontology C3 vs C2.GSE33382", 
       x = "Clusters")



ggplot(c3_vs_c2_hallmark, aes(x = colnames(c3_vs_c2_hallmark), y = rownames(c3_vs_c2_hallmark), fill = c3_vs_c2)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Spectral") +
  theme_classic() + 
  scale_x_discrete(labels = "C3 vs C2") +
  labs(title = "GSEA between clusters Hallmarks of cancer C3 vs C2. GSE33382", 
       x = "Clusters")








