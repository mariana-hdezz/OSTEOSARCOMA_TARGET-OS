library(oligo)
library(GEOquery)
library(dplyr)
library(limma)
library(hgu133a.db)

# getGEOSuppFiles(GEO = "GSE21257", makeDirectory = TRUE, baseDir = "d:/")



# getGEOSuppFiles(GEO = "GSE33382", makeDirectory = TRUE, baseDir = "d:/")

untar("D:/GSE33382/GSE33382_RAW.tar", exdir = "d:/GSE33382/")

gpl <- getGEO("GPL10295", AnnotGPL = TRUE)

annot <- Table(gpl) %>% 
  dplyr::select(ID,
                Symbol,
                Illumina_ProbeID,
                Illumina_Gene)

gset <- getGEO("GSE21257", GSEMatrix = TRUE)

pheno <- pData(gset[[1]])

unique(pheno$data_processing)

counts_data_norm <- exprs(gset[[1]]) %>% 
  as.data.frame()

counts_data_norm <- counts_data_norm[rownames(counts_data_norm) %in% annot$ID, ]

annot <- annot[annot$ID %in% rownames(counts_data_norm), ]

counts_data <- counts_data_norm %>% 
  rownames_to_column("ID") %>% 
  left_join(annot, by = "ID") %>% 
  select(- ID) 

counts_data$var <- apply(counts_data[, !(colnames(counts_data) %in% "Symbol")], 1, function(i){var(i, na.rm = TRUE)})

counts_data <- counts_data[!is.na(counts_data$Symbol) & counts_data$Symbol != "", ]

counts_data <- 
  counts_data %>% 
  group_by(Symbol) %>% 
  slice_max(order_by = var, n = 1, with_ties = FALSE) %>% 
  ungroup() %>% 
  select(- var) %>% 
  column_to_rownames("Symbol")


dist_counts <- get_dist(t(counts_data), method = "manhattan")

# Clustering

hc_counts <- hclust(dist_counts, method = "ward.D2")


# Clustering characteristics ------------------------------------------

# Dendrogram

plot(hc_counts,
     labels = FALSE,
     hang = -1,
     main = "CLUSTERING JERÁRQUICO TARGET-OS")

rect.hclust(hc_counts, k = 3, border = "purple") # bottom up approach

# Tree of clusters

clusters <- cutree(hc_counts, k = 2)

sil <- silhouette(clusters, dist_counts)

mean(sil[, "sil_width"])
