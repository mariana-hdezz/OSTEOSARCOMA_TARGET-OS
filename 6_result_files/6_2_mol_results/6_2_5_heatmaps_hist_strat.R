

diffex_telan <- list.files(path = "results/diffex_gsea_gse/", pattern = "^telangiectatic", full.names = TRUE)
diffex_osteo <- list.files(path = "results/diffex_gsea_gse/", pattern = "^osteoblastic", full.names = TRUE)
diffex_chond <- list.files(path = "results/diffex_gsea_gse/", pattern = "^chondroblastic", full.names = TRUE)
gsea_telan <- list.files(path = "results/diffex_gsea_gse/", pattern = "telangiectatic.csv$", full.names = TRUE)
gsea_osteo <- list.files(path = "results/diffex_gsea_gse/", pattern = "osteoblastic.csv$", full.names = TRUE)
gsea_chond <- list.files(path = "results/diffex_gsea_gse/", pattern = "chondroblastic.csv$", full.names = TRUE)


files_diffex <- 
  c(
  diffex_telan,
  diffex_osteo,
  diffex_chond
)

files_gsea <- 
  c(
    gsea_telan, 
    gsea_osteo, 
    gsea_chond 
  )

names(files_diffex) <-  gsub("^results/diffex_gsea_gse/(.*)_\\.csv$", "\\1", files_diffex)

names(files_gsea) <-  gsub("^results/diffex_gsea_gse/(.*)\\.csv$", "\\1", files_gsea)


diffex <- lapply(files_diffex, function(i){
  
  read.csv(i) %>% 
    column_to_rownames("X")

})


gsea <- lapply(files_gsea, function(i){
  
  read.csv(i) %>% 
    column_to_rownames("X")
  
})


for (i in seq_along(diffex)){
  name <- names(diffex)[[i]]
  obj <- diffex[[i]]
  assign(name, obj)
  }


for (i in seq_along(gsea)){
  name <- names(gsea)[[i]]
  obj <- gsea[[i]]
  assign(name, obj)
}



c3vc1_telangiectatic <-
  data.frame(path = gsea_c3_vs_c1_GO_telangiectatic$Description,
             c3vc1_telangiectatic = gsea_c3_vs_c1_GO_telangiectatic$NES)

c1vc2_osteoblastic <-
  data.frame(path = gsea_c1_vs_c2_GO_osteoblastic$Description,
             c1vc2_osteoblastic = gsea_c1_vs_c2_GO_osteoblastic$NES)

c3vc1_osteoblastic <-
  data.frame(path = gsea_c3_vs_c1_GO_osteoblastic$Description,
             c3vc1_osteoblastic = gsea_c3_vs_c1_GO_osteoblastic$NES)

c3vc2_osteoblastic <-
  data.frame(path = gsea_c3_vs_c2_GO_osteoblastic$Description,
             c3vc2_osteoblastic = gsea_c3_vs_c2_GO_osteoblastic$NES)

c1vc2_chondroblastic <-
  data.frame(path = gsea_c1_vs_c2_GO_chondroblastic$Description,
             c1vc2_chondroblastic = gsea_c1_vs_c2_GO_chondroblastic$NES)

c3vc1_chondroblastic <-
  data.frame(path = gsea_c3_vs_c1_GO_chondroblastic$Description,
             c3vc1_chondroblastic = gsea_c3_vs_c1_GO_chondroblastic$NES)

c3vc2_chondroblastic <-
  data.frame(path = gsea_c3_vs_c2_GO_chondroblastic$Description,
             c3vc2_chondroblastic = gsea_c3_vs_c2_GO_chondroblastic$NES)


object_name <- c("c1vc2_osteoblastic", "c3vc1_osteoblastic", "c3vc2_osteoblastic", "c1vc2_chondroblastic", "c3vc1_chondroblastic", "c3vc2_chondroblastic")

merge1 <- c3vc1_telangiectatic

for (i in seq_along(object_name)) {
  
  obj <- get(object_name[i])
  
  merge1 <- merge(merge1, obj, by = "path", all = TRUE)
  
  
}

merge1[is.na(merge1)] <- 0




# Calculate distances for dendogram

row_go <- hclust(dist(merge1 %>% 
                        column_to_rownames("path")))
col_go <- hclust(dist(t(merge1 %>% 
                          column_to_rownames("path"))))

# Dendograms

tree_right_go <- ggtree(row_go) + 
  scale_x_reverse() + 
  scale_y_continuous(expand = c(0, 0))

tree_top_go <- ggtree(col_go, hang = -1) + 
  layout_dendrogram() + 
  scale_y_reverse(expand = c(0, 0))


heat_obj_hist <- merge1 %>% 
  pivot_longer(
    cols = c("c3vc1_telangiectatic", all_of(object_name)),
    names_to = "comparison", 
    values_to = "NES"
  )


heatmap_hist <- 
  heat_obj_hist %>% 
  ggplot(aes(x = comparison, y = path, fill = NES)) +
  geom_tile() +
  theme_classic() + 
  scale_fill_distiller(palette = "Spectral")  
  
  
heatmap_hist%>% 
  insert_right(tree_right_go, width = 0.1) %>% 
  insert_top(tree_top_go, height = 0.1)






c3_vs_c1_hm_telangiectatic <-
  data.frame(path = gsea_c3_vs_c1_hm_telangiectatic$Description,
             c3_vs_c1_hm_telangiectatic = gsea_c3_vs_c1_hm_telangiectatic$NES)

c1_vs_c2_hm_osteoblastic <-
  data.frame(path = gsea_c1_vs_c2_hm_osteoblastic$Description,
             c1_vs_c2_hm_osteoblastic = gsea_c1_vs_c2_hm_osteoblastic$NES)

c3_vs_c1_hm_osteoblastic <-
  data.frame(path = gsea_c3_vs_c1_hm_osteoblastic$Description,
             c3_vs_c1_hm_osteoblastic = gsea_c3_vs_c1_hm_osteoblastic$NES)

c3_vs_c2_hm_osteoblastic <-
  data.frame(path = gsea_c3_vs_c2_hm_osteoblastic$Description,
             c3_vs_c2_hm_osteoblastic = gsea_c3_vs_c2_hm_osteoblastic$NES)

c1_vs_c2_hm_chondroblastic <-
  data.frame(path = gsea_c1_vs_c2_hm_chondroblastic$Description,
             c1_vs_c2_hm_chondroblastic = gsea_c1_vs_c2_hm_chondroblastic$NES)

c3_vs_c1_hm_chondroblastic <-
  data.frame(path = gsea_c3_vs_c1_hm_chondroblastic$Description,
             c3_vs_c1_hm_chondroblastic = gsea_c3_vs_c1_hm_chondroblastic$NES)

c3_vs_c2_hm_chondroblastic <-
  data.frame(path = gsea_c3_vs_c2_hm_chondroblastic$Description,
             c3_vs_c2_hm_chondroblastic = gsea_c3_vs_c2_hm_chondroblastic$NES)




object_name2 <- c("c1_vs_c2_hm_osteoblastic", "c3_vs_c1_hm_osteoblastic", "c3_vs_c2_hm_osteoblastic", "c1_vs_c2_hm_chondroblastic", "c3_vs_c1_hm_chondroblastic", "c3_vs_c2_hm_chondroblastic")

merge2 <- c3_vs_c1_hm_telangiectatic

for (i in seq_along(object_name2)) {
  
  obj <- get(object_name2[i])
  
  merge2 <- merge(merge2, obj, by = "path", all = TRUE)
  
  
}

merge2[is.na(merge2)] <- 0




# Calculate distances for dendogram

row_go <- hclust(dist(merge2 %>% 
                        column_to_rownames("path")))
col_go <- hclust(dist(t(merge2 %>% 
                          column_to_rownames("path"))))

# Dendograms

tree_right_go <- ggtree(row_go) + 
  scale_x_reverse() + 
  scale_y_continuous(expand = c(0, 0))

tree_top_go <- ggtree(col_go, hang = -1) + 
  layout_dendrogram() + 
  scale_y_reverse(expand = c(0, 0))


heat_obj_hist <- merge2 %>% 
  pivot_longer(
    cols = c("c3_vs_c1_hm_telangiectatic", all_of(object_name2)),
    names_to = "comparison", 
    values_to = "NES"
  )


heatmap_hist <- 
  heat_obj_hist %>% 
  ggplot(aes(x = comparison, y = path, fill = NES)) +
  geom_tile() +
  theme_classic() + 
  scale_fill_distiller(palette = "Spectral")  


heatmap_hist%>% 
  insert_right(tree_right_go, width = 0.1) %>% 
  insert_top(tree_top_go, height = 0.1)






rm(list = ls())
gc()

