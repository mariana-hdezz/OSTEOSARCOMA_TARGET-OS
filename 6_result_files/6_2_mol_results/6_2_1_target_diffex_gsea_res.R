
#------------------- TARGET-OS RESULTS EASY ACCESS -----------------------

# ------------- 1-DATA LOADING FOR DIFFERENTIAL EXPRESSION RESULTS ------------

# Load 37 gene signature
gene_signature <- scan("./output_data/gene_signature.csv", sep = ",", what = character())

# Load CSV with results of the differential expression

res_1_vs_2 <-  read.csv("./results/diffex_gsea_target/res_1_vs_2.csv")
res_3_vs_2 <-   read.csv("./results/diffex_gsea_target/res_3_vs_2.csv")
res_3_vs_1 <-   read.csv("./results/diffex_gsea_target/res_3_vs_1.csv")

# Load CSV with results of the differential expression filtered by p-values

res_sig_1v2  <-  read.csv("./results/diffex_gsea_target/res_sig_1v2.csv")
res_sig_3v2  <-  read.csv("./results/diffex_gsea_target/res_sig_3v2.csv")
res_sig_3v1  <-  read.csv("./results/diffex_gsea_target/res_sig_3v1.csv")


# ANALYSIS 

# Differentialy expressed genes in cluster 3 utilizing cluster 1 as reference
dim(res_sig_3v1) # View the total of diffex genes

res_sig_3v1 %>% 
  filter(log2FoldChange >= 0)

res_sig_3v1$X[res_sig_3v1$log2FoldChange > 0]


sort(res_sig_3v1$log2FoldChange[res_sig_3v1$log2FoldChange > 0], decreasing = TRUE)

res_view <- res_sig_3v1 %>% 
  arrange(desc(by = log2FoldChange))

res_sig_3v1$X[res_sig_3v1$log2FoldChange < 0]


genes_inter <- intersect(res_sig_3v1$X, gene_signature )

cat(genes_inter3, sep = ", ")

genes_inter2 <- res_sig_3v1[res_sig_3v1$X %in% gene_signature, ]

genes_inter3 <- genes_inter2$X[genes_inter2$log2FoldChange < 0]




# --------- 2-DATA LOADING FOR CLUSTER COMPARISSON GSEA GO RESULTS ------------

gsea_df_GO3.vs.1 <-  read_csv("./results/diffex_gsea_target/gsea_df_GO3_vs_1.csv")
gsea_df_GO3.vs.2 <-  read_csv("./results/diffex_gsea_target/gsea_df_GO3_vs_2.csv")
gsea_df_GO1.vs.2 <-  read_csv("./results/diffex_gsea_target/gsea_df_GO1_vs_2.csv")


# ------- 3-DATA LOADING FOR CLUSTER COMPARISSON GSEA HALLMARK RESULTS ---------

gsea_df_HM3.vs.1 <-  read_csv("./results/diffex_gsea_target/gsea_df_HM3_vs_1.csv")
gsea_df_HM3.vs.2 <-  read_csv("./results/diffex_gsea_target/gsea_df_HM3_vs_2.csv")
gsea_df_HM1.vs.2 <-  read_csv("./results/diffex_gsea_target/gsea_df_HM1_vs_2.csv")


# ------- 4-DATA LOADING FOR "MEAN-BASED" GSEA GO RESULTS ------------

c1_cent_GO <-  read_csv("./results/diffex_gsea_target/c1_cent_GO.csv")
c2_cent_GO <-  read_csv("./results/diffex_gsea_target/c2_cent_GO.csv")
c3_cent_GO <-  read_csv("./results/diffex_gsea_target/c3_cent_GO.csv")


# ------- 5-DATA LOADING FOR "MEAN-BASED" GSEA HM RESULTS ------------

c1_cent_hallmark <-  read_csv("./results/diffex_gsea_target/c1_cent_hallmark.csv")
c2_cent_hallmark <-  read_csv("./results/diffex_gsea_target/c2_cent_hallmark.csv")
c3_cent_hallmark <-  read_csv("./results/diffex_gsea_target/c3_cent_hallmark.csv")



# ANALYSIS

c1_cent_GO %>%
  rownames_to_column("pathway") %>%
  filter(c1 > 0) %>%
  arrange(desc(c1)) %>%
  head(50)

c2_cent_GO %>%
  rownames_to_column("pathway") %>%
  filter(c2 > 0) %>%
  arrange(desc(c2)) %>%
  head(50)

c3_cent_GO %>%
  rownames_to_column("pathway") %>%
  filter(c3 > 0) %>%
  arrange(desc(c3)) %>%
  head(50)





# ANALYSIS

# View top 30 (interchangeable) enriched genes for each cluster

c1_cent_hallmark %>%
  rownames_to_column("pathway") %>%
  filter(c1 > 0) %>%
  arrange(desc(c1)) %>%
  head(30)

c2_cent_hallmark %>%
  rownames_to_column("pathway") %>%
  filter(c2 > 0) %>%
  arrange(desc(c2)) %>%
  head(30)

c3_cent_hallmark %>%
  rownames_to_column("pathway") %>%
  filter(c3 > 0) %>%
  arrange(desc(c3)) %>%
  head(30)