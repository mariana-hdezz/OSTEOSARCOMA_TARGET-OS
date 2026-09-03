
#############################################################################
#> Script to a fast view the results from differential expression and GSEA of TCGA-OS
#> 
#> This correspond to the script used for our comprehension of the 
#>  results with what we thought was needed at the time, but further analysis can 
#>  be conducted for whom presents interest with the following objects that are saved 
#>  in previous scripts.
#> 
#> 
#> Inputs:
#> ##> # This corresponds to the 37 gene signature:
##> 
##> gene_signature.csv
##> 
#>
#>  
#>  This corresponds to the overall results of the diffex in each comparison:
#>  
#>  res_1_vs_2.csv
#>  res_3_vs_2.csv
#>  res_3_vs_1.csv
##> 
##> 
##> This corresponds to the results after filtering by LFC and adj p value:
##> 
##> res_sig_1v2.csv
##> res_sig_3v2.csv
##> res_sig_3v1.csv
##> 
##> 
##> This are the results of the GSEA for each cluster comparison for both GO and 
##> Hallmark terms
##> 
##> gsea_df_GO3.vs.1.csv
##> gsea_df_GO3.vs.2.csv
##> gsea_df_GO1.vs.2.csv
##> gsea_df_HM3.vs.1.csv 
##> gsea_df_HM3.vs.2.csv 
##> gsea_df_HM1.vs.2.csv 
##> 
##> This are the results of the mean-based GSEA for each cluster for both GO
##> and Hallmark terms:
##> 
##> c1_cent_GO.csv
##> c2_cent_GO.csv
##> c3_cent_GO.csv
##> c1_cent_hallmark.csv
##> c2_cent_hallmark.csv
##> c3_cent_hallmark.csv
##> 
##> 
############################################################################


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


# ------ Analysis ------

# Total significant differentialy expressed genes in each cluster comparisson 
dim(res_sig_1v2)
dim(res_sig_3v2)
dim(res_sig_3v1)

# Significant genes for cluster 1 (logfoldchange  >0) when comparing vs cluster 2
res_sig_1v2$X[res_sig_1v2$log2FoldChange > 0]
sort(res_sig_1v2$log2FoldChange[res_sig_1v2$log2FoldChange > 0], decreasing = TRUE) 

# Significant genes for cluster 2 (logfoldchange <0) when comparing vs cluster 1
res_sig_1v2$X[res_sig_1v2$log2FoldChange < 0]
sort(res_sig_1v2$log2FoldChange[res_sig_1v2$log2FoldChange > 0], decreasing = TRUE) 

# Significant genes for cluster 3 (logfoldchange  >0) when comparing vs cluster 2
res_sig_3v2$X[res_sig_3v2$log2FoldChange > 0]

# Significant genes for cluster 2 (logfoldchange <0) when comparing vs cluster 3
res_sig_3v2$X[res_sig_3v2$log2FoldChange < 0]

# Significant genes for cluster 3 (logfoldchange >0) when comparing vs cluster 3
res_sig_3v1$X[res_sig_3v1$log2FoldChange > 0]

# Significant genes for cluster 1 (logfoldchange  >0) when comparing vs cluster 3
res_sig_3v1$X[res_sig_3v1$log2FoldChange > 0]


# View which of the selected genes in each comparison
intersect_1v2 <- intersect(res_sig_1v2$X, gene_signature)
intersect_3v2 <- intersect(res_sig_3v2$X, gene_signature)
intersect_3v1 <-intersect(res_sig_3v1$X, gene_signature)

# List
cat(intersect_1v2, sep = ", ")
cat(intersect_3v2, sep = ", ")
cat(intersect_3v1, sep = ", ")



# --------- 2-DATA LOADING FOR CLUSTER COMPARISSON GSEA GO RESULTS ------------

gsea_df_GO3.vs.1 <-  read_csv("./results/diffex_gsea_target/gsea_df_GO3_vs_1.csv")
gsea_df_GO3.vs.2 <-  read_csv("./results/diffex_gsea_target/gsea_df_GO3_vs_2.csv")
gsea_df_GO1.vs.2 <-  read_csv("./results/diffex_gsea_target/gsea_df_GO1_vs_2.csv")

# ------ Analysis ------
gsea_df_GO3.vs.1 %>%
  arrange(desc(NES)) %>%
  dplyr::select(Description, NES, p.adjust) %>%
  head(50)

gsea_df_GO3.vs.2 %>%
  arrange(desc(NES)) %>%
  dplyr::select(Description, NES, p.adjust) %>%
  head(50)

gsea_df_GO1.vs.2 %>%
  arrange(desc(NES)) %>%
  dplyr::select(Description, NES, p.adjust) %>%
  head(50)



# ------- 3-DATA LOADING FOR CLUSTER COMPARISSON GSEA HALLMARK RESULTS ---------

gsea_df_HM3.vs.1 <-  read_csv("./results/diffex_gsea_target/gsea_df_HM3_vs_1.csv")
gsea_df_HM3.vs.2 <-  read_csv("./results/diffex_gsea_target/gsea_df_HM3_vs_2.csv")
gsea_df_HM1.vs.2 <-  read_csv("./results/diffex_gsea_target/gsea_df_HM1_vs_2.csv")

# ------ Analysis ------
gsea_df_HM3.vs.1 %>%
  arrange(desc(NES)) %>%
  dplyr::select(Description, NES, p.adjust)

gsea_df_HM3.vs.2 %>%
  arrange(desc(NES)) %>%
  dplyr::select(Description, NES, p.adjust) 

gsea_df_HM1.vs.2 %>%
  arrange(desc(NES)) %>%
  dplyr::select(Description, NES, p.adjust)



# ------- 4-DATA LOADING FOR "MEAN-BASED" GSEA GO RESULTS ------------

c1_cent_GO <-  read_csv("./results/diffex_gsea_target/c1_cent_GO.csv")
c2_cent_GO <-  read_csv("./results/diffex_gsea_target/c2_cent_GO.csv")
c3_cent_GO <-  read_csv("./results/diffex_gsea_target/c3_cent_GO.csv")

# ------ Analysis ------
c1_cent_GO %>%
  rownames_to_column("pathway") %>%
  filter(c1 > 0) %>%
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


# ------- 5-DATA LOADING FOR "MEAN-BASED" GSEA HM RESULTS ------------

c1_cent_hallmark <-  read_csv("./results/diffex_gsea_target/c1_cent_hallmark.csv")
c2_cent_hallmark <-  read_csv("./results/diffex_gsea_target/c2_cent_hallmark.csv")
c3_cent_hallmark <-  read_csv("./results/diffex_gsea_target/c3_cent_hallmark.csv")

# ------ Analysis ------
# View top 30 (interchangeable) enriched genes for each cluster

c1_cent_hallmark %>%
  rownames_to_column("pathway") %>%
  filter(c1 > 0) %>%
  arrange(desc(c1)) 

c2_cent_hallmark %>%
  rownames_to_column("pathway") %>%
  filter(c2 > 0) %>%
  arrange(desc(c2)) 

c3_cent_hallmark %>%
  rownames_to_column("pathway") %>%
  filter(c3 > 0) %>%
  arrange(desc(c3)) 



rm(list = ls())
gc()