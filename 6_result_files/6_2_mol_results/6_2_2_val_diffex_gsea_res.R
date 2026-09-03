#############################################################################
#> Script to observe results from differential expression and GSEA of the validation sets
#> 
#> Inputs:
#> 
#>  This corresponds to the full diff expression result. It is not
#>  further analysed on thsi script but we leave it read for anybody interested in
#>  analyzing it
#> 
##> gse21257_res_1v2_.csv
##> gse21257_res_3v1_.csv
##> gse21257_res_3v2_.csv
##> gse33382_res_1v2_.csv
##> gse33382_res_3v1_.csv
##> gse33382_res_3v2_.csv
##> 
##> This corresponds to the results after filtering by LFC and adj p value
##> 
##> gse21257_res_sig_1v2_.csv
##> gse21257_res_sig_3v1_.csv
##> gse21257_res_sig_3v2_.csv
##> gse33382_res_sig_1v2_.csv
##> gse33382_res_sig_3v1_.csv
##> gse33382_res_sig_3v2_.csv
##> 
##> This are the results of the GSEA fot each comparison for both GO and 
##> Hallmark terms
##> 
##> gsea_c1_vs_c2_GO_gse21257.csv
##> gsea_c1_vs_c2_GO_gse33382.csv
##> gsea_c1_vs_c2_hm_gse21257.csv
##> gsea_c1_vs_c2_hm_gse33382.csv
##> gsea_c3_vs_c1_GO_gse21257.csv
##> gsea_c3_vs_c1_GO_gse33382.csv
##> gsea_c3_vs_c1_hm_gse21257.csv
##> gsea_c3_vs_c1_hm_gse33382.csv
##> gsea_c3_vs_c2_GO_gse21257.csv
##> gsea_c3_vs_c2_GO_gse33382.csv
##> gsea_c3_vs_c2_hm_gse21257.csv
##> gsea_c3_vs_c2_hm_gse33382.csv
##> 
##> 
##> gene_signature_gse
############################################################################

library(dplyr)

# Load data

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



res_c1_vs_c2_gse21257 <- read.csv("results/diffex_gsea_gse/gse21257_res_1v2_.csv") %>% 
  column_to_rownames("X")
res_c1_vs_c2_gse33382 <- read.csv("results/diffex_gsea_gse/gse33382_res_1v2_.csv") %>% 
  column_to_rownames("X")
res_c3_vs_c2_gse21257 <- read.csv("results/diffex_gsea_gse/gse21257_res_3v2_.csv") %>% 
  column_to_rownames("X")
res_c3_vs_c2_gse33382 <- read.csv("results/diffex_gsea_gse/gse33382_res_3v2_.csv") %>% 
  column_to_rownames("X")
res_c3_vs_c1_gse21257 <- read.csv("results/diffex_gsea_gse/gse21257_res_3v1_.csv") %>% 
  column_to_rownames("X")
res_c3_vs_c1_gse33382 <- read.csv("results/diffex_gsea_gse/gse33382_res_3v1_.csv") %>% 
  column_to_rownames("X")



res_sig_c1_vs_c2_gse21257 <- read.csv("results/diffex_gsea_gse/gse21257_res_sig_1v2_.csv") %>% 
  column_to_rownames("X")
res_sig_c1_vs_c2_gse33382 <- read.csv("results/diffex_gsea_gse/gse33382_res_sig_1v2_.csv") %>% 
  column_to_rownames("X")
res_sig_c3_vs_c2_gse21257 <- read.csv("results/diffex_gsea_gse/gse21257_res_sig_3v2_.csv") %>% 
  column_to_rownames("X")
res_sig_c3_vs_c2_gse33382 <- read.csv("results/diffex_gsea_gse/gse33382_res_sig_3v2_.csv") %>% 
  column_to_rownames("X")
res_sig_c3_vs_c1_gse21257 <- read.csv("results/diffex_gsea_gse/gse21257_res_sig_3v1_.csv") %>% 
  column_to_rownames("X")
res_sig_c3_vs_c1_gse33382 <- read.csv("results/diffex_gsea_gse/gse33382_res_sig_3v1_.csv") %>% 
  column_to_rownames("X")


gene_signature_gse <- scan("output_data/gene_signature_gse.csv", sep = ",", what = character()) 

# Genes that where found as differentially expressed in both cohorts

common_1v2 <- intersect(rownames(res_sig_c1_vs_c2_gse21257), rownames(res_sig_c1_vs_c2_gse33382 ))

common_3v2 <- intersect(rownames(res_sig_c3_vs_c2_gse21257), rownames(res_sig_c3_vs_c2_gse33382 ))

common_3v1 <- intersect(rownames(res_sig_c3_vs_c1_gse21257), rownames(res_sig_c3_vs_c1_gse33382 ))

# Results of those genes for GSE21257, only 1 is shown since the results are similar to GSE33382
 
## C1 VS C2 (C2 is negative LFC)

res_sig_c1_vs_c2_gse21257[common_1v2, ]
res_sig_c1_vs_c2_gse33382[common_1v2, ]
common_1v2[common_1v2 %in% gene_signature_gse]

## C3 VS C2 (C2 is negative LFC)

res_sig_c3_vs_c2_gse21257[common_3v2, ]
res_sig_c3_vs_c2_gse33382[common_3v2, ]
common_3v2[common_3v2 %in% gene_signature_gse]

## C3 VS 1C (C1 is negative LFC)

res_sig_c3_vs_c1_gse21257[common_3v1, ] %>% 
  arrange(desc(logFC))

res_sig_c3_vs_c1_gse33382[common_3v1, ]%>% 
  arrange(desc(logFC))

common_3v1[common_3v1 %in% gene_signature_gse]

## C1 VS C2 (C2 is negative NES)

gsea_c1_vs_c2_GO_gse21257[gsea_c1_vs_c2_GO_gse21257$ID %in% intersect(gsea_c1_vs_c2_GO_gse21257$ID, gsea_c1_vs_c2_GO_gse33382$ID), ]
gsea_c1_vs_c2_hm_gse21257[gsea_c1_vs_c2_hm_gse21257$ID %in% intersect(gsea_c1_vs_c2_hm_gse21257$ID, gsea_c1_vs_c2_hm_gse33382$ID), ]

## C3 VS C2 (C2 is negative NES)

gsea_c3_vs_c2_GO_gse21257[gsea_c3_vs_c2_GO_gse21257$ID %in% intersect(gsea_c3_vs_c2_GO_gse21257$ID, gsea_c3_vs_c2_GO_gse33382$ID), ]
gsea_c3_vs_c2_hm_gse21257[gsea_c3_vs_c2_hm_gse21257$ID %in% intersect(gsea_c3_vs_c2_hm_gse21257$ID, gsea_c3_vs_c2_hm_gse33382$ID), ]

## C3 VS 1C (C1 is negative NES)

gsea_c3_vs_c1_GO_gse21257[gsea_c3_vs_c1_GO_gse21257$ID %in% intersect(gsea_c3_vs_c1_GO_gse21257$ID, gsea_c3_vs_c1_GO_gse33382$ID), ]
gsea_c3_vs_c1_hm_gse21257[gsea_c3_vs_c1_hm_gse21257$ID %in% intersect(gsea_c3_vs_c1_hm_gse21257$ID, gsea_c3_vs_c1_hm_gse33382$ID), ]
