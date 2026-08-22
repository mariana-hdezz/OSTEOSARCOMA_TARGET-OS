
#############################################################################
#> Since not all genes in GSE21257 and GSE33382 are present in owr signature as 
#> symbol (which could be outdated). This script searches for the true genes 
#> in common (even if with another symbol version) and create the gene 
#> signature that is compatible with those GSE. 
#> 
#> Note that 2 genes are added manually since there was no compatibility in probes 
#> But where preesen tin the original symbol comparison. Or the previous symbol name
#> was manually found (In supplementary we have the evidence that it is the same gene)
#> 
#> Inputs: annot, gene_signature
#> 
#> Outputs: gene_signature_gse
#> 
#############################################################################


library(biomaRt)
library(lumi)
library(illuminaHumanv2.db)
library(lumiHumanIDMapping)

gene_signature <- scan("output_data/gene_signature.csv", sep = ",", what = character())
annot <- readRDS("output_data/annot.RDS")

# 1. Convert nuID to ILMN Probe ID

annot$ILMN_ID <- nuID2IlluminaID(as.character(annot$ID), chipVersion = "Human-6 v2")

# The following lines correspond to mapping to ENSEMBL

ensembl_annot <- data.frame(Gene = unlist(mget(
  x = as.character(annot$ILMN_ID[, 5]) , envir = illuminaHumanv2ENSEMBL
)))


mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl") 


myannot <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filters = "ensembl_gene_id", 
                 values =  ensembl_annot$Gene,   
                 mart = mart)

common <- intersect(gene_signature, myannot$hgnc_symbol) # Keep available genes in signature form ensembl


y <- intersect(gene_signature, annot$Symbol) # Original available genes

"TMEM49" %in% annot$Symbol # Hand searcehd gene equivalent to VMP1 

true_common <- unique(c(common, y, "TMEM49")) 

probes <- mapIds( # VMP1 and SLC45A4 are not mapped to probe ids, still SLC45A4 appears in the GLP object, and VMP1 appears as TNEN49
  illuminaHumanv2.db,
  keys = true_common,
  "PROBEID",
  keytype = "SYMBOL",
  multiVals = "first"
)

gene_signature_gse <- c(annot$Symbol[annot$ILMN_ID[,5] %in% probes], "TMEM49")

write.table(matrix(gene_signature_gse, nrow = 1), file = "output_data/gene_signature_gse.csv", sep = ",", row.names = FALSE, col.names = FALSE)

rm(list = ls())
