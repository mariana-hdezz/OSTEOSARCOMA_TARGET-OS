library(biomaRt)
library(lumi)
library(illuminaHumanv2.db)
library(lumiHumanIDMapping)

# 1. Convert nuID to ILMN Probe ID

annot$ILMN_ID <- nuID2IlluminaID(as.character(annot$ID), chipVersion = "Human-6 v2")

# The following lines correspond to mapping to ENSEMBL, ENTREZ and SYMBOL from the previous step


ensembl_annot <- data.frame(Gene = unlist(mget(
  x = as.character(annot$ILMN_ID[, 5]) , envir = illuminaHumanv2ENSEMBL
)))


mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl") 

#We create myannot, with GC content, biotype, info for length & names per transcript

myannot <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filters = "ensembl_gene_id", 
                 values =  ensembl_annot$Gene,  #annotate the genes in the count matrix 
                 mart = mart)

common <- intersect(pucky, myannot$hgnc_symbol)


y <- intersect(pucky, annot$Symbol)

"TMEM49" %in% annot$Symbol

true_common <- unique(c(common, y, "TMEM49"))

probes <- mapIds(
  illuminaHumanv2.db,
  keys = true_common,
  "PROBEID",
  keytype = "SYMBOL",
  multiVals = "first"
)

pucky_gse <- c(annot$Symbol[annot$ILMN_ID[,5] %in% probes], "TMEM49", "SLC45A4")











