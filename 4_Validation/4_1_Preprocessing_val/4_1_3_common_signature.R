library(biomaRt)
library(lumi)
library(illuminaHumanv2.db)
library(lumiHumanIDMapping)

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

common <- intersect(pucky, myannot$hgnc_symbol) # Keep available genes in signature form ensembl


y <- intersect(pucky, annot$Symbol) # Original available genes

"TMEM49" %in% annot$Symbol # Hand searcehd gene equivalent to VMP1 

true_common <- unique(c(common, y, "TMEM49")) 

probes <- mapIds( # VMP1 and SLC45A4 are not mapped to probe ids, still SLC45A4 appears in the GLP object, and VMP1 appears as TNEN49
  illuminaHumanv2.db,
  keys = true_common,
  "PROBEID",
  keytype = "SYMBOL",
  multiVals = "first"
)

pucky_gse <- c(annot$Symbol[annot$ILMN_ID[,5] %in% probes], "TMEM49", "SLC45A4")











