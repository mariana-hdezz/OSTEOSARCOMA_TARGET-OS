
library(lumi)
library(illuminaHumanv2.db)

# 1. Convert nuID to ILMN Probe ID

annot$ILMN_ID <- nuID2IlluminaID(as.character(annot$ID), chipVersion = "Human-6 v2")

# The following lines correspond to mapping to ENSEMBL, ENTREZ and SYMBOL from the previous step


ensembl_annot <- data.frame(Gene = unlist(mget(
  x = as.character(annot$ILMN_ID[, 5]) , envir = illuminaHumanv2ENSEMBL
)))


symbol_anot <- data.frame(Gene = unlist(mget(
  x = as.character(annot$ILMN_ID[, 5]) , envir = illuminaHumanv2SYMBOL
)))

accnum_anot <- data.frame(Gene = unlist(mget(
  x = as.character(annot$ILMN_ID[, 5]) , envir = illuminaHumanv2ACCNUM
)))

# Convert signature to its ENSEMBL counterpart

pucky_ensembl <- gene_dist$ensembl[gene_dist$gene_name %in% pucky]

# Signature in symbol anottation from probe IDs

pucky[pucky %in% symbol_anot$Gene]

# Signature from GLP

pucky[pucky %in% annot$Symbol]

 # Obtaiun gene name in ensembl version of the signature that can be found in thge direct ensembl annot

pucky_ensembl[pucky_ensembl %in% ensembl_annot$Gene] # signature in ensembl annotation

sym_in_EnPucky <- gene_dist$gene_name[gene_dist$ensembl %in% pucky_ensembl[pucky_ensembl %in% ensembl_annot$Gene]] # Obtain symbol form available ensembl


common_pucky <- unique(c(sym_in_EnPucky, pucky[pucky %in% annot$Symbol], pucky[pucky %in% symbol_anot$Gene]))


pucky_probe <- mapIds(
  illuminaHumanv2.db,
  keys = common_pucky,
  column = "PROBEID",
  keytype = "SYMBOL",
  multiVals = "first"
)




pucky_probe %in% annot$ILMN_ID[,5]
probes <- as.character(annot$ILMN_ID[,1][annot$ILMN_ID[,5] %in% pucky_probe])

pucky_gse <- unique(c(annot$Symbol[annot$Illumina_Gene %in% probes], pucky[pucky %in% annot$Symbol]))

pucky[pucky%in%pucky_gse]


pucky[!(pucky%in%pucky_gse)]
pucky_gse[!(pucky_gse %in% pucky)]


rm(list = ls()[!(ls() %in% c("vst_counts", "pucky", "counts_data", "metadata_os", "pucky_gse", "annot"))])
