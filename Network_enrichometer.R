# # Enrichmentometer (Fixed & Robust)

library(clusterProfiler)
library(org.Hs.eg.db)
library(igraph)

# 1. Map unique network symbols to ENSEMBL IDs ONCE
gene_conv <- bitr(
  V(g)$name, 
  fromType = "SYMBOL", 
  toType   = "ENSEMBL", 
  OrgDb    = org.Hs.eg.db
)

# Clean/format mapping table

gene_conv$SYMBOL <- toupper(trimws(gene_conv$SYMBOL))

# 2. Extract module ENSEMBL IDs directly from the mapping table
# Note: Using V(g)$name and lc$membership (or V(g)$member) from the same graph 'g'

module_genes <- lapply(seq_along(lc), function(i) {

  mod_symbols <- toupper(trimws(V(g)$name[V(g)$member == i]))
  
  mapped_ensembl <- gene_conv$ENSEMBL[gene_conv$SYMBOL %in% mod_symbols]
  
  return(unique(mapped_ensembl))
})

names(module_genes) <- paste0("Module_", seq_along(lc))

# 3. Run vectorized enrichment across all modules
comp_go <- compareCluster(
  geneClusters  = module_genes,
  fun           = "enrichGO",
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",  
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE,
  minGSSize = 30
)


# 4. View and visualize results
enrich_mod <- as.data.frame(comp_go)


enrich_mod %>% 
  dplyr::count(Description) %>% 
  filter(n > (0.05 * length(seq_along(lc))))


comp_go_simplified <- clusterProfiler::simplify(
  comp_go, 
  cutoff      = 0.7,      # Similarity threshold (0.7 is standard; lower = more strict)
  by          = "p.adjust", # Keep the term with the best p-value in each cluster
  select_fun  = min,
  measure     = "Wang"
)

# Check how many terms were reduced
dim(as.data.frame(comp_go))            # Original row count
dim(as.data.frame(comp_go_simplified)) # Reduced row count


simplif_df <- as.data.frame(comp_go_simplified)

simplif_df %>% 
  dplyr::count(Description) %>% 
  filter(n > 7)

