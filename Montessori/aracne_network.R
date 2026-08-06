library(igraph)
library(dplyr)
library(ggraph)
library(RColorBrewer)
library(dplyr)
library(igraph)


aracne <- read.delim(
  "results/network.txt",
  header = TRUE,
  sep = "\t"
)


g <- graph_from_data_frame(aracne)


# 1. Initial attributes

# Num edges

gsize(g) 


# Num nodes

gorder(g)


# 2. Parameter

# 2.1 Degree

V(g)$degree <- degree(g, mode = c("All")) 

V(g)$degree_in <- degree(g, mode = c("in")) 

V(g)$degree_out <- degree(g, mode = c("out")) 

# Max degree

max(V(g)$degree)
V(g)$name[which.max(V(g)$degree)]

# Mean

(2 * gsize(g)) / gorder(g)

# Eigenvector

V(g)$eigen <- eigen_centrality(g)$vector

# Max eigenvector

max(V(g)$eigen)
V(g)$name[which.max(V(g)$eigen)]

# Betweenes centrality

V(g)$bc <- betweenness(g)

# Max betweennes centrality 

max(V(g)$bc)
V(g)$name[which.max(V(g)$bc)]

# Edge betweenness

E(g)$eb <- edge_betweenness(g)

# Max edge betweenness

max(E(g)$eb)
E(g)[which.max(E(g)$eb)]

# Network density

edge_density(g) 

transitivity(g, type = "local")

# Degree dist

plot(degree_distribution(g, cumulative = FALSE))

# Dist distribution

dist_data <- distance_table(g, directed = TRUE)

df <- data.frame(
  Path_Length = 1:length(dist_data$res),
  Count = dist_data$res
)

ggplot(df, aes(x = factor(Path_Length), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "white") 


# Community
set.seed(124)
und_g <- as_undirected(g)

set.seed(123)

lc <- cluster_louvain(und_g)

communities(lc)

V(g)$member <- membership(lc)

n_of_members <- data.frame(names = names(table(V(g)$member)),
                           values = as.numeric(as.character(table(V(g)$member))))

ggplot(n_of_members, aes(x = names, y = values)) +
  geom_col()

module_list <- list()

for (i in 1:length(communities(lc))) {
  if(length(V(g)$name[V(g)$member == i]) > 50){
    und_g <- subgraph(g, V(g)$name[V(g)$member == i])
    
    if(median(E(und_g)$eb) > 2 & median(degree(und_g, V(und_g), "all")) >= 1){
      
      cat("#############\n")
      print(paste("Module:", i))
      cat("-------------")
      print(summary(degree(und_g, V(und_g), "all"))) 
      print(summary(V(und_g)$bc))
      print(summary(V(und_g)$eigen))
      print(summary(E(und_g)$eb))
      module_list[[i]] <- i
      
      cat("//////////////\n")
    }else{
      
    }
  }else{
    
  }
}

library(clusterProfiler)
library(org.Hs.eg.db)

module_list <- as.character(Filter(Negate(is.null), module_list))

module_list


for (i in module_list) {
  cat("\n##############\n", i, "\n")
  genes_to_test <- V(g)$name[V(g)$member == i]
  
  genes_to_test <- toupper(trimws(genes_to_test))
  
  # Convert Symbols to Entrez IDs
  gene_conv <- bitr(genes_to_test, 
                    fromType = "SYMBOL", 
                    toType   = "ENSEMBL", 
                    OrgDb    = org.Hs.eg.db)
  
  
  # Run enrichment using the Entrez IDs
  go_results <- enrichGO(gene = gene_conv$ENSEMBL,
                         OrgDb = org.Hs.eg.db,
                         keyType = 'ENSEMBL', 
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         readable = TRUE) 
  # 3. View the results
  print(head(go_results, n = 10))
}

sig_genes <- list()

for (i in c(3)) {
  
  sig_genes[[i]] <- V(g)$name[V(g)$member == i & V(g)$degree > 5]
  
}


sig_g <- as.character(unlist(sig_genes))


module_list

# Find all target nodes and their 2-step neighbors
target_nodes <- V(g)[member %in% c(12, 17)]
all_ego_nodes <- ego(g, order = 2, nodes = target_nodes, mode = "all")

# Flatten the list of vertex sequences into a single unique set of nodes
combined_nodes <- unique(do.call(c, all_ego_nodes))

# Build one combined subgraph
combined_subgraph <- induced_subgraph(g, combined_nodes)


x <- as.undirected(combined_subgraph)


# --- Step 1: Filter Labels ---
d <- degree(x, mode = "all")
threshold <- 5  # Adjust this value based on your network's degree distribution

v_labels <- V(x)$name  # Get node names
v_labels[d <= threshold] <- NA  # Hide labels for low-degree nodes

# --- Step 2: Expand Layout ---
l <- layout_with_kk(x)
l <- l * 2.5  # Increase multiplier to disperse nodes more

# --- Step 3: Plot ---
plot(
  x,
  layout = l,
  vertex.label = v_labels,
  vertex.label.dist = 0.5,
  vertex.label.cex = 0.6,
  vertex.color = V(x)$member,
  vertex.size = sqrt(V(x)$degree) * 10,
  edge.color = adjustcolor("black", alpha.f = 0.1),
  rescale = FALSE,
  xlim = range(l[,1]), 
  ylim = range(l[,2])
)

E(g)[E(g)$eb == max(E(g)$eb[V(g)$member == 1])]

module_list
list_1 <- list()
for (i in gene_list_vst_sur){
  print(i)
  list_1[[i]] <- (V(g)$member[V(g)$name == i])
  cat("-------\n")
}
table(unlist(list_1))
