
library(paletteer)

col_scheme <- c("#088BBEFF", "#C75DABFF", "#F5D000FF", "black", "steelblue")

title_size <- 
  theme(plot.title = element_text(size = 5))
  

line_width <- 8

base_size <- 4

text_size <- 3

p_size <- 2

p_coord <- c(130, 0.5)

ph_coord <- function(i){annotate(
  geom = "text", 
  x = 150,          
  y = 0.20,         
  label = paste0("PH assumption ", round(i$table[1,3], 2)), 
  color = "black", 
  fontface = "bold"
)
}

fill_col <- c(
  "low" = "#43d6ff",
  "mid" = "#eddd8a",
  "high" = "#e953c3"
)

legend <- theme(
  legend.title = element_text(size = 3),     
  legend.text = element_text(size = 3),      
  legend.key.size = unit(1, "mm"),            
  legend.key.height = unit(3, "mm"),        
  legend.key.width = unit(3, "mm")
) 
