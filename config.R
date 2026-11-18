
library(paletteer)

col_scheme <- c("#088BBEFF", "#C75DABFF", "#F5D000FF", "black", "steelblue")

fill <- paletteer::scale_fill_paletteer_c("grDevices::PuOr", direction = -1)

line_width <- 8

base_size <- 6

text_size <- 3

fill <- c(
  "low" = "#088BBEFF",
  "mid" = "#C75DABFF",
  "high" = "#F5D000FF"
)
