
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

boruta_list <- readRDS("results/boruta/boruta_conf.RDS")
boruta_tent <- readRDS("results/boruta/boruta_tent.RDS")

confirmed <- as.data.frame(table(unlist(boruta_list))) %>% 
  rename(Confirmed = "Freq")
tentative_confirmed <- as.data.frame(table(unlist(boruta_tent))) %>%  
  rename(Tentative = "Freq") 


gene_counts_long <- confirmed %>%
  full_join(tentative_confirmed, by = "Var1") %>%
  mutate(
    Confirmed = as.numeric(Confirmed),
    Tentative = as.numeric(Tentative)
  ) %>%
  pivot_longer(
    cols = c(Confirmed, Tentative), 
    names_to = "Status", 
    values_to = "Count"
  ) %>%
  filter(!is.na(Count)) %>% 
  mutate(Var1 = fct_reorder(Var1, Count, .fun = sum, .desc = TRUE))

ggplot(gene_counts_long, aes(x = Var1, y = Count, fill = Status)) +
  geom_col(position = "stack") +
  labs(
    x = "Gene (Var1)",
    y = "Total Appearances",
    fill = "Category",
    title = "Gene Counts: Confirmed vs. Tentative"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

gene_counts_long %>% 
  group_by(Var1) %>% 
  summarise(sum = sum(Count)) %>% 
  filter(sum > 75)
