
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(patchwork)

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

p1 <- ggplot(gene_counts_long, aes(x = Var1, y = Count, fill = Status)) +
  geom_col(position = "stack") +
  labs(
    x = "Gene (Var1)",
    y = "Total Appearances",
    fill = "Category",
    title = "Gene Counts: Confirmed vs. Tentative"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p1

gene_counts_long %>% 
  group_by(Var1) %>% 
  summarise(sum = sum(Count)) %>% 
  filter(sum > median(sum))

total_counts <- gene_counts_long %>% 
  group_by(Var1) %>% 
  summarise(sum = sum(Count))


boruta_list_pert <- readRDS("results/boruta/boruta_conf_pert.RDS")
boruta_tent_pert <- readRDS("results/boruta/boruta_tent_pert.RDS")

confirmed_pert <- as.data.frame(table(unlist(boruta_list_pert))) %>% 
  rename(Confirmed = "Freq")
tentative_confirmed_pert <- as.data.frame(table(unlist(boruta_tent_pert))) %>%  
  rename(Tentative = "Freq") 


gene_counts_long_pert <- confirmed_pert %>%
  full_join(tentative_confirmed_pert, by = "Var1") %>%
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

ggplot(gene_counts_long_pert, aes(x = Var1, y = Count, fill = Status)) +
  geom_col(position = "stack") +
  labs(
    x = "Gene (Var1)",
    y = "Total Appearances",
    fill = "Category",
    title = "Gene Counts: Confirmed vs. Tentative"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

gene_counts_long_pert %>% 
  group_by(Var1) %>% 
  summarise(sum = sum(Count)) %>% 
  filter(sum > 75)


total_counts_pert <- gene_counts_long_pert %>% 
  group_by(Var1) %>% 
  summarise(sum = sum(Count))




confirmed_pert_5 <- as.data.frame(table(unlist(boruta_list_pert[1:24]))) %>% 
  rename(Confirmed = "Freq")

confirmed_pert_10 <- as.data.frame(table(unlist(boruta_list_pert[25:49]))) %>% 
  rename(Confirmed = "Freq")

confirmed_pert_15 <- as.data.frame(table(unlist(boruta_list_pert[50:74]))) %>% 
  rename(Confirmed = "Freq")

confirmed_pert_20 <- as.data.frame(table(unlist(boruta_list_pert[75:100]))) %>% 
  rename(Confirmed = "Freq")


total_counts_appear_all <- total_counts[total_counts$Var1 %in% confirmed_pert_5$Var1 & total_counts$Var1 %in% confirmed_pert_10$Var1 & total_counts$Var1 %in% confirmed_pert_15$Var1 & total_counts$Var1 %in% confirmed_pert_20$Var1, ]


p2 <- ggplot(gene_counts_long[gene_counts_long$Var1 %in% total_counts_appear_all$Var1, ], aes(x = Var1, y = Count, fill = Status)) +
  geom_col(position = "stack") +
  labs(
    x = "Gene (Var1)",
    y = "Total Appearances",
    fill = "Category",
    title = "Gene Counts: Confirmed vs. Tentative"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


p3 <- ggplot(gene_counts_long, aes(x = Var1, y = Count, fill = Var1 %in% total_counts_appear_all$Var1 )) +
  geom_col(position = "stack") +
  labs(
    x = "Gene (Var1)",
    y = "Total Appearances",
    fill = "Category",
    title = "Gene Counts: Confirmed vs. Tentative"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p1 / p3

total_counts_pert[total_counts_pert$Var1 %in% total_counts_appear_all$Var1, ]


total_counts_appear_all$rank_all <- seq_along(total_counts_appear_all$Var1)

total_counts_pert$rank_pert <- seq_along(total_counts_pert$Var1)


total_counts_appear_all %>% 
  left_join(total_counts_pert, by = "Var1") %>% 
  pivot_longer(cols = c(rank_all, rank_pert),
               names_to = "rank",
               values_to = "value") %>% 
  ggplot(aes(x = Var1, y = value)) +
  geom_col(position = "stack") +
  facet_grid(~ rank)


total_counts_appear_all %>% 
  left_join(total_counts_pert, by = "Var1") %>% 
  ggplot(aes(y = reorder(Var1, -rank_all))) +
  geom_segment(aes(x = rank_all, xend = rank_pert, yend = Var1), color = "grey70") +
  geom_point(aes(x = rank_all, color = "All"), size = 3) +
  geom_point(aes(x = rank_pert, color = "Perturbed"), size = 3) +
  scale_color_manual(values = c("All" = "#1F77B4", "Perturbed" = "#FF7F0E")) +
  labs(x = "Rank", y = "Gene (Var1)", color = "Condition") +
  theme_minimal()


