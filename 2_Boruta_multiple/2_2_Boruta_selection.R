
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(patchwork)
library(tidymodels)
library(survminer)
library(broom)
library(survival)
library(Boruta)




boruta_list <- readRDS("results/boruta/boruta_conf.RDS")
boruta_tent <- readRDS("results/boruta/boruta_tent.RDS")
boruta_sign <- readRDS("results/boruta/boruta_signature.RDS")

# This part is for demosntrating that the mean, median and mode belong to 37 as gene set size
# it is commented because it is heavy

# set_size <- list()
# 
# for (i in 1:100) {
#   
#   x <- boruta_sign[[i]]
#   
#   x <- TentativeRoughFix(x)
#   
#   set_size[[i]] <- length(x$finalDecision[x$finalDecision == "Confirmed"])
#   
# }
# summary(unlist(set_size))

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


total_counts <- gene_counts_long %>% 
  group_by(Var1) %>% 
  summarise(Count = sum(Count))

gene_signature <- as.character(total_counts$Var1[1:37])


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


# Perturbation analysis ---------------------------------------------------


boruta_list_pert <- readRDS("results/boruta/boruta_conf_pert.RDS")
boruta_tent_pert <- readRDS("results/boruta/boruta_tent_pert.RDS")

total_counts_pert <-as.data.frame(table(unlist(boruta_list_pert))) %>% 
    rename(Confirmed = "Freq") %>% 
  arrange(desc(Confirmed)) 

confirmed_pert_5 <- as.data.frame(table(unlist(boruta_list_pert[1:25]))) %>% 
    rename(Confirmed = "Freq") %>% 
  mutate(pert = "5%")

confirmed_pert_10 <- as.data.frame(table(unlist(boruta_list_pert[25:49]))) %>% 
    rename(Confirmed = "Freq") %>% 
  mutate(pert = "10%")

confirmed_pert_15 <- as.data.frame(table(unlist(boruta_list_pert[50:74]))) %>% 
    rename(Confirmed = "Freq") %>% 
  mutate(pert = "15%")
 
confirmed_pert_20 <- as.data.frame(table(unlist(boruta_list_pert[75:100]))) %>% 
    rename(Confirmed = "Freq") %>% 
  mutate(pert = "20%")


total_counts_appear_all <- total_counts[total_counts$Var1 %in% confirmed_pert_5$Var1 & total_counts$Var1 %in% confirmed_pert_10$Var1 & total_counts$Var1 %in% confirmed_pert_15$Var1 & total_counts$Var1 %in% confirmed_pert_20$Var1, ]



total_counts_appear_all$rank_all <- seq_along(total_counts_appear_all$Var1)

total_counts$rank_all <- seq_along(total_counts$Var1)

total_counts_pert$rank_pert <- seq_along(total_counts_pert$Var1)



# Visualization and gene signature creation ----------------------------------

total_counts_appear_all %>% 
  left_join(total_counts_pert, by = "Var1") %>% 
  filter(Var1 %in% gene_signature) %>% 
  ggplot(aes(y = reorder(Var1, -rank_all))) +
  geom_segment(aes(x = rank_all, xend = rank_pert, yend = Var1), color = "grey70") +
  geom_point(aes(x = rank_all, color = "All"), size = 3) +
  geom_point(aes(x = rank_pert, color = "Perturbed"), size = 3) +
  scale_color_manual(values = c("All" = "#1F77B4", "Perturbed" = "#FF7F0E")) +
  labs(x = "Rank", y = "Gene (Var1)", color = "Condition") +
  theme_minimal()



gene_signature <- as.character(total_counts$Var1[1:37])


p2 <- ggplot(gene_counts_long, aes(x = Var1, y = Count, fill = Var1 %in% gene_signature )) +
  geom_col(position = "stack") +
  labs(
    x = "Gene (Var1)",
    y = "Total Appearances",
    fill = "Category",
    title = "Gene Counts: Confirmed vs. Tentative"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  geom_hline(yintercept = 35)

p1 / p2


write.table(matrix(gene_signature, nrow = 1), file = "output_data/gene_signature.csv", sep = ",", row.names = FALSE, col.names = FALSE)

rm(list = ls())
gc()