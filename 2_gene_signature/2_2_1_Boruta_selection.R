
#############################################################################
#> Gene signature selection based on the 100 Boruta iterations, 100
#> perturbed Boruta iterations and Lasso penalized logistic regression
#> 
#> 
#> Inputs: The various objects obtained at each Boruta script
#> boruta_conf = The genes htat where confirmed even before tentative rough fix
#> boruta_tent = Only includes the genes confirmed after tentative decision force
#> boruta_signature = Full object of each iteration for more complete analysis
#> 
#> Outputs: gene_signature_gse
#> 
#############################################################################

# Libraries

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

# Load data

# 100 Iteration Boruta


boruta_list <- readRDS("results/boruta/boruta_conf.RDS")
boruta_tent <- readRDS("results/boruta/boruta_tent.RDS")
boruta_sign <- readRDS("results/boruta/boruta_signature.RDS")

# Perturbation boruta

boruta_list_pert <- readRDS("results/boruta/boruta_conf_pert.RDS")
boruta_tent_pert <- readRDS("results/boruta/boruta_tent_pert.RDS")

# This part is for demosntrating that the mean, median and mode belong to 37 as gene set size

set_size <- list()

for (i in 1:100) { # For loop that obtains each signatue full object, calculates the length and then assigns it to list named by number of iteration

  x <- boruta_sign[[i]]

  set_size[[i]] <- length(x$finalDecision[x$finalDecision == "Confirmed"])

}

# Summary

summary(unlist(set_size))

# Mode

table(unlist(set_size))[table(unlist(set_size)) == max(table(unlist(set_size)))]

# Table with name and freuency of apperance for pure confirmed

confirmed <- as.data.frame(table(unlist(boruta_list))) %>% 
  dplyr::rename(Confirmed = "Freq")

# Table with name and freuency of apperance for genes confirmed after forcing tenttaive decision

tentative_confirmed <- as.data.frame(table(unlist(boruta_tent))) %>%  
  dplyr::rename(Tentative = "Freq") 

# Join previous objects with values for each gene of mode of being selected by each iteration

pre_gene_counts_long <- confirmed %>%
  full_join(tentative_confirmed, by = "Var1") %>%
  mutate(
    Confirmed = as.numeric(Confirmed),
    Tentative = as.numeric(Tentative)
  ) 

# Pivot long for plot

gene_counts_long <- pre_gene_counts_long %>%
  pivot_longer(
    cols = c(Confirmed, Tentative), 
    names_to = "Status", 
    values_to = "Count"
  ) %>%
  filter(!is.na(Count)) %>% 
  mutate(Var1 = fct_reorder(Var1, Count, .fun = sum, .desc = TRUE))

# Sum the counts independentof if initially tentative or initially confirmed

total_counts <- gene_counts_long %>% 
  group_by(Var1) %>% 
  summarise(Count = sum(Count))

# Keep top 37 like in the paper

gene_signature <- as.character(total_counts$Var1[1:37])


# Perturbation analysis ---------------------------------------------------

# Here we only focuse on those confirmed directly not on those confirmed after tentative decision

# We create a generallobject as well as individual objects for each pert analysis
# (5, 1, 15 and 25%). Since each analysis ran for 25

total_counts_pert <-as.data.frame(table(unlist(boruta_list_pert))) %>% 
    dplyr::rename(Confirmed = "Freq") %>% 
  arrange(desc(Confirmed)) 

confirmed_pert_5 <- as.data.frame(table(unlist(boruta_list_pert[1:25]))) %>% 
    dplyr::rename(Confirmed = "Freq") %>% 
  mutate(pert = "5%")

confirmed_pert_10 <- as.data.frame(table(unlist(boruta_list_pert[26:50]))) %>% 
    dplyr::rename(Confirmed = "Freq") %>% 
  mutate(pert = "10%")

confirmed_pert_15 <- as.data.frame(table(unlist(boruta_list_pert[51:75]))) %>% 
    dplyr::rename(Confirmed = "Freq") %>% 
  mutate(pert = "15%")
 
confirmed_pert_20 <- as.data.frame(table(unlist(boruta_list_pert[76:100]))) %>% 
    dplyr::rename(Confirmed = "Freq") %>% 
  mutate(pert = "20%")

# Keep genes from original 100 iterations that appeared at least once in each perturbation analysis

total_counts_appear_all <- total_counts[total_counts$Var1 %in% confirmed_pert_5$Var1 & total_counts$Var1 %in% confirmed_pert_10$Var1 & total_counts$Var1 %in% confirmed_pert_15$Var1 & total_counts$Var1 %in% confirmed_pert_20$Var1, ]

# Rank by counts

total_counts_appear_all$rank_all <- seq_along(total_counts_appear_all$Var1)

total_counts$rank_all <- seq_along(total_counts$Var1)

# Assign values of rank oin pert

total_counts_pert$rank_pert <- seq_along(total_counts_pert$Var1)


# Visualization 

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


# Lasso

source("./7_isolated_functions/2_2_2_lasso_elasticNet.R")


# Characteristics of gene signature oif chosen by Lasso


total_counts %>% 
  filter(Var1 %in% lasso_sign) %>% 
  summarise(
    max_c = max(Count),
    min_r = min(rank_all),
    mean_c = mean(Count),
    mean_r = mean(rank_all),
    medi_c = median(Count),
    medi_r = median(rank_all),
    min_c = min(Count),
    max_r = max(rank_all),
    sd_c = sd(Count)
            
            )

# Characteristics of gene signature oif chosen by top 37 genes


total_counts %>% 
  filter(Var1 %in% gene_signature) %>% 
  summarise(
    max_c = max(Count),
    min_r = min(rank_all),
    mean_c = mean(Count),
    mean_r = mean(rank_all),
    medi_c = median(Count),
    medi_r = median(rank_all),
    min_c = min(Count),
    max_r = max(rank_all),
    sd_c = sd(Count)
    
  )


write.table(matrix(gene_signature, nrow = 1), file = "output_data/gene_signature.csv", sep = ",", row.names = FALSE, col.names = FALSE)


# Visualizte the genes selected by distinct methods and combinations

pre_gene_counts_long %>%
  mutate(appear = case_when(
    Var1 %in% lasso_sign & 
      Var1 %in% gene_signature & 
       Var1 %in% total_counts_appear_all$Var1 ~ "Appears in all",
    !(Var1 %in% lasso_sign | 
      Var1 %in% gene_signature | 
      Var1 %in% total_counts_appear_all$Var1) ~ "Appears in none",
    Var1 %in% lasso_sign & 
      !(Var1 %in% gene_signature | 
      Var1 %in% total_counts_appear_all$Var1) ~ "Appears only in Lasso",
    Var1 %in% gene_signature & 
      !(Var1 %in% lasso_sign | 
          Var1 %in% total_counts_appear_all$Var1) ~ "Appears only in Top 37 genes",
    Var1 %in% total_counts_appear_all$Var1 & 
      !(Var1 %in% lasso_sign | 
          Var1 %in% gene_signature) ~ "Appears only in perturbation",
    !Var1 %in% total_counts_appear_all$Var1 & 
      (Var1 %in% lasso_sign | 
          Var1 %in% gene_signature) ~ "Appears in Lasso and Top 37 but not pert",
    !Var1 %in% lasso_sign & 
      (Var1 %in% gene_signature | 
          Var1 %in% total_counts_appear_all$Var1) ~ "Appears in Top 37 and perturbations but not Lasoo",
    !Var1 %in% gene_signature & 
      (Var1 %in% lasso_sign | 
          Var1 %in% total_counts_appear_all$Var1) ~ "Appears in perturbations and Lasso but not Top 37",

    )
         ) %>%
  pivot_longer(
    cols = c(Confirmed, Tentative), 
    names_to = "Status", 
    values_to = "Count"
  ) %>%
  filter(!is.na(Count)) %>% 
  mutate(Var1 = fct_reorder(Var1, Count, .fun = sum, .desc = TRUE)) %>%
ggplot(aes(x = Var1, y = Count, fill = appear)) +
  geom_col(position = "stack") +
  labs(
    x = "Gene (Var1)",
    y = "Total Appearances",
    fill = "Category",
    title = "Genes that appear in the different conditions"
  ) +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_discrete(palette = "Accent")



rm(list = ls())
gc()
