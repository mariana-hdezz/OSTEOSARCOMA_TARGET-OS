library(tidymodels)
library(dplyr)
library(tidyr)
library(tibble)


vst_counts <- as.data.frame(readRDS("output_data/vst_counts.RDS"))
metadata_os <- readRDS("output_data/metadata_os.RDS")

metadata_surv <- 
  metadata_os %>% 
  dplyr::select(sample,
         survival_stat) %>% 
  drop_na(survival_stat)

vst_counts <- as.data.frame(t(vst_counts[as.character(total_counts$Var1), metadata_surv$sample]))

all(rownames(vst_counts) == metadata_surv$sample)

vst_counts$surv <- metadata_surv$survival_stat

# 2.- Preparing recipe and model ------------------------------------------

# 2.1 Recipe

lr_rec <- recipe(surv ~ ., data = vst_counts) %>% # Recurrence object created in preprocessing for linear regression
  step_dummy(all_nominal_predictors()) %>% 
  step_zv(all_predictors()) %>% # Eliminates variables with a single value
  step_nzv(all_predictors()) # Eliminates highly sparse variables


# 2.2 Model

lr_mod <- linear_reg(
  penalty = tune(),    # lambda establishes the severith of the penalty
  mixture = tune()     # alpha establishes the type, 1 being lasso, 0 being ridge, and 0.5 being elasticnet
) %>%
  set_engine("glmnet") # Engine that permits penalizing by elasticnet, ridge, and lasso

# 2.3 Workflow

lr_wf <- workflow() %>%
  add_model(lr_mod) %>%
  add_recipe(lr_rec)


# 3.- Selecting best penalizing parameters ------------------------------------

  set.seed(123)
  
  # 3.1 Folds for evaluating with resamples
  
  folds <- vfold_cv(
    vst_counts,
    v = 10,
    strata = surv
  )
  
  # 3.2 Grid for penalizing range
  
  grid <- grid_regular(
    penalty(range = c( - 4, 1)),   
    mixture(range = c(0.5, 1)),
    levels = 10
  )
  
  # 3.3 Running the different penalization methods
  
  res_ml <- tune_grid(
    lr_wf,
    resamples = folds,
    grid = grid,
    metrics = metric_set(rmse), 
    control = control_grid(save_pred = TRUE)
  )
  
  # 3.3.1 Observe metrics
  
  collect_metrics(res_ml)
  
  # 3.3.2 Object with best parameters for penalizing
  
  best_params <- select_best(res_ml, metric = "rmse")


# 4.- Actual training -----------------------------------------------------

# 4.1 Final workflow with the selecting the best parameter tested previously

final_wf <- finalize_workflow(lr_wf, best_params)

# 4.2 Final fit with training data

final_fit <- fit(final_wf, data = vst_counts)



# 4.2.2 Observing genes that are maintained after penalziation


coef_tbl <- tidy(final_fit) %>%
  filter(estimate != 0) %>%
  arrange(desc(abs(estimate)))

cat(coef_tbl$term, sep = ", ")

lasso_sign <- coef_tbl$term[-1]
