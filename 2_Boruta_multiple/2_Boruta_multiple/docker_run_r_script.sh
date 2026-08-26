docker compose up --build

docker compose exec r-analysis Rscript 2_Boruta_multiple/2_1_1_Boruta_surv_bin.R

docker compose exec r-analysis Rscript 2_Boruta_multiple/2_1_2_Boruta_surv_bin_pert.R
