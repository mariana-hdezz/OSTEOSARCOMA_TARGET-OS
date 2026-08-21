docker compose up --build

docker compose exec r-boruta-analysis Rscript 2_Boruta_multiple/Boruta_surv_bin.R

docker compose exec r-boruta-analysis Rscript 2_Boruta_multiple/Boruta_surv_bin_perm.R
