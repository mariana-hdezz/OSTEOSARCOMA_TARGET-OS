docker compose up
docker compose exec r-analysis Rscript 5_molecular_analysis/5_1_target_mol/5_1_1_diffEx_gsea_target.R 
docker compose exec r-analysis Rscript 5_molecular_analysis/5_1_target_mol/5_1_2_mean_ranked_gsea.R 
docker compose exec r-analysis Rscript 5_molecular_analysis/5_2_val_mol/5_2_1_diffEx_val.R 
docker compose exec r-analysis Rscript 5_molecular_analysis/5_2_val_mol/5_2_3_mean_ranked_gsea_val.R
docker compose exec r-analysis Rscript 5_molecular_analysis/5_2_val_mol/5_2_4_hist_strat_diffex_gsea.R
docker compose down