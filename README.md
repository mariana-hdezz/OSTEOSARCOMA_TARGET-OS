# OSTEOSARCOMA_TARGET-OS

**Titulo del articulo**

Computational pipeline to identify gene seignature and cluster TARGET-OS, GSE21257 and GSE33382 patients as well as to characterize clinical and moleculat characteristics

---

## Overview

1. Download and preprocess TARGET-OS RNA-Seq data; GSE212257 and GSE33382 microarray
2. Obtain gene signature for TARGET OS and adapt comp<atibility with both validation datasets
3. Hierarchcical clustering on TARGE-TOS patients / nearest centroids to classify the validation sets patients
4. Clinical analysis of the clusters on all dataset
5. Molecular analysis of the clusters on all dataset
6. Observe the results and figures

---

## Directory and Scripts

| **Directory** | *Subdirectory* | `Script` | Description |
|-----------|--------------|--------|-------------|
| **1_preprocessing_data** | *1_1_training_preprocess* | `1_1_1_training_preprocess.R` | Download TARGET-OS RNA-seq (TCGAbiolinks); pre-process metadata, filter protein-coding and duplicated genes, obtain both raw and RPKM data, normalize raw data with VST |
| **1_preprocessing_data** | *1_2_validation_preprocess* | `1_2_1_gse33382_gse21257_preprocess.R` | Download GSE21257 and GSE33382 (GEOquery); pre-process metadata, filter duplicated genes, batch correction (Limma) |
| **1_preprocessing_data** | *1_2_validation_preprocess* | `1_2_1_gse33382_gse21257_preprocess.R` | Download GSE21257 and GSE33382 (GEOquery); download GPL for both; pre-process metadata, filter duplicated genes, batch correction (Limma) |
| **1_preprocessing_data** | *1_2_validation_preprocess* | `1_2_2_gse39055_preprocess.R` | Download GSE39055(GEOquery); pre-process metadata, filter duplicated genes (This script has less impact on the paper and subsequent analysis |
| **2_gene_signature** | *none* | `2_0_1_docker_composer_script.sh`  `2_0_2_docker_run_boruta.sh` `2_0_3_docker_boruta_pert.sh`  | Shell script to run the docker for the Boruta analysis ordered by the compose, the 100 iteration Boruta and the 100 iteration perturbated Boruta |
| **2_gene_signature** | *none* | `2_1_1_Boruta_surv_bin.R` `2_1_2_Boruta_surv_bin_pert.R`| Perform the Boruta algorithm with (`2_1_2_Boruta_surv_bin_pert.R`) and without  (`2_1_1_Boruta_surv_bin.R`) perturbation. Preferentially run inside the docker with the previous scripts and not directly |
| **2_gene_signature** | *none* | `2_2_1_Boruta_selection.R` | After having obtained the results from boruta, this script creates the gene signature from the results. It automatically obtaines the top 37 rnked genes but also performs the Lasso (via sourcing of `2_2_2_lasso_elasticNet.R`. Plots supplementary Figure 1 |
| **2_gene_signature** | *none* | `2_3_common_signature.R` | Adapt the signature created on `2_2_1_Boruta_selection.R` utilizing the GPL obtained on `1_2_1_gse33382_gse21257_preprocess.R` to adapt the signature to the validation sets. Notice that at one point TMEM49 is added manually and that is because it is equivalent to VMP1 but the pipeline did not map VMP1 to TMEM49 |
| **3_clustering_assignment** | *3_1_cluster_train* | `3_1_hierarchical_clustering.R` | Hierarchcial clustering utilizing the gene signature to classify TARGET-OS patients |
| **3_clustering_assignment** | *3_1_cluster_train* | `3_2_cluster_validation` | Obtain nearest centroid weights on the clusters of TARGET-OS and then classify the patients in the validation sets |
| **4_clinical_analysis** | *4_1_target_clin_anal* | `4_1_1_survival_analys_target.R` | Log rank test for both survival and recurrence between clusters, creates kaplan meier curve, performs chi squared test between clusters based on Huvos grade on TARGET-OS |
| **4_clinical_analysis** | *4_2_validation_clin_anal* | `4_2_validation_clin_anal` | Log rank test for both survival and recurrence between clusters, creates kaplan meier curve, performs chi squared test between clusters based on Huvos grade and histologic subtype on validation sets |
| **5_molecular_analysis** | *5_1_target_mol* | `5_1_1_diffEx_gsea_target.R` | Differential expression and log fold ranked GSEA between clusters in TARGET-OS |
| **5_molecular_analysis** | *5_1_target_mol* | `5_1_1_diffEx_gsea_target.R` | Mean ranked GSEA for each cluster in TARGET-OS |
| **5_molecular_analysis** | *5_2_val_mol* | `5_2_1_diffEx_val.R` | Differential expression for both GSE21257 and GSE33382 for each clsuter compariuson, it also creates the results for log fold ranekd GSEA by sourcing `5_2_2_gsea_val.R` |
| **5_molecular_analysis** | *5_2_val_mol* | `5_2_3_mean_ranked_gsea_val.R` | Mean ranked GSEA for each cluster in validation sets |
| **5_molecular_analysis** | *5_2_val_mol* | `5_2_4_hist_strat_diffex_gsea.R` | Differential expression for the validation sets merged counts for each clsuter comparison stratefied by histologic subtype, it also creates the results for log fold ranekd GSEA by sourcing `5_2_2_gsea_val.R` |
| **6_result_files** | *6_1_clinical_res* | `6_1_1_surv_anal_res.R.R` | Creates Table 1 comparing the log rank tests between clusters on both TARGET-OS and GSE21257. Showcases the results for the chi squared test between clusters for both Huvos grade (for all datasets) and histologic subtype (on validation datasets). Plots Kaplan meier curves |
| **6_result_files** | *6_2_mol_results* | `6_2_1_target_diffex_gsea_res.R` | Simplified results from differential expression and GSEA on TARGET-OS. For complete tables in **results** directory. |
| **6_result_files** | *6_2_mol_results* | `6_2_2_val_diffex_gsea_res.R` | Simplified results from differential expression and GSEA on Validation sets. For complete tables in **results** directory. |
| **6_result_files** | *6_2_mol_results* | `6_2_3_heatmaps_target.R` `6_2_4_heatmaps_gse.R` `6_2_5_heatmaps_hist_strat.R` `6_2_6_concatenated_heatmaps.R` | Plotting all the heatmaps from the paper. The firt 3 scripts create the heatmaps and `6_2_6_concatenated_heatmaps.R` plots the actuall figure used in the paper |  
| **6_result_files** | *6_2_mol_results* | `6_2_7_pca.R` | Principal component analysis based on the signature for both TARGET-OS and the validation sets |
| **7_isolated_functions** | *none* | `2_2_2_lasso_elasticNet.R` | Script sourced for `2_2_1_Boruta_selection.R` to apply logistic regression penalized by Elastic net/lasso on the genes seelected by `2_1_1_Boruta_surv_bin.R` with survival as outcome |
| **7_isolated_functions** | *none* | `5_2_2_gsea_val.R` | Script sourced for `5_2_1_diffEx_val.R` and `5_2_4_hist_strat_diffex_gsea.R` to apply GSEA and save the results with the names corresponding to the comparison |

---

## Dependencies

R version::

renv available

---

## Manual interventions

Only on `1_1_1_training_preprocess.R` is a manuall input asaked for and it is to establish the directory of download for those data 

---

## Execution order

The whole pipeline is numbered

Numbered execution suggested since that assures that all necessary inputs on a certain script have already been obtained

The only exceptions lie in directoy **7_isolated_functions** which are not directly run, they are sourced to other scripts.

---
No manuall directory creation is necessary except for the desired directory to download the TARGET-OS ddata
Script 1 creates the output_data and resuts directory which will have diufferent objects. Other subdirectories for results are created on the script that needs it. 
After running the full script the OSTEOSARCOMA_TARGET-OS direcetory will contain about 308 mb without taking into account the TARGET-OS data. It does not create temporary files for most data
---
