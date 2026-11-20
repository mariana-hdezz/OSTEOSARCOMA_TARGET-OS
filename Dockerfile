FROM rocker/r-base:4.5.0

# 1. Install required system dependencies

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    curl \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libglpk-dev \
    libpng-dev \
    libssl-dev \
    libx11-dev \
    libxml2 \
    libxml2-dev \
    pandoc \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# 2. Copy DESCRIPTION file
COPY DESCRIPTION DESCRIPTION

# Optional: Set GitHub PAT to avoid API rate limits during package download
ARG GITHUB_PAT
ENV GITHUB_PAT=$GITHUB_PAT

# 3. Install remotes and install ONLY the dependencies declared in DESCRIPTION
RUN R -e 'install.packages("remotes", repos = "https://cloud.r-project.org")' \
 && R -e 'remotes::install_deps(dependencies = TRUE)'

# Copy lockfile and activation scripts

COPY renv.lock renv.lock
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json

# Explicitly pass only the package names you need

RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos = 'https://cloud.r-project.org')"

RUN R --no-save << 'EOF'
lock <- renv::lockfile_read("renv.lock")
pkgs <- c("DESeq2", "dplyr", "Boruta", "ranger", "tidyr", "tibble",
          "apeglm", "clusterProfiler", "msigdbr", "AnnotationDbi", "org.Hs.eg.db", "fgsea")

if (is.null(lock$Packages$limma)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install("limma", version = '3.68.5', ask = FALSE, update = FALSE)
  renv::restore(packages = pkgs, prompt = FALSE)
}else{
 renv::restore(packages = c(pkgs, "limma"), prompt = FALSE) 
}

EOF

# 4. Set up output directories and copy project files

RUN mkdir -p /app/results/boruta && chmod 777 /app/results /app/results/boruta
RUN mkdir -p /app/results/diffex_gsea_gse && chmod 777 /app/results /app/results/diffex_gsea_gse
RUN mkdir -p /app/results/diffex_gsea_target && chmod 777 /app/results /app/results/diffex_gsea_target
RUN mkdir -p /app/output_data && chmod 777 /app/output_data
COPY 7_isolated_functions/2_1_1_Boruta_surv_bin.R /app/7_isolated_functions/2_1_1_Boruta_surv_bin.R
COPY 7_isolated_functions/2_1_2_Boruta_surv_bin_pert.R /app/7_isolated_functions/2_1_2_Boruta_surv_bin_pert.R 
COPY 7_isolated_functions/5_2_2_gsea_val.R /app/7_isolated_functions/5_2_2_gsea_val.R
COPY 5_molecular_analysis /app/5_molecular_analysis
COPY output_data /app/output_data

CMD ["sleep", "infinity"]
