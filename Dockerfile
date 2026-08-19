FROM rocker/r-base:4.5.0

# 1. Install required system dependencies

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libglpk-dev \
    libpng-dev \
    libssl-dev \
    libx11-dev \
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

# 4. Set up output directories and copy project files
RUN mkdir -p /app/results/boruta && chmod 777 /app/results /app/results/boruta
COPY 2_Boruta_multiple/Boruta_surv_bin.R /app/2_Boruta_multiple/Boruta_surv_bin.R

CMD ["Rscript", "2_Boruta_multiple/Boruta_surv_bin.R"]
