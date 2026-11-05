# R/Bioconductor + read-annotation tools for smrna-nextflow
# Bioconductor 3.20 (R 4.4.x)

FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL org.opencontainers.image.description="smrna-nextflow: R/Bioconductor, bowtie and bedtools"
LABEL org.opencontainers.image.source="https://github.com/damouzo/smrna-nextflow"
LABEL org.opencontainers.image.licenses="MIT"

RUN apt-get update && apt-get install -y --no-install-recommends \
    bowtie \
    bedtools \
    && rm -rf /var/lib/apt/lists/*

RUN R -e " \
    BiocManager::install(c( \
        'DESeq2', \
        'ashr', \
        'clusterProfiler', \
        'org.Hs.eg.db', \
        'org.Mm.eg.db' \
    ), ask=FALSE, update=FALSE) \
"

RUN R -e " \
    install.packages(c( \
        'ggplot2', 'ggrepel', 'pheatmap', 'data.table', \
        'optparse', 'jsonlite', 'yaml', 'rmarkdown' \
    ), repos='https://cloud.r-project.org', Ncpus=4) \
"

# Fail the build if a core package does not load
RUN R -e " \
    library(DESeq2); library(ashr); library(clusterProfiler); \
    library(ggplot2); library(data.table); library(optparse); \
    library(jsonlite); library(yaml); library(rmarkdown); \
    sessionInfo() \
" && bowtie --version | head -1 && bedtools --version

WORKDIR /data
