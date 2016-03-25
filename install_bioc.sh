#!/bin/bash

Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite("BiocInstaller"); biocLite("BSgenome.Hsapiens.UCSC.hg19"); biocLite("Biostrings"); install.packages("devtools"); biocLite("GenomicRanges"); devtools::install_github("mskilab/gUtils")'
