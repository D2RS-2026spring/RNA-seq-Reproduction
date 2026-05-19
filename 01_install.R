if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos="https://cran.rstudio.com/")
BiocManager::install(c("DESeq2","ggplot2","pheatmap","fastqcr"))