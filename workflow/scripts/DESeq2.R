options(warn = -1)

# 1.install packages ------------------------------------------------------

if(!require("tidyverse",quietly = TRUE)){
  install.packages("tidyverse",dependencies = TRUE,
                   repos = "http://cran.us.r-project.org")
}

library(tidyverse,quietly = TRUE) %>% suppressMessages()

# make sure you haved install the c compiler, such as gcc in linux
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repo = "https://lib.ugent.be/CRAN")
library(BiocManager) %>% suppressMessages()

if (!require("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")
library(DESeq2,quietly = TRUE) %>% suppressMessages()

# 2.Import files ----------------------------------------------------------
# import count matrix
transcript_count <- read.csv(snakemake@input[["transcript_count"]], row.names=1) %>% 
  filter(rowSums(., na.rm = TRUE)>0)
gene_count <- read.csv(snakemake@input[["gene_count"]], row.names=1) %>% 
  filter(rowSums(., na.rm = TRUE)>0)

# import group file
group <- read.table(snakemake@params[["sample_group"]],row.names=1,header = T)
names(group[,1]) <- "group"


# 3.Deseq2 analysis -------------------------------------------------------

# import data as deseq2 object
dds <- DESeqDataSetFromMatrix(countData = gene_count, colData = group, design = ~ group)
# deseq2 analysis
dds <- DESeq(dds)
# deseq2 result
res <- results(dds)
# extract diff gene list
diff_gene_deseq2 <- res[order(res$padj),] %>% as.data.frame
# export as .csv
write.csv(diff_gene_deseq2,file="results/quant/stringtie/diff_gene_deseq2.csv")