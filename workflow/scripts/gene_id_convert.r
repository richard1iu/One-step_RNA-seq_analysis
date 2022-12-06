library(dplyr)
library(rtracklayer)

# input gene_count matrix
expr.gene <- read.csv("C:/Users/grasslab/Desktop/snakemake/One-step_RNA-seq_analysis/results/quant/stringtie/gene_count_matrix.csv")

# input annotated.gtf
ensembl_anno <- rtracklayer::import("C:/Users/grasslab/Desktop/snakemake/One-step_RNA-seq_analysis/results/quant/stringtie/gffcompare/gffcompare.annotated.gtf") %>% as.data.frame()

# combine two file
anno_result <- dplyr::left_join(expr.gene,ensembl_anno[,(11:14)],by ="gene_id")

# remove duplicate gene id
anno_result <- anno_result[!duplicated(anno_result$gene_id),]



# Deseq2 -----------------------------------------------------------------
library(DESeq2)

database <- as.matrix(read.csv("transcript_count_matrix.csv", row.names="transcript_id"))
condition <- factor(c("control","control","KD","KD"))
coldata <- data.frame(row.names = colnames(database), condition)
countData <- countData[, rownames(colData)]
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
resordered <- res[order(res$padj),]
summary(res)
write.csv(as.data.frame(resordered),file="results.csv")
