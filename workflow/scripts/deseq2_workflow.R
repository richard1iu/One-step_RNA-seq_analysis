options(warn = -1)

# 1.install packages ------------------------------------------------------

if (!require("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse",
    dependencies = TRUE,
    repos = "http://cran.us.r-project.org"
  )
}
library(tidyverse, quietly = TRUE)

# make sure you haved install the c compiler, such as gcc in linux
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repo = "https://lib.ugent.be/CRAN")
}
library(BiocManager)

if (!require("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}
library(DESeq2, quietly = TRUE)

if (!require("apeglm", quietly = TRUE)) {
  BiocManager::install("apeglm")
}
library(apeglm, quietly = TRUE)
library(pheatmap)
library(RColorBrewer)
BiocManager::install("pathview")
# 2.Import files ----------------------------------------------------------
# import count matrix
transcript_count <- read.csv("D:/Database/snakemake_rna-seq/results/quant/stringtie/gene_count_matrix.csv", row.names = 1) %>%
  filter(rowSums(., na.rm = TRUE) > 0)
gene_count <- read.csv("D:/Database/snakemake_rna-seq/results/quant/stringtie/gene_count_matrix.csv", row.names = 1) %>%
  filter(rowSums(., na.rm = TRUE) > 0)

# import group file
group <- read.table("D:/Database/snakemake_rna-seq/data/sample_group.txt", row.names = 1, header = T)
names(group[, 1]) <- "group"

# 3.Deseq2 analysis -------------------------------------------------------

# import data as deseq2 object
dds <- DESeqDataSetFromMatrix(countData = gene_count, colData = group, design = ~group)
# deseq2 analysis
dds <- DESeq(dds)

# deseq2 result
res <- results(dds, pAdjustMethod = "fdr", alpha = 0.05)
summary(res)

# extract diff gene list
diff_gene_list <- res %>% as.data.frame() %>% 
  merge(gene_count, by = "row.names") %>% 
  mutate(sig = case_when(log2FoldChange > 1 & pvalue < 0.1 ~"up",
                         log2FoldChange < -1 & pvalue < 0.1~"down",
                         TRUE~"insig"))
# export as .csv
write.csv(diff_gene_list, file = "results/quant/stringtie/diff_gene_list.csv")

# (optional) shrinkage if lot's of low abundance gene
# it will reduce foldchange but don't modify padj
shf_count <- lfcShrink(dds, coef = resultsNames(dds)[2], type = "apeglm")

# (optional) normalization
normal_Count <- counts(dds, normalized = TRUE) %>% as.data.frame()

# (optional) log2 transfer: samples >30 using vst
log2_count <- rlog(dds)


# simple plot -------------------------------------------------------------

# PCA plot
DESeq2::plotPCA(log2_count, intgroup = "group", ntop = 500)+
  theme_bw() +
  geom_jitter(size = 3) + 
  scale_y_continuous(limits = c(-5, 5)) +
  ggtitle(label = "Principal Component Analysis (PCA)", 
          subtitle = "Top 500 most variable genes") 

# mean-difference plot
DESeq2::plotMA(res, ylim = c(-10, 10), alpha = 0.1, main = "MA plot: ")

#
DESeq2::plotDispEsts(dds)


# volcano plot ------------------------------------------------------------
ggplot(data = diff_gene_list, 
       aes(x = log2FoldChange, y = -1 * log10(pvalue))) +
  geom_point(aes(color = sig)) +
  scale_color_manual(values = c("#546de5", "grey", "red")) +
  labs(title = "Volcano Plot: ",
       x = expression(log[2](FC)),
       y = expression(-log[10](padj))) +
  geom_hline(yintercept = 0.975, linetype = 4) + 
  geom_vline(xintercept = c(-1, 1), linetype = 4) +
  theme_bw() +
  theme(panel.grid = element_blank())


# heatmap -----------------------------------------------------------------

diff_sig_gene <- diff_gene_list %>% filter(sig != "insig")
heatmap_data <- log2_count %>% assay() %>% as.data.frame %>% .[diff_sig_gene$Row.names,]

## group
ann_group = data.frame(
  Group = factor(colData(log2_count)$group),
  row.names = colData(log2_count) %>% row.names
)

## color
ann_colors = list(
  Group = c(CK = "lightblue", Treat = "darkorange")
)

## heatmap
pheatmap(mat = heatmap_data, 
         color = colorRampPalette(brewer.pal(9, "YlOrBr"))(255), 
         scale = "row",
         annotation_col = ann_group, 
         annotation_colors = ann_colors,
         fontsize = 6.5, 
         cellwidth = 55,
         show_colnames = F)


# Ensembl id transfer (optional) ---------------------------------------------------------------
readCount <- gene_count %>%
  dplyr::mutate(gene_id = rownames(gene_count)) %>%
  tidyr::separate(col = gene_id, into = c("ensembl_gene_id", "gene_name"), sep = "\\|", remove = TRUE) %>%
  dplyr::select(ensembl_gene_id, rownames(sampleGroup)) %>%
  dplyr::mutate(ensembl_gene_id = map_chr(ensembl_gene_id, geneID)) %>%
  as.data.frame() %>%
  column_to_rownames("ensembl_gene_id")

#
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
geneMap <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
  filters = "ensembl_gene_id",
  values = rownames(readCount), mart = ensembl
) %>%
  dplyr::filter(!(is.na(hgnc_symbol) & is.na(entrezgene_id))) %>%
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)

# add genome database
library(org.Mm.eg.db) 

# add gene name
res$description <- mapIds(x = org.Mm.eg.db,
                          keys = row.names(res),
                          column = "GENENAME",
                          keytype = "SYMBOL",
                          multiVals = "first")

# add gene symbol
res$symbol <- row.names(res)

# add ENTREZ ID
res$entrez <- mapIds(x = org.Mm.eg.db,
                     keys = row.names(res),
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")

# add ENSEMBL
res$ensembl <- mapIds(x = org.Mm.eg.db,
                      keys = row.names(res),
                      column = "ENSEMBL",
                      keytype = "SYMBOL",
                      multiVals = "first")


# function analysis -------------------------------------------------------
sig_entrez <- subset(res, is.na(entrez) == FALSE)
gene_matrix <- sig_entrez$log2FoldChange
names(gene_matrix) <- sig_entrez$entrez

# kegg enrich
kegg_enrich <- enrichKEGG(gene = names(gene_matrix),
                          organism = 'mouse',
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.10)

barplot(kegg_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "KEGG Enrichment Pathways",
        font.size = 8)

# go enrich 
go_enrich <- enrichGO(gene = names(gene_matrix),
                      OrgDb = 'org.Mm.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)

# pathview
## pathway.id: KEGG pathway identifier
pathview(gene.data = gene_matrix, 
         pathway.id = "04070", 
         species = "mouse")
