options(warn = -1)

# 1.install packages ------------------------------------------------------

if(!require("tidyverse", quietly = TRUE)){
  install.packages("tidyverse",dependencies = TRUE,
                   repos = "http://cran.us.r-project.org")
}
library(tidyverse, quietly = TRUE)

if (!require("data.table", quietly = TRUE))
  BiocManager::install("data.table")
library(data.table, quietly = TRUE) %>% suppressMessages()

# make sure you haved install the c compiler, such as gcc in linux
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repo = "https://lib.ugent.be/CRAN") %>% 
  suppressMessages()
library(BiocManager) %>% suppressMessages() 

if (!require("apeglm", quietly = TRUE))
  BiocManager::install("apeglm") %>% suppressMessages()
library(apeglm, quietly = TRUE) %>% suppressMessages()

if (!require("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2") %>% suppressMessages()
library(DESeq2, quietly = TRUE) %>% suppressMessages()

if (!require("Cairo", quietly = TRUE))
  BiocManager::install("Cairo") %>% suppressMessages()
library(Cairo, quietly = TRUE) %>% suppressMessages()

# 2.Import files ----------------------------------------------------------
# import count matrix
transcript_count <- fread(snakemake@input[["transcript_count"]]) %>%
  as_tibble() %>%
  column_to_rownames("transcript_id") %>%
  filter(rowSums(., na.rm = TRUE) > 0)

gene_count <- fread(snakemake@input[["gene_count"]]) %>%
  as_tibble() %>%
  column_to_rownames("gene_id") %>%
  filter(rowSums(., na.rm = TRUE) > 0)

# import group file
group <- read.table(snakemake@params[["sample_group"]], row.names = 1, header = T)
names(group[, 1]) <- "group"


# 3.Deseq2 analysis -------------------------------------------------------

# import data as deseq2 object
dds <- DESeqDataSetFromMatrix(countData = gene_count, colData = group, design = ~ group)
# deseq2 analysis
dds <- DESeq(dds)
# deseq2 result
res <- results(dds)
# extract diff gene list
all_gene_deseq2 <- res[order(res$padj), ] %>% as.data.frame()
diff_gene_deseq2 <- all_gene_deseq2 %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>% 
  merge(gene_count, by = "row.names") %>% 
  mutate(sig = case_when(log2FoldChange > 2  & padj < 0.05 ~ "up",
                         log2FoldChange < -2 & padj < 0.05 ~ "down",
                         TRUE ~ "insig"))

# export as .csv
write.csv(all_gene_deseq2, file = "results/quant/stringtie/all_gene_deseq2.csv")
write.csv(diff_gene_deseq2, file = "results/quant/stringtie/diff_gene_deseq2.csv")


# (optional) shrinkage if lot's of low abundance gene
# it will reduce foldchange but don't modify padj
shf_count <- lfcShrink(dds, coef = resultsNames(dds)[2], type = "apeglm")

# (optional) normalization
normal_count <- counts(dds, normalized = TRUE) %>% as.data.frame()

# (optional) log2 transfer: samples >30 using vst
log2_count <- rlog(dds)


# simple plot -------------------------------------------------------------

# PCA plot
p1 <- DESeq2::plotPCA(log2_count, intgroup = "group", ntop = 500) +
  theme_bw() +
  geom_jitter(size = 3) +
  scale_y_continuous(limits = c(-5, 5)) +
  ggtitle(label = "Principal Component Analysis (PCA)",
          subtitle = "Top 500 most variable genes") 
ggsave("DESeq2_PCA.pdf",p1)

# mean-difference plot
p2 <- DESeq2::plotMA(res, ylim = c(-10, 10), alpha = 0.1, main = "MA plot: ")
ggsave("DESeq2_plotMA.pdf",p2)

#
p3 <- DESeq2::plotDispEsts(dds)
ggsave("DESeq2_DispEsts.pdf",p3)

# volcano plot ------------------------------------------------------------
p4 <- ggplot(data = diff_gene_deseq2, 
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

ggsave("DESeq2_volcano.pdf",p4)


# heatmap -----------------------------------------------------------------

diff_sig_gene <- diff_gene_deseq2 %>% filter(sig != "insig")
heatmap_data <- log2_count %>% assay() %>% as.data.frame %>% .[diff_sig_gene$Row.names,]

## group
ann_group <- data.frame(
  Group = factor(colData(log2_count)$group),
  row.names = colData(log2_count) %>% row.names
)

## color
ann_colors <- list(
  Group = c(CK = "lightblue", Treat = "darkorange")
)

# 1. creat a cairo
Cairo::CairoPDF( 
  filename = "DESeq2_heatmap.pdf",
  width = 12,
  height = 12,
  units = "in",
  dpi = 300)     

## heatmap plot
pheatmap(mat = heatmap_data,
         color = colorRampPalette(brewer.pal(9, "YlOrBr"))(255),
         scale = "row",
         annotation_col = ann_group,
         annotation_colors = ann_colors,
         fontsize = 6.5,
         cellwidth = 55,
         show_colnames = F)

# close
dev.off() 

