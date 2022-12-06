library(edgeR)
library(RColorBrewer)
library(limma)


## set group 
levels <- group %>% distinct(group) %>% as_vector()
col.group <- factor(group$group,levels = levels)


# transform into edgeR object
deg <- DGEList(gene_count)

# TMM normalization
deg <- calcNormFactors(deg)

# low-expression threshold
cutoff <- 12 # set manually
cut <- which(apply(cpm(deg), 1, max) > cutoff)
deg <- deg[cut,] 

# Multidimensional scaling (MDS) plot -------------------------------------
## set color
color <- factor(group$group,levels = levels)
levels(color) <- brewer.pal(nlevels(color), "Set1") 
color <- as.character(color)
## plot MDS
plotMDS(deg, labels=col.group, col=color) 
title(main="Sample MDS")

# Find diff gene using voom -----------------------------------------------

## voom transform and variance weight
## counts to log2 CPM (counts per million reads)£¬
## per million reads are determined by the norm.factors of calcNormFactors
## fit modol of log2CPM and estimate the residual
## mean expression (red line) sqrt(residual standard deviation)£»

## generate group model
design <- model.matrix(~0 + col.group)
colnames(design) <- c("CK","Treat")

## diff compare matrix
contr.group <- makeContrasts(CK-Treat, levels=design)

### interactions of diverse group
# design <- model.matrix(~col.group * col.group2)
### batch effects
# design <- model.matrix(~0 + col.group + col.batch)

## normalization
dge_norm <- voom(deg, design, plot = T)
fit <- lmFit(dge_norm, design)
coef(fit) %>% head()

## group compare of each gene
contr.fit <- contrasts.fit(fit, contr.group)
## shrinks standard errors using empirical bayes
contr.fit.bayes <- eBayes(contr.fit)

## diff gene list
diff_gene <- topTable(contr.fit.bayes, sort.by = "P", n = Inf)
diff_gene %>% filter(adj.P.Val < 0.05 & abs(logFC)>2) %>% dim()

diff_gene_limma <- diff_gene %>% as.data.frame() %>% 
  merge(gene_count, by = "row.names") %>% 
  mutate(sig = case_when(logFC > 1 & adj.P.Val < 0.1 ~"up",
                         logFC < -1 & adj.P.Val < 0.1~"down",
                         TRUE~"insig"))

# export as .csv
write.csv(diff_gene_limma, file = "results/quant/stringtie/diff_gene_limma.csv")

# Plot --------------------------------------------------------------------

ggplot(data = diff_gene_limma, 
       aes(x = logFC, y = -1 * log10(adj.P.Val))) +
  geom_point(aes(color = sig)) +
  scale_color_manual(values = c("#546de5", "#d2dae2","#ff4757")) +
  labs(title = "Volcano Plot: ",
       x = expression(log[2](FC)),
       y = expression(-log[10](padj))) +
  geom_hline(yintercept = 0.975, linetype = 4) + 
  geom_vline(xintercept = c(-1, 1), linetype = 4) +
  theme_bw() +
  theme(panel.grid = element_blank())

## 
diff_sig_gene <- diff_gene_limma %>% filter(sig != "insig")
heatmap_data <- gene_count %>% .[diff_sig_gene$Row.names,]

## group
ann_group = data.frame(
  Group = factor(group$group),
  row.names = group %>% row.names
)

## color
ann_colors = list(
  Group = c(CK = "lightblue", Treat = "darkorange")
)

## heatmap
pheatmap(mat = heatmap_data[1:20,], 
         color = colorRampPalette(brewer.pal(9, "YlOrBr"))(255), 
         scale = "row",
         annotation_col = ann_group, 
         annotation_colors = ann_colors,
         fontsize = 6.5, 
         cellwidth = 55,
         show_colnames = F)
