#
library(limma)
library(edgeR)
library(DESeq2)
library(biomaRt)
library(VennDiagram)

# 2.Import files ----------------------------------------------------------
# import count matrix
transcript_count <- read.csv("D:/Database/snakemake_rna-seq/results/quant/stringtie/gene_count_matrix.csv", row.names = 1) %>%
  filter(rowSums(., na.rm = TRUE) > 0)
gene_count <- read.csv("D:/Database/snakemake_rna-seq/results/quant/stringtie/gene_count_matrix.csv", row.names = 1) %>%
  filter(rowSums(., na.rm = TRUE) > 0)

## import group file
group <- read.table("D:/Database/snakemake_rna-seq/data/sample_group.txt", row.names = 1, header = T)
names(group[, 1]) <- "group"
## set group 
levels <- group %>% distinct(group) %>% as_vector()
col.group <- factor(group$group,levels = levels)

# filter ------------------------------------------------------------------

#filter_gene <- function(gene_count, col.group){
  min.samples <- min(sapply(levels(col.group), function(x){length(which(col.group %in% x))}))
  rpm <- colSums(gene_count)/1000000
  filter_ind <- gene_count > rpm
  filter_ind_rowsums <- apply(filter_ind, 1, sum)
  return(gene_count[filter_ind_rowsums > min.samples,])
}
#gene_count <- filter_gene(gene_count, col.group)


#  diff gene analysis -----------------------------------------------------
## function
deseq2_quant <- function(gene_count, group.file, group){
  ## require(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = gene_count, colData = group.file, design = ~ group)
  # deseq2 analysis
  dds <- DESeq(dds)
  # deseq2 result
  res <- results(dds)
  res$disp <- dispersions(dds)
  diff_gene_deseq2 <- res
  #
  return(diff_gene_deseq2)
}
edgeR_quant <- function(gene_count, col.group){
  # transform into edgeR object
  deg <- DGEList(gene_count,group = col.group)
  
  ## group matrix
  design <- model.matrix(~ 0 + col.group)
  colnames(design) <- c("CK","Treat")
  ## diff compare matrix
  contrast.group <- makeContrasts(CK-Treat, levels=design)
  
  # TMM normalization
  deg <- calcNormFactors(deg)
  keep <- filterByExpr(deg, design)
  keep <- rowSums(cpm(deg) > 1) >= 1
  deg <- deg[keep,]
  
  ## the dispersal results combine:
  ## estimateGLMTagwiseDisp(), estimateGLMCommonDisp(), and estimateGLMTrendedDisp()
  deg <- estimateDisp(deg, design, robust = TRUE) 
  #plotBCV(deg)
  
  ## negative general binomal model
  fit <- glmFit(deg, design, robust = TRUE)
  lrt <- glmLRT(fit)
  diff_gene_edgeR <- topTags(lrt,n = nrow(deg),sort.by = "none")
  
  ## negative general log binomal model (more strict)
  #fit <- glmQLFit(deg, design, robust = TRUE)       
  #lrt <- glmQLFTest(fit)    
  #topTags(lrt)
  
  ## 
  #diff_gene_edgeR <- decideTestsDGE(lrt, adjust.method = 'fdr', p.value = 0.05)  
  
  return(diff_gene_edgeR)
}
limma_quant <- function(gene_count, col.group){
  # transform into edgeR object
  deg <- DGEList(gene_count)
  
  # TMM normalization
  deg <- calcNormFactors(deg)
  
  # low-expression threshold
  #cutoff <- 12 # set manually
  #cut <- which(apply(cpm(deg), 1, max) > cutoff)
  #deg <- deg[cut,] 
  
  ## generate group model
  design <- model.matrix(~0 + col.group)
  colnames(design) <- c("CK","Treat")
  
  ## diff compare matrix
  contr.group <- makeContrasts(CK-Treat, levels=design)

  ## normalization
  dge_norm <- voom(deg, design, plot = T)
  fit <- lmFit(dge_norm, design)
  # coef(fit) %>% head()
  
  ## group compare of each gene
  contr.fit <- contrasts.fit(fit, contr.group)
  ## shrinks standard errors using empirical bayes
  contr.fit.bayes <- eBayes(contr.fit)
  
  ## diff gene list
  diff_gene_limma <- topTable(contr.fit.bayes, sort.by = "none", n = Inf)
  return(diff_gene_limma)
}

## deseq2
res_deseq2 <- deseq2_quant(gene_count,group,group)

## edgeR 
res_edgeR <- edgeR_quant(gene_count,col.group)

## limma
res_limma <- limma_quant(gene_count,col.group)

all(rownames(res_deseq2)==rownames(res_limma))


# benchmark  --------------------------------------------------------------
res_summary <- data.frame(deseq2_padj=res_deseq2$padj,
                          edgeR_padj=res_edgeR$table$FDR,
                          limma_padj=res_limma$adj.P.Val,
                          deseq2_logfc=res_deseq2$log2FoldChange,
                          edgeR_logfc=res_edgeR$table$logFC,
                          limma_logfc=res_limma$logFC,
                          mean_exprs=log10(gene_count+1) %>% rowMeans()
)

## for plotting purposes
breaks <- gene_count %>% as.matrix() %>% as.numeric() %>% 
  max() %>% log10() %>% seq(0,., 0.01)

## cumulative significant genes vs. log10 mean expression
min.fdr <- 1
cgn_deseq2 <- sapply(breaks,
                     function(x) nrow(subset(res_summary, mean_exprs < x & deseq2_padj < min.fdr)))

cgn_edgeR <- sapply(breaks,
                    function(x) nrow(subset(res_summary, mean_exprs < x & edgeR_padj < min.fdr)))

cgn_limma <- sapply(breaks,
                    function(x) nrow(subset(res_summary, mean_exprs < x & limma_padj < min.fdr)))
#

cum_numbers <- data.frame(breaks,cgn_deseq2,cgn_edgeR,cgn_limma)
cum_proportions <- cum_numbers %>% 
  transmute(breaks = breaks,
            cgp_deseq2 = cgn_deseq2/max(cgn_deseq2),
            cgp_edgeR = cgn_edgeR/max(cgn_edgeR),
            cgp_limma = cgn_limma/max(cgn_limma)
            )

##
cum_numbers <- pivot_longer(cum_numbers,cols = 2:4,
                            names_to = "methods",
                            values_to = "cum_numbers")

cum_proportions <- pivot_longer(cum_proportions,cols = 2:4,
                         names_to = "methods",
                         values_to = "cum_proportions")

##

ggplot(cum_numbers,aes(breaks,cum_numbers,color = methods))+geom_line()+
  labs(x = "log10 expression",
       y = paste0("Cumulative significant genes (FDR < ",min.fdr,")"))+
  theme_bw()

ggplot(cum_proportions,aes(breaks,cum_proportions,color = methods))+geom_line()+
  labs(x = "log10 expression",
       y = paste0("Cumulative significant genes (FDR < ",min.fdr,")"))+
  theme_bw()

##
hist(subset(res_summary, limma_padj < min.fdr)[, "mean_exprs"], 
     freq=FALSE,
     col=rgb(1,0,0,1/4), xlab = "log10 expression",
     ylab = "number of sig. genes (fdr<0.25)",
     main = "Proportion of number of significant genes by log10 count", xaxt = "n")

hist(subset(res_summary, edgeR_padj < min.fdr)[, "mean_exprs"],freq=FALSE,
     col=rgb(0,0,1,1/4), add = T)
labels <- sapply( 1:5,function(i) as.expression(bquote(10^ .(i))) )
axis(1,at=1:5,labels=labels)
box()
legend("topright", c("DEseq2", "edgeR"),
       col=c(rgb(1,0,0,1/4), rgb(0,0,1,1/4), rgb(1,1,0,1/4)), lwd=10)

labels <- sapply(1:5,function(i) as.expression(bquote(10^ .(i))))


