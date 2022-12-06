library(edgeR)


## set group 
levels <- group %>% distinct(group) %>% as_vector()
col.group <- factor(group$group,levels = levels)

# transform into edgeR object
deg <- DGEList(gene_count,group = col.group)

# TMM normalization
deg <- calcNormFactors(deg)
keep <- filterByExpr(deg, design)
keep <- rowSums(cpm(deg) > 1) >= 1
deg <- deg[keep,]

## group matrix
design <- model.matrix(~ 0 + col.group)
names(design) <- c("CK","Treat")
## diff compare matrix
contrast.group <- makeContrasts(CK-Treat, levels=design)

## the dispersal results combine:
## estimateGLMTagwiseDisp(), estimateGLMCommonDisp(), and estimateGLMTrendedDisp()
deg <- estimateDisp(deg, design, robust = TRUE) 
plotBCV(deg)

## negative general binomal model
fit <- glmFit(deg, design, robust = TRUE)
lrt <- glmLRT(fit)
topTags(lrt)

## negative general log binomal model (more strict)
#fit <- glmQLFit(deg, design, robust = TRUE)       
#lrt <- glmQLFTest(fit)    
#topTags(lrt)

## 
diff_gene_edgeR <- decideTestsDGE(lrt, adjust.method = 'fdr', p.value = 0.05)  

summary(diff_gene_edgeR)
?decideTestsDGE


