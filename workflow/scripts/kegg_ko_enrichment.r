# https://www.cnblogs.com/zhanmaomao/p/11529589.html
# https://www.jianshu.com/p/133257d344d3
#
library(tidyverse)
library(clusterProfiler) # GSEA
library(topGO) # GO plot
library(AnnotationHub) #
library(pathview) # KEGG pathway plot
library(gage) # # KEGG pathway
library(gageData)
library(DOSE)
library(BiocManager)
# load results of DEG analysis (DESeq2, limma etc.)
diff_gene_deseq2 <- read_csv("C:/Users/grasslab/Desktop/snakemake/One-step_RNA-seq_analysis/results/quant/stringtie/diff_gene_deseq2.csv")
## gene name of sig differential gene
gene_id <- diff_gene_deseq2[, 1] %>% as.matrix() %>% as.character()

# 1.get org.db ------------------------------------------------------
## 1.1 search your species in the org.db
# change the cache dir if needed
# rappdirs::user_cache_dir(appname="AnnotationHub") # old version cache
# tools::R_user_dir("AnnotationHub", which="cache") # new version cache

hub <- AnnotationHub()  # search all hub
hub$dataprovider %>% unique() # database sourece
query(hub, "Apis_cerana")  # interest species
tar_org <- hub[["AH86625"]] # download org of interest species

## 1.2.install the org.db you made manually
install.packages("F:/RNA-seq_reads/paried/alfafa_zm4_reference_genome/org.Msativa.eg.db/",repos = NULL, type = "source") 
library(org.Msativa.eg.db)
tar_org <- org.Msativa.eg.db

### glimpse all key
keytypes(tar_org)

### glimpse different keytypes
keys(tar_org, keytype = "GID") %>% head(10) #

# search corresponding database id according to current gene id
goAnno <- select(tar_org, keys = gene_id,keytype = "GID", columns = c("GO","Ko","Pathway"))

# 2.convert gene id to ENTREZID (optional for model species) ----------
## 2.1.mapIDs
# DEG.entrez_id <- mapIds(x = tar_org,          # species org
#                        keys = gene_id,       # diff gene id
#                        keytype = "SYMBOL",   # gene id type
#                        column = "ENTREZID")  # convert to ENTREZID

# DEG.entrez_id <- na.omit(DEG.entrez_id)

# ## 2.2.bitr
# gene.df <- bitr(gene_id, fromType = "ENSEMBL",
#               toType = c("SYMBOL", "ENTREZID"),
#               OrgDb = tar_org)

# 3.GO analysis -------------------------------------------------------
## 3.1.biological process
erich.go.BP <- enrichGO(gene = gene_id,       # diff gene id
                        OrgDb = tar_org,            # org db
                        keyType = "GID",       # gene id type
                        ont = "BP",                 # BP,MF,CC
                        pvalueCutoff = 0.5,         # fisher
                        qvalueCutoff = 0.5)         # p-adj

## 3.2.plot
dotplot(erich.go.BP)
barplot(erich.go.BP)

## 3.3 tree plot are very large
pdf(file="./enrich.go.bp.tree.pdf",width = 10,height = 15)
plotGOgraph(erich.go.BP)
dev.off()

## 3.4 network graph
### network graph
enrichMap(erich.go.BP, vertex.label.cex=1.2, layout = igraph::layout.kamada.kawai)

### GO plot
plotGOgraph(erich.go.BP)

# 4.KEGG pathway -------------------------------------------------------
# ## for model species using ENTREZID
# enrich.KEGG.BP <- enrichKEGG(gene = gene_id,        # diff gene id (ENTREZID)
#                              keyType = "kegg",      # key type
#                              organism = "soe",      # species
#                              pvalueCutoff = 0.05,
#                              pAdjustMethod = "BH",
#                              qvalueCutoff = 0.1)

## non model species using your KEGG files
### 1.import kegg files
pathway2gene <- read.table("F:/RNA-seq_reads/paried/alfafa_zm4_reference_genome/alfalfa_pathway2gene", header = T, sep = "\t")
pathway2name <- read.table("F:/RNA-seq_reads/paried/alfafa_zm4_reference_genome/alfalfa_pathway2name", header = T, sep = "\t")

### 2.enricher
enrich_kegg <- enricher(gene_id,
                TERM2GENE = pathway2gene,
                TERM2NAME = pathway2name,
                pvalueCutoff = 1,
                qvalueCutoff = 1,
                pAdjustMethod = "BH",
                minGSSize = 1)
dotplot(enrich_kegg)

## plot
barplot(enrich_kegg, showCategory = 10)
dotplot(enrich_kegg, showCategory = 10)
cnetplot(enrich_kegg, showCategory = 5)
cnetplot(enrich_kegg, circular = T, colorEdge = T) 

## Enrichment Map
pairwise_ter
emapplot(enrich_kegg)

## path plot

## browse KEGG
browseKEGG(enrich_kegg,"soe00564")    #plot specific pathway

## pathview
pathview(gene.data = gene_id,     # gene id (Entrez_ID)
         pathway.id = "soe00564", # KEGG id
         species = "soe",
         kegg.native = TRUE,      # TRUE: integrate pathway; False: only input gene pathway
)

# 5.Gene Set Enrichment Analysis (GSEA) -------------------------------------------------------
## 5.1.gene list 
genelist <- sig_diff_geneset$log2FoldChange
names(genelist) <- sig_diff_geneset$ENSEMBL

## 5.2.sort 
genelist <- sort(genelist, decreasing = TRUE)

## 5.3.GSEA 
gsemf <- gseGO(genelist,
               OrgDb = org.Mm.eg.db,
               keyType = "ENSEMBL",
               ont="BP"
)

## 5.4.Plot
gseaplot(gsemf, geneSetID="GO:0001819")

## 5.5.convert id
gene_df <- bitr(gene,
  fromType = "ENSEMBL", toType = c("SYMBOL", "ENTREZID"),
  OrgDb = tar_org
)

foldchanges <- sig_diff_genese$log2FoldChange
names(foldchanges) <- gene_df$ENTREZID


# 6.Pathway analysis --------------------------------------------------------
##
keggres <- gage(foldchanges, gsets = kegg.sets.mm, same.dir = TRUE)

## Look at both up (greater), down (less), and stat.
lapply(keggres, head)

##
kegg_path = data.frame(id = rownames(keggres$greater), keggres$greater) %>%
  filter(row_number() <= 10) %>%
  .$id %>%
  as.character()

## Get the IDs.
keggresids <- substr(keggrespathways, start = 1, stop = 8)

## plot function
plot_pathway = function(pid) {
  pathview(gene.data = foldchanges, pathway.id = pid, species = "mmu", new.signature = FALSE)
}
## batch plot
batch_plot <- sapply(keggresids, plot_pathway)
