# https://www.cnblogs.com/zhanmaomao/p/11529589.html
# https://www.jianshu.com/p/133257d344d3
#
library(clusterProfiler) # GSEA
library(topGO) # GO plot
library(AnnotationHub) #
library(pathview) # KEGG pathway plot
library(gage) # # KEGG pathway
library(gageData)
BiocManager::install("gageData")

#
# change the cache dir if needed
# rappdirs::user_cache_dir(appname="AnnotationHub") # old version cache
# tools::R_user_dir("AnnotationHub", which="cache") # new version cache

# 2.get orgDb

##
hub <- AnnotationHub()  # search all hub
hub$dataprovider %>% unique() # database sourece
query(hub, "Apis_cerana")  # interest species
tar_org <- hub[["AH86625"]] # download org of interest species

##
# glimpse all key
keytypes(tar_org)

# glimpse gene id in specific key
keys(tar_org,keytype = "SYMBOL") %>% head(10) #

# search corresponding database id according to current gene id
select(tar_org, keys= "ATP6",keytype = "SYMBOL", columns=c("ENTREZID","GO"))

# batch search
kk <- keys(tar_org, keytype = "SYMBOL")
goAnno <- select(tar_org, keys = kk, keytype = "SYMBOL", 
                 columns = c("GOALL","ENTREZID", "ONTOLOGYALL"))

# gene name of sig differential gene
gene_id<-sig_diff_geneset[,1] 

# convert gene id to ENTREZID
## 1.mapIDs
DEG.entrez_id = mapIds(x = tar_org,          # species org
                       keys = gene_id,       # diff gene id
                       keytype = "SYMBOL",   # gene id type
                       column = "ENTREZID")  # convert to ENTREZID

DEG.entrez_id = na.omit(DEG.entrez_id)

## bitr
gene.df<-bitr(gene_id, fromType = "ENSEMBL", 
              toType = c("SYMBOL","ENTREZID"),
              OrgDb = tar_org)

# GO 
## biological process
erich.go.BP = enrichGO(gene = DEG.entrez_id,      ###差异基因ID
                       OrgDb = tar_org,           ###数据库
                       keyType = "ENTREZID",       ##基因ID类型
                       ont = "BP",                 #:BP,MF,CC
                       pvalueCutoff = 0.5,         ###fisher检验对p值
                       qvalueCutoff = 0.5)         ###对p值进行校对对q值，一般大于p值

## plot
dotplot(erich.go.BP)
barplot(erich.go.CC)+  scale_size(range=c(2, 12))+
  scale_x_discrete(labels=function(ego_bp) str_wrap(ego_bp,width = 25))
# tree plot are very large
pdf(file="./enrich.go.bp.tree.pdf",width = 10,height = 15)
plotGOgraph(erich.go.BP)
dev.off()


# KEGG pathway
enrich.KEGG.BP <- enrichKEGG(gene = diff_gene_id,   # diff gene id (ENTREZID)
                             keyType = "kegg",     # key type
                             organism = "soe",     # species 
                             pvalueCutoff = 0.05,                
                             pAdjustMethod = "BH",
                             qvalueCutoff = 0.1,)

# plot
barplot(enrich.KEGG.BP, showCategory = 10)
dotplot(enrich.KEGG.BP, showCategory = 10)
cnetplot(enrich.KEGG.BP, showCategory = 5)
cnetplot(enrich.KEGG.BP,circular=T,colorEdge=T) 

# Enrichment Map
emapplot(enrich.KEGG.BP)

# path plot

## browse KEGG
browseKEGG(enrich.KEGG.BP,"soe00564")    ###画出某一特定pathway的图
## pathview
pathview(gene.data = gene_id,     # gene id (Entrez_ID)
         pathway.id = "soe00564", # KEGG id
         species = "soe",
         kegg.native = TRUE,      # TRUE: integrate pathway; False: only input gene pathway
)

# Gene Set Enrichment Analysis (GSEA)
# gene list 
genelist <- sig_diff_geneset$log2FoldChange
names(genelist) <- sig_diff_geneset$ENSEMBL

# sort 
genelist <- sort(genelist, decreasing = TRUE)

# GSEA
gsemf <- gseGO(genelist,
               OrgDb = org.Mm.eg.db,
               keyType = "ENSEMBL",
               ont="BP"
)

# Plot
gseaplot(gsemf, geneSetID="GO:0001819")

#
gene.df<-bitr(gene, fromType = "ENSEMBL", 
              toType = c("SYMBOL","ENTREZID"),
              OrgDb = tar_org)

foldchanges = sig_diff_genese$log2FoldChange
names(foldchanges)= gene.df$ENTREZID


# pathway analysis --------------------------------------------------------
#
keggres = gage(foldchanges, gsets = kegg.sets.mm, same.dir = TRUE)

# Look at both up (greater), down (less), and stat.
lapply(keggres, head)

#
kegg_path = data.frame(id=rownames(keggres$greater),
                             keggres$greater) %>% 
  tbl_df() %>% 
  filter(row_number() <= 10) %>% 
  .$id %>%
  as.character()
  
# Get the IDs.
keggresids = substr(keggrespathways, start=1, stop=8)

# plot function
plot_pathway = function(pid) {
  pathview(gene.data=foldchanges, pathway.id=pid, species="mmu", new.signature=FALSE)
}
# batch plot
batch_plot = sapply(keggresids, plot_pathway)
