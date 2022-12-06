# https://www.jianshu.com/p/f2e4dbaae719

# set
options(stringsAsFactors = F)
options(warn = -1)

# load packages
library(dplyr)
library(stringr)
library(jsonlite)
library(AnnotationForge)

# import protein annotation file from output of eggnog or other mapper software
emapper <- read.table("alfafa.annotations", header=TRUE, sep = "\t",quote = "")
emapper[emapper==""]<-NA

# 1.Extract GO and KEGG info -----------------------------------------------

# split single row to multiple row
gos_list <- function(x){
  the_gos <- str_split(x[2], ",", simplify = FALSE)[[1]]
  df_temp <- data.frame(GID = rep(x[1], length(the_gos)),
                        GO = the_gos,
                        EVIDENCE = rep("IEA", length(the_gos))
  )
  return(df_temp)
}

# extract GO
gene2gol <- emapper %>% dplyr::select(query, GOs) %>% na.omit() %>% as.matrix() %>% apply(.,1,gos_list)
gene2gol_df <- do.call(rbind.data.frame,gene2gol)
gene2go <- gene2gol_df
gene2go$GO[gene2go$GO=="-"]<-NA
gene2go<-na.omit(gene2go)

# extract KEGG
gene2ko <- emapper %>% dplyr::select(GID = query, Ko = KEGG_ko)
gene2ko$Ko[gene2ko$Ko=="-"]<-NA
gene2kol<- gene2ko %>% na.omit() %>% as.matrix() %>% apply(.,1,gos_list)

gene2kol_df <- do.call(rbind.data.frame, gene2kol)
gene2ko <- gene2kol_df[,1:2]
colnames(gene2ko) <- c("GID","Ko")
gene2ko$Ko <- gsub("ko:","",gene2ko$Ko)

# 2.kegg function -----------------------------------------------------------

# import KEGG json file (https://www.genome.jp/kegg-bin/get_htext?ko00001)
kegg_json <- fromJSON("F:/Database/ko00001.json")

# kegg function
update_kegg <- function(kegg_json = kegg_json) {
  # organize data set
  pathway2name <- tibble(Pathway = character(), Name = character())
  ko2pathway <- tibble(Ko = character(), Pathway = character())
  
  #
  kegg <- kegg_json
  
  #
  for (a in seq_along(kegg[["children"]][["children"]])) {
    A <- kegg[["children"]][["name"]][[a]]
    
    for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
      B <- kegg[["children"]][["children"]][[a]][["name"]][[b]]
      
      for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
        
        pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
        pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
        pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
        pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
        kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
        kos <- str_match(kos_info, "K[0-9]*")[,1]
        ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))
        }
      }
    }
  # output kegg as .RData
  save(pathway2name, ko2pathway, file = "kegg_info.RData")
}

# update kegg database
update_kegg()

# load kegg_info.RData
load(file = "kegg_info.RData")

# combine ko and pathway
gene2pathway <- gene2ko %>% 
  left_join(ko2pathway, by = "Ko") %>% 
  dplyr::select(GID, Pathway) %>% 
  na.omit()

# distinct your gene2pathway
# gene2go <- unique(gene2go)
# gene2go <- gene2go[!duplicated(gene2go),]
# gene2ko <- gene2ko[!duplicated(gene2ko),]
# gene2pathway <- gene2pathway[!duplicated(gene2pathway),]
# gene_info <- gene_info[!duplicated(gene_info),]

# 3.make org ----------------------------------------------------------------

# org info: search the name in ncbi(https://www.ncbi.nlm.nih.gov/taxonomy)
# or you can define it as you wish
taxid = 3879
genus = "Medicago" 
species = "sativa"

# gene info
gene_info <- emapper %>% dplyr::select(GID = query, GENENAME = Preferred_name) %>% na.omit()

# make org
makeOrgPackage(gene_info=gene_info,
               go=gene2go,
               ko=gene2ko,
               pathway=gene2pathway,
               version="1.34.1",     # version of this package
               maintainer = "Richard1iu <ly159632874@163.com>",  # change it as your name and email
               author = "Richard1iu <ly159632874@163.com>",      # change it as your name and email
               tax_id=tax_id,    # species id or anything you like
               genus=genus,      # species genus or anything you like
               species=species,  # species name or anything you like
               goTable="go")

# 4.org use -----------------------------------------------------------------

# install the package you have made
install.packages("./org.Medicago.eg.db", repos=NULL)
library(org.Medicago.eg.db)

# glimpse the org
columns(org.Medicago.eg.db)
keys(org.Medicago.eg.db)
select(org.Medicago.eg.db, keys = "CW07G09620", columns = c("GO"))

# format it
pathway2name$Name <- gsub(" \\[BR:ko[0-9]{5}\\]", "",pathway2name$Name)
pathway2name<- na.omit(pathway2name)
pathway2gene <-gene2pathway[, c("Pathway","GID")]

# output
write.table(pathway2name,file = "./pathway2name", sep = "\t", quote = F, row.names = F)
write.table(pathway2gene,file = "./pathway2gene", sep = "\t", quote = F, row.names = F)


# 5.enrichment analysis -----------------------------------------------------

# just import the files below, no need to load org
pathway2gene <- read.table("./pathway2gene",header = T,sep = "\t")
pathway2name <- read.table("./pathway2name",header = T,sep = "\t")

# import the diff gene
gene <- read.csv("/root/total_diff_gene.csv")
gene_list <- gene[,1]

# KEGG pathway enrichment
ekp <- enricher(gene_list, 
                TERM2GENE = pathway2gene, 
                TERM2NAME = pathway2name, 
                pvalueCutoff = 1,  # 表示全部保留，可以设为0.05作为阈值
                qvalueCutoff = 1, # 表示全部保留，可以设为0.05作为阈值
                pAdjustMethod = "BH",
                minGSSize = 1)
dotplot(ekp)

# GO enrichment
library(org.Medicago.eg.db)
ego <- enrichGO(gene=gene_list,
                OrgDb=org.Medicago.eg.db,
                keyType="GID",
                ont="ALL",   #CC/BP/MF可选
                qvalueCutoff = 0.05,
                pvalueCutoff =0.05)
dotplot(ego)

# output
ego_results<-as.data.frame(ego)
write.table(ego_results, file = "./ego_results.txt", quote = F)

