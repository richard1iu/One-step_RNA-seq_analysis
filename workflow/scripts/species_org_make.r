# https://blog.csdn.net/u012110870/article/details/102804467
#
library(tidyverse)
library(clusterProfiler)
library(GO.db)

##
species_annotation <- read.delim("F:/RNA-seq_reads/paried//alfafa_zm4_reference_genome/ATH_GO_GOSLIM.txt",header=FALSE)
species_annotation <- species_annotation[,c(1,6,8,10)]
colnames(species_annotation) <- c("geneID","GOTerm","Ont","Source")

# BP
species_GO_BP <- species_annotation %>% 
  dplyr::filter(Ont == "P") %>% 
  select(GOTerm,geneID) %>% head(500) %>% 
  buildGOmap()

# MF
species_GO_MF <- species_annotation %>% 
  dplyr::filter(Ont == "F") %>% 
  select(GOTerm,geneID) %>% 
  buildGOmap()

# CC
species_GO_BP <- species_annotation %>% 
  dplyr::filter(Ont == "C") %>% 
  select(GOTerm,geneID) %>% 
  buildGOmap()

#
goname_BP <- AnnotationDbi::select(x=GO.db, keys = species_GO_BP$GO,  keytype = "GOID",columns = "TERM" )

# enrich analysis
diff_gene_id <- species_GO_BP$Gene %>% head(50)
ego <- enricher(diff_gene_id,TERM2GENE = species_GO_BP, TERM2NAME=goname_BP)

#
ego@ontology <- "BP" # 按照需求改成，CC或MF
simplify(ego)

#
list(ego= ego@result, ego_new=ego@result) %>% merge_result() %>% dotplot()

