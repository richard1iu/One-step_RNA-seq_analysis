# https://www.jianshu.com/p/f2e4dbaae719

# set
options(stringsAsFactors = F)
options(warn = -1)

# load packages
library(tidyverse)
library(jsonlite)
library(AnnotationForge)

# import protein annotation file from output of eggnog or other mapper software
emapper <- read.table("F:/RNA-seq_reads/paried/alfafa_zm4_reference_genome/alfalfa_emapper_annotations_R_input.txt", header=TRUE, sep = "\t",quote = "")
emapper[emapper==""] <- NA


# 1.Split single row to multiple row and extract GO info ----------------------
##  1.1 one-step using tidyr::separate_rows
gene_go <- emapper %>% dplyr::select(query, GOs) %>% na.omit() %>%
                tidyr::separate_rows(GOs, sep = ",") %>%
                dplyr::mutate(EVIDENCE = "IEA") %>%
                dplyr::rename(GO = GOs) %>% 
                dplyr::filter(GO != "-")

# ##  1.2 function using str_split and apply
# ### functions
# gos_list <- function(x){
#   the_gos <- str_split(x[2], pattern = ",", simplify = FALSE)[[1]]
#   df_temp <- data.frame(GID = rep(x[1], length(the_gos)),
#                         GO = the_gos,
#                         EVIDENCE = rep("IEA", length(the_gos))
#   )
#   return(df_temp)
# }
# ### extract GO
# gene2gol <- emapper %>% dplyr::select(query, GOs) %>% na.omit() %>% apply(.,1,gos_list)
# gene2go <- do.call(rbind.data.frame, gene2gol) %>% filter(GO != "-")

# ## 1.3 using str_split and sapply
# the_gos <- str_split(emapper$GOs, ",", simplify = FALSE)
# gene2go <- data.frame(GID = rep(emapper$query_name,
#                                 times = sapply(the_gos, length)
#                                 ),
#                       GO = unlist(the_gos),
#                       EVIDENCE = "IEA")

# 2.Extract KEGG info -----------------------------------------------
# extract KEGG
gene_ko <- emapper %>% dplyr::select(GID = query, Ko = KEGG_ko) %>%
                dplyr::filter(Ko != "-") %>%
                tidyr::separate_rows(Ko, sep = ",") %>%
                dplyr::mutate(Ko = str_replace_all(Ko, "ko:", ""))

# 3.Extract kegg functions -------------------------------------------

## import KEGG json file (https://www.genome.jp/kegg-bin/get_htext?ko00001)
kegg_json <- fromJSON("F:/Database/Kegg/ko00001.json")

## kegg function
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

## update kegg database
update_kegg(kegg_json)

## load kegg_info.RData
load(file = "F:/Database/Kegg/kegg_info.RData")

## combine ko and pathway
gene_pathway <- gene_ko %>%
  left_join(ko2pathway, by = "Ko") %>%
  dplyr::select(GID, Pathway) %>%
  na.omit()

## gene info
gene_info <- emapper %>% dplyr::select(GID = query, GENENAME = Preferred_name) %>% na.omit()

## remove duplicates
gene_go <- unique(gene_go)
gene_go <- gene_go[!duplicated(gene_go), ]
gene_ko <- gene_ko[!duplicated(gene_ko), ]
gene_pathway <- gene_pathway[!duplicated(gene_pathway), ]
gene_info <- gene_info[!duplicated(gene_info), ]

# 4.Make org ----------------------------------------------------------------

# org info: search the name in ncbi(https://www.ncbi.nlm.nih.gov/taxonomy)
# or you can define it as you wish
taxid = 3879
genus = "Medicago"
species = "sativa"

# make org
makeOrgPackage(
  gene_info = gene_info,
  go = gene_go,
  ko = gene_ko,
  pathway = gene2pathway,
  version = "1.34.1", # version of this package
  maintainer = "Richard1iu <ly159632874@163.com>", # change it as your name and email
  author = "Richard1iu <ly159632874@163.com>", # change it as your name and email
  tax_id = taxid, # species id or anything you like
  genus = genus, # species genus or anything you like
  species = species, # species name or anything you like
  goTable = "go"
)

# 5.org use -----------------------------------------------------------------

## install the package you have made
install.packages("F:/RNA-seq_reads/paried/alfafa_zm4_reference_genome/org.Msativa.eg.db/",repos = NULL, type = "source")
library(org.Msativa.eg.db)

## glimpse the org
columns(org.Msativa.eg.db)
keys(org.Msativa.eg.db,keytype = "GID") %>% head()
select(org.Msativa.eg.db, keys = "Msa0022050", columns = c("GO"))

## format it
pathway2name$Name <- gsub(" \\[BR:ko[0-9]{5}\\]", "", pathway2name$Name)
pathway2name <- na.omit(pathway2name)
pathway2gene <- gene2pathway[, c("Pathway", "GID")]

## output
write.table(pathway2name,file = "F:/RNA-seq_reads/paried/alfafa_zm4_reference_genome/alfalfa_pathway2name", sep = "\t", quote = F, row.names = F)
write.table(pathway2gene,file = "F:/RNA-seq_reads/paried/alfafa_zm4_reference_genome/alfalfa_pathway2gene", sep = "\t", quote = F, row.names = F)


# 6.enrichment analysis -----------------------------------------------------

## import only the files below, no need to load org
pathway2gene <- read.table("F:/RNA-seq_reads/paried/alfafa_zm4_reference_genome/alfalfa_pathway2gene", header = T, sep = "\t")
pathway2name <- read.table("F:/RNA-seq_reads/paried/alfafa_zm4_reference_genome/alfalfa_pathway2name", header = T, sep = "\t")

## import the diff gene
gene <- read.csv("/root/total_diff_gene.csv")
gene_list <- gene[, 1]

## KEGG pathway enrichment
ekp <- enricher(gene_list,
                TERM2GENE = pathway2gene,
                TERM2NAME = pathway2name,
                pvalueCutoff = 1,
                qvalueCutoff = 1,
                pAdjustMethod = "BH",
                minGSSize = 1)
dotplot(ekp)

## GO enrichment
library(org.Msativa.eg.db)
ego <- enrichGO(gene = gene_list,
                OrgDb = org.Msativa.eg.db,
                keyType = "GID",
                ont = "ALL",   # CC/BP/MF
                qvalueCutoff = 0.05,
                pvalueCutoff = 0.05)
dotplot(ego)

## output
ego_results <- as.data.frame(ego)
write.table(ego_results, file = "./ego_results.txt", quote = F)