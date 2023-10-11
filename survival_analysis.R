require(cluster)
require(factoextra)
require(readxl)

rm(list=ls())

    setwd("/path_to_gene_signature...")
    firma_up = read_xlsx("DEGs_DZ_vs_LZ.xlsx", sheet = 1)
    firma_down = read_xlsx("DEGs_DZ_vs_LZ.xlsx", sheet = 2)
    
    ## specify the signature name (for graphs)
    signature_name = "DZ_vs_LZ"
    g_up<- "DZ-like"
    g_down<- "LZ-like"


##  weights for the ordering
firma_cor<- rbind.data.frame(firma_up, firma_down)
firma_cor$logpval<- ifelse(firma_cor$logFC>0 , -log(firma_cor$P.Value, 10), log(firma_cor$P.Value, 10))

geniF = c(firma_cor$geni) ; length(geniF)

dati_cor = firma_cor[, c("geni",   "logpval")] 


## load gene expression data - matrix format
## rownames are genes,  colnames are patient IDs
setwd("/...path")
dati = readRDS( "/path.../dataset_name.rds")


### clinical information - 
# first column patient ID, second column fullow-up state (1 dead, 0 censor), 
# third column follow-up time, fourth column COO class (ABC, GCB, Unclassified)
info = read_xlsx("...dataset_name.xlsx")

names(info) = c("geo_accession", "fu_state", "follow_up", "class")


## function
source("/path.../function_survival_analysis.R")

##### specify dataset name 
dataset= "dataset_name"

#specify  percentile for group separation
th= 0.3333

## printing dicerctory
dir = paste0("/path...")

funz_sopravvivenza(dati=dati , geniF= geniF, info=info,dati_cor = dati_cor, th=th, dir=dir, g_up=g_up, g_down=g_down, 
                   graph_title=paste0(gsub("_", " ", signature_name)," - ", dataset,  " et al."))


