
rm(list=ls())
require(gplots)
require(readxl)
require(limma) 
require(writexl)

## load the data
setwd("/.../")
info<- readRDS( "Info_data_CD68.rds")
dati<- readRDS("Harmonized_data_CD68.rds")

identical(info$lab, colnames(dati))

head(info)

## exclusion of IF distant AOIs
info$Tonsil_type<- ifelse(info$Interfollicular_type=="Distant", "NA",info$Tonsil_type )

##
colnames(dati) <- paste(colnames(dati), info$Tonsil_type, sep="-")
colnames(dati) <-gsub("CD20CD3CD68\\.\\.\\.", "",  colnames(dati) )
colnames(dati) <-gsub("CD3CD20CD68\\.\\.\\.", "",  colnames(dati) )
colnames(dati) <-gsub("CD20CD68CD3\\.\\.\\.", "",  colnames(dati) )
colnames(dati) <-gsub("CD20CD68CD3\\.\\.\\.", "",  colnames(dati) )
colnames(dati) <-gsub("CD20CD68", "",  colnames(dati) )

########## 
ROI_one = "DZ"
ROI_base = "LZ"

col_one =grep(ROI_one,colnames(dati))
col_base = c(grep(ROI_base ,colnames(dati)) )   

one_DE = dati[, col_one] ;ncol(one_DE)
base_DE = dati[, col_base] ;ncol(base_DE)

A = cbind(base_DE, one_DE); max(A) 
type= c(rep(paste0("aaa", ROI_base), ncol(base_DE)),rep(ROI_one, ncol(one_DE) ))

design <- model.matrix(~ type)
fit <- lmFit(A, design)
fit <- eBayes(fit)
ris1 =   topTable(fit, number=nrow(dati), adjust="BH" )

ris = data.frame(geni=rownames(ris1), ris1 )

th<-0.58
UP_genes = subset(ris, logFC>th & adj.P.Val<0.05 )  ; nrow(UP_genes)
LOW_genes = subset(ris, logFC<(-th) & adj.P.Val<0.05 )  ; nrow(LOW_genes)
