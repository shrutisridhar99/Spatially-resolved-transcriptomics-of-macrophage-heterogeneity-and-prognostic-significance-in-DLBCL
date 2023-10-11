
rm(list=ls())
require(gplots)
require(readxl)
require(limma)
require(writexl)

## load the data
setwd("/.../")
info<- readRDS( "Info_data_CD68.rds")
dati<- readRDS("Harmonized_data_CD68.rds")

##
info$label<- "NA"
info$label[ grep("DLBCL",  info$Tissue_type)] <- "DLBCL"

## remove the not-pretreatment and not-CHOP.R patients
info$label<- ifelse(info$pretreatment.or.Relapse.sample..0.pretreatment.1.first.relapse.2.second.relapse.!="0", "NA", info$label)
info$label <- ifelse(info$treatment..0.CHOP.R..1.others.!="0", "NA", info$label)

## remove outlier samples
out<- c("Susan.DLBCL.TMA.B2.tonsil.17.12.right...005...CD3" ,"Susan.DLBCL.TMA.B2.tonsil.17.12.right...005...CD20","Susan.DLBCL.TMA.B2.tonsil.17.12.right...010...CD68")
info$label[which(info$lab %in% out)]<- "NA"

info$label<- ifelse(info$Tonsil_type=="GC", "GC", info$label)

take<- which(info$label=="DLBCL" | info$label=="GC" )
info<- info[ take , ]

A<- dati[ , take ]


#####
ROI_one = "DLBCL"
ROI_base = "GC"

### run the DEA
design <- model.matrix(~0+info$label)
colnames(design) <- c("DLBCL", "GC")
corfit <- duplicateCorrelation(A,design,block=info$Study_number)
corfit$consensus
fit <- lmFit(A,design,block=info$Study_number,correlation=corfit$consensus)

cm <- makeContrasts(
              comp = DLBCL-GC,
     levels=design)
  
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
ris0<-   topTable(fit2, coef="comp", number=nrow(dati), adjust="BH" )

ris<- data.frame(geni=rownames(ris0), ris0 )  

th<- 0.58
UP_genes = subset(ris, logFC>th & adj.P.Val<0.05 )
LOW_genes = subset(ris, logFC<(-th) & adj.P.Val<0.05 )

