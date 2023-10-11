#https://davislaboratory.github.io/GeoMXAnalysisWorkflow/articles/GeoMXAnalysisWorkflow.html

rm(list=ls())
require(readxl)
require(limma)
require(edgeR)

setwd("/...")
dati0<- data.frame(read_xlsx("Q3 NORMALIZATION newest 2 final version for giorgio.xlsx", sheet = 3))

rownames(dati0)<- dati0$TargetName
dati<- dati0[, -1]

###info
info <- data.frame(read_xlsx("Q3 NORMALIZATION newest 2 final version for giorgio.xlsx", sheet = 1))

data_dlbcl_tonsil<- c("DLBCL A3", "DLBCL B2", "DLBCL C1", "Tonsil CD20CD68NGFR")
data_tonsil<- c("Tonsil 17-4B+TCP CD20CD3CD68", "Tonsil 21-2B+TCP CD20CD68CD3", "Tonsil1+2 CD3CD20CD68")

info$dataset<- info$ScanLabel

info$dataset[which(info$ScanLabel %in% data_tonsil ) ]<- "Tonsil data 1"
info$dataset[which(info$ScanLabel %in% data_dlbcl_tonsil  ) ]<- "NUH DLBCL tonsil data"

info$lab<- gsub(" ", ".", info$SegmentDisplayName)
info$lab<- gsub("\\|", ".", info$lab)
info$lab<- gsub("\\+", ".", info$lab)
info$lab<- gsub("-", ".", info$lab)
info$lab<- gsub("4tonsil", "X4tonsil", info$lab)

info<- info[ match(names(dati),info$lab ) , ]

identical(info$lab, names(dati))

###
Mask="CD68"
take<- which(info$SegmentLabel==Mask)

info<- info[take, ]
dati<- dati[ , take]

label = info$dataset
label_sample = data.frame(sample= names(dati) ,label)

dati_input<- log(as.matrix(dati), 2) 

tipo<- info$Tonsil_type
tipo[which(info$Interfollicular_type=="Distant")]<- "IF.distant"

tipo[which(info$DFS_status=="1")]=  "DLBCL_relapse"
tipo[which(info$DFS_status=="0")]= "DLBCL_no.relapse"

tipo[which(tipo=="NA" & info$Tissue_type=="Susan DLBCL")]<-  "NA.DLBCL"
tipo[which(tipo=="NA" & info$Tissue_type=="NUH DLBCL")]<-  "NA.DLBCL"

## remove batch effect
design <- model.matrix(~ tipo)
dati_noBatch = removeBatchEffect(dati_input, batch=label, design = design)
