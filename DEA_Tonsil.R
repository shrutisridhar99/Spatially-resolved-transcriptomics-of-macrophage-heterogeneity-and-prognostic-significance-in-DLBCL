rm(list=ls())
require(gplots)
require(readxl)
require(limma)
require(writexl)

##
setwd("/path...")
dati0<- data.frame(read_xlsx("Q3 normalization seq 50 filter 5% target.xlsx", sheet=3))

M<- as.matrix(dati0[, -1])
rownames(M)<- dati0$TargetName

info<- data.frame(read_xlsx("Q3 normalization seq 50 filter 5% target.xlsx", sheet=1))

info$lab <- gsub(" ", "\\.", info$SegmentDisplayName)
info$lab <- gsub("\\|", "\\.", info$lab)
info$lab <- gsub("\\+", "\\.", info$lab)
info$lab <- gsub("-", "\\.", info$lab)

identical(info$lab, colnames(M))

## inserisco le ROI label nei nomi di colonna
colnames(M) <- paste(colnames(M), info$Region_type_new, sep="-")

colnames(M) <-gsub("CD20CD3CD68\\.\\.\\.", "",  colnames(M) )
colnames(M) <-gsub("CD3CD20CD68\\.\\.\\.", "",  colnames(M) )
colnames(M) <-gsub("CD20CD68CD3\\.\\.\\.", "",  colnames(M) )

## select a mask of interest
mask <- "CD68"

dati<- M[ , grep(paste0("\\.\\.\\.", Mask), colnames(M))]


## Select the ROIs of interest
ROI_one = "DZ"
ROI_base = "LZ"

## Select the colors for the graphs
color_one = "red"
color_base = "deepskyblue"

## start of the code for differential expression analysis
comparison = dir = paste0(ROI_one, "_vs_", ROI_base)

col_one =grep(ROI_one,colnames(dati))
col_base = c(grep(ROI_base ,colnames(dati)) )  

one_DE = dati[, col_one] ;ncol(one_DE)
base_DE = dati[, col_base] ;ncol(base_DE)


A = log(cbind(base_DE, one_DE),2)
type= c(rep(paste0("aaa", ROI_base), ncol(base_DE)),rep(ROI_one, ncol(one_DE) ))

design <- model.matrix(~ type)
fit <- lmFit(A, design)
fit <- eBayes(fit)
ris1 =   topTable(fit, number=nrow(dati), adjust="BH" )

ris = data.frame(geni=rownames(ris1), ris1 )

UP_bonf = subset(ris, logFC>0 & adj.P.Val<0.05 )[, c(1,2,4,5,6)]  ;nrow(UP_bonf)
LOW_bonf = subset(ris, logFC<0 & adj.P.Val<0.05 )[, c(1,2,4,5,6)] ; nrow(LOW_bonf)

## save res
print.ris = ris[, c(1,2,4,5,6)]
write.table(print.ris, paste0("DE_Allgenes_",comparison, ".txt" ) )

print.ris = subset(print.ris,adj.P.Val<0.05)

### scatter plot
pdf(paste0("scatterplot_",ROI_one,"_vs_",ROI_base,"_bonf.pdf"))
plot(ris$logFC, -log(ris$P.Value), ylab="-log(p-value)", xlab="log(FC)")
points (UP_bonf$logFC, -log(UP_bonf$P.Value), col=color_one)
points (LOW_bonf$logFC, -log(LOW_bonf$P.Value), col=color_base)
legend("topleft",       
       legend = c( paste(nrow(UP_bonf),"UP genes in",ROI_one,"ROIs" ), paste(nrow(LOW_bonf),"UP genes in",ROI_base,"ROIs" )),
       col = c(color_one, color_base),bg="white",
       pch = 1)
dev.off()



######
########## heatmap
sign = subset(ris , adj.P.Val<0.05 ) ; dim(sign)
geni_sig = sign$geni

M=  as.matrix(A[ which(rownames(A) %in% geni_sig  ),  ])

my_palette <- colorRampPalette(c( "blue4","blue",  "white", "red",  "darkred"))(n = 200)

library(pheatmap)
my_sample_col <- data.frame(ROI = rep( c(ROI_base, ROI_one),  c(ncol(base_DE),ncol( one_DE)) ))
row.names(my_sample_col) <- colnames(M)

annotation_colors = list(    ROI = c( ROI_base=color_base, ROI_one=color_one ))
names(annotation_colors$ROI)<-c(ROI_base, ROI_one)

pdf(paste0("pheatmap_DEgene_",ROI_one,"_vs_",ROI_base,".pdf"),width=9,height=7)
pheatmap(M ,scale="row", fontsize_row = 1,  fontsize_col = 0.01, 
         cluster_rows = T,
         cluster_cols=T,
         row.names=T,
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean",
         clustering_method = "complete",
         annotation_col = my_sample_col, color=my_palette, annotation_colors=annotation_colors  )
dev.off()
