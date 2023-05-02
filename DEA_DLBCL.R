rm(list=ls())
require(gplots)
require(readxl)
require(limma)
require(writexl)

##
setwd("/path...")
M<- readRDS("Q3_filter_5_gene_DLBCL_and_tonsil_TMA_noReplicates_new.rds")

## select a mask of interest
mask <- "CD68"

dati<- M[, grep(mask, colnames(M))]


## Select the ROIs of interest
ROI_one = "DZ"
ROI_base = "LZ"

## Select the colors for the graphs
color_one = "red"
color_base = "deepskyblue"


## start of the code for differential expression analysis
comparison = dir = paste0(ROI_one, "_vs_", ROI_base)

col_one =grep(ROI_one,colnames(dati))
col_base = grep(ROI_base ,colnames(dati))   

one_DE = dati[, col_one] ; dim(one_DE)
base_DE = dati[, col_base] ; dim(base_DE)

A = log(cbind(base_DE, one_DE),2)
type= c(rep(paste0("aaabase", ROI_base), ncol(base_DE)),rep(ROI_one, ncol(one_DE) ))

design <- model.matrix(~ type)
fit <- lmFit(A, design)
fit <- eBayes(fit)
ris1 =   topTable(fit, number=nrow(dati), adjust="BH" )

ris = data.frame(geni=rownames(ris1), ris1 )

## select a logFC threshold 
th<- 0

UP_bonf = subset(ris, logFC>th & adj.P.Val<0.05 )[, c(1,2,4,5,6)] 
LOW_bonf = subset(ris, logFC<(-th) & adj.P.Val<0.05 )[, c(1,2,4,5,6)] 


## save res
print.ris = ris[, c(1,2,4,5,6)]
write.table(print.ris, paste0("DE_Allgenes_",comparison, ".txt" ) )

##
if(ROI_base=="GC-"){ROI_base="GC"}

### scatter plot
pdf(paste0("scatterplot_",ROI_one,"_vs_",ROI_base,"_th-",th,".pdf"))
plot(ris$logFC, -log(ris$P.Value), ylab="-log(p-value)", xlab="log(FC)")
points (UP_bonf$logFC, -log(UP_bonf$P.Value), col=color_one)
points (LOW_bonf$logFC, -log(LOW_bonf$P.Value), col=color_base)
legend("topright",  bg="white",     
       legend = c( paste(nrow(UP_bonf),"UP genes in",ROI_one,"ROIs" ), paste(nrow(LOW_bonf),"UP genes in",ROI_base,"ROIs" )),
       col = c(color_one, color_base),
       pch = 1)
dev.off()


######
########## heatmap
sign = ris[ which(ris$geni %in% c(UP_bonf$geni, LOW_bonf$geni)) ,  ] ; dim(sign)
geni_sig = sign$geni

M3=  as.matrix(A[ which(rownames(A) %in% geni_sig  ),  ])

library(pheatmap)
my_sample_col <- data.frame(ROI =colnames(A) )
row.names(my_sample_col) <- colnames(A)

if(ROI_one=="DLBCL"){
  my_sample_col$ROI[grep("DLBCL", my_sample_col$ROI) ]= "DLBCL"
}else{
my_sample_col$ROI[grep("DLBCL-relapse", my_sample_col$ROI) ]= "DLBCL-relapse"
my_sample_col$ROI[grep("DLBCL-not.relapse", my_sample_col$ROI) ]= "DLBCL-not.relapse"
}

my_sample_col$ROI[grep("GC", my_sample_col$ROI) ]= "GC"

my_sample_col$ROI[grep("DZ", my_sample_col$ROI) ]= "DZ"
my_sample_col$ROI[grep("LZ", my_sample_col$ROI) ]= "LZ"


my_palette <- colorRampPalette(c( "blue4","blue",  "white", "red",  "darkred"))(n = 200)

annotation_colors = list( ROI=c("x"=color_base, "y"=color_one))
names(annotation_colors$ROI)<- c(ROI_base, ROI_one)

setwd(paste0("/Volumes/Bertolaz/DSP_data_Min/DEA_DLBCL_noReplicates/",mask, "/", dir ))
pdf(paste0("pheatmap_DEgene_",dir,".pdf"),width=9,height=7)
pheatmap(M3 ,scale="row", fontsize_row = 1, 
         cluster_rows = T,
         cluster_cols=T,
         row.names=T,
         col.names=F, 
         fontsize_col = 0.001,
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean",
         clustering_method = "complete",
         annotation_col = my_sample_col,  annotation_colors = annotation_colors,
         color=my_palette, 
         treeheight_row=0) 
dev.off()
