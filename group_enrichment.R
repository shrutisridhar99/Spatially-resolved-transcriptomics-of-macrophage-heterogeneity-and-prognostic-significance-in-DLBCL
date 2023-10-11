
rm(list=ls())

require(writexl)

## dataset that contains patients classified based on percentiles (e.g., DZ-like and LZ-like patients)
setwd(paste0("/path..."))
clust<-read.table("cluster_data.txt", h=T)

     ### clinical information - 
    # geo_accession contains column patient ID, 
    # class has to be a categorical variable for the enrichment analysis (e.g., COO, genetic subtype)
    info = read.table("info_example.txt", h=T)
    

    ## merge the information
join <- merge(clust, info, by.x = "sample",by.y="geo_accession" )

cluster= names(table(join$cluster))
cluster= cluster[-which(cluster=="intermediate")]

class<- unique(join$class)

comb<- expand.grid(cluster, class)

comb<- comb[ order(as.vector(comb$Var1), as.vector(comb$Var2)), ]

names(comb)<- c("cluster", "class")

## for cycle to calculate Fisher p-values
  m<- kn <- q<-pval <- numeric()

for(i in 1:nrow(comb)){

N<- nrow(join)
m[i] = nrow(subset(join,class==comb$class[i]) )
n = N - m[i]

q[i] = nrow(subset(join,cluster==comb$cluster[i] & class==comb$class[i] ) )

kn[i] = nrow(subset(join, cluster ==comb$cluster[i] ) )

pval[i] = phyper(q[i] -1, m[i], n, kn[i], lower.tail = FALSE, log.p = FALSE)

} # fine ciclo in i

res <- data.frame(comparison= "signature_name" , 
                  dataset="dataset_name",
                  comb , n.group= kn,  n.class=m ,   overlap = q,  pval = signif(pval, 2) )
