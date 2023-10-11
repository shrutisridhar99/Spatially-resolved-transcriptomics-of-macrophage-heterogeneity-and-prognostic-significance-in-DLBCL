
require(readxl)

funz_enrich_gsea <- function(input_genes, correction="BH",
                            dir,    pool=NULL ){
  
  
  N<- length(pool)
  m<- length(input_genes)  
  n<- N-m      

    setwd(dir)
    hallmark<-readLines("Hallmark_2022.txt")
    
    k=q= pval= Genes = Term=   numeric()
    
    for(i in 1:length(hallmark) ){  
      
      ## pathway gsea - intersecato nel pool 
      gsea = as.vector(unlist(strsplit( hallmark[i], split="\t")))
      
      geni = gsea[-c(1:2)]
      geni<- intersect(geni, pool)
      
      Term[i]<- gsea[1]
      
      k[i]<- length( geni)
      
      if(length(geni)>0 & k[i]>1){
        
        intersezione<- intersect(geni, input_genes)
        q[i]<- length(intersezione)
        
        pval[i] = phyper(q[i] -1, m, n, k[i], lower.tail = FALSE, log.p = FALSE)  
        
        Genes[i]<- paste(intersezione, collapse = ";")
        
      }else{
        q[i] =pval[i]=Genes[i]= NA
      }
      
    } #end for in i
    
    
    ris<- data.frame(library="GSEA_hallmark",Term, adj.pval=pval, found=q, size=k , Genes )  
    
    if(length(which(is.na(ris$adj.pval)))>0){
      ris<- ris[-which(is.na(ris$adj.pval)), ]
    }
    
    ris<- ris[ order(ris$adj.pval) , ]  
    ris$adj.pval<- p.adjust( ris$adj.pval, correction)
    
    ris_bonf<- subset(ris, adj.pval<0.05)
    

  print(paste0("number of sig. terms: ", nrow(ris_bonf)))
  
  return(ris_bonf)
  
} # end function
