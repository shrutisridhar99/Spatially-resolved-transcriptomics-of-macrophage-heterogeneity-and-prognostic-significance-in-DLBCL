#### survival analysis between groups obtained from quantile separation (tertile separation by default)

funz_sopravvivenza = function(dati , geniF, info, dir, th=0.3333, dati_cor , subset="all", clustMethod="ward.D2", PheatScale="row",
                              g_up="high-rank", g_down="low-rank",  sumCox=FALSE, 
                              graph_title="", surv_name="OS" ,Pheatmap=TRUE , 
                              colore_down="deepskyblue" , colore_intermediate="pink", colore_up="forestgreen"){
  
  g_up_gene<- gsub("like", "UP",  g_up)
  g_down_gene<- gsub("like", "UP",g_down)
  
  setwd(paste(dir))
  
  if(subset=="ABC"){
    info =     subset(info,class=="ABC") }
  
  if(subset=="GCB"){
    info =     subset(info,class=="GCB") }
  
  if(subset=="unknown"){
    info =     subset(info,class!="GCB" & class!="ABC" & class!="TypeIII" ) }
  
  datiF1 = dati[ rownames(dati) %in% geniF, ]
  datiF = datiF1[ , which(colnames(datiF1) %in% info$geo_accession) ]
  
  ## ordering
  datiF = datiF[ order(rownames(datiF)), ]
  
  print(paste0("signature = " ,length(geniF) , " - ", "intersect = ", length(intersect(rownames(dati) , geniF)) ))
  print(paste0("sample = " ,ncol(dati) , " - ", "info.patient = ", ncol(datiF)))
  
  ###  create a dataset for weighed sum
  names(dati_cor)=c("geni","corr.vec")
  dati_cor =  dati_cor[ which( dati_cor$geni %in% rownames(datiF)) , ]
  dati_cor =  dati_cor[ order( dati_cor$geni), ]
  
  print( identical(rownames(datiF),  dati_cor$geni))
  
  funz_pesi  =function(x){ 
    y= sum(x* dati_cor$corr.vec )
    return(y) }
  
  ss = apply(datiF, 2, funz_pesi)
  
  q1 = quantile(ss, th)
  q2 = quantile(ss, 1-th)
  
  cluster= ifelse(ss<=q1, g_down, 
                  ifelse(ss>=q2, g_up, "intermediate"))
  
  color_clust_down = colore_down
  color_clust_up = colore_up
  
  makeTransparent<-function(someColor, alpha=100)
  {
    newColor<-col2rgb(someColor)
    apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                                blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
  }
  
  pdf("density_plot.pdf")
  plot(density(ss), main=" ", xlab = "weighed total expression per patient")
  abline(h=0)
  a = density(ss)
  polygon(c(a$x[a$x>=q2], max(a$x), q2), c(a$y[a$x>=q2], 0, 0), col=makeTransparent( color_clust_up))
  polygon(c(min(a$x), a$x[a$x<=q1], q1), c(0, a$y[a$x<=q1], 0), col=makeTransparent( color_clust_down))
  polygon( c(q1, q1, a$x[a$x>=q1 & a$x<=q2], q2), c(0, a$y[a$x>q1][1], a$y[a$x>=q1 & a$x<=q2], a$y[a$x<q2][1]), col=makeTransparent(colore_intermediate))
  dev.off()  
  
  require(gridExtra)
  freq = data.frame(table(cluster))
  pdf("freq_groups.pdf", width = 3, height = 2)
  grid.table(freq)
  dev.off()
  
##
  dati_clust = data.frame(sample = colnames(datiF),  cluster, ss   )
  datiF = datiF[ , order(dati_clust$ss)  ]
  
  # save groups 
    write.table(dati_clust, "cluster_data.txt", row.names = F)

  
  ### pheatmap
  if( subset=="all"){
    ## my sample col with cell of origin
    my_sample_col0 <- data.frame( cluster   )
    my_sample_col0$geo_accession = rownames(my_sample_col0)
    my_sample_col1 = merge(my_sample_col0, info, by.x = "geo_accession",  by.y = "geo_accession")
    rownames(my_sample_col1) = my_sample_col1$geo_accession  
    names(my_sample_col1)[5]="cell of origin"
    names(my_sample_col1)[2]="group"
    my_sample_col1$`cell of origin`=  ifelse( my_sample_col1$`cell of origin`=="UNC", "Unclassified", my_sample_col1$`cell of origin`)
    my_sample_col1$`cell of origin`=  ifelse( my_sample_col1$`cell of origin`=="TypeIII", "Unclassified", my_sample_col1$`cell of origin`)
    my_sample_col = my_sample_col1[ , c(2,5)]
    
    annotation_colors = list(    group = c(g_down=color_clust_down,g_up=color_clust_up, "intermediate"=colore_intermediate ),
                                 `cell of origin`=c("ABC"="brown4", "GCB"="cyan4", "Unclassified"="darkorange" )  ,  #"TypeIII"="darkorange",
                                 gene=c("Negative"="violet", "Positive"="orange"))
  
    
    if( length(grep("MHG",  my_sample_col1$`cell of origin`))>0){
      annotation_colors = list(    group = c(g_down=color_clust_down,g_up=color_clust_up, "intermediate"=colore_intermediate ),
                                   `cell of origin`=c("ABC"="brown4", "GCB"="cyan4", "Unclassified"="darkorange", "MHG"="violet" )  ,  #"TypeIII"="darkorange",
                                   gene=c("Negative"="violet", "Positive"="orange"))    }
    
  }else{
    group= cluster
    my_sample_col <- data.frame( group   )
    annotation_colors = list(    group = c("down"=color_clust_down,"up"=color_clust_up, "intermediate"=colore_intermediate),
                                 gene=c( "Negative"="violet", "Positive"="orange"))
  }
  
  names(annotation_colors$group)[1:2]<-c(g_down, g_up)
  names(annotation_colors$gene)[1:2]<-c(g_down_gene, g_up_gene)
  
  my_sample_row = data.frame(gene= dati_cor$corr.vec)
  rownames(  my_sample_row) = dati_cor$geni
  my_sample_row$gene = ifelse(  my_sample_row$gene<0, g_down_gene, g_up_gene )
  
  
  my_palette <- colorRampPalette(c("royalblue3","blue" ,"ivory","red","red2"))(n = 100)
 
  if(nrow(datiF)<= 10){ fontsize = 10} 
  if(nrow(datiF)<= 20 & nrow(datiF)>10){ fontsize = 6}
  if(nrow(datiF)<= 50 & nrow(datiF)>20 ){ fontsize = 4}
  if(nrow(datiF)<= 100 & nrow(datiF)>50 ){ fontsize = 3}
  if(nrow(datiF)<= 150 & nrow(datiF)>100 ){ fontsize = 2}
  if( nrow(datiF)>100 ){ fontsize = 1}
  
  
ddd<-  apply(datiF, 1, sum)  #remove null expression
if(length(  which(ddd==0))>0){
  datiF<- datiF[-which(ddd==0), ]
  ddd<-  apply(datiF, 1, sum) }

if(max(table(ddd))>1){  #remove duplicates
datiF<-   datiF[!duplicated(datiF) , ] }

sss<-  apply(datiF, 1, sd) # remove rows without variability
if(length(  which(sss==0))>0){
  datiF<- datiF[-which(sss==0), ] }

  if(Pheatmap==TRUE){
    require(pheatmap)
    pdf("heatmap_rowScale_new.pdf",width=9,height=7)
    pheatmap(datiF,  scale= PheatScale, 
             fontsize_col = 0.001,   fontsize_row = fontsize,
             cluster_rows = T,
             cluster_cols=F,
             row.names=T,
             col.names=rep("", ncol(datiF)),
             labels_col=FALSE,
             clustering_distance_rows="euclidean",
             clustering_distance_cols="euclidean",
             clustering_method =   clustMethod,
             annotation_col = my_sample_col , annotation_row = my_sample_row ,
             annotation_colors = annotation_colors,
             color=my_palette, main = graph_title )
    dev.off()
  }
  

############################################################
  ##### survival analysis
  gruppi = data.frame(camp=names(cluster), group=cluster )
  gruppi = subset(gruppi, group != "intermediate")
  
  join = merge(info, gruppi, by.x = "geo_accession", by.y = "camp")
  

  library(survival)
  library(survminer)
  library(dplyr)
  
  surv_object <- Surv(time = join$follow_up, event = join$fu_state)
  fit1 <- surv_fit(surv_object ~ group , data = join)
  
  #### get survival values
  summ = surv_summary(fit1, join)
  
  summ$time_rounded = trunc(summ$time)
  summ$upper= round(  summ$upper, 2)
  summ$lower= round(  summ$lower, 2)
  
  sum1 = subset(summ,group==g_down)
  sum3 = subset(summ,group==g_up)
  
  n_risk1 =   round( c(sum1$n.risk[ which(sum1$time_rounded>=0)[1] ],
                       sum1$n.risk[ which(sum1$time_rounded>=12)[1] ],
                       sum1$n.risk[ which(sum1$time_rounded>=24)[1] ],
                       sum1$n.risk[ which(sum1$time_rounded>=36)[1] ],
                       sum1$n.risk[ which(sum1$time_rounded>=48)[1] ],
                       sum1$n.risk[ which(sum1$time_rounded>=60)[1] ]) , 2)
  
  cluster1 = round( c(1, sum1$surv[ which(sum1$time_rounded>=12)[1] ],
                      sum1$surv[ which(sum1$time_rounded>=24)[1] ],
                      sum1$surv[ which(sum1$time_rounded>=36)[1] ],
                      sum1$surv[ which(sum1$time_rounded>=48)[1] ],
                      sum1$surv[ which(sum1$time_rounded>=60)[1] ]) , 2)
  
  IC1 =  c("-", paste(sum1$lower[ which(sum1$time_rounded>=12)[1] ],sum1$upper[ which(sum1$time_rounded>=12)[1] ], sep="-"), 
           paste(sum1$lower[ which(sum1$time_rounded>=24)[1] ],sum1$upper[ which(sum1$time_rounded>=24)[1] ], sep="-"),
           paste(sum1$lower[ which(sum1$time_rounded>=36)[1] ],sum1$upper[ which(sum1$time_rounded>=36)[1] ], sep="-"),
           paste(sum1$lower[ which(sum1$time_rounded>=48)[1] ],sum1$upper[ which(sum1$time_rounded>=48)[1] ], sep="-"),
           paste(sum1$lower[ which(sum1$time_rounded>=60)[1] ],sum1$upper[ which(sum1$time_rounded>=60)[1] ], sep="-"))
  
 ##
  n_risk3 =round( c(sum3$n.risk[ which(sum3$time_rounded>=0)[1] ],
                    sum3$n.risk[ which(sum3$time_rounded>=12)[1] ],
                    sum3$n.risk[ which(sum3$time_rounded>=24)[1] ],
                    sum3$n.risk[ which(sum3$time_rounded>=36)[1] ],
                    sum3$n.risk[ which(sum3$time_rounded>=48)[1] ],
                    sum3$n.risk[ which(sum3$time_rounded>=60)[1] ]) , 2)
  
  cluster3 =round( c( 1, sum3$surv[ which(sum3$time_rounded>=12)[1] ],
                      sum3$surv[ which(sum3$time_rounded>=24)[1] ],
                      sum3$surv[ which(sum3$time_rounded>=36)[1] ],
                      sum3$surv[ which(sum3$time_rounded>=48)[1] ],
                      sum3$surv[ which(sum3$time_rounded>=60)[1] ]) , 2)
  
  IC3 =  c("-", paste(sum3$lower[ which(sum3$time_rounded>=12)[1] ],sum3$upper[ which(sum3$time_rounded>=12)[1] ], sep="-"), 
           paste(sum3$lower[ which(sum3$time_rounded>=24)[1] ],sum3$upper[ which(sum3$time_rounded>=24)[1] ], sep="-"),
           paste(sum3$lower[ which(sum3$time_rounded>=36)[1] ],sum3$upper[ which(sum3$time_rounded>=36)[1] ], sep="-"),
           paste(sum3$lower[ which(sum3$time_rounded>=48)[1] ],sum3$upper[ which(sum3$time_rounded>=48)[1] ], sep="-"),
           paste(sum3$lower[ which(sum3$time_rounded>=60)[1] ],sum3$upper[ which(sum3$time_rounded>=60)[1] ], sep="-"))
  
  follow_up = c("risk pop.", "survival", "CI","risk pop.", "survival", "CI")
  gr = c("", g_down, "", "", g_up, "")
  S = rbind(n_risk1, cluster1,IC1, n_risk3,  cluster3 , IC3)  
  colnames(S) = c("t0","t12","t24","t36","t48","t60")     
  
  S =   data.frame(gr, follow_up, S)
  
  require(gridExtra)
  pdf("survival_values.pdf", width = 9, height = 4)
  grid.table(S)
  dev.off()
  
  write.table(S, "survival_values.txt", row.names = F, quote=F, sep=" & ")
  
  ##############################
  ######
  
  suv_lab<- data.frame(g=c(g_up,  g_down), color= c( color_clust_up,color_clust_down))
  suv_lab<-  suv_lab[order( suv_lab$g), ]
  
  pdf("survival_plot.pdf", width = 8, onefile=F)
  h <-ggsurvplot(fit1, data = join, pval = T, pval.method =T , palette = suv_lab$color,  
                 legend = c(0.79,0.1) , 
                 title=paste0(graph_title), 
                 legend.labs=suv_lab$g ) #
  print(h)
  dev.off()
  
  
  ### cox model
  join$group[which(join$group==g_down)]= paste0("aaa",g_down )
  
  cox <- coxph(surv_object  ~ group , data = join)
  
  ris = summary(cox)
  ris2 = data.frame(group = rownames(ris$coefficients)  , cbind(round(ris$coefficients[, 1],2),
                                                                round(ris$conf.int,2), ris$coefficients[, 5]))
  names(ris2) = c("group", "coef", "RR", "1/RR", "IC.RR.low", "IC.RR.up", "pval")
  
  ### proportional hazard assumption test
  test_chi =cox.zph(cox)
  
  ris2$zph<- signif( test_chi$table[1, 3], 2)

  write.csv(ris2, "Cox_ris.csv")
  

}  ## end function
