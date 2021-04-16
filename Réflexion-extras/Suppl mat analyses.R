library(purrr)  # for the "compact" function (removing empty objects from a list)
library(stringr)

Comparison=function(condition="Hom_bal",param="bias",fct="meandiff",chemin){

  Filestxt=list.files("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/",pattern = ".*good_mes.txt")

  for (i in seq_len(length(Filestxt))){
    
    File=read.table(paste0("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/",Filestxt[i]),sep=";",dec=",")  
    # Supprimer tous les "Inf" du fichier
    
    # Selected row
    if (condition=="Hom_bal"){
      id <- as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1==File$n2)&(File$sd1.sd2 == 1),]))
    } else if (condition=="Hom_rnull"){
      id_Hom_nN=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 == 1),]))
      id_Hom_Nn=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 == 1),]))
      id <- c(id_Hom_nN,id_Hom_Nn)
    } else if (condition=="Het_bal") {
      id_sdSD_bal=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1==File$n2)&(File$sd1.sd2 < 1),]))
      id_SDsd_bal=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1==File$n2)&(File$sd1.sd2 > 1),]))
      id <-c(id_sdSD_bal,id_SDsd_bal)
    } else if (condition=="Het_rpos"){
      id_sdSD_nN=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 < 1),]))
      id_SDsd_Nn=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 > 1),]))
      id <-c(id_sdSD_nN,id_SDsd_Nn)
    } else if(condition=="Het_rneg") {
      id_SDsd_nN=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 > 1),]))
      id_sdSD_Nn=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 < 1),]))
      id <-c(id_SDsd_nN,id_sdSD_Nn)
    }

    # Selected columns 
    if (param=="bias"){
      col = c(1:5,6:11)
      main1 = "bias, as a function of "
      ylabel= "bias"
    } else if (param == "relbias"){
      col = c(1:5,12:17)
      main1 = "relative bias, as a function of "
      ylabel= "relative bias"
    } else if (param == "eff") {
      col = c(1:5,18:23)
      main1 = "Variance, as a function of "
      ylabel= "Variance"
    } else if (param == "releff"){
      col = c(1:5,24:29)
      main1 = "relative variance, as a function of "
      ylabel= "relative variance"
    }
    
    Sel <- File[id,col]

    # Do we want to check evoluations as a function of mean differences, as a function of sd-ratio, or both?        
    if (fct=="meandiff"){
      discr <- list(n1=Sel$n1,n2=Sel$n2,sdratio=Sel$sd1.sd2) # all parameters but the mean
      main2 <- "the mean difference"
      xlabel <- "Mean difference"
      xlimit <- c(1,4)
    } else if (fct=="sdratio"){
      discr <- list(n1=Sel$n1,n2=Sel$n2,meandiff=Sel$m1.m2) # all parameters but the sd-ratio
      main2 <- "the sd ratio"
      xlabel <- "SD-ratio"
      xlimit <- c(.1,10)
    } 

    # Split values as a function of all parameters but the one of interest    
    X1 <- split(round(Sel[,6],3),discr) # Cohen
    X1 <- compact(X1)


    X2 <- split(round(Sel[,7],3),discr) # Hedges
    X2 <- compact(X2)

    X3 <- split(round(Sel[,8],3),discr) # Glass1
    X3 <- compact(X3)

    X4 <- split(round(Sel[,9],3),discr) # Glass2
    X4 <- compact(X4)

    X5 <- split(round(Sel[,10],3),discr) # Shieh
    X5 <- compact(X5)

    X6 <- split(round(Sel[,11],3),discr) # Transformed shieh
    X6 <- compact(X6)
    
    
    allval <- unlist(c(X1,X2,X3,X4,X5,X6))
    
    for (k in seq_len(length(X1))){
      
      Param <- str_extract_all(Filestxt[i], "[[:digit:]]+\\.*[[:digit:]]*")
      G1 <- as.numeric(Param[[1]][2])
      G2 <- as.numeric(Param[[1]][4])
      
      N1=as.numeric(str_extract_all(names(X1)[k],"[[:digit:]]+\\.")[[1]][1])
      N2=as.numeric(str_extract_all(names(X1)[k],"[[:digit:]]+\\.")[[1]][2])
      
      if (fct=="meandiff"){
        params<-paste0("sd-ratio=",str_extract_all(names(X1)[k], "[[:digit:]]+\\.*[[:digit:]]*")[[1]][2])
      } else {params<-paste0("m1-m2=",str_extract_all(names(X1)[k], "[[:digit:]]+\\.*[[:digit:]]*")[[1]][2])}
      main3=paste0("n1=",N1,";n2=",N2,";",params)
      
      setwd(chemin)
      png(file=paste0("plot",i, "&",k,".png"),width=2000,height=2000, units = "px", res = 300)  
       
      plot(0,ylim=c(min(allval),max(allval)),xlim=xlimit,main=paste0(main1,main2,"\n when G1 = ",G1, " and G2 = ", G2,"\n ",main3),col="white",ylab=ylabel,xaxt="n",bty="n",xlab=xlabel)
      if (fct=="meandiff"){axis(1,at=1:4,labels=1:4)
      } else {
        
        if (condition=="Hom_bal"){
               axis(1,at=1,labels=1)
        } else if (condition=="Hom_rnull"){
               axis(1,at=1,labels=1)
        } else if (condition=="Het_bal") {
               axis(1,at=1:6,labels=c(.1,.25,.5,2,4,10))
        } else if (condition=="Het_rpos"){
               axis(1,at=1:6,labels=c(.1,.25,.5,2,4,10))
        } else if(condition=="Het_rneg") {
               axis(1,at=1:6,labels=c(.1,.25,.5,2,4,10))
        }
        
      }

      lines(X2[[k]],col="darkgrey") # Hedges 
      points(X2[[k]],col="darkgrey",pch=19,cex=.5) # Hedges
      text(X2[[k]], labels=X2[[k]], cex= 0.7, offset = 10)
      
      lines(X1[[k]],col="black",lty=2) #Cohen  
      points(X1[[k]],col="black",pch=19,cex=.5) #Cohen 
      text(X1[[k]], labels=X1[[k]], cex= 0.7, offset = 10)
      
      lines(X4[[k]],col="pink") # Glass sd2
      points(X4[[k]],col="pink",pch=19,cex=.5) # Glass sd2
      text(X4[[k]], labels=X4[[k]], cex= 0.7, offset = 10)
      
      lines(X3[[k]],col="violet",lty=2) # Glass sd1
      points(X3[[k]],col="violet",pch=19,cex=.5) # Glass sd1
      text(X3[[k]], labels=X3[[k]], cex= 0.7, offset = 10)
      
      lines(X6[[k]],col="lightblue") # Transformed shieh
      points(X6[[k]],col="lightblue",pch=19,cex=.5) # Transformed shieh
      text(X6[[k]], labels=X6[[k]], cex= 0.7, offset = 10)
      
      lines(X5[[k]],col="darkblue",lty=2) # Shieh
      points(X5[[k]],col="darkblue",pch=19,cex=.5) # Shieh
      text(X5[[k]], labels=X5[[k]], cex= 0.7, offset = 10)

      diff=c(max(X1[k][[1]])-min(X1[k][[1]]),
             max(X2[k][[1]])-min(X2[k][[1]]),
             max(X3[k][[1]])-min(X3[k][[1]]),
             max(X4[k][[1]])-min(X4[k][[1]]),
             max(X5[k][[1]])-min(X5[k][[1]]),
             max(X6[k][[1]])-min(X6[k][[1]]))

      ratio=c(max(X1[k][[1]])/min(X1[k][[1]]),
              max(X2[k][[1]])/min(X2[k][[1]]),
              max(X3[k][[1]])/min(X3[k][[1]]),
              max(X4[k][[1]])/min(X4[k][[1]]),
              max(X5[k][[1]])/min(X5[k][[1]]),
              max(X6[k][[1]])/min(X6[k][[1]]))
      
      if(param=="bias" || param=="relbias"){
            if(max(diff)==diff[1]){
              Legend=paste0("Max diff = Cohen (=",round(max(X1[k][[1]])-min(X1[k][[1]]),3),")")
            } else if (max(diff)==diff[2]){
              Legend=paste0("Max diff = Hedges (=",round(max(X2[k][[1]])-min(X2[k][[1]]),3),")")
            } else if (max(diff)==diff[3]){
              Legend=paste0("Max diff = Glass1 (=",round(max(X3[k][[1]])-min(X3[k][[1]]),3),")")
            } else if (max(diff)==diff[4]){
              Legend=paste0("Max diff = Glass2 (=",round(max(X4[k][[1]])-min(X4[k][[1]]),3),")")
            } else if (max(diff)==diff[5]){
              Legend=paste0("Max diff = Shieh (=",round(max(X5[k][[1]])-min(X5[k][[1]]),3),")")
            } else if (max(diff)==diff[6]){
              Legend=paste0("Max diff = TrShieh (=",round(max(X6[k][[1]])-min(X6[k][[1]]),3),")")
            }
      } else if (param=="eff"||param=="releff"){
             if(max(ratio)==ratio[1]){
             Legend=paste0("Max ratio = Cohen (=",round(max(X1[k][[1]])/min(X1[k][[1]]),3),")")
             } else if (max(ratio)==ratio[2]){
             Legend=paste0("Max ratio = Hedges (=",round(max(X2[k][[1]])/min(X2[k][[1]]),3),")")
             } else if (max(ratio)==ratio[3]){
             Legend=paste0("Max ratio = Glass1 (=",round(max(X3[k][[1]])/min(X3[k][[1]]),3),")")
             } else if (max(ratio)==ratio[4]){
             Legend=paste0("Max ratio = Glass2 (=",round(max(X4[k][[1]])/min(X4[k][[1]]),3),")")
             } else if (max(ratio)==ratio[5]){
             Legend=paste0("Max ratio = Shieh (=",round(max(X5[k][[1]])/min(X5[k][[1]]),3),")")
             } else if (max(ratio)==ratio[6]){
             Legend=paste0("Max diff = TrShieh (=",round(max(X6[k][[1]])/min(X6[k][[1]]),3),")")
        }
        
      }

      legend("top",legend=Legend,col="black",pch=19,bty="n")
      
      dev.off()
      }

  }

}

### REM: la partie "max diff" dans mon script est pertinente pour meandiff, mais pas pr sdratio ###

## As a function of the sd-ratio

# condition Het_bal
Comparison(condition="Het_bal",param="bias",fct="sdratio",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/sdratio/G1=0/Het_bal/bias")  
Comparison(condition="Het_bal",param="relbias",fct="sdratio",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/sdratio/G1=0/Het_bal/relative bias")  
Comparison(condition="Het_bal",param="eff",fct="sdratio",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/sdratio/G1=0/Het_bal/variance")  
Comparison(condition="Het_bal",param="releff",fct="sdratio",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/sdratio/G1=0/Het_bal/variance relative")  

# condition Het_rneg
Comparison(condition="Het_rneg",param="bias",fct="sdratio",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/sdratio/G1=0/Het_rneg/bias")  
Comparison(condition="Het_rneg",param="relbias",fct="sdratio",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/sdratio/G1=0/Het_rneg/relative bias")  
Comparison(condition="Het_rneg",param="eff",fct="sdratio",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/sdratio/G1=0/Het_rneg/variance")  
Comparison(condition="Het_rneg",param="releff",fct="sdratio",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/sdratio/G1=0/Het_rneg/variance relative")  

# condition Het_rpos
Comparison(condition="Het_rpos",param="bias",fct="sdratio",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/sdratio/G1=0/Het_rpos/bias")  
Comparison(condition="Het_rpos",param="relbias",fct="sdratio",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/sdratio/G1=0/Het_rpos/relative bias")  
Comparison(condition="Het_rpos",param="eff",fct="sdratio",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/sdratio/G1=0/Het_rpos/variance")  
Comparison(condition="Het_rpos",param="releff",fct="sdratio",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/sdratio/G1=0/Het_rpos/variance relative")  




# Variations as a function of mean diff=

# Condition Hom_bal
Comparison(condition="Hom_bal",param="bias",fct="meandiff",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/Diff de moyenne/G1=0/Hom_bal/bias")  
Comparison(condition="Hom_bal",param="relbias",fct="meandiff",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/Diff de moyenne/G1=0/Hom_bal/relative bias")  
Comparison(condition="Hom_bal",param="eff",fct="meandiff",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/Diff de moyenne/G1=0/Hom_bal/variance")  
Comparison(condition="Hom_bal",param="releff",fct="meandiff",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/Diff de moyenne/G1=0/Hom_bal/variance relative")  

# Condition Hom_rnull
Comparison(condition="Hom_rnull",param="bias",fct="meandiff",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/Diff de moyenne/G1=0/Hom_rnul/bias") 
Comparison(condition="Hom_rnull",param="relbias",fct="meandiff",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/Diff de moyenne/G1=0/Hom_rnul/relative bias")  
Comparison(condition="Hom_rnull",param="eff",fct="meandiff",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/Diff de moyenne/G1=0/Hom_rnul/variance")  
Comparison(condition="Hom_rnull",param="releff",fct="meandiff",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/Diff de moyenne/G1=0/Hom_rnul/variance relative")  

# condition Het_bal
Comparison(condition="Het_bal",param="bias",fct="meandiff",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/Diff de moyenne/G1=0/Het_bal/bias")  
Comparison(condition="Het_bal",param="relbias",fct="meandiff",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/Diff de moyenne/G1=0/Het_bal/relative bias")  
Comparison(condition="Het_bal",param="eff",fct="meandiff",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/Diff de moyenne/G1=0/Het_bal/variance")  
Comparison(condition="Het_bal",param="releff",fct="meandiff",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/Diff de moyenne/G1=0/Het_bal/variance relative")  

# condition Het_rpos
Comparison(condition="Het_rpos",param="bias",fct="meandiff",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/Diff de moyenne/G1=0/Het_rpos/bias")  
Comparison(condition="Het_rpos",param="relbias",fct="meandiff",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/Diff de moyenne/G1=0/Het_rpos/relative bias")  
Comparison(condition="Het_rpos",param="eff",fct="meandiff",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/Diff de moyenne/G1=0/Het_rpos/variance")  
Comparison(condition="Het_rpos",param="releff",fct="meandiff",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/Diff de moyenne/G1=0/Het_rpos/variance relative")  

# condition Het_rneg
Comparison(condition="Het_rneg",param="bias",fct="meandiff",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/Diff de moyenne/G1=0/Het_rneg/bias")  
Comparison(condition="Het_rneg",param="relbias",fct="meandiff",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/Diff de moyenne/G1=0/Het_rneg/relative bias")  
Comparison(condition="Het_rneg",param="eff",fct="meandiff",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/Diff de moyenne/G1=0/Het_rneg/variance")  
Comparison(condition="Het_rneg",param="releff",fct="meandiff",chemin="C:/Users/Marie/Dropbox/graphs analyses goodness/Diff de moyenne/G1=0/Het_rneg/variance relative")  
