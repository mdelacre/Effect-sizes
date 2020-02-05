library(stringr)

# First of all, I'm going to make different graph for each category


# Codes for subcategories

# bal vs. (nN + Nn) 
##### bal = equal n across groups 
########## nN = the first group has the smallest sample size 
########## Nn = the first group has the biggest sample size

# Hom vs. Het (sdSD + SDsd) 
##### Hom = homoscedasticity
##### Het = heteroscedasticity
########## sdSD = the first group has the smallest sd 
########## SDsd = the first group has the biggest sd

## When there is heteroscedasticity:
# If sdSD & nN or SDsd & Nn: positive correlation between n and sd (rpos)
# If sdSD & Nn or SDsd & nN: negative correlation between n and sd (rneg)
# If sdSD & bal or SDsd & bal: no correlation between n and sd (rnull)


  for (j in seq_len(length(list.files(Path)))){
   Path <-  "C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of measures/"
   setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Graphs")

    File=read.table(paste0(Path,list.files(Path)[j]),header=T,sep=";",dec=",")  
    
    # Subdivise file into five categories
    # "Hom_bal", "Hom_unbal", "Het_bal","Het_rpos","Het_rneg"
    
    ### Conditions id 
    
    # Homoscedasticity, balanced designs ("Hom_bal")
    id_Hom_bal=as.numeric(rownames(File[(File$n1==File$n2)&(File$sd1.sd2 == 1),]))
    # Homoscedasticity, unbalanced designs ("Hom_unbal")  
    id_Hom_nN=as.numeric(rownames(File[(File$n1<File$n2)&(File$sd1.sd2 == 1),]))
    id_Hom_Nn=as.numeric(rownames(File[(File$n1>File$n2)&(File$sd1.sd2 == 1),]))
    # Heteroscedasticity, balanced designs ("Het_bal")
    id_sdSD_bal=as.numeric(rownames(File[(File$n1==File$n2)&(File$sd1.sd2 < 1),]))
    id_SDsd_bal=as.numeric(rownames(File[(File$n1==File$n2)&(File$sd1.sd2 > 1),]))
    # Heteroscedasticity, positive correlation between n and sd ("Het_rpos")
    id_sdSD_nN=as.numeric(rownames(File[(File$n1<File$n2)&(File$sd1.sd2 < 1),]))
    id_SDsd_Nn=as.numeric(rownames(File[(File$n1>File$n2)&(File$sd1.sd2 > 1),]))
    # Heteroscedasticity, negative correlation between n and sd ("Het_rneg")
    id_SDsd_nN=as.numeric(rownames(File[(File$n1<File$n2)&(File$sd1.sd2 > 1),]))
    id_sdSD_Nn=as.numeric(rownames(File[(File$n1>File$n2)&(File$sd1.sd2 < 1),]))
    
    # Study each condition separately
    Conditions_id <- list(id_Hom_bal=id_Hom_bal,
                          id_Hom_rnull=c(id_Hom_nN,id_Hom_Nn),
                          id_Het_bal=c(id_sdSD_bal,id_SDsd_bal), 
                          id_Het_rpos=c(id_sdSD_nN,id_SDsd_Nn),
                          id_Het_rneg=c(id_SDsd_nN,id_sdSD_Nn))
    
    for (i in seq_len(length(Conditions_id))){ 
      
      Sel <- File[Conditions_id[[i]],]
      
      n2val <- as.numeric(levels(factor(Sel$n2)))
      n1val<- as.numeric(levels(factor(Sel$n1)))
      combi <- expand.grid(n2val, n1val)

      # Matrix containing biases information
      res <- matrix(0,9,8)  
      names<-expand.grid("MSE_",c("Cohen","Hedge","Glass1","Glass2","Shieh","Shieh_corr"))
      colnames(res) <- c("n1","n2",paste0(names[,1],names[,2]))
      res[,1:2] <- cbind(combi[,1],combi[,2])
      res[,3] <- tapply(Sel$MSE_Cohen,list(Sel$n1,Sel$n2),mean)[1:9]
      res[,4] <- tapply(Sel$MSE_Hedge,list(Sel$n1,Sel$n2),mean)[1:9]
      res[,5] <- tapply(Sel$MSE_Glass1,list(Sel$n1,Sel$n2),mean)[1:9]
      res[,6] <- tapply(Sel$MSE_Glass2,list(Sel$n1,Sel$n2),mean)[1:9]
      res[,7] <- tapply(Sel$MSE_Shieh,list(Sel$n1,Sel$n2),mean)[1:9]
      res[,8] <- tapply(Sel$MSE_Shieh_corr,list(Sel$n1,Sel$n2),mean)[1:9]
      # Select only rows with no "NA"  
      res <- subset(res,res[,3] != "NA") 
      # Select only bias columns
      res.MSE <- t(res[,3:8])
      colnames(res.MSE) <- paste0(res[,1],":",res[,2])
      
      param <- str_extract_all(list.files(Path)[j], "[[:digit:]]+\\.*[[:digit:]]*")
      G1 <- as.numeric(param[[1]][2])
      G2 <- as.numeric(param[[1]][4])
      
      png(file=paste0("MSE,G1=",G1, " & G2=",G2,";",names(Conditions_id)[i], ".png"),width=2500,height=1700, units = "px", res = 300)  
      
    
      par(xpd = F,mar = c(4,5,0.5,0),mfrow = c(2,1))   
      
      # legend
      plot(1,1,bty="n",xaxt="n",yaxt="n",ylim=c(.62,.67),main="",xlab="",ylab="",pch=19,type="o")
      legend("center", 
             legend=c("Cohen's d","Hedge's g","Glass's delta (delta = sd1)","Glass's delta (delta = sd2)","Shieh's d","Corrected Shieh's d"),
             fill=c("black","grey25","grey50","grey70","grey90","white"),
             bty="n"
            )
      

      # plot for the variance
      barplot(res.MSE, 
              col = c("black","grey25","grey50","grey70","grey90","white"),
              beside = TRUE,
              ylab = expression(paste( "MSE(",hat(delta), ") =",hat(delta)[efficiency]+hat(delta)[bias]^2)),
              xlab = "n1 : n2",
              args.legend = list(
                x = length(res.bias)*1.2,
                y = max(res.bias)+.5,
                bty="n"
              ))
      
     dev.off()
      
    }
    
  }




