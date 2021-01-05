library(stringr)

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------
#     GRAPH: PLOTTING bias, relative bias, and EFFICIENCY
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

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


setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/when mu1=mu2/")

png(file="legend.png",width=1500,height=1000, units = "px", res = 300)  

plot(1,1,bty="n",xaxt="n",yaxt="n",ylim=c(.62,.67),main="",xlab="",ylab="",pch=19,type="o")
legend("center", 
       legend=c("Cohen's d","Hedge's g","Glass's delta (delta = sd1)","Glass's delta (delta = sd2)","Shieh's d","Corrected Shieh's d"),
       fill=c("black","grey25","grey50","grey70","grey90","white"),
       bty="n"
)

dev.off()

Path <-  "C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/"
for (j in seq_len(length(list.files(Path,pattern = ".*good_mes.txt")))){
  

  File=read.table(paste0(Path,list.files(Path,pattern = ".*good_mes.txt")[j]),header=T,sep=";",dec=",")  
  
  # Subdivise file into five categories
  # "Hom_bal", "Hom_unbal", "Het_bal","Het_rpos","Het_rneg"
  
  ### Conditions id 
  
  # Homoscedasticity, balanced designs ("Hom_bal")
  id_Hom_bal=as.numeric(rownames(File[(File$m1.m2==0)&(File$n1==File$n2)&(File$sd1.sd2 == 1),]))
  # Homoscedasticity, unbalanced designs ("Hom_unbal")  
  id_Hom_nN=as.numeric(rownames(File[(File$m1.m2==0)&(File$n1<File$n2)&(File$sd1.sd2 == 1),]))
  id_Hom_Nn=as.numeric(rownames(File[(File$m1.m2==0)&(File$n1>File$n2)&(File$sd1.sd2 == 1),]))
  # Heteroscedasticity, balanced designs ("Het_bal")
  id_sdSD_bal=as.numeric(rownames(File[(File$m1.m2==0)&(File$n1==File$n2)&(File$sd1.sd2 < 1),]))
  id_SDsd_bal=as.numeric(rownames(File[(File$m1.m2==0)&(File$n1==File$n2)&(File$sd1.sd2 > 1),]))
  # Heteroscedasticity, positive correlation between n and sd ("Het_rpos")
  id_sdSD_nN=as.numeric(rownames(File[(File$m1.m2==0)&(File$n1<File$n2)&(File$sd1.sd2 < 1),]))
  id_SDsd_Nn=as.numeric(rownames(File[(File$m1.m2==0)&(File$n1>File$n2)&(File$sd1.sd2 > 1),]))
  # Heteroscedasticity, negative correlation between n and sd ("Het_rneg")
  id_SDsd_nN=as.numeric(rownames(File[(File$m1.m2==0)&(File$n1<File$n2)&(File$sd1.sd2 > 1),]))
  id_sdSD_Nn=as.numeric(rownames(File[(File$m1.m2==0)&(File$n1>File$n2)&(File$sd1.sd2 < 1),]))
  
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
    names<-expand.grid("bias_",c("Cohen","Hedge","Glass1","Glass2","Shieh","Shieh_corr"))
    colnames(res) <- c("n1","n2",paste0(names[,1],names[,2]))
    res[,1:2] <- cbind(combi[,1],combi[,2])
    res[,3] <- tapply(Sel$bias_Cohen,list(Sel$n1,Sel$n2),mean)[1:9]
    res[,4] <- tapply(Sel$bias_Hedge,list(Sel$n1,Sel$n2),mean)[1:9]
    res[,5] <- tapply(Sel$bias_Glass1,list(Sel$n1,Sel$n2),mean)[1:9]
    res[,6] <- tapply(Sel$bias_Glass2,list(Sel$n1,Sel$n2),mean)[1:9]
    res[,7] <- tapply(Sel$bias_Shieh,list(Sel$n1,Sel$n2),mean)[1:9]
    res[,8] <- tapply(Sel$bias_Shieh_corr,list(Sel$n1,Sel$n2),mean)[1:9]
    # Select only rows with no "NA"  
    res <- subset(res,res[,3] != "NA") 
    # Select only bias columns
    res.bias <- t(res[,3:8])
    colnames(res.bias) <- paste0(res[,1],":",res[,2])
    
    param <- str_extract_all(list.files(Path,pattern=".*good_mes.txt")[j], "[[:digit:]]+\\.*[[:digit:]]*")
    if(param[[1]][2]==2.08){
      G1 <- -as.numeric(param[[1]][2])
    } else (G1 <- as.numeric(param[[1]][2]))
    G2 <- as.numeric(param[[1]][4])
    
    
    # Matrix containing variances information
    res3 <- matrix(0,9,8)  
    names<-expand.grid("var_",c("Cohen","Hedge","Glass1","Glass2","Shieh","Shieh_corr"))
    colnames(res3) <- c("n1","n2",paste0(names[,1],names[,2]))
    res3[,1:2] <- cbind(combi[,1],combi[,2]) 
    res3[,3] <- tapply(Sel$eff_Cohen,list(Sel$n1,Sel$n2),mean)[1:9]
    res3[,4] <- tapply(Sel$eff_Hedge,list(Sel$n1,Sel$n2),mean)[1:9]
    res3[,5] <- tapply(Sel$eff_Glass1,list(Sel$n1,Sel$n2),mean)[1:9]
    res3[,6] <- tapply(Sel$eff_Glass2,list(Sel$n1,Sel$n2),mean)[1:9]
    res3[,7] <- tapply(Sel$eff_Shieh,list(Sel$n1,Sel$n2),mean)[1:9]
    res3[,8] <- tapply(Sel$eff_Shieh_corr,list(Sel$n1,Sel$n2),mean)[1:9]
    # Select only rows with no "NA"  
    res3 <- subset(res3,res3[,3] != "NA") 
    # Select only bias columns
    res.eff <- t(res3[,3:8])
    colnames(res.eff) <- paste0(res3[,1],":",res3[,2])
    
    setwd(paste0("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/when mu1=mu2/",names(Conditions_id)[i]))
    #dir.create(names(Conditions_id)[i])
    png(file=paste0("bias_eff,G1=",G1, " & G2=",G2,";",names(Conditions_id)[i], ".png"),width=900,height=1700, units = "px", res = 300)  
    
    par(xpd = FALSE,mar = c(4,5,1.5,0),mfrow = c(2,1))   
    
    
    # plot for the bias
    
    if (j==1){ylabelbias=expression(paste("E(" , hat(delta) , ") -",delta ))
    } else {ylabelbias=""}
    
    barplot(res.bias, 
            col = c("black","grey25","grey50","grey70","grey90","white"),
            beside = TRUE,
            main=paste0("G1=",G1,"; G2=",G2),
            xaxt="n",
            cex.lab=1.5,
            cex.main=1.5,
            ylab=ylabelbias
    )
    
    # plot for the relative bias
    
    if (j==1){ylabelbias=expression(paste("(E(" , hat(delta) , ") -",delta,")/",delta ))
    } else {ylabelbias=""}
    
    barplot(res.eff, 
            col = c("black","grey25","grey50","grey70","grey90","white"),
            beside = TRUE,
            ylab = ylabeleff,
            cex.lab=1.5,
            xaxt="n",
            args.legend = list(
              x = length(res.eff)*1.2,
              y = max(res.eff)+.5,
              bty="n"
            ))
    

    # plot the the relative variance
    if (j==1){ylabeleff=expression(paste("Var(" , hat(delta) , ")/",delta^2))
    } else {ylabeleff=""}
    
    dev.off()
    
  }
  
}
