library(stringr)

### Ajout de glass 1 et 2 non biaisés dans mes bdd
### Idem pour Shieh's d
### Multiplicaiton de transformed Shieh's d par 2, pour utiliser les benchmarks du d de Cohen
### Ajout de notre transformation débiaisée (mais là suis pas sûre à 100%...)

Mainfolder="C:/Users/Marie/Documents/ES MEASURES/"
subfolder=list.files(Mainfolder)
Folder=paste0(Mainfolder,subfolder)

for (i in 4){
  
  for (j in seq_len(length(list.files(Folder[i])))){ 
    
    filepath = paste0(Folder[i],"/",list.files(Folder[i])[j])
    file=readRDS(filepath)
    param <- str_extract_all(list.files(Folder[i])[j], "[[:digit:]]+\\.*[[:digit:]]*")
    if(as.numeric(param[[1]][2])==2.08){
      G1=-2.08
    } else {G1=as.numeric(param[[1]][2])}
    G2 <- as.numeric(param[[1]][4])
    n1 <- as.numeric(param[[1]][5])
    n2 <- as.numeric(param[[1]][6])
    N <- n1+n2
    m1 <- as.numeric(param[[1]][7])
    m2 <- as.numeric(param[[1]][8])
    sd1 <- as.numeric(param[[1]][9])
    sd2 <- as.numeric(param[[1]][10])
    
    # Calculer les glass's non biaisés
    unbiased_glass_sd1 <- file[,11]*gamma((n1-1)/2)/(sqrt((n1-1)/2)*gamma((n1-2)/2)) # Hedge's correction
    unbiased_glass_sd2 <- file[,12]*gamma((n2-1)/2)/(sqrt((n2-1)/2)*gamma((n2-2)/2)) # Hedge's correction

    # calculer le Shieh non biaisés
    # Faire l'équivalent avec la correction du biais de Hedges:
    df <- (file[,3]^2/n1+file[,4]^2/n2)^2/((file[,3]^2/n1)^2/(n1-1)+(file[,4]^2/n2)^2/(n2-1))
    correction <- gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2))
    unbiased_shieh <- file[,13]*correction

    # multiplier par 2 notre estimation du Shieh (pour les benchmark de Cohen)
    transformed_shieh <- file[,14]*2
    
    # Faire l'équivalent avec la correction du biais de Hedges:
    df <- (file[,3]^2/n1+file[,4]^2/n2)^2/((file[,3]^2/n1)^2/(n1-1)+(file[,4]^2/n2)^2/(n2-1))
    corr <- gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2))
    unbiased_transformed_shieh <- file[,17]*corr
    
    newfile<- cbind(file[,1:12],unbiased_glass_sd1,unbiased_glass_sd2,file[,13],unbiased_shieh,transformed_shieh,unbiased_transformed_shieh)
    colnames(newfile)[15]="Shieh" 
    
    # Extraction of the ES matrix
    chem <- "C:/Users/Marie/Documents/ES MEASURES/"
    sschem<- paste0("G1=",G1,",G2=",G2) 
    setwd(dir=paste0(chem,sschem)) # destination file  
    
    fname=paste0("G1=",G1,", G2=",G2,", n=[",n1,",",n2,"], means=[",m1,",",m2,"], sds=[",sd1,",",sd2,"].rds")
    saveRDS(newfile, file = fname)
    
  }
  
}



# Note: multiplier notre transformation par 2 n'en modifie absolument pas le biais et la variance relatifs.


#--------------------------------------------------------------------------------------------------------------------------------------------------------------------
#     GRAPH: PLOTTING bias, relative bias, and variance
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

library(stringr)

Path <-  "C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/"
for (j in seq_len(length(list.files(Path,pattern = ".*shiehvsshiehbis.txt")))){
  

  File=read.table(paste0(Path,list.files(Path,pattern = ".*shiehvsshiehbis.txt")[j]),header=T,sep=";",dec=",")  
  # Subdivise file into five categories
  # "Hom_bal", "Hom_unbal", "Het_bal","Het_rpos","Het_rneg"
  
  ### Conditions id 
  
  # Homoscedasticity, balanced designs ("Hom_bal")
  id_Hom_bal=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1==File$n2)&(File$sd1.sd2 == 1),]))
  # Homoscedasticity, unbalanced designs ("Hom_unbal")  
  id_Hom_nN=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 == 1),]))
  id_Hom_Nn=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 == 1),]))
  # Heteroscedasticity, balanced designs ("Het_bal")
  id_sdSD_bal=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1==File$n2)&(File$sd1.sd2 < 1),]))
  id_SDsd_bal=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1==File$n2)&(File$sd1.sd2 > 1),]))
  # Heteroscedasticity, positive correlation between n and sd ("Het_rpos")
  id_sdSD_nN=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 < 1),]))
  id_SDsd_Nn=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 > 1),]))
  # Heteroscedasticity, negative correlation between n and sd ("Het_rneg")
  id_SDsd_nN=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 > 1),]))
  id_sdSD_Nn=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 < 1),]))
  
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
    res <- matrix(0,9,4)  
    names<-expand.grid("bias_",c("Shieh_corr","Shieh_corrbis"))
    colnames(res) <- c("n1","n2",paste0(names[,1],names[,2]))
    res[,1:2] <- cbind(combi[,1],combi[,2])
    res[,3] <- tapply(Sel$bias_Shieh_corr,list(Sel$n1,Sel$n2),mean)[1:9]
    res[,4] <- tapply(Sel$bias_Shieh_corrbis,list(Sel$n1,Sel$n2),mean)[1:9]
    
    # Select only rows with no "NA"  
    res <- subset(res,res[,3] != "NA") 
    # Select only bias columns
    res.bias <- t(res[,3:4])
    colnames(res.bias) <- paste0(res[,1],":",res[,2])
    
    # Matrix containing relative biases information
    res2 <- matrix(0,9,4)  
    names<-expand.grid("relbias_",c("Shieh_corr","Shieh_corrbis"))
    colnames(res2) <- c("n1","n2",paste0(names[,1],names[,2]))
    res2[,1:2] <- cbind(combi[,1],combi[,2])
    res2[,3] <- tapply(Sel$relbias_Shieh_corr,list(Sel$n1,Sel$n2),mean)[1:9]
    res2[,4] <- tapply(Sel$relbias_Shieh_corrbis,list(Sel$n1,Sel$n2),mean)[1:9]
    
    # Select only rows with no "NA"  
    res2 <- subset(res2,res2[,3] != "NA") 
    # Select only bias columns
    res.relbias <- t(res2[,3:4])
    colnames(res.relbias) <- paste0(res2[,1],":",res2[,2])
    
    param <- str_extract_all(list.files(Path,pattern=".*shiehvsshiehbis.txt")[j], "[[:digit:]]+\\.*[[:digit:]]*")
    if(param[[1]][2]==2.08){
      G1 <- -as.numeric(param[[1]][2])
    } else (G1 <- as.numeric(param[[1]][2]))
    G2 <- as.numeric(param[[1]][4])
    
    # Matrix containing variances information
    res3 <- matrix(0,9,4)  
    names<-expand.grid("var_",c("Shieh_corr","Shieh_corrbis"))
    colnames(res3) <- c("n1","n2",paste0(names[,1],names[,2]))
    res3[,1:2] <- cbind(combi[,1],combi[,2]) 
    res3[,3] <- tapply(Sel$eff_Shieh_corr,list(Sel$n1,Sel$n2),mean)[1:9]
    res3[,4] <- tapply(Sel$eff_Shieh_corrbis,list(Sel$n1,Sel$n2),mean)[1:9]
    
    # Select only rows with no "NA"  
    res3 <- subset(res3,res3[,3] != "NA") 
    # Select only bias columns
    res.eff <- t(res3[,3:4])
    colnames(res.eff) <- paste0(res3[,1],":",res3[,2])
    
    # Matrix containing MSE information
    res4 <- matrix(0,9,4)  
    names<-expand.grid("releff_",c("Shieh_corr","Shieh_corrbis"))
    colnames(res4) <- c("n1","n2",paste0(names[,1],names[,2]))
    res4[,1:2] <- cbind(combi[,1],combi[,2])
    res4[,3] <- tapply(Sel$releff_Shieh_corr,list(Sel$n1,Sel$n2),mean)[1:9]
    res4[,4] <- tapply(Sel$releff_Shieh_corrbis,list(Sel$n1,Sel$n2),mean)[1:9]
    
    # Select only rows with no "NA"  
    res4 <- subset(res4,res4[,3] != "NA") 
    # Select only bias columns
    res.releff <- t(res4[,3:4])
    colnames(res.releff) <- paste0(res4[,1],":",res4[,2])
    
    setwd(paste0("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Réflexion interprétation simulations/Plots/",names(Conditions_id)[i]))
    png(file=paste0("bias_eff,G1=",G1, " & G2=",G2,";",names(Conditions_id)[i], ".png"),width=900,height=1700, units = "px", res = 300)  
    
    par(mar = c(4,5,1.5,0),mfrow = c(4,1))   
    
    
    # plot for the bias
    
    if (j==1){ylabelbias=expression(paste("E(" , hat(delta) , ") -",delta ))
    } else {ylabelbias=""}
    
    barplot(res.bias, 
            col = c("black","grey70"),
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
    
    barplot(res.relbias, 
            col = c("black","grey70"),
            beside = TRUE,
            xaxt="n",
            cex.lab=1.5,
            cex.main=1.5,
            ylab=ylabelbias
    )
    # plot for the variance
    
    if (j==1){ylabeleff=expression(paste("Var(" , hat(delta) , ")"))
    } else {ylabeleff=""}
    
    barplot(res.eff, 
            col = c("black","grey70"),
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
    
    barplot(res.releff, 
            col = c("black","grey70"),
            beside = TRUE,
            ylab = ylabeleff,
            cex.lab=1.5,
            cex.names=.8,
            args.legend = list(
              x = length(res.eff)*1.2,
              y = max(res.eff)+.5,
              bty="n"
              
            ))
    
    dev.off()
    
  }
  
}

