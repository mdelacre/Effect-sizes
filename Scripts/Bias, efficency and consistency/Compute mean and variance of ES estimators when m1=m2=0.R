library(stringr)

setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of measures")
dir.create("When m1=m2")

Mainfolder="C:/Users/Marie/Documents/ES MEASURES/"
subfolder=list.files(Mainfolder)
Folder=paste0(Mainfolder,subfolder)

for (i in seq_len(length(Folder))){
  
  # Extract files where m1=m2=0
  Files <- list.files(Folder[i])[str_detect(list.files(Folder[i]), "0,0")==TRUE]
  
  # set up empty container for all estimated parameters
  average_ES <-matrix(0,length(Files),19)
  ES_indic <- c("mean_","var_")
  estimator <- c("Cohen","Hedge","Glass1","Glass2","Shieh", "Shieh_corr")
  columns_names=expand.grid(estimator,ES_indic)
  colnames(average_ES) <- c("n1","n2","n1/n2","m1-m2","sd1/sd2","G1","G2",paste0(columns_names[,2],columns_names[,1]))
  
  for (j in seq_len(length(Files))){

    filepath = paste0(Folder[i],"/",Files[j])
    file=readRDS(filepath)
    param <- str_extract_all(Files[j], "[[:digit:]]+\\.*[[:digit:]]*")
    G1 <- as.numeric(param[[1]][2])
    G2 <- as.numeric(param[[1]][4])
    n1 <- as.numeric(param[[1]][5])
    n2 <- as.numeric(param[[1]][6])
    m1 <- as.numeric(param[[1]][7])
    m2 <- as.numeric(param[[1]][8])
    sd1 <- as.numeric(param[[1]][9])
    sd2 <- as.numeric(param[[1]][10])

    # Compute ES means
    mean_cohen <- mean(file[,9])
    mean_hedge <- mean(file[,10])
    mean_glass1 <- mean(file[,11]) 
    mean_glass2 <- mean(file[,12])
    mean_shieh <- mean(file[,13])
    mean_shieh_corr <- mean(file[,14])
    
    # Compute Efficiency
    var_cohen <- var(file[,9])
    var_hedge <- var(file[,10])
    var_glass1 <- var(file[,11])
    var_glass2 <- var(file[,12])
    var_shieh <- var(file[,13])
    var_shieh_corr <- var(file[,14])
    
    average_ES[j,1] <- n1
    average_ES[j,2] <- n2
    average_ES[j,3] <- n1/n2
    average_ES[j,4] <- m1-m2
    average_ES[j,5] <- sd1/sd2
    average_ES[j,6] <- G1
    average_ES[j,7] <- G2
    
    average_ES[j,8] <- mean_cohen
    average_ES[j,9] <- mean_hedge
    average_ES[j,10] <- mean_glass1
    average_ES[j,11] <- mean_glass2
    average_ES[j,12] <- mean_shieh
    average_ES[j,13] <- mean_shieh_corr
    
    average_ES[j,14] <- var_cohen
    average_ES[j,15] <- var_hedge
    average_ES[j,16] <- var_glass1
    average_ES[j,17] <- var_glass2
    average_ES[j,18] <- var_shieh
    average_ES[j,19] <- var_shieh_corr
    
  }

  setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of measures/When m1=m2")
  
  sub<-paste0("G1=",G1,",G2=",G2) # Sign is not extracted but it's not a big deal
  write.table(average_ES,paste0(sub,",average_ES.txt"),sep=";",dec=",")
  
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
Path <-  "C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of measures/When m1=m2/"
setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Graphs")

for (j in seq_len(length(list.files(Path)))){

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
    
    # Matrix containing ES means information
    res.mean <- matrix(0,1,6)  
    names<-expand.grid("average_",c("Cohen","Hedge","Glass1","Glass2","Shieh","Shieh_corr"))
    colnames(res.mean) <- paste0(names[,1],names[,2])
    res.mean[,1] <- mean(Sel$mean_Cohen)
    res.mean[,2] <- mean(Sel$mean_Hedge)
    res.mean[,3] <- mean(Sel$mean_Glass1)
    res.mean[,4] <- mean(Sel$mean_Glass2)
    res.mean[,5] <- mean(Sel$mean_Shieh)
    res.mean[,6] <- mean(Sel$mean_Shieh_corr)
    # Select only rows with no "NA"
    
    # Matrix containing ES variances information
    res.var <- matrix(0,1,6)  
    names<-expand.grid("variance_",c("Cohen","Hedge","Glass1","Glass2","Shieh","Shieh_corr"))
    colnames(res.var) <- paste0(names[,1],names[,2])
    res.var[,1] <- var(Sel$mean_Cohen)
    res.var[,2] <- var(Sel$mean_Hedge)
    res.var[,3] <- var(Sel$mean_Glass1)
    res.var[,4] <- var(Sel$mean_Glass2)
    res.var[,5] <- var(Sel$mean_Shieh)
    res.var[,6] <- var(Sel$mean_Shieh_corr)
    # Select only rows with no "NA"


    param <- str_extract_all(list.files(Path)[j], "[[:digit:]]+\\.*[[:digit:]]*")
    G1 <- as.numeric(param[[1]][2])
    G2 <- as.numeric(param[[1]][4])
    
    #setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Graphs/")
    #dir.create("when m1=m2=0")
    setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Graphs/when m1=m2=0")
    
    png(file=paste0("mean,G1=",G1, " & G2=",G2,";",names(Conditions_id)[i], ".png"),width=2500,height=1700, units = "px", res = 300)  
    
    
    par(xpd = F,mar = c(4,5,0.5,0),mfrow = c(3,1))   
    
    # legend
    plot(1,1,bty="n",xaxt="n",yaxt="n",ylim=c(.62,.67),main="",xlab="",ylab="",pch=19,type="o")
    legend("center", 
           legend=c("Cohen's d","Hedge's g","Glass's delta (delta = sd1)","Glass's delta (delta = sd2)","Shieh's d","Corrected Shieh's d"),
           fill=c("black","grey25","grey50","grey70","grey90","white"),
           bty="n"
    )
    
    
    # plot for the average ES measure
    barplot(res.mean, 
            col = c("black","grey25","grey50","grey70","grey90","white"),
            beside = TRUE,
            ylab = "average estimator value",
            ylim = c(-.005,.005),
            args.legend = list(
              x = length(res.bias)*1.2,
              y = max(res.bias)+.5,
              bty="n"
            ))
    abline(h=0,lwd=2)
    
    # plot for the variance of ES measures
    barplot(res.var, 
            col = c("black","grey25","grey50","grey70","grey90","white"),
            beside = TRUE,
            ylab = "average variance of ES estimates",
            ylim = c(-.005,.005),
            args.legend = list(
              x = length(res.bias)*1.2,
              y = max(res.bias)+.5,
              bty="n"
            ))
    abline(h=0,lwd=2)
    
        dev.off()
    
  }
  
}

getwd()
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################

library(stringr)

setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of measures")
dir.create("When m1=m2")

Mainfolder="C:/Users/Marie/Documents/ES MEASURES/"
subfolder=list.files(Mainfolder)
Folder=paste0(Mainfolder,subfolder)

for (i in seq_len(length(Folder))){
  
  
  Files <- list.files(Folder[i])[str_detect(list.files(Folder[i]), "0,0")==TRUE]
  
  # set up empty container for all estimated parameters
  good_mes <-matrix(0,length(Files),5+3*6)
  goodness_indic <- c("bias_","eff_","MSE_")
  estimator <- c("Cohen","Hedge","Glass1","Glass2","Shieh", "Shieh_corr")
  columns_names=expand.grid(estimator,goodness_indic)
  colnames(good_mes) <- c("n1","n2","n1/n2","m1-m2","sd1/sd2",paste0(columns_names[,2],columns_names[,1]))
  
  for (j in seq_len(length(Files))){
    
    filepath = paste0(Folder[i],"/",Files[j])
    file=readRDS(filepath)
    param <- str_extract_all(Files[j], "[[:digit:]]+\\.*[[:digit:]]*")
    G1 <- as.numeric(param[[1]][2])
    G2 <- as.numeric(param[[1]][4])
    n1 <- as.numeric(param[[1]][5])
    n2 <- as.numeric(param[[1]][6])
    m1 <- as.numeric(param[[1]][7])
    m2 <- as.numeric(param[[1]][8])
    sd1 <- as.numeric(param[[1]][9])
    sd2 <- as.numeric(param[[1]][10])
    
    cohen_delta <- (m1-m2)/sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
    hedge_delta <- cohen_delta* (1-3/(4*(n1+n2-9)))
    glass_delta1 <- (m1-m2)/sd1
    glass_delta2 <- (m1-m2)/sd2
    q1 <- n1/(n1+n2)
    q2 <- n2/(n1+n2)
    shieh_delta <- (m1-m2)/sqrt(sd1^2/q1+sd2^2/q2)       
    nratio=n1/n2
    shieh_delta_corr <- shieh_delta*((nratio+1)/sqrt(nratio))       
    
    # Compute bias
    bias_cohen <- mean(file[,9]) - cohen_delta
    bias_hedge <- mean(file[,10]) - hedge_delta
    bias_glass1 <- mean(file[,11]) - glass_delta1 
    bias_glass2 <- mean(file[,12]) - glass_delta2
    bias_shieh <- mean(file[,13]) - shieh_delta
    bias_shieh_corr <- mean(file[,14]) - shieh_delta_corr
    
    # Compute Efficiency
    eff_cohen <- var(file[,9])
    eff_hedge <- var(file[,10])
    eff_glass1 <- var(file[,11])
    eff_glass2 <- var(file[,12])
    eff_shieh <- var(file[,13])
    eff_shieh_corr <- var(file[,14])
    
    # MSE
    MSE_cohen <- eff_cohen + bias_cohen^2
    MSE_hedge <- eff_hedge + bias_hedge^2
    MSE_glass1 <- eff_glass1 + bias_glass1^2
    MSE_glass2 <- eff_glass2 + bias_glass2^2
    MSE_shieh <- eff_shieh + bias_shieh^2
    MSE_shieh_corr <- eff_shieh_corr + bias_shieh_corr^2
    
    good_mes[j,1] <- n1
    good_mes[j,2] <- n2
    good_mes[j,3] <- n1/n2
    good_mes[j,4] <- m1-m2
    good_mes[j,5] <- sd1/sd2
    
    good_mes[j,6] <- bias_cohen
    good_mes[j,7] <- bias_hedge
    good_mes[j,8] <- bias_glass1
    good_mes[j,9] <- bias_glass2
    good_mes[j,10] <- bias_shieh
    good_mes[j,11] <- bias_shieh_corr
    
    good_mes[j,12] <- eff_cohen
    good_mes[j,13] <- eff_hedge
    good_mes[j,14] <- eff_glass1
    good_mes[j,15] <- eff_glass2
    good_mes[j,16] <- eff_shieh
    good_mes[j,17] <- eff_shieh_corr
    
    good_mes[j,18] <- MSE_cohen
    good_mes[j,19] <- MSE_hedge
    good_mes[j,20] <- MSE_glass1
    good_mes[j,21] <- MSE_glass2
    good_mes[j,22] <- MSE_shieh
    good_mes[j,23] <- MSE_shieh_corr
    
  }
  
  setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of measures/When m1=m2")
  
  sub<-paste0("G1=",G1,",G2=",G2) # Sign is not extracted but it's not a big deal
  write.table(good_mes,paste0(sub,",good_mes.txt"),sep=";",dec=",")
  
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------
#     GRAPH: PLOTTING bias AND EFFICIENCY
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


for (j in seq_len(length(list.files(Path)))){
  Path <-  "C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of measures/When m1=m2/"
  setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Graphs/When m1=m2")
  
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
    
    # Matrix containing variances information
    res2 <- matrix(0,9,8)  
    names<-expand.grid("var_",c("Cohen","Hedge","Glass1","Glass2","Shieh","Shieh_corr"))
    colnames(res2) <- c("n1","n2",paste0(names[,1],names[,2]))
    res2[,1:2] <- cbind(combi[,1],combi[,2]) 
    res2[,3] <- tapply(Sel$eff_Cohen,list(Sel$n1,Sel$n2),mean)[1:9]
    res2[,4] <- tapply(Sel$eff_Hedge,list(Sel$n1,Sel$n2),mean)[1:9]
    res2[,5] <- tapply(Sel$eff_Glass1,list(Sel$n1,Sel$n2),mean)[1:9]
    res2[,6] <- tapply(Sel$eff_Glass2,list(Sel$n1,Sel$n2),mean)[1:9]
    res2[,7] <- tapply(Sel$eff_Shieh,list(Sel$n1,Sel$n2),mean)[1:9]
    res2[,8] <- tapply(Sel$eff_Shieh_corr,list(Sel$n1,Sel$n2),mean)[1:9]
    # Select only rows with no "NA"  
    res2 <- subset(res2,res2[,3] != "NA") 
    # Select only bias columns
    res.eff <- t(res2[,3:8])
    colnames(res.eff) <- paste0(res2[,1],":",res2[,2])
    
    param <- str_extract_all(list.files(Path)[j], "[[:digit:]]+\\.*[[:digit:]]*")
    G1 <- as.numeric(param[[1]][2])
    G2 <- as.numeric(param[[1]][4])
    
    png(file=paste0("bias_eff,G1=",G1, " & G2=",G2,";",names(Conditions_id)[i], ".png"),width=2500,height=1700, units = "px", res = 300)  
    
    
    par(xpd = F,mar = c(4,5,0.5,0),mfrow = c(3,1))   
    
    # legend
    plot(1,1,bty="n",xaxt="n",yaxt="n",ylim=c(.62,.67),main="",xlab="",ylab="",pch=19,type="o")
    legend("center", 
           legend=c("Cohen's d","Hedge's g","Glass's delta (delta = sd1)","Glass's delta (delta = sd2)","Shieh's d","Corrected Shieh's d"),
           fill=c("black","grey25","grey50","grey70","grey90","white"),
           bty="n"
    )
    
    
    # plot for the bias
    barplot(res.bias, 
            col = c("black","grey25","grey50","grey70","grey90","white"),
            beside = TRUE,
            ylab = expression(paste("E(" , hat(delta) , ") -",delta )),
            xlab = "n1 : n2")
    
    # plot for the variance
    barplot(res.eff, 
            col = c("black","grey25","grey50","grey70","grey90","white"),
            beside = TRUE,
            ylab = expression(paste("Var(" , hat(delta) , ")")),
            xlab = "n1 : n2",
            args.legend = list(
              x = length(res.bias)*1.2,
              y = max(res.bias)+.5,
              bty="n"
            ))
    
    dev.off()
    
  }
  
}


#--------------------------------------------------------------------------------------------------------------------------------------------------------------------
#     GRAPH: PLOTTING MSE
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

for (j in seq_len(length(list.files(Path)))){
  Path <-  "C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of measures/When m1=m2/"
  setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Graphs/When m1=m2")

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
    
    # Matrix containing MSE information
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
