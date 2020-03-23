library(stringr)

setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs")
dir.create("Confidence intervals")

Mainfolder="C:/Users/Marie/Documents/CI measures/"
subfolder=list.files(Mainfolder)

for (k in seq_len(length(subfolder))){
  subfolder2=list.files(paste0(Mainfolder,subfolder[k],"/"))  
  Folder=paste0(Mainfolder,subfolder[k],subfolder2)

  
  list.files(Mainfolder)
  
for (i in seq_len(length(Folder))){

  # set up empty container for all estimated parameters
  good_mes <-matrix(0,length(list.files(Folder[i])),6+5)
       # 5 columns dedicated to the simulated coverage probability
       # 6 columns for the simulation design
  
  colnames(good_mes) <- c("Confidence.level","n1","n2","m1-m2","sd1/sd2",
                          "NRH0_student","Cohen_SCP","Hedges_SCP",
                          "NRH0_welch","Shieh_SCP","transfoshieh_SCP")

    for (j in seq_len(length(list.files(Folder[i])))){
    
    filepath = paste0(Folder[i],"/",list.files(Folder[i])[j])
    file=readRDS(filepath) 
    param <- str_extract_all(list.files(Folder[i])[j], "[[:digit:]]+\\.*[[:digit:]]*")
    conf.level <- param[[1]][1]
    n1 <- as.numeric(param[[1]][6])
    n2 <- as.numeric(param[[1]][7])
    m1 <- as.numeric(param[[1]][8])
    m2 <- as.numeric(param[[1]][9])
    sd1 <- as.numeric(param[[1]][10])
    sd2 <- as.numeric(param[[1]][11])
    
    cohen_delta <- (m1-m2)/sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
    Hedges_delta <- cohen_delta* (1-3/(4*(n1+n2-9)))
    q1 <- n1/(n1+n2)
    q2 <- n2/(n1+n2)
    shieh_delta <- (m1-m2)/sqrt(sd1^2/q1+sd2^2/q2)       
    nratio=n1/n2
    sigma_unbal <- sqrt((1-n1/(n1+n2))*sd1^2 + (1-n2/(n1+n2))*sd2^2)    
    sigma_bal <-   sqrt((sd1^2+sd2^2)/2)
    shieh_delta_corr <- shieh_delta*(((nratio+1)*sigma_unbal)/(2*sigma_bal*sqrt(nratio)))  
    
    # Compute the proportion of tests conducting to not rejecting the null hypothesis
    # In order to check consistency with confidence limits
    Student_NRH0 <- sum(file[,1]>= (1-conf.level))/length(A[,1])
    Welch_NRH0 <- sum(file[,6]>= (1-conf.level))/length(A[,6])
 
    # Compute the Simulated coverage probability
    Cohen_SCP <- sum(file[,4]<= cohen_delta & file[,5]>=cohen_delta)/length(file[,4]) # Simulated coverage probability of theoretical Hedge's gs value when computing Cohen's limits
    Hedges_SCP <-  sum(file[,4]<= Hedges_delta & file[,5]>=Hedges_delta)/length(file[,4]) # Simulated coverage probability of theoretical Cohen's ds value when computing Cohen's limits
    Shieh_SCP <-  sum(file[,8]<= shieh_delta & file[,9]>=shieh_delta)/length(file[,8]) # Simulated coverage probability of theoretical Shieh's ds value when computing Shieh's limits
    transfoshieh_SCP <- sum(file[,11]<= shieh_delta_corr & file[,12]>=shieh_delta_corr)/length(file[,11]) # Simulated coverage probability of theoretical transformed Shieh's ds value when computing transformed Shieh's limits

    colnames(good_mes) <- c("Confidence.level","n1","n2","m1-m2","sd1/sd2",
                            "NRH0_student","Cohen_SCP","Hedges_SCP",
                            "NRH0_welch","Shieh_SCP","transfoshieh_SCP")
    
    
    good_mes[j,1] <- conf.level
    good_mes[j,2] <- n1
    good_mes[j,3] <- n2
    good_mes[j,4] <- m1-m2
    good_mes[j,5] <- sd1/sd2
    
    good_mes[j,6] <- Student_NRH0
    good_mes[j,7] <- Cohen_SCP
    good_mes[j,8] <- Hedges_SCP 
    good_mes[j,9] <- Welch_NRH0
    good_mes[j,10] <- Shieh_SCP 
    good_mes[j,11] <- transfoshieh_SCP

  }
  
  setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Confidence intervals")
  shapeparam<-str_extract_all(Folder[i], "[[:digit:]]+\\.*[[:digit:]]*")
  sub<-paste0(subfolder[k],"G1=",shapeparam[[1]][3],",G2=",shapeparam[[1]][5]) # Sign is not extracted but it's not a big deal
  write.table(good_mes,paste0(sub,",good_mes.txt"),sep=";",dec=",")
  
}
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


setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Graphs")

for (j in seq_len(length(list.files(Path)))){
  Path <-  "C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of measures/"

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


