library(stringr)

setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs")
dir.create("Quality of ES measures")

Mainfolder="C:/Users/Marie/Documents/ES MEASURES/"
subfolder=list.files(Mainfolder)
Folder=paste0(Mainfolder,subfolder)

for (i in seq_len(length(Folder))){

  # set up empty container for all estimated parameters
  good_mes <-matrix(0,length(list.files(Folder[i])),5+3*6)
  goodness_indic <- c("bias_","eff_","MSE_")
  estimator <- c("Cohen","Hedge","Glass1","Glass2","Shieh", "Shieh_corr")
  columns_names=expand.grid(estimator,goodness_indic)
  colnames(good_mes) <- c("n1","n2","n1/n2","m1-m2","sd1/sd2",paste0(columns_names[,2],columns_names[,1]))

    for (j in seq_len(length(list.files(Folder[i])))){
    
    filepath = paste0(Folder[i],"/",list.files(Folder[i])[j])
    file=readRDS(filepath)
    param <- str_extract_all(list.files(Folder[i])[j], "[[:digit:]]+\\.*[[:digit:]]*")
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
    sigma_unbal <- sqrt((1-n1/(n1+n2))*sd1^2 + (1-n2/(n1+n2))*sd2^2)    
    sigma_bal <-   sqrt((sd1^2+sd2^2)/2)
    shieh_delta_corr <- shieh_delta*(((nratio+1)*sigma_unbal)/(2*sigma_bal*sqrt(nratio)))  
    
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
  
  setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures")
  
  shapeparam<-str_extract_all(Folder[i], "[[:digit:]]+\\.*[[:digit:]]*")
  sub<-paste0("G1=",shapeparam[[1]][2],",G2=",shapeparam[[1]][4]) # Sign is not extracted but it's not a big deal
  write.table(good_mes,paste0(sub,",good_mes.txt"),sep=";",dec=",")
  
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------
#     GRAPH: PLOTTING bias, EFFICIENCY and MSE
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


setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/")
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
    
    # Matrix containing MSE information
    res3 <- matrix(0,9,8)  
    names<-expand.grid("MSE_",c("Cohen","Hedge","Glass1","Glass2","Shieh","Shieh_corr"))
    colnames(res3) <- c("n1","n2",paste0(names[,1],names[,2]))
    res3[,1:2] <- cbind(combi[,1],combi[,2])
    res3[,3] <- tapply(Sel$MSE_Cohen,list(Sel$n1,Sel$n2),mean)[1:9]
    res3[,4] <- tapply(Sel$MSE_Hedge,list(Sel$n1,Sel$n2),mean)[1:9]
    res3[,5] <- tapply(Sel$MSE_Glass1,list(Sel$n1,Sel$n2),mean)[1:9]
    res3[,6] <- tapply(Sel$MSE_Glass2,list(Sel$n1,Sel$n2),mean)[1:9]
    res3[,7] <- tapply(Sel$MSE_Shieh,list(Sel$n1,Sel$n2),mean)[1:9]
    res3[,8] <- tapply(Sel$MSE_Shieh_corr,list(Sel$n1,Sel$n2),mean)[1:9]
    # Select only rows with no "NA"  
    res3 <- subset(res3,res3[,3] != "NA") 
    # Select only bias columns
    res.MSE <- t(res3[,3:8])
    colnames(res.MSE) <- paste0(res3[,1],":",res3[,2])
    
    param <- str_extract_all(list.files(Path,pattern=".*.txt")[j], "[[:digit:]]+\\.*[[:digit:]]*")
    G1 <- as.numeric(param[[1]][2])
    G2 <- as.numeric(param[[1]][4])
    
     setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs")
    #dir.create(names(Conditions_id)[i])
    
    setwd(paste0("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/",names(Conditions_id)[i]))
    png(file=paste0("bias_eff,G1=",G1, " & G2=",G2,";",names(Conditions_id)[i], ".png"),width=1000,height=1700, units = "px", res = 300)  
    
    par(xpd = F,mar = c(4,5,1.5,0),mfrow = c(3,1))   
    

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


    
    # plot for the variance
    
    if (j==1){ylabeleff=expression(paste("Var(" , hat(delta) , ")"))
    } else {ylabeleff=""}
    
    barplot(res.eff, 
            col = c("black","grey25","grey50","grey70","grey90","white"),
            beside = TRUE,
            ylab = ylabeleff,
            xaxt="n",
            cex.lab=1.5,
            args.legend = list(
              x = length(res.eff)*1.2,
              y = max(res.eff)+.5,
              bty="n"
            ))
    
    # plot for the MSE
    
    if (j==1){ylabelMSE=expression(paste( "MSE(",hat(delta), ")"))
    } else {ylabelMSE=""}
    
    barplot(res.MSE, 
            col = c("black","grey25","grey50","grey70","grey90","white"),
            beside = TRUE,
            ylab = ylabelMSE,
            xlab = "n1 : n2",
            cex.lab=1.5,
            args.legend = list(
              x = length(res.MSE)*1.2,
              y = max(res.MSE)+.5,
              bty="n"
            ))
    
    
    dev.off()
    
  }
  
}

