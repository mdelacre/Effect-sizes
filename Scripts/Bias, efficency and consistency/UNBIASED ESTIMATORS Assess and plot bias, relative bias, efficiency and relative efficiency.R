library(stringr)

#setwd("C:/Users/Admin/Documents/Github_projects/Effect-sizes/Scripts outputs")
#dir.create("Quality of ES measures")

#Mainfolder="C:/Users/Admin/Documents/ES MEASURES/"
#subfolder=list.files(Mainfolder)
#Folder=paste0(Mainfolder,subfolder)

#for (i in seq_len(length(Folder))){

  # set up empty container for all estimated parameters
#  good_mes <-matrix(0,length(list.files(Folder[i])),5+4*5)
#  goodness_indic <- c("bias_","relbias_","eff_","releff_")
#  estimator <- c("Hedges","Glass1","Glass2","Shieh", "cohen_delta_prime")
#  columns_names=expand.grid(estimator,goodness_indic)
#  colnames(good_mes) <- c("n1","n2","n1/n2","m1-m2","sd1/sd2",paste0(columns_names[,2],columns_names[,1]))

  #    for (j in seq_len(length(list.files(Folder[i])))){
    
      #    filepath = paste0(Folder[i],"/",list.files(Folder[i])[j])
    #    file=readRDS(filepath) 
    
    #    param <- str_extract_all(list.files(Folder[i])[j], "[[:digit:]]+\\.*[[:digit:]]*")
    #    n1 <- as.numeric(param[[1]][5])
    #    n2 <- as.numeric(param[[1]][6])
    #    m1 <- as.numeric(param[[1]][7])
    #    m2 <- as.numeric(param[[1]][8])
    #    sd1 <- as.numeric(param[[1]][9])
    #    sd2 <- as.numeric(param[[1]][10])
    
    #    cohen_delta <- (m1-m2)/sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
    #    glass_delta1 <- (m1-m2)/sd1
    #    glass_delta2 <- (m1-m2)/sd2
    #    q1 <- n1/(n1+n2)
    #    q2 <- n2/(n1+n2)
    #    shieh_delta <- (m1-m2)/sqrt(sd1^2/q1+sd2^2/q2)       
    #    nratio=n1/n2
    #    cohen_delta_prime <- (m1-m2)/sqrt((sd1^2+sd2^2)/2)  
    
# Compute bias
#bias_cohen <- mean(file[,9]) - cohen_delta
#bias_glass1 <- mean(file[,11]) - glass_delta1 
#bias_glass2 <- mean(file[,12]) - glass_delta2
#bias_shieh <- mean(file[,15]) - shieh_delta
#bias_cohen_delta_prime <- mean(file[,17]) - cohen_delta_prime

# Compute relative bias
#relbias_cohen <- (mean(file[,9]) - cohen_delta)/cohen_delta
#relbias_glass1 <- (mean(file[,11]) - glass_delta1)/glass_delta1 
#relbias_glass2 <- (mean(file[,12]) - glass_delta2)/glass_delta2
#relbias_shieh <- (mean(file[,15]) - shieh_delta)/shieh_delta
#relbias_cohen_delta_prime <- (mean(file[,17]) - cohen_delta_prime)/cohen_delta_prime

# Compute variance
#eff_cohen <- var(file[,9])
#eff_glass1 <- var(file[,11])
#eff_glass2 <- var(file[,12])
#eff_shieh <- var(file[,15])
#eff_cohen_delta_prime <- var(file[,17])

# Compute relative variance
#releff_cohen <- var(file[,9])/cohen_delta^2
#releff_glass1 <- var(file[,11])/glass_delta1^2 
#releff_glass2 <- var(file[,12])/glass_delta2^2
#releff_shieh <- var(file[,15])/shieh_delta^2
#releff_cohen_delta_prime <- var(file[,17])/cohen_delta_prime^2
    #    good_mes[j,1] <- n1
    #    good_mes[j,2] <- n2
    #    good_mes[j,3] <- n1/n2
    #    good_mes[j,4] <- m1-m2
    #    good_mes[j,5] <- sd1/sd2
    
    #    good_mes[j,6] <- bias_cohen
    #    good_mes[j,7] <- bias_glass1
    #    good_mes[j,8] <- bias_glass2
    #    good_mes[j,9] <- bias_shieh
    #    good_mes[j,10] <- bias_cohen_delta_prime

    #   good_mes[j,11] <- relbias_cohen
    #    good_mes[j,12] <- relbias_glass1
    #    good_mes[j,13] <- relbias_glass2
    #    good_mes[j,14] <- relbias_shieh
    #    good_mes[j,15] <- relbias_cohen_delta_prime
    
    #    good_mes[j,16] <- eff_cohen
    #    good_mes[j,17] <- eff_glass1
    #    good_mes[j,18] <- eff_glass2
    #    good_mes[j,19] <- eff_shieh
    #    good_mes[j,20] <- eff_cohen_delta_prime
    
    #    good_mes[j,21] <- releff_cohen
    #    good_mes[j,22] <- releff_glass1
    #    good_mes[j,23] <- releff_glass2
    #    good_mes[j,24] <- releff_shieh
    #    good_mes[j,25] <- releff_cohen_delta_prime
    
    #  }
  
  #  setwd("C:/Users/Admin/Documents/Github projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators")

  #  shapeparam<-str_extract_all(Folder[i], "[[:digit:]]+\\.*[[:digit:]]*")
  #  sub<-paste0("G1=",shapeparam[[1]][2],",G2=",shapeparam[[1]][4]) # Sign is not extracted but it's not a big deal
  #  write.table(good_mes,paste0(sub,"_unbiasedest_properties.txt"),sep=";",dec=",")
  
#}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------
#     GRAPH: PLOTTING raw bias and variances 
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

#     Conditions a, b and c, as a function of the n-ratio
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

setwd("C:/Users/Admin/Documents/Github projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators")

# png(file="legend.png",width=1500,height=1000, units = "px", res = 300)  
# plot(1,1,bty="n",xaxt="n",yaxt="n",ylim=c(.62,.67),main="",xlab="",ylab="",pch=19,type="o")
# legend("center", 
#       legend=c(expression(paste("Hedges' ",g[s])),expression(paste("Glass's ",g[s],"(",sigma," =",S[1],")")),expression(paste("Glass's ",g[s],"(",sigma," =",S[2],")")),expression(paste("Shieh's ",g[s])),
#                expression(paste("Hedges' ",g[s],"'"))),
#       fill=c("black","grey40","grey60","grey80","white"),
#       bty="n"
# )
# dev.off()
getwd()
png(file="legend.png",width=8000,height=800, units = "px", res = 300)  

plot(1,1,bty="n",xaxt="n",yaxt="n",ylim=c(.62,.67),main="",xlab="",ylab="",pch=19,type="o")
legend("center", 
       legend=c(expression(paste("Hedges' ",g)),expression(paste("Glass's ",g,"(",sigma," =",S[1],")")),expression(paste("Glass's ",g,"(",sigma," =",S[2],")")),expression(paste("Shieh's ",g)),
                expression(paste("Hedges' ",g,"*"))),
       fill=c("black","grey40","grey60","grey80","white"),
       bty="n",horiz=T,cex=2
)

dev.off()

Path <-  "C:/Users/Admin/Documents/Github projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators"

for (j in seq_len(length(list.files(Path)))){
  

  File=read.table(paste0(Path,list.files(Path)[j]),header=T,sep=";",dec=",")  
  
  # Extract three subcategories from File: "Hom_bal", "Hom_unbal", "Het_bal"
  
  ### Conditions id 
  
  # Homoscedasticity, balanced designs ("Hom_bal")
  id_Hom_bal=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1==File$n2)&(File$sd1.sd2 == 1),]))
  # Homoscedasticity, unbalanced designs ("Hom_unbal")  
  id_Hom_nN=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 == 1),]))
  id_Hom_Nn=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 == 1),]))
  # Heteroscedasticity, balanced designs ("Het_bal")
  id_sdSD_bal=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1==File$n2)&(File$sd1.sd2 < 1),]))
  id_SDsd_bal=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1==File$n2)&(File$sd1.sd2 > 1),]))
  
  # Study each condition separately
  Conditions_id <- list(id_Hom_bal=id_Hom_bal,
                        id_Hom_unbal=c(id_Hom_nN,id_Hom_Nn),
                        id_Het_bal=c(id_sdSD_bal,id_SDsd_bal))
  
  for (i in seq_len(length(Conditions_id))){ 
    
    Sel <- File[Conditions_id[[i]],]
    
    n2val <- as.numeric(levels(factor(Sel$n2)))
    n1val<- as.numeric(levels(factor(Sel$n1)))
    combi <- expand.grid(n1val, n2val)
    
    # Matrix containing biases information
    res <- matrix(0,9,7)  
    names<-expand.grid("bias_",c("Hedges","Glass1","Glass2","Shieh","cohen_delta_prime"))
    colnames(res) <- c("n1","n2",paste0(names[,1],names[,2]))
    res[,1:2] <- cbind(combi[,1],combi[,2])
    res[,3] <- tapply(Sel$bias_Hedges,list(Sel$n1,Sel$n2),mean)[1:9]
    res[,4] <- tapply(Sel$bias_Glass1,list(Sel$n1,Sel$n2),mean)[1:9]
    res[,5] <- tapply(Sel$bias_Glass2,list(Sel$n1,Sel$n2),mean)[1:9]
    res[,6] <- tapply(Sel$bias_Shieh,list(Sel$n1,Sel$n2),mean)[1:9]
    res[,7] <- tapply(Sel$bias_cohen_delta_prime,list(Sel$n1,Sel$n2),mean)[1:9]
    # Select only rows with no "NA"  
    res <- subset(res,res[,3] != "NA") 
    # Select only bias columns
    res.bias <- t(res[,3:7])
    colnames(res.bias) <- paste0(res[,1],":",res[,2])
    
    param <- str_extract_all(list.files(Path)[j], "[[:digit:]]+\\.*[[:digit:]]*")
    if(param[[1]][2]==2.08){
      G1 <- -as.numeric(param[[1]][2])
    } else (G1 <- as.numeric(param[[1]][2]))
    G2 <- as.numeric(param[[1]][4])
    
    # Matrix containing variances information
    res3 <- matrix(0,9,7)  
    names<-expand.grid("var_",c("Hedges","Glass1","Glass2","Shieh","cohen_delta_prime"))
    colnames(res3) <- c("n1","n2",paste0(names[,1],names[,2]))
    res3[,1:2] <- cbind(combi[,1],combi[,2]) 
    res3[,3] <- tapply(Sel$eff_Hedges,list(Sel$n1,Sel$n2),mean)[1:9]
    res3[,4] <- tapply(Sel$eff_Glass1,list(Sel$n1,Sel$n2),mean)[1:9]
    res3[,5] <- tapply(Sel$eff_Glass2,list(Sel$n1,Sel$n2),mean)[1:9]
    res3[,6] <- tapply(Sel$eff_Shieh,list(Sel$n1,Sel$n2),mean)[1:9]
    res3[,7] <- tapply(Sel$eff_cohen_delta_prime,list(Sel$n1,Sel$n2),mean)[1:9]
    # Select only rows with no "NA"  
    res3 <- subset(res3,res3[,3] != "NA") 
    # Select only bias columns
    res.eff <- t(res3[,3:7])
    colnames(res.eff) <- paste0(res3[,1],":",res3[,2])
    
    setwd("C:/Users/Admin/Documents/Github projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/Raw estimators of goodness/")
    #dir.create(names(Conditions_id)[i])

    setwd(paste0("C:/Users/Admin/Documents/Github projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/Raw estimators of goodness/",names(Conditions_id)[i]))
    if(G1==-2.08){
      g1=2.08
    } else {g1=G1}
    png(file=paste0("bias_eff,G1=",g1, " & G2=",G2,";",names(Conditions_id)[i], ".png"),width=1400,height=1700, units = "px", res = 300)  
    
    # plot for the bias
    
    par(mar = c(4,5,1.5,0),mfrow = c(2,1))   
    
    
    if (j==1){ylabelbias=expression(paste("E(" , hat(delta) , ") -",delta ))
    } else {ylabelbias=""}
    
    barplot(res.bias, 
            col = c("black","grey40","grey60","grey80","white"),
            beside = TRUE,
            main=paste0("G1=",G1,"; G2=",G2),
            xaxt="n",
            cex.lab=1.5,
            cex.main=1.5,
            ylab=ylabelbias
    )
    
    # plot for the variance (need to be modified later)
    
    if (j==1){ylabeleff=expression(paste("Var(" , hat(delta) , ")"))
    } else {ylabeleff=""}
    
    barplot(res.eff, 
            col = c("black","grey40","grey60","grey80","white"),
            beside = TRUE,
            ylab = ylabeleff,
            cex.lab=1.5,
            xaxt="n",
            args.legend = list(
              x = length(res.eff)*1.2,
              y = max(res.eff)+.5,
              bty="n"
            ))
    

    dev.off()
    
  }
  
}

#     Condition c, as a function of the SD-ratio
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Because the comparison pattern is very similar whatever n = 20, 50 or 100
# We will only plot the results as a function of the SD-ratio when n = 20

Path <-  "C:/Users/Admin/Documents/Github projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/"

for (j in seq_len(length(list.files(Path)))){
  
  
  File=read.table(paste0(Path,list.files(Path)[j]),header=T,sep=";",dec=",")  

  # Extract "Het_bal" subcategory (when n1=n2=20)
  
  id_sdSD_bal=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1==File$n2)&(File$n1==20)&(File$sd1.sd2 < 1),]))
  id_SDsd_bal=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1==File$n2)&(File$n1==20)&(File$sd1.sd2 > 1),]))
  id_Het_bal <- c(id_sdSD_bal,id_SDsd_bal)
  
  Sel <- File[id_Het_bal,]
  
  # Possible SD-ratios
  sdval <- as.numeric(levels(factor(Sel$sd1.sd2)))
  
  # Matrix containing biases information
  res <- matrix(0,6,6)  
  names<-expand.grid("bias_",c("Hedges","Glass1","Glass2","Shieh","cohen_delta_prime"))
  colnames(res) <- c("sd-ratio",paste0(names[,1],names[,2]))
  res[,1] <- sdval
  res[,2] <- tapply(Sel$bias_Hedges,list(Sel$sd1.sd2),mean)[1:6]
  res[,3] <- tapply(Sel$bias_Glass1,list(Sel$sd1.sd2),mean)[1:6]
  res[,4] <- tapply(Sel$bias_Glass2,list(Sel$sd1.sd2),mean)[1:6]
  res[,5] <- tapply(Sel$bias_Shieh,list(Sel$sd1.sd2),mean)[1:6]
  res[,6] <- tapply(Sel$bias_cohen_delta_prime,list(Sel$sd1.sd2),mean)[1:6]
  # Select only bias columns
  res.bias <- t(res[,2:6])
  colnames(res.bias) <- paste0(1,":",1/sdval)
  
  param <- str_extract_all(list.files(Path)[j], "[[:digit:]]+\\.*[[:digit:]]*")
  if(param[[1]][2]==2.08){
    G1 <- -as.numeric(param[[1]][2])
  } else (G1 <- as.numeric(param[[1]][2]))
  G2 <- as.numeric(param[[1]][4])
  
  # Matrix containing variances information
  res3 <- matrix(0,6,6)  
  names<-expand.grid("var_",c("Hedges","Glass1","Glass2","Shieh","cohen_delta_prime"))
  colnames(res3) <- c("sd-ratio",paste0(names[,1],names[,2]))
  res3[,1] <- sdval 
  res3[,2] <- tapply(Sel$eff_Hedges,list(Sel$sd1.sd2),mean)[1:6]
  res3[,3] <- tapply(Sel$eff_Glass1,list(Sel$sd1.sd2),mean)[1:6]
  res3[,4] <- tapply(Sel$eff_Glass2,list(Sel$sd1.sd2),mean)[1:6]
  res3[,5] <- tapply(Sel$eff_Shieh,list(Sel$sd1.sd2),mean)[1:6]
  res3[,6] <- tapply(Sel$eff_cohen_delta_prime,list(Sel$sd1.sd2),mean)[1:6]
  # Select only bias columns
  res.eff <- t(res3[,2:6])
  colnames(res.eff) <- paste0(1,":",1/sdval)
  
  setwd(paste0("C:/Users/Admin/Documents/Github projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/Raw estimators of goodness/id_Het_bal/sd-ratio"))
  if(G1==-2.08){
    g1=2.08
  } else {g1=G1}
  png(file=paste0("bias_eff,G1=",g1, " & G2=",G2,";",names(Conditions_id)[i], ".png"),width=1400,height=1700, units = "px", res = 300)  
  
  par(mar = c(4,5,1.5,0),mfrow = c(2,1))   
  
  # plot for the bias
  
  if (j==1){ylabelbias=expression(paste("E(" , hat(delta) , ") -",delta ))
  } else {ylabelbias=""}
  
  barplot(res.bias, 
          col = c("black","grey40","grey60","grey80","white"),
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
          col = c("black","grey40","grey60","grey80","white"),
          beside = TRUE,
          ylab = ylabeleff,
          cex.lab=1.5,
          args.legend = list(
            x = length(res.eff)*1.2,
            y = max(res.eff)+.5,
            bty="n"
          ))
  
  
  dev.off()
  
}

#     Condition d and e
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

# If sdSD & nN or SDsd & Nn: positive correlation between n and sd (rpos; condition d)
# If sdSD & Nn or SDsd & nN: negative correlation between n and sd (rneg; condition e)

# Note: There are too many conditions to represent them all in a single Figure. 
# We will therefore generate 3 figures per condition, as a function of the Total sample size (n1+n2)

Path <-  "C:/Users/Admin/Documents/Github projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/"

plot_hetr <- function(totalN){for (j in seq_len(length(list.files(Path)))){
  
  File=read.table(paste0(Path,list.files(Path)[j]),header=T,sep=";",dec=",")  
  
  # Extract "Het_bal" subcategory (for a specific total sample size)
  
  if (totalN == 70){
    
    # Heteroscedasticity, when the first group is the smaller one
    id_sdSD_nN=as.numeric(rownames(File[(File$n1+File$n2==70) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 < 1),]))
    id_SDsd_nN=as.numeric(rownames(File[(File$n1+File$n2==70) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 > 1),]))
    # Heteroscedasticity, when the second group is the larger one
    id_sdSD_Nn=as.numeric(rownames(File[(File$n1+File$n2==70) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 < 1),]))
    id_SDsd_Nn=as.numeric(rownames(File[(File$n1+File$n2==70) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 > 1),]))
    
    Conditions_id <- list(id_Het_firstsmaller=c(id_sdSD_nN,id_SDsd_nN),id_Het_firstlarger=c(id_sdSD_Nn,id_SDsd_Nn))
    
  } else if (totalN == 120){
    
    # Heteroscedasticity, positive correlation between n and sd ("Het_rpos")
    id_sdSD_nN=as.numeric(rownames(File[(File$n1+File$n2==120) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 < 1),]))
    id_SDsd_nN=as.numeric(rownames(File[(File$n1+File$n2==120) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 > 1),]))
    # Heteroscedasticity, negative correlation between n and sd ("Het_rneg")
    id_sdSD_Nn=as.numeric(rownames(File[(File$n1+File$n2==120) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 < 1),]))
    id_SDsd_Nn=as.numeric(rownames(File[(File$n1+File$n2==120) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 > 1),]))
    
    Conditions_id <- list(id_Het_firstsmaller=c(id_sdSD_nN,id_SDsd_nN),id_Het_firstlarger=c(id_sdSD_Nn,id_SDsd_Nn))
    
  } else if (totalN == 150){
    
    # Heteroscedasticity, positive correlation between n and sd ("Het_rpos")
    id_sdSD_nN=as.numeric(rownames(File[(File$n1+File$n2==150) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 < 1),]))
    id_SDsd_nN=as.numeric(rownames(File[(File$n1+File$n2==150) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 > 1),]))
    # Heteroscedasticity, negative correlation between n and sd ("Het_rneg")
    id_sdSD_Nn=as.numeric(rownames(File[(File$n1+File$n2==150) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 < 1),]))
    id_SDsd_Nn=as.numeric(rownames(File[(File$n1+File$n2==150) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 > 1),]))
    
    Conditions_id <- list(id_Het_firstsmaller=c(id_sdSD_nN,id_SDsd_nN),id_Het_firstlarger=c(id_sdSD_Nn,id_SDsd_Nn))
    
  }  
  
  for (i in seq_len(length(Conditions_id))){ 
    
  Sel <- File[Conditions_id[[i]],]

  # Possible SD-ratios and nratios
  sdval <- as.numeric(levels(factor(Sel$sd1.sd2)))
  n2val <- as.numeric(levels(factor(Sel$n2)))
  n1val<- as.numeric(levels(factor(Sel$n1)))
  combi <- expand.grid(sdval, n1val, n2val)
  
  # Matrix containing biases information
  res <- matrix(0,6,8)  
  names<-expand.grid("bias_",c("Hedges","Glass1","Glass2","Shieh","cohen_delta_prime"))
  colnames(res) <- c("n1","n2","sd-ratio",paste0(names[,1],names[,2]))
  res[,1:3] <- cbind(combi[,2],combi[,3],combi[,1])
  res[,4] <- tapply(Sel$bias_Hedges,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:6]
  res[,5] <- tapply(Sel$bias_Glass1,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:6]
  res[,6] <- tapply(Sel$bias_Glass2,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:6]
  res[,7] <- tapply(Sel$bias_Shieh,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:6]
  res[,8] <- tapply(Sel$bias_cohen_delta_prime,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:6]
  
  # Select only bias columns
  res.bias <- t(res[,4:8])
  colnames(res.bias) <- paste0(res[,1],":",res[,2],"\n 1:",1/res[,3])
  
  param <- str_extract_all(list.files(Path)[j], "[[:digit:]]+\\.*[[:digit:]]*")
  if(param[[1]][2]==2.08){
    G1 <- -as.numeric(param[[1]][2])
  } else (G1 <- as.numeric(param[[1]][2]))
  G2 <- as.numeric(param[[1]][4])
  
  # Matrix containing variances information
  res3 <- matrix(0,6,8)  
  names<-expand.grid("bias_",c("Hedges","Glass1","Glass2","Shieh","cohen_delta_prime"))
  colnames(res3) <- c("n1","n2","sd-ratio",paste0(names[,1],names[,2]))
  res3[,1:3] <- cbind(combi[,2],combi[,3],combi[,1])
  res3[,4] <- tapply(Sel$eff_Hedges,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:6]
  res3[,5] <- tapply(Sel$eff_Glass1,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:6]
  res3[,6] <- tapply(Sel$eff_Glass2,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:6]
  res3[,7] <- tapply(Sel$eff_Shieh,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:6]
  res3[,8] <- tapply(Sel$eff_cohen_delta_prime,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:6]
  
  # Select only variances columns
  res.eff <- t(res3[,4:8])
  colnames(res.eff) <- paste0(res3[,1],":",res3[,2],"\n (1:",1/res3[,3],")")

  setwd(paste0("C:/Users/Admin/Documents/Github projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/Raw estimators of goodness/",names(Conditions_id)[i]))
  if(G1==-2.08){
    g1=2.08
  } else {g1=G1}
  png(file=paste0("G1=",g1, " & G2=",G2,";",names(Conditions_id)[i],"; N=",totalN,".png"),width=1400,height=1700, units = "px", res = 300)  
  
  par(mar = c(4,5,1.5,0),mfrow = c(2,1))   
  
  # plot for the bias
  
  if (j==1){ylabelbias=expression(paste("E(" , hat(delta) , ") -",delta ))
  } else {ylabelbias=""}
  
  barplot(res.bias, 
          col = c("black","grey40","grey60","grey80","white"),
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
          col = c("black","grey40","grey60","grey80","white"),
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
}

plot_hetr(totalN=70)
plot_hetr(totalN=120)
plot_hetr(totalN=150)



#--------------------------------------------------------------------------------------------------------------------------------------------------------------------
#     GRAPH: PLOTTING relative bias and variance 
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

#     Conditions a, b and c, as a function of the n-ratio
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

setwd("C:/Users/Admin/Documents/Github projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/")

Path <-  "C:/Users/Admin/Documents/Github projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/"

for (j in seq_len(length(list.files(Path)))){
  
  File=read.table(paste0(Path,list.files(Path)[j]),header=T,sep=";",dec=",")  
  
  # Extract three subcategories from File: "Hom_bal", "Hom_unbal", "Het_bal"
  
  ### Conditions id 
  
  # Homoscedasticity, balanced designs ("Hom_bal")
  id_Hom_bal=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1==File$n2)&(File$sd1.sd2 == 1),]))
  # Homoscedasticity, unbalanced designs ("Hom_unbal")  
  id_Hom_nN=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 == 1),]))
  id_Hom_Nn=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 == 1),]))
  # Heteroscedasticity, balanced designs ("Het_bal")
  id_sdSD_bal=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1==File$n2)&(File$sd1.sd2 < 1),]))
  id_SDsd_bal=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1==File$n2)&(File$sd1.sd2 > 1),]))
  
  # Study each condition separately
  Conditions_id <- list(id_Hom_bal=id_Hom_bal,
                        id_Hom_unbal=c(id_Hom_nN,id_Hom_Nn),
                        id_Het_bal=c(id_sdSD_bal,id_SDsd_bal))
  
  for (i in seq_len(length(Conditions_id))){ 
    
    Sel <- File[Conditions_id[[i]],]
    
    n2val <- as.numeric(levels(factor(Sel$n2)))
    n1val<- as.numeric(levels(factor(Sel$n1)))
    combi <- expand.grid(n1val, n2val)
    
    # Matrix containing biases information
    res2 <- matrix(0,9,7)  
    names<-expand.grid("bias_",c("Hedges","Glass1","Glass2","Shieh","cohen_delta_prime"))
    colnames(res2) <- c("n1","n2",paste0(names[,1],names[,2]))
    res2[,1:2] <- cbind(combi[,1],combi[,2])
    res2[,3] <- tapply(Sel$relbias_Hedges,list(Sel$n1,Sel$n2),mean)[1:9]
    res2[,4] <- tapply(Sel$relbias_Glass1,list(Sel$n1,Sel$n2),mean)[1:9]
    res2[,5] <- tapply(Sel$relbias_Glass2,list(Sel$n1,Sel$n2),mean)[1:9]
    res2[,6] <- tapply(Sel$relbias_Shieh,list(Sel$n1,Sel$n2),mean)[1:9]
    res2[,7] <- tapply(Sel$relbias_cohen_delta_prime,list(Sel$n1,Sel$n2),mean)[1:9]
    # Select only rows with no "NA"  
    res2 <- subset(res2,res2[,3] != "NA") 
    # Select only bias columns
    res.relbias <- t(res2[,3:7])
    colnames(res.relbias) <- paste0(res2[,1],":",res2[,2])
    
    param <- str_extract_all(list.files(Path)[j], "[[:digit:]]+\\.*[[:digit:]]*")
    if(param[[1]][2]==2.08){
      G1 <- -as.numeric(param[[1]][2])
    } else (G1 <- as.numeric(param[[1]][2]))
    G2 <- as.numeric(param[[1]][4])
    
    # Matrix containing variances information
    res4 <- matrix(0,9,7)  
    names<-expand.grid("var_",c("Hedges","Glass1","Glass2","Shieh","cohen_delta_prime"))
    colnames(res4) <- c("n1","n2",paste0(names[,1],names[,2]))
    res4[,1:2] <- cbind(combi[,1],combi[,2]) 
    res4[,3] <- tapply(Sel$releff_Hedges,list(Sel$n1,Sel$n2),mean)[1:9]
    res4[,4] <- tapply(Sel$releff_Glass1,list(Sel$n1,Sel$n2),mean)[1:9]
    res4[,5] <- tapply(Sel$releff_Glass2,list(Sel$n1,Sel$n2),mean)[1:9]
    res4[,6] <- tapply(Sel$releff_Shieh,list(Sel$n1,Sel$n2),mean)[1:9]
    res4[,7] <- tapply(Sel$releff_cohen_delta_prime,list(Sel$n1,Sel$n2),mean)[1:9]
    # Select only rows with no "NA"  
    res4 <- subset(res4,res4[,3] != "NA") 
    # Select only bias columns
    res.releff <- t(res4[,3:7])
    colnames(res.releff) <- paste0(res4[,1],":",res4[,2])
    
    setwd("C:/Users/Admin/Documents/Github projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/Relative estimators of goodness/")
    #dir.create(names(Conditions_id)[i])
    
    setwd(paste0("C:/Users/Admin/Documents/Github projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/Relative estimators of goodness/",names(Conditions_id)[i]))
    if(G1==-2.08){
      g1=2.08
    } else {g1=G1}
    png(file=paste0("bias_eff,G1=",g1, " & G2=",G2,";",names(Conditions_id)[i], ".png"),width=1400,height=1700, units = "px", res = 300)  
    
    par(mar = c(4,5,1.5,0),mfrow = c(2,1))   
    
    # plot for the relative bias
    
    if (j==1){ylabelbias=expression(paste("[E(" , hat(delta) , ")-",delta,"]/",delta))
    } else {ylabelbias=""}
    
    if (j==4){ylim=c(0,1.6)
    } else {ylim=c(0,.12)}
    
    barplot(res.relbias, 
            col = c("black","grey40","grey60","grey80","white"),
            beside = TRUE,
            xaxt="n",
            cex.lab=1.2,
            las=1,
            main=paste0("G1=",G1,"; G2=",G2),
            cex.main=1.5,
            ylim=ylim,
            ylab=ylabelbias
    )

    # plot for the relative variance
    
    if (j==1){ylabeleff=expression(paste("Var(" , hat(delta) , ")/",delta^2))
    } else {ylabeleff=""}
    
    if (j==4){ylim=c(0,4)
    } else {ylim=c(0,.6)}
    
    barplot(res.releff, 
            col = c("black","grey40","grey60","grey80","white"),
            beside = TRUE,
            ylab = ylabeleff,
            las=1,
            ylim=ylim,
            cex.lab=1.5,
            cex.names=.8,
            args.legend = list(
              x = length(res.releff)*1.2,
              y = max(res.releff)+.5,
              bty="n"
            ))
    
    par(xpd=TRUE)
    text(-2,-(ylim[2]-ylim[1])/5.5,labels=expression(bold(paste(n["1"],":",n["2"],":"))))
    
    dev.off()
    
  }
  
}

#     Condition c, as a function of the SD-ratio
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Because the comparison pattern is very similar whatever n = 20, 50 or 100par
# We will only plot the results as a function of the SD-ratio when n = 100
# Reason: we don't entirely remove the effect of delta, when dividing the variance by delta²
# But the effect of delta² is much reduced with larger sample size.
# We therefore chose n=100 in order to reduce a result that is an artefact of the way we computed
# The relative variance.

Path <-  "C:/Users/Admin/Documents/Github projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/"

for (j in seq_len(length(list.files(Path)))){
  
  
  File=read.table(paste0(Path,list.files(Path)[j]),header=T,sep=";",dec=",")  
  
  # Extract "Het_bal" subcategory (when n1=n2=20)
  
  id_sdSD_bal=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1==File$n2)&(File$n1==20)&(File$sd1.sd2 < 1),]))
  id_SDsd_bal=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1==File$n2)&(File$n1==20)&(File$sd1.sd2 > 1),]))
  id_Het_bal <- c(id_sdSD_bal,id_SDsd_bal)
  
  Sel <- File[id_Het_bal,]
  
  # Possible SD-ratios
  sdval <- as.numeric(levels(factor(Sel$sd1.sd2)))
  
  # Matrix containing biases information
  res <- matrix(0,6,6)  
  names<-expand.grid("bias_",c("Hedges","Glass1","Glass2","Shieh","cohen_delta_prime"))
  colnames(res) <- c("sd-ratio",paste0(names[,1],names[,2]))
  res[,1] <- sdval
  res[,2] <- tapply(Sel$relbias_Hedges,list(Sel$sd1.sd2),mean)[1:6]
  res[,3] <- tapply(Sel$relbias_Glass1,list(Sel$sd1.sd2),mean)[1:6]
  res[,4] <- tapply(Sel$relbias_Glass2,list(Sel$sd1.sd2),mean)[1:6]
  res[,5] <- tapply(Sel$relbias_Shieh,list(Sel$sd1.sd2),mean)[1:6]
  res[,6] <- tapply(Sel$relbias_cohen_delta_prime,list(Sel$sd1.sd2),mean)[1:6]
  # Select only bias columns
  res.relbias <- t(res[,2:6])
  colnames(res.relbias) <- paste0(1,":",1/sdval)
  
  param <- str_extract_all(list.files(Path)[j], "[[:digit:]]+\\.*[[:digit:]]*")
  if(param[[1]][2]==2.08){
    G1 <- -as.numeric(param[[1]][2])
  } else (G1 <- as.numeric(param[[1]][2]))
  G2 <- as.numeric(param[[1]][4])
  
  # Matrix containing variances information
  res3 <- matrix(0,6,6)  
  names<-expand.grid("var_",c("Hedges","Glass1","Glass2","Shieh","cohen_delta_prime"))
  colnames(res3) <- c("sd-ratio",paste0(names[,1],names[,2]))
  res3[,1] <- sdval 
  res3[,2] <- tapply(Sel$releff_Hedges,list(Sel$sd1.sd2),mean)[1:6]
  res3[,3] <- tapply(Sel$releff_Glass1,list(Sel$sd1.sd2),mean)[1:6]
  res3[,4] <- tapply(Sel$releff_Glass2,list(Sel$sd1.sd2),mean)[1:6]
  res3[,5] <- tapply(Sel$releff_Shieh,list(Sel$sd1.sd2),mean)[1:6]
  res3[,6] <- tapply(Sel$releff_cohen_delta_prime,list(Sel$sd1.sd2),mean)[1:6]
  # Select only bias columns
  res.releff <- t(res3[,2:6])
  colnames(res.releff) <- paste0(1,":",1/sdval)
  
  setwd(paste0("C:/Users/Admin/Documents/Github projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/Relative estimators of goodness/id_Het_bal/sd-ratio"))
  if(G1==-2.08){
    g1=2.08
  } else {g1=G1}
  
  png(file=paste0("bias_eff,G1=",g1, " & G2=",G2,";","id_Het_bal.png"),width=1400,height=1700, units = "px", res = 300)  
  
  par(mar = c(4,5,1.5,0),mfrow = c(2,1))   
  
  # plot for the relative bias
  
  if (j==1){ylabelbias=expression(paste("[E(" , hat(delta) , ")-",delta,"]/",delta))
  } else {ylabelbias=""}
  
  if (j==4){ylim=c(0,1.6)
  } else {ylim=c(0,.12)}
  
  barplot(res.relbias, 
          col = c("black","grey40","grey60","grey80","white"),
          beside = TRUE,
          main=paste0("G1=",G1,"; G2=",G2),
          xaxt="n",
          ylim=ylim,
          cex.lab=1.5,
          cex.main=1.5,
          ylab=ylabelbias
  )
  
  # plot for the relative variance
  
  if (j==1){ylabeleff=expression(paste("Var(" , hat(delta) , ")/",delta^2))
  } else {ylabeleff=""}

  if (j==4){ylim=c(0,4) # ylim=c(0,13) 
  } else {ylim=c(0,.4)} # ylim=c(0,3)
  
  barplot(res.releff, 
          col = c("black","grey40","grey60","grey80","white"),
          beside = TRUE,
          ylab = ylabeleff,
          ylim=ylim,
          cex.lab=1.5,
          args.legend = list(
            x = length(res.releff)*1.2,
            y = max(res.releff)+.5,
            bty="n"
          ))
  
  par(xpd=TRUE)
  text(-2,-(ylim[2]-ylim[1])/5.5,labels=expression(bold(paste("\u03c3"["1"],":","\u03c3"["2"],":"))))
  

  dev.off()
  
}


#     Condition d and e
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

# If sdSD & nN or SDsd & Nn: positive correlation between n and sd (rpos; condition d)
# If sdSD & Nn or SDsd & nN: negative correlation between n and sd (rneg; condition e)

# Note: There are too many conditions to represent them all in a single Figure. 
# We will therefore generate 3 figures per condition, as a function of the Total sample size (n1+n2)

### For a constant N, effect of the SD

Path <-  "C:/Users/Admin/Documents/Github projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/"

plot_hetr <- function(totalN){for (j in seq_len(length(list.files(Path)))){
  
  File=read.table(paste0(Path,list.files(Path)[j]),header=T,sep=";",dec=",")  
  
  # Extract "Het_unbal" subcategory (for a specific total sample size)
  
  if (totalN == 70){
    
    # Heteroscedasticity, positive correlation between n and sd ("Het_rpos")
    id_sdSD_nN=as.numeric(rownames(File[(File$n1+File$n2==70) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 < 1),]))
    id_SDsd_nN=as.numeric(rownames(File[(File$n1+File$n2==70) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 > 1),]))
    # Heteroscedasticity, negative correlation between n and sd ("Het_rneg")
    id_sdSD_Nn=as.numeric(rownames(File[(File$n1+File$n2==70) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 < 1),]))
    id_SDsd_Nn=as.numeric(rownames(File[(File$n1+File$n2==70) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 > 1),]))
    
    Conditions_id <- list(id_Het_firstsmaller=c(id_sdSD_nN,id_SDsd_nN),id_Het_firstlarger=c(id_sdSD_Nn,id_SDsd_Nn))
    
  } else if (totalN == 120){
    
    # Heteroscedasticity, positive correlation between n and sd ("Het_rpos")
    id_sdSD_nN=as.numeric(rownames(File[(File$n1+File$n2==120) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 < 1),]))
    id_SDsd_nN=as.numeric(rownames(File[(File$n1+File$n2==120) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 > 1),]))
    # Heteroscedasticity, negative correlation between n and sd ("Het_rneg")
    id_sdSD_Nn=as.numeric(rownames(File[(File$n1+File$n2==120) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 < 1),]))
    id_SDsd_Nn=as.numeric(rownames(File[(File$n1+File$n2==120) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 > 1),]))
    
    Conditions_id <- list(id_Het_firstsmaller=c(id_sdSD_nN,id_SDsd_nN),id_Het_firstlarger=c(id_sdSD_Nn,id_SDsd_Nn))
    
  } else if (totalN == 150){
    
    # Heteroscedasticity, positive correlation between n and sd ("Het_rpos")
    id_sdSD_nN=as.numeric(rownames(File[(File$n1+File$n2==150) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 < 1),]))
    id_SDsd_nN=as.numeric(rownames(File[(File$n1+File$n2==150) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 > 1),]))
    # Heteroscedasticity, negative correlation between n and sd ("Het_rneg")
    id_sdSD_Nn=as.numeric(rownames(File[(File$n1+File$n2==150) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 < 1),]))
    id_SDsd_Nn=as.numeric(rownames(File[(File$n1+File$n2==150) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 > 1),]))
    
    Conditions_id <- list(id_Het_firstsmaller=c(id_sdSD_nN,id_SDsd_nN),id_Het_firstlarger=c(id_sdSD_Nn,id_SDsd_Nn))
    
  }  
  
  
  for (i in seq_len(length(Conditions_id))){ 
    
    Sel <- File[Conditions_id[[i]],]
    
    # Possible SD-ratios and nratios
    sdval <- as.numeric(levels(factor(Sel$sd1.sd2)))
    n2val <- as.numeric(levels(factor(Sel$n2)))
    n1val<- as.numeric(levels(factor(Sel$n1)))
    combi <- expand.grid(sdval, n1val, n2val)
    
    # Matrix containing biases information
    res <- matrix(0,6,8)  
    names<-expand.grid("bias_",c("Hedges","Glass1","Glass2","Shieh","cohen_delta_prime"))
    colnames(res) <- c("n1","n2","sd-ratio",paste0(names[,1],names[,2]))
    res[,1:3] <- cbind(combi[,2],combi[,3],combi[,1])
    res[,4] <- tapply(Sel$relbias_Hedges,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:6]
    res[,5] <- tapply(Sel$relbias_Glass1,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:6]
    res[,6] <- tapply(Sel$relbias_Glass2,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:6]
    res[,7] <- tapply(Sel$relbias_Shieh,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:6]
    res[,8] <- tapply(Sel$relbias_cohen_delta_prime,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:6]
    
    # Select only bias columns
    res.relbias <- t(res[,4:8])
    colnames(res.relbias) <- paste0(res[,1],":",res[,2],"\n 1:",1/res[,3])
    
    param <- str_extract_all(list.files(Path)[j], "[[:digit:]]+\\.*[[:digit:]]*")
    if(param[[1]][2]==2.08){
      G1 <- -as.numeric(param[[1]][2])
    } else (G1 <- as.numeric(param[[1]][2]))
    G2 <- as.numeric(param[[1]][4])
    
    # Matrix containing variances information
    res3 <- matrix(0,6,8)  
    names<-expand.grid("eff_",c("Hedges","Glass1","Glass2","Shieh","cohen_delta_prime"))
    colnames(res3) <- c("n1","n2","sd-ratio",paste0(names[,1],names[,2]))
    res3[,1:3] <- cbind(combi[,2],combi[,3],combi[,1])
    res3[,4] <- tapply(Sel$releff_Hedges,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:6]
    res3[,5] <- tapply(Sel$releff_Glass1,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:6]
    res3[,6] <- tapply(Sel$releff_Glass2,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:6]
    res3[,7] <- tapply(Sel$releff_Shieh,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:6]
    res3[,8] <- tapply(Sel$releff_cohen_delta_prime,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:6]
    
    # Select only variances columns
    res.releff <- t(res3[,4:8])
    colnames(res.releff) <- paste0(res3[,1],":",res3[,2],"\n (1:",1/res3[,3],")")
    
    setwd(paste0("C:/Users/Admin/Documents/Github projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/Relative estimators of goodness/",names(Conditions_id)[i]))
    if(G1==-2.08){
      g1=2.08
    } else {g1=G1}
    png(file=paste0("G1=",g1, " & G2=",G2,";",names(Conditions_id)[i],"; N=",totalN,".png"),width=1400,height=1700, units = "px", res = 300)  
    
    par(mar = c(4,5,1.5,0),mfrow = c(2,1))   
    
    # plot for the relative bias
    
    if (j==1){ylabelbias=expression(paste("[E(" , hat(delta) , ")-",delta,"]/",delta))
    } else {ylabelbias=""}

    if (j==4){ylim=c(0,1.6)
    } else {ylim=c(0,.12)}
    
    barplot(res.relbias, 
            col = c("black","grey40","grey60","grey80","white"),
            beside = TRUE,
            main=paste0("G1=",G1,"; G2=",G2),
            xaxt="n",
            ylim=ylim,
            cex.lab=1.5,
            cex.main=1.5,
            ylab=ylabelbias
    )
    
    # plot for the variance
    
    if (j==1){ylabeleff=expression(paste("Var(" , hat(delta) , ")/",delta^2))
    } else {ylabeleff=""}

    if (j==4){ylim=c(0,4) #ylim=c(0,11)
    } else {ylim=c(0,.5)} #ylim=c(0,2) 
    
    barplot(res.releff, 
            col = c("black","grey40","grey60","grey80","white"),
            beside = TRUE,
            ylab = ylabeleff,
            cex.lab=1.5,
            cex.names=.8,
            ylim=ylim,
            args.legend = list(
              x = length(res.releff)*1.2,
              y = max(res.releff)+.5,
              bty="n"
            ))
    
    par(xpd=TRUE)
    text(-2,-(ylim[2]-ylim[1])/10.5,labels=expression(bold(paste(n["1"],":",n["2"],":"))))
    text(-2,-(ylim[2]-ylim[1])/5,labels=expression(bold(paste("\u03c3"["1"],":","\u03c3"["2"],":"))))
    
    dev.off()
    
  }
}
}

plot_hetr(totalN=70)
plot_hetr(totalN=120)
plot_hetr(totalN=150)


# Graphs are really hard to read, because due to the way we computed the relative variance, 

### For a constant SD, effect of N
Path <-  "C:/Users/Admin/Documents/Github projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/"

#Ratiolarger/smaller
plot_hetr <- function(ratio){for (j in seq_len(length(list.files(Path)))){
  
  File=read.table(paste0(Path,list.files(Path)[j]),header=T,sep=";",dec=",")  
  
  # Extract "Het_unbal" subcategory (for a specific total sample size)
  
  if (ratio == 10){
    
    # Heteroscedasticity
    id_sdSD_nN=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 == .1),]))
    id_sdSD_Nn=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 == .1),]))
    id_SDsd_Nn=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 == 10),]))
    id_SDsd_nN=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 == 10),]))
    
    Conditions_id <- list(id_Het_firstsmaller=c(id_sdSD_nN,id_sdSD_Nn),id_Het_firstlarger=c(id_SDsd_nN,id_SDsd_Nn))
    
  } else if (ratio == 4){
    
    id_sdSD_nN=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 == .25),]))
    id_sdSD_Nn=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 == .25),]))
    id_SDsd_Nn=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 == 4),]))
    id_SDsd_nN=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 == 4),]))
    
    Conditions_id <- list(id_Het_firstsmaller=c(id_sdSD_nN,id_sdSD_Nn),id_Het_firstlarger=c(id_SDsd_nN,id_SDsd_Nn))
    
  } else if (ratio == 2){
    
    id_sdSD_nN=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 == .5),]))
    id_sdSD_Nn=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 == .5),]))
    id_SDsd_Nn=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 == 2),]))
    id_SDsd_nN=as.numeric(rownames(File[(File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 == 2),]))
    
    Conditions_id <- list(id_Het_firstsmaller=c(id_sdSD_nN,id_sdSD_Nn),id_Het_firstlarger=c(id_SDsd_nN,id_SDsd_Nn))
    
  }  
  

  for (i in seq_len(length(Conditions_id))){ 
    
    Sel <- File[Conditions_id[[i]],]
    
    # Possible SD-ratios and nratios
    sd1val <- 1
    sd2val <- 1/as.numeric(levels(factor(Sel$sd1.sd2)))
    nval <- as.numeric(levels(factor(Sel$n1.n2)))
    combi <- expand.grid(nval, sd1val, sd2val)
    
    # Matrix containing biases information
    res <- matrix(0,6,9)  
    names<-expand.grid("bias_",c("Hedges","Glass1","Glass2","Shieh","cohen_delta_prime"))
    colnames(res) <- c("sd1","sd2","n1","n2",paste0(names[,1],names[,2]))
    res[,1:2] <- cbind(combi[,2],combi[,3])
    for (k in 1:length(combi[,1])){
      if(combi[k,1]==.2){
        res[k,3]<-20
        res[k,4] <- 100
      } else if (combi[k,1]==.4){
        res[k,3]<-20
        res[k,4] <- 50 
      } else if(combi[k,1]==.5){
        res[k,3]<-50 
        res[k,4] <- 100
      } else if(combi[k,1]==2){
        res[k,3]<-100 
        res[k,4] <- 50
      } else if(combi[k,1]==2.5){
        res[k,3]<-50 
        res[k,4] <- 20
      } else if(combi[k,1]==5){
        res[k,3]<-100 
        res[k,4] <- 20}
    }
    res[,5] <- tapply(Sel$relbias_Hedges,list(Sel$n1.n2,Sel$sd1.sd2),mean)[1:6]
    res[,6] <- tapply(Sel$relbias_Glass1,list(Sel$n1.n2,Sel$sd1.sd2),mean)[1:6]
    res[,7] <- tapply(Sel$relbias_Glass2,list(Sel$n1.n2,Sel$sd1.sd2),mean)[1:6]
    res[,8] <- tapply(Sel$relbias_Shieh,list(Sel$n1.n2,Sel$sd1.sd2),mean)[1:6]
    res[,9] <- tapply(Sel$relbias_Hedges,list(Sel$n1.n2,Sel$sd1.sd2),mean)[1:6]
    
    # Select only bias columns
    res.relbias <- t(res[,5:9])
    colnames(res.relbias) <- paste0("1:",res[,2],"\n (",res[,3],":",res[,4],")")

    setwd("C:/Users/Admin/Documents/Github_projects/Effect-sizes/Réflexion-extras/interaction_n_sd_shieh&hedges")
    write.table(res.relbias,"res.relbias.txt",sep=";",dec=",")
    
    param <- str_extract_all(list.files(Path)[j], "[[:digit:]]+\\.*[[:digit:]]*")
    if(param[[1]][2]==2.08){
      G1 <- -as.numeric(param[[1]][2])
    } else (G1 <- as.numeric(param[[1]][2]))
    G2 <- as.numeric(param[[1]][4])
    
    # Matrix containing variances information
    res3 <- matrix(0,6,9)  
    names<-expand.grid("eff_",c("Hedges","Glass1","Glass2","Shieh","cohen_delta_prime"))
    colnames(res3) <- c("n1","n2","sd1","sd2",paste0(names[,1],names[,2]))
    res3[,1:2] <- cbind(combi[,2],combi[,3])
    for (k in 1:length(combi[,1])){
      if(combi[k,1]==.2){
        res3[k,3]<-20
        res3[k,4] <- 100
      } else if (combi[k,1]==.4){
        res3[k,3]<-20
        res3[k,4] <- 50 
      } else if(combi[k,1]==.5){
        res3[k,3]<-50 
        res3[k,4] <- 100
      } else if(combi[k,1]==2){
        res3[k,3]<-100 
        res3[k,4] <- 50
      } else if(combi[k,1]==2.5){
        res3[k,3]<-50 
        res3[k,4] <- 20
      } else if(combi[k,1]==5){
        res3[k,3]<-100 
        res3[k,4] <- 20}
    }
    res3[,5] <- tapply(Sel$releff_Hedges,list(Sel$n1.n2,Sel$sd1.sd2),mean)[1:6]
    res3[,6] <- tapply(Sel$releff_Glass1,list(Sel$n1.n2,Sel$sd1.sd2),mean)[1:6]
    res3[,7] <- tapply(Sel$releff_Glass2,list(Sel$n1.n2,Sel$sd1.sd2),mean)[1:6]
    res3[,8] <- tapply(Sel$releff_Shieh,list(Sel$n1.n2,Sel$sd1.sd2),mean)[1:6]
    res3[,9] <- tapply(Sel$releff_cohen_delta_prime,list(Sel$n1.n2,Sel$sd1.sd2),mean)[1:6]
    
    # Select only variances columns
    res.releff <- t(res3[,5:9])
    colnames(res.releff) <- paste0("1:",res[,2],"\n (",res[,3],":",res[,4],")")

    setwd("C:/Users/Admin/Documents/Github_projects/Effect-sizes/Réflexion-extras/interaction_n_sd_shieh&hedges")
    write.table(res.releff,"res.releff.txt",sep=";",dec=",")
    
    setwd(paste0("C:/Users/Admin/Documents/Github projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/Relative estimators of goodness/",names(Conditions_id)[i]))
    if(G1==-2.08){
      g1=2.08
    } else {g1=G1}
    png(file=paste0("G1=",g1, " & G2=",G2,";",names(Conditions_id)[i],"; SDR=",ratio,".png"),width=1400,height=1700, units = "px", res = 300)  
    
    par(mar = c(4,5,1.5,0),mfrow = c(2,1))   
    
    # plot for the relative bias
    
    if (j==1){ylabelbias=expression(paste("[E(" , hat(delta) , ")-",delta,"]/",delta))
    } else {ylabelbias=""}
    
    if (j==4){ylim=c(0,1.6)
    } else {ylim=c(0,.16)}
    
    barplot(res.relbias, 
            col = c("black","grey40","grey60","grey80","white"),
            beside = TRUE,
            main=paste0("G1=",G1,"; G2=",G2),
            xaxt="n",
            ylim=ylim,
            cex.lab=1.5,
            cex.name=.6,
            cex.main=1.5,
            ylab=ylabelbias
    )
    
    # plot for the variance
    
    if (j==1){ylabeleff=expression(paste("Var(" , hat(delta) , ")/",delta^2))
    } else {ylabeleff=""}

    if (i==1){
      if(ratio==10){    
          if (j==4){ylim=c(0,11)
          } else {ylim=c(0,2)}  
      } else if(ratio==4) {
          if (j==4){ylim=c(0,4) 
          } else {ylim=c(0,.4)}  
      } else if (ratio==2){
          if (j==4){ylim=c(0,2) 
      } else {ylim=c(0,.2)}}
    }  else if (i==2) {
      if(ratio==10){    
        if (j==4){ylim=c(0,1.6)
        } else {ylim=c(0,.16)} 
      } else if(ratio==4) {
        if (j==4){ylim=c(0,1.3) 
        } else {ylim=c(0,.13)}  
      } else if (ratio==2){
        if (j==4){ylim=c(0,1.3) 
        } else {ylim=c(0,.13)}}
    }
    
    barplot(res.releff, 
            col = c("black","grey40","grey60","grey80","white"),
            beside = TRUE,
            ylab = ylabeleff,
            cex.lab=1.5,
            cex.names=.6,
            ylim=ylim,
            args.legend = list(
              x = length(res.releff)*1.2,
              y = max(res.releff)+.5,
              bty="n"
            ))
    
    par(xpd=TRUE)
    text(-2,-(ylim[2]-ylim[1])/10.5,labels=expression(bold(paste("\u03c3"["1"],":","\u03c3"["2"],":"))),cex=.8)
    text(-2,-(ylim[2]-ylim[1])/5,labels=expression(bold(paste(n["1"],":",n["2"],":"))),cex=.8)
    
    dev.off()
    
  }
}
}

plot_hetr(ratio=10)
plot_hetr(ratio=4)
plot_hetr(ratio=2)




