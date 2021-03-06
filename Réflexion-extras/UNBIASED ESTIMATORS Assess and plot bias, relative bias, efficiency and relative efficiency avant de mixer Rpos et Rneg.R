library(stringr)

#setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs")
#dir.create("Quality of ES measures")

#Mainfolder="C:/Users/Marie/Documents/ES MEASURES/"
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
    
    # Compute bias View(file)
    #    bias_Hedges <- mean(file[,10]) - cohen_delta # Hedges' g_s is a corrected estimate of cohen's delta!
    #    bias_glass1 <- mean(file[,13]) - glass_delta1 
    #    bias_glass2 <- mean(file[,14]) - glass_delta2
    #    bias_shieh <- mean(file[,16]) - shieh_delta
    #    bias_cohen_delta_prime <- mean(file[,18]) - cohen_delta_prime
    
    # Compute relative bias
    #    relbias_Hedges <- (mean(file[,10]) - cohen_delta)/cohen_delta
    #    relbias_glass1 <- (mean(file[,13]) - glass_delta1)/glass_delta1 
    #    relbias_glass2 <- (mean(file[,14]) - glass_delta2)/glass_delta2
    #    relbias_shieh <- (mean(file[,16]) - shieh_delta)/shieh_delta
    #    relbias_cohen_delta_prime <- (mean(file[,18]) - cohen_delta_prime)/cohen_delta_prime
    
    # Compute variance
    #    eff_Hedges <- var(file[,10])
    #    eff_glass1 <- var(file[,13])
    #    eff_glass2 <- var(file[,14])
    #    eff_shieh <- var(file[,16])
    #    eff_cohen_delta_prime <- var(file[,18])
    
    # Compute relative variance
    #    releff_Hedges <- var(file[,10])/cohen_delta^2
    #    releff_glass1 <- var(file[,13])/glass_delta1^2 
    #    releff_glass2 <- var(file[,14])/glass_delta2^2
    #    releff_shieh <- var(file[,16])/shieh_delta^2
    #    releff_cohen_delta_prime <- var(file[,18])/cohen_delta_prime^2
    
    #    good_mes[j,1] <- n1
    #    good_mes[j,2] <- n2
    #    good_mes[j,3] <- n1/n2
    #    good_mes[j,4] <- m1-m2
    #    good_mes[j,5] <- sd1/sd2
    
    #    good_mes[j,6] <- bias_Hedges
    #    good_mes[j,7] <- bias_glass1
    #    good_mes[j,8] <- bias_glass2
    #    good_mes[j,9] <- bias_shieh
    #    good_mes[j,10] <- bias_cohen_delta_prime
    
    #   good_mes[j,11] <- relbias_Hedges
    #    good_mes[j,12] <- relbias_glass1
    #    good_mes[j,13] <- relbias_glass2
    #    good_mes[j,14] <- relbias_shieh
    #    good_mes[j,15] <- relbias_cohen_delta_prime
    
    #    good_mes[j,16] <- eff_Hedges
    #    good_mes[j,17] <- eff_glass1
    #    good_mes[j,18] <- eff_glass2
    #    good_mes[j,19] <- eff_shieh
    #    good_mes[j,20] <- eff_cohen_delta_prime
    
    #    good_mes[j,21] <- releff_Hedges
    #    good_mes[j,22] <- releff_glass1
    #    good_mes[j,23] <- releff_glass2
    #    good_mes[j,24] <- releff_shieh
    #    good_mes[j,25] <- releff_cohen_delta_prime
    
    #  }
  
  #  setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/Data summary/Unbiased estimators")

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

setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/")

# png(file="legend.png",width=1500,height=1000, units = "px", res = 300)  
# plot(1,1,bty="n",xaxt="n",yaxt="n",ylim=c(.62,.67),main="",xlab="",ylab="",pch=19,type="o")
# legend("center", 
#       legend=c(expression(paste("Hedges' ",g[s])),expression(paste("Glass's ",g[s],"(",sigma," =",S[1],")")),expression(paste("Glass's ",g[s],"(",sigma," =",S[2],")")),expression(paste("Shieh's ",g[s])),
#                expression(paste("Hedges' ",g[s],"'"))),
#       fill=c("black","grey40","grey60","grey80","white"),
#       bty="n"
# )
# dev.off()

png(file="legend.png",width=8000,height=800, units = "px", res = 300)  

plot(1,1,bty="n",xaxt="n",yaxt="n",ylim=c(.62,.67),main="",xlab="",ylab="",pch=19,type="o")
legend("center", 
       legend=c(expression(paste("Hedges' ",g[s])),expression(paste("Glass's ",g[s],"(",sigma," =",S[1],")")),expression(paste("Glass's ",g[s],"(",sigma," =",S[2],")")),expression(paste("Shieh's ",g[s])),
                expression(paste("Hedges' ",g[s],"'"))),
       fill=c("black","grey40","grey60","grey80","white"),
       bty="n",horiz=T,cex=2
)

dev.off()

Path <-  "C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/Data summary/Unbiased estimators/"

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
    combi <- expand.grid(n2val, n1val)
    
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
    
    setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/Raw estimators of goodness/")
    #dir.create(names(Conditions_id)[i])

    setwd(paste0("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/Raw estimators of goodness/",names(Conditions_id)[i]))
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

Path <-  "C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/Data summary/Unbiased estimators/"

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
  
  setwd(paste0("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/Raw estimators of goodness/id_Het_bal/sd-ratio"))
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

Path <-  "C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/Data summary/Unbiased estimators/"

plot_hetr <- function(totalN){for (j in seq_len(length(list.files(Path)))){
  
  File=read.table(paste0(Path,list.files(Path)[j]),header=T,sep=";",dec=",")  
  
  # Extract "Het_bal" subcategory (for a specific total sample size)
  
  if (totalN == 70){
    
    # Heteroscedasticity, positive correlation between n and sd ("Het_rpos")
    id_sdSD_nN=as.numeric(rownames(File[(File$n1+File$n2==70) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 < 1),]))
    id_SDsd_Nn=as.numeric(rownames(File[(File$n1+File$n2==70) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 > 1),]))
    # Heteroscedasticity, negative correlation between n and sd ("Het_rneg")
    id_SDsd_nN=as.numeric(rownames(File[(File$n1+File$n2==70) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 > 1),]))
    id_sdSD_Nn=as.numeric(rownames(File[(File$n1+File$n2==70) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 < 1),]))
    
    Conditions_id <- list(id_Het_rpos=c(id_sdSD_nN,id_SDsd_Nn),id_Het_rneg=c(id_SDsd_nN,id_sdSD_Nn))
    
  } else if (totalN == 120){

    # Heteroscedasticity, positive correlation between n and sd ("Het_rpos")
    id_sdSD_nN=as.numeric(rownames(File[(File$n1+File$n2==120) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 < 1),]))
    id_SDsd_Nn=as.numeric(rownames(File[(File$n1+File$n2==120) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 > 1),]))
    # Heteroscedasticity, negative correlation between n and sd ("Het_rneg")
    id_SDsd_nN=as.numeric(rownames(File[(File$n1+File$n2==120) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 > 1),]))
    id_sdSD_Nn=as.numeric(rownames(File[(File$n1+File$n2==120) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 < 1),]))
    
    Conditions_id <- list(id_Het_rpos=c(id_sdSD_nN,id_SDsd_Nn),id_Het_rneg=c(id_SDsd_nN,id_sdSD_Nn))
    
  } else if (totalN == 150){

    # Heteroscedasticity, positive correlation between n and sd ("Het_rpos")
    id_sdSD_nN=as.numeric(rownames(File[(File$n1+File$n2==150) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 < 1),]))
    id_SDsd_Nn=as.numeric(rownames(File[(File$n1+File$n2==150) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 > 1),]))
    # Heteroscedasticity, negative correlation between n and sd ("Het_rneg")
    id_SDsd_nN=as.numeric(rownames(File[(File$n1+File$n2==150) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 > 1),]))
    id_sdSD_Nn=as.numeric(rownames(File[(File$n1+File$n2==150) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 < 1),]))
    
    Conditions_id <- list(id_Het_rpos=c(id_sdSD_nN,id_SDsd_Nn),id_Het_rneg=c(id_SDsd_nN,id_sdSD_Nn))

  }  


  for (i in seq_len(length(Conditions_id))){ 
    
  Sel <- File[Conditions_id[[i]],]

  # Possible SD-ratios and nratios
  sdval <- as.numeric(levels(factor(Sel$sd1.sd2)))
  n2val <- as.numeric(levels(factor(Sel$n2)))
  n1val<- as.numeric(levels(factor(Sel$n1)))
  combi <- expand.grid(sdval, n2val, n1val)
  
  # Matrix containing biases information
  res <- matrix(0,24,8)  
  names<-expand.grid("bias_",c("Hedges","Glass1","Glass2","Shieh","cohen_delta_prime"))
  colnames(res) <- c("n1","n2","sd-ratio",paste0(names[,1],names[,2]))
  res[,1:3] <- cbind(combi[,2],combi[,3],combi[,1])
  res[,4] <- tapply(Sel$bias_Hedges,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:24]
  res[,5] <- tapply(Sel$bias_Glass1,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:24]
  res[,6] <- tapply(Sel$bias_Glass2,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:24]
  res[,7] <- tapply(Sel$bias_Shieh,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:24]
  res[,8] <- tapply(Sel$bias_cohen_delta_prime,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:24]
  
  # Select only rows with no "NA"  
  res <- subset(res,res[,4] != "NA") 
  
  # Select only bias columns
  res.bias <- t(res[,4:8])
  colnames(res.bias) <- paste0(res[,1],":",res[,2],"\n 1:",1/res[,3])
  
  param <- str_extract_all(list.files(Path)[j], "[[:digit:]]+\\.*[[:digit:]]*")
  if(param[[1]][2]==2.08){
    G1 <- -as.numeric(param[[1]][2])
  } else (G1 <- as.numeric(param[[1]][2]))
  G2 <- as.numeric(param[[1]][4])
  
  # Matrix containing variances information
  res3 <- matrix(0,24,8)  
  names<-expand.grid("bias_",c("Hedges","Glass1","Glass2","Shieh","cohen_delta_prime"))
  colnames(res3) <- c("n1","n2","sd-ratio",paste0(names[,1],names[,2]))
  res3[,1:3] <- cbind(combi[,2],combi[,3],combi[,1])
  res3[,4] <- tapply(Sel$eff_Hedges,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:24]
  res3[,5] <- tapply(Sel$eff_Glass1,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:24]
  res3[,6] <- tapply(Sel$eff_Glass2,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:24]
  res3[,7] <- tapply(Sel$eff_Shieh,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:24]
  res3[,8] <- tapply(Sel$eff_cohen_delta_prime,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:24]
  
  # Select only rows with no "NA"  
  res3 <- subset(res3,res3[,4] != "NA") 
  
  # Select only variances columns
  res.eff <- t(res3[,4:8])
  colnames(res.eff) <- paste0(res3[,1],":",res3[,2],"\n (1:",1/res3[,3],")")

  setwd(paste0("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/Raw estimators of goodness/",names(Conditions_id)[i]))
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

setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/")

Path <-  "C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/Data summary/Unbiased estimators/"

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
    combi <- expand.grid(n2val, n1val)
    
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
    
    setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/Relative estimators of goodness/")
    #dir.create(names(Conditions_id)[i])
    
    setwd(paste0("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/Relative estimators of goodness/",names(Conditions_id)[i]))
    if(G1==-2.08){
      g1=2.08
    } else {g1=G1}
    png(file=paste0("bias_eff,G1=",g1, " & G2=",G2,";",names(Conditions_id)[i], ".png"),width=1400,height=1700, units = "px", res = 300)  
    
    par(mar = c(4,5,1.5,0),mfrow = c(2,1))   
    
    # plot for the relative bias
    
    if (j==1){ylabelbias=expression(paste("(E(" , hat(delta) , ") -",delta,")/",delta ))
    } else {ylabelbias=""}
    
    if (j==4){ylim=c(0,1.5)
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
    
    # plot for the variance (need to be modified later)
    
    if (j==1){ylabeleff=expression(paste("Var(" , hat(delta) , ")"))
    } else {ylabeleff=""}
    
    if (j==4){ylim=c(0,6)
    } else {ylim=c(0,1)}
    
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
    
    
    dev.off()
    
  }
  
}

#     Condition c, as a function of the SD-ratio
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Because the comparison pattern is very similar whatever n = 20, 50 or 100
# We will only plot the results as a function of the SD-ratio when n = 20

Path <-  "C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/Data summary/Unbiased estimators/"

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
  
  setwd(paste0("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/Relative estimators of goodness/id_Het_bal/sd-ratio"))
  if(G1==-2.08){
    g1=2.08
  } else {g1=G1}
  png(file=paste0("bias_eff,G1=",g1, " & G2=",G2,";",names(Conditions_id)[i], ".png"),width=1400,height=1700, units = "px", res = 300)  
  
  par(mar = c(4,5,1.5,0),mfrow = c(2,1))   
  
  # plot for the bias
  
  if (j==1){ylabelbias=expression(paste("E(" , hat(delta) , ") -",delta ))
  } else {ylabelbias=""}
  
  if (j==4){ylim=c(0,1.4)
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
  
  if (j==1){ylabeleff=expression(paste("Var(" , hat(delta) , ")"))
  } else {ylabeleff=""}
  
  barplot(res.releff, 
          col = c("black","grey40","grey60","grey80","white"),
          beside = TRUE,
          ylab = ylabeleff,
          cex.lab=1.5,
          args.legend = list(
            x = length(res.releff)*1.2,
            y = max(res.releff)+.5,
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

Path <-  "C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/Data summary/Unbiased estimators/"

plot_hetr <- function(totalN){for (j in seq_len(length(list.files(Path)))){
  
  File=read.table(paste0(Path,list.files(Path)[j]),header=T,sep=";",dec=",")  
  
  # Extract "Het_bal" subcategory (for a specific total sample size)
  
  if (totalN == 70){
    
    # Heteroscedasticity, positive correlation between n and sd ("Het_rpos")
    id_sdSD_nN=as.numeric(rownames(File[(File$n1+File$n2==70) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 < 1),]))
    id_SDsd_Nn=as.numeric(rownames(File[(File$n1+File$n2==70) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 > 1),]))
    # Heteroscedasticity, negative correlation between n and sd ("Het_rneg")
    id_SDsd_nN=as.numeric(rownames(File[(File$n1+File$n2==70) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 > 1),]))
    id_sdSD_Nn=as.numeric(rownames(File[(File$n1+File$n2==70) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 < 1),]))
    
    Conditions_id <- list(id_Het_rpos=c(id_sdSD_nN,id_SDsd_Nn),id_Het_rneg=c(id_SDsd_nN,id_sdSD_Nn))
    
  } else if (totalN == 120){
    
    # Heteroscedasticity, positive correlation between n and sd ("Het_rpos")
    id_sdSD_nN=as.numeric(rownames(File[(File$n1+File$n2==120) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 < 1),]))
    id_SDsd_Nn=as.numeric(rownames(File[(File$n1+File$n2==120) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 > 1),]))
    # Heteroscedasticity, negative correlation between n and sd ("Het_rneg")
    id_SDsd_nN=as.numeric(rownames(File[(File$n1+File$n2==120) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 > 1),]))
    id_sdSD_Nn=as.numeric(rownames(File[(File$n1+File$n2==120) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 < 1),]))
    
    Conditions_id <- list(id_Het_rpos=c(id_sdSD_nN,id_SDsd_Nn),id_Het_rneg=c(id_SDsd_nN,id_sdSD_Nn))
    
  } else if (totalN == 150){
    
    # Heteroscedasticity, positive correlation between n and sd ("Het_rpos")
    id_sdSD_nN=as.numeric(rownames(File[(File$n1+File$n2==150) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 < 1),]))
    id_SDsd_Nn=as.numeric(rownames(File[(File$n1+File$n2==150) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 > 1),]))
    # Heteroscedasticity, negative correlation between n and sd ("Het_rneg")
    id_SDsd_nN=as.numeric(rownames(File[(File$n1+File$n2==150) & (File$m1.m2!=0)&(File$n1<File$n2)&(File$sd1.sd2 > 1),]))
    id_sdSD_Nn=as.numeric(rownames(File[(File$n1+File$n2==150) & (File$m1.m2!=0)&(File$n1>File$n2)&(File$sd1.sd2 < 1),]))
    
    Conditions_id <- list(id_Het_rpos=c(id_sdSD_nN,id_SDsd_Nn),id_Het_rneg=c(id_SDsd_nN,id_sdSD_Nn))
    
  }  
  
  
  for (i in seq_len(length(Conditions_id))){ 
    
    Sel <- File[Conditions_id[[i]],]
    
    # Possible SD-ratios and nratios
    sdval <- as.numeric(levels(factor(Sel$sd1.sd2)))
    n2val <- as.numeric(levels(factor(Sel$n2)))
    n1val<- as.numeric(levels(factor(Sel$n1)))
    combi <- expand.grid(sdval, n2val, n1val)
    
    # Matrix containing biases information
    res <- matrix(0,24,8)  
    names<-expand.grid("bias_",c("Hedges","Glass1","Glass2","Shieh","cohen_delta_prime"))
    colnames(res) <- c("n1","n2","sd-ratio",paste0(names[,1],names[,2]))
    res[,1:3] <- cbind(combi[,2],combi[,3],combi[,1])
    res[,4] <- tapply(Sel$relbias_Hedges,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:24]
    res[,5] <- tapply(Sel$relbias_Glass1,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:24]
    res[,6] <- tapply(Sel$relbias_Glass2,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:24]
    res[,7] <- tapply(Sel$relbias_Shieh,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:24]
    res[,8] <- tapply(Sel$relbias_cohen_delta_prime,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:24]
    
    # Select only rows with no "NA"  
    res <- subset(res,res[,4] != "NA") 
    
    # Select only bias columns
    res.relbias <- t(res[,4:8])
    colnames(res.relbias) <- paste0(res[,1],":",res[,2],"\n 1:",1/res[,3])
    
    param <- str_extract_all(list.files(Path)[j], "[[:digit:]]+\\.*[[:digit:]]*")
    if(param[[1]][2]==2.08){
      G1 <- -as.numeric(param[[1]][2])
    } else (G1 <- as.numeric(param[[1]][2]))
    G2 <- as.numeric(param[[1]][4])
    
    # Matrix containing variances information
    res3 <- matrix(0,24,8)  
    names<-expand.grid("bias_",c("Hedges","Glass1","Glass2","Shieh","cohen_delta_prime"))
    colnames(res3) <- c("n1","n2","sd-ratio",paste0(names[,1],names[,2]))
    res3[,1:3] <- cbind(combi[,2],combi[,3],combi[,1])
    res3[,4] <- tapply(Sel$releff_Hedges,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:24]
    res3[,5] <- tapply(Sel$releff_Glass1,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:24]
    res3[,6] <- tapply(Sel$releff_Glass2,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:24]
    res3[,7] <- tapply(Sel$releff_Shieh,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:24]
    res3[,8] <- tapply(Sel$releff_cohen_delta_prime,list(Sel$sd1.sd2,Sel$n1,Sel$n2),mean)[1:24]
    
    # Select only rows with no "NA"  
    res3 <- subset(res3,res3[,4] != "NA") 
    
    # Select only variances columns
    res.releff <- t(res3[,4:8])
    colnames(res.releff) <- paste0(res3[,1],":",res3[,2],"\n (1:",1/res3[,3],")")
    
    setwd(paste0("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Quality of ES measures/Graphs/Unbiased estimators/Relative estimators of goodness/",names(Conditions_id)[i]))
    if(G1==-2.08){
      g1=2.08
    } else {g1=G1}
    png(file=paste0("G1=",g1, " & G2=",G2,";",names(Conditions_id)[i],"; N=",totalN,".png"),width=1400,height=1700, units = "px", res = 300)  
    
    par(mar = c(4,5,1.5,0),mfrow = c(2,1))   
    
    # plot for the bias
    
    if (j==1){ylabelbias=expression(paste("E(" , hat(delta) , ") -",delta ))
    } else {ylabelbias=""}

    if (j==4){ylim=c(0,1.4)
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
    
    if (j==1){ylabeleff=expression(paste("Var(" , hat(delta) , ")"))
    } else {ylabeleff=""}
    
    barplot(res.releff, 
            col = c("black","grey40","grey60","grey80","white"),
            beside = TRUE,
            ylab = ylabeleff,
            cex.lab=1.5,
            cex.names=.8,
            args.legend = list(
              x = length(res.releff)*1.2,
              y = max(res.releff)+.5,
              bty="n"
            ))
    
    
    dev.off()
    
  }
}
}

plot_hetr(totalN=70)
plot_hetr(totalN=120)
plot_hetr(totalN=150)





