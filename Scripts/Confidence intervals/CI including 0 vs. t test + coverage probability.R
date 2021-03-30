# We determine if the real population effect size is inside the confidence interval (in order to asses the coverage probability)
# For Supplemental Material 4, we also computed pvalues (based on t-values) and determined if 0 is 
# out the C.I around a point estimate using the same quantity as the t statistic.

library(stringr)

Mainfolder="D:/Documents/CI.MEASURES/"
subfolder=list.files(Mainfolder)
Folder=paste0(Mainfolder,subfolder)

for (i in seq_len(length(Folder))){ 

  for (j in seq_len(length(list.files(Folder[i])))){

    filepath = paste0(Folder[i],"/",list.files(Folder[i])[j])
    file=readRDS(filepath)
    
    ttests <- c("tcohen","tglass1","tglass2","tcohenprime","tshieh")
    col.ptests <- do.call(paste, c(expand.grid("pvalue",ttests), sep = "_"))
    parameters <- c("nullincluded","ESincluded") 
    est <- c("Bcohen","Bglass1","Bglass2","Bcohen d'","Bshieh","Uhedge","Uglass1","Uglass2","Ucohen d'","Ushieh")
    col.CI <- do.call(paste, c(expand.grid(parameters,est), sep = "_"))
    res<-matrix(0,length(file[,1]),length(c(col.ptests,col.CI))) 
    colnames(res) <- c(col.ptests,col.CI)

    # Extracting population parameters values and sample sizes from file names
    param <- str_extract_all(list.files(Folder[i])[j], "[[:digit:]]+\\.*[[:digit:]]*")
    n1 <- as.numeric(param[[1]][5])
    n2 <- as.numeric(param[[1]][6])
    m1 <- as.numeric(param[[1]][7])
    m2 <- as.numeric(param[[1]][8])
    sd1 <- as.numeric(param[[1]][9])
    sd2 <- as.numeric(param[[1]][10])

    #Cohen's d
    cohen_delta <- (m1-m2)/sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
    res[,1] <- 2*(1-pt(abs(file[,1]),df=file[,2]))
    #Glass's d sd1 
    glass_delta1 <- (m1-m2)/sd1
    res[,2] <- 2*(1-pt(abs(file[,4]),df=file[,5]))
    #Glass's d sd2 
    glass_delta2 <- (m2-m1)/sd2
    res[,3] <- 2*(1-pt(abs(file[,7]),df=file[,8]))
    #Cohen's d' 
    cohenprime_delta <- (m1-m2)/sqrt((sd1^2+sd2^2)/2)
    res[,4] <- 2*(1-pt(abs(file[,10]),df=file[,11]))
    #Shieh's d 
    q1 <- n1/(n1+n2)
    q2 <- n2/(n1+n2)
    Shieh_delta <- (m1-m2)/sqrt(sd1^2/q1+sd2^2/q2)       
    res[,5] <- 2*(1-pt(abs(file[,13]),df=file[,14]))
    
    
    #CI Cohen's d
    res[,6] <- (0>=file[,17] & 0<=file[,18])  
    res[,7] <- (cohen_delta>=file[,17] & cohen_delta<=file[,18])  
    #CI Glass's d1    
    res[,8] <- (0>=file[,20] & 0<=file[,21])  
    res[,9] <- (glass_delta1>=file[,20] & glass_delta1<=file[,21])  
    #CI Glass's d2
    res[,10] <- (0>=file[,23] & 0<=file[,24])  
    res[,11] <- (glass_delta2>=file[,23] & glass_delta2<=file[,24])  
    #CI Cohen's dprime
    res[,12] <- (0>=file[,26] & 0<=file[,27])  
    res[,13] <- (cohenprime_delta>=file[,26] & cohenprime_delta<=file[,27])  
    #CI Shieh's d
    res[,14] <- (0>=file[,29] & 0<=file[,30])  
    res[,15] <- (Shieh_delta>=file[,29] & Shieh_delta<=file[,30])  
    
    #CI Hedges's g
    res[,16] <- (0>=file[,32] & 0<=file[,33])  
    res[,17] <- (cohen_delta>=file[,32] & cohen_delta<=file[,33])  
    #Glass's g sd1 
    res[,18] <- (0>=file[,35] & 0<=file[,36])  
    res[,19] <- (glass_delta1>=file[,35] & glass_delta1<=file[,36])  
    #Glass's g sd2 
    res[,20] <- (0>=file[,38] & 0<=file[,39])  
    res[,21] <- (glass_delta2>=file[,38] & glass_delta2<=file[,39])  
    #Cohen's g' 
    res[,22] <- (0>=file[,41] & 0<=file[,42])  
    res[,23] <- (cohenprime_delta>=file[,41] & cohenprime_delta<=file[,42])  
    #Shieh's g 
    res[,24] <- (0>=file[,44] & 0<=file[,45])  
    res[,25] <- (Shieh_delta>=file[,44] & Shieh_delta<=file[,45])  
   
    chem <- "D:/Documents/CIvs.test/"
    if(param[[1]][2]==2.08){
      G1 <- -2.08  
    } else {G1 <- param[[1]][2]}
    G2 = param[[1]][4]
    
    sschem<- paste0("G1=",G1,",G2=",G2) 
    setwd(dir=paste0(chem,sschem)) # destination file  
    
    fname=paste0("G1=",G1,", G2=",G2,", n=[",n1,",",n2,"], means=[",m1,",",m2,"], sds=[",sd1,",",sd2,"].rds")
    saveRDS(res, file = fname)
    
  }
  
}

#---------------------------------------------------------------------------------------------------------------------------

# For each file, we compute the coverage probability

Mainfolder="D:/Documents/CIvs.test/"
subfolder=list.files(Mainfolder)
Folder=paste0(Mainfolder,subfolder)

for (i in seq_len(length(Folder))){ 
  
  # set up empty container for all estimated parameters
  parameter <- "C.P_"
  est <- c("cohen.d","glass1.d","glass2.d","cohen.dprime","shieh.d","hedges.g","glass1.g","glass2.g","hedges.gprime","shieh.g")
  col.CI <- do.call(paste0, c(expand.grid(parameter,est)))
  
  finalres <-matrix(0,length(list.files(Folder[i])),length(c("n1","n2","m1","m2","sd1","sd2",col.CI)))
  colnames(finalres) <- c("n1","n2","m1","m2","sd1","sd2",col.CI)
  
  for (j in seq_len(length(list.files(Folder[i])))){
    
    filepath = paste0(Folder[i],"/",list.files(Folder[i])[j])
    file=readRDS(filepath)
    
    # Extracting population parameters values and sample sizes from file names
    param <- str_extract_all(list.files(Folder[i])[j], "[[:digit:]]+\\.*[[:digit:]]*")
    n1 <- as.numeric(param[[1]][5])
    n2 <- as.numeric(param[[1]][6])
    m1 <- as.numeric(param[[1]][7])
    m2 <- as.numeric(param[[1]][8])
    sd1 <- as.numeric(param[[1]][9])
    sd2 <- as.numeric(param[[1]][10])
    
    finalres[j,1:6] <- c(n1,n2,m1,m2,sd1,sd2)    

    # Coverage probability = proportion of C.I including the population effect size
    
    CIincludingdelta_Bcohen <- sum(file[,7])/length(file[,7])
    finalres[j,7] <- CIincludingdelta_Bcohen
    CIincludingdelta_Bglass1 <- sum(file[,9])/length(file[,9])
    finalres[j,8] <- CIincludingdelta_Bglass1
    CIincludingdelta_Bglass2 <- sum(file[,11])/length(file[,11])
    finalres[j,9] <- CIincludingdelta_Bglass2
    CIincludingdelta_Bcohenprime <- sum(file[,13])/length(file[,13])
    finalres[j,10] <- CIincludingdelta_Bcohenprime
    CIincludingdelta_Bshieh <- sum(file[,15])/length(file[,15])
    finalres[j,11] <- CIincludingdelta_Bshieh
    
    CIincludingdelta_Ucohen <- sum(file[,17])/length(file[,17])
    finalres[j,12] <- CIincludingdelta_Ucohen
    CIincludingdelta_Uglass1 <- sum(file[,19])/length(file[,19])
    finalres[j,13] <- CIincludingdelta_Uglass1
    CIincludingdelta_Uglass2 <- sum(file[,21])/length(file[,21])
    finalres[j,14] <- CIincludingdelta_Uglass2
    CIincludingdelta_Ucohenprime <- sum(file[,23])/length(file[,23])
    finalres[j,15] <- CIincludingdelta_Ucohenprime
    CIincludingdelta_Ushieh <- sum(file[,25])/length(file[,25])
    finalres[j,16] <- CIincludingdelta_Ushieh
    
  }
  
  setwd("D:/Documents/Github_projects/Effect-sizes/Scripts outputs/Confidence intervals/Data summary/.txt files") 
  shape <- str_extract_all(Folder[i], "[[:digit:]]+\\.*[[:digit:]]*")
  if(shape[[1]][2]==2.08){
    G1 <- -2.08  
  } else {G1 <- shape[[1]][2]}
  G2 = shape[[1]][4]
  
  fname=paste0("CovProb_G1=",G1,", G2=",G2,".txt")
  write.table(finalres,fname,sep=";",dec=",")
  
}

#-------------------------------------------------------------------------------------------------------

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

setwd("D:/Documents/Github_projects/Effect-sizes/Scripts outputs/Confidence intervals/Graphs/Unbiased estimators/")

Path <-  "D:/Documents/Github_projects/Effect-sizes/Scripts outputs/Confidence intervals/Data summary/.txt files/"

for (j in seq_len(length(list.files(Path)))){
  
  File=read.table(paste0(Path,list.files(Path)[j]),header=T,sep=";",dec=",")  
  
  # Extract three subcategories from File: "Hom_bal", "Hom_unbal", "Het_bal"
  
  ### Conditions id 
  
  # Homoscedasticity, balanced designs ("Hom_bal")
  id_Hom_bal = as.numeric(rownames(File[(File$n1==File$n2)&(File$sd1 == File$sd2),]))
  # Homoscedasticity, unbalanced designs ("Hom_unbal")  
  id_Hom_nN = as.numeric(rownames(File[(File$n1 < File$n2)&(File$sd1 == File$sd2),]))
  id_Hom_Nn = as.numeric(rownames(File[(File$n1 > File$n2)&(File$sd1 == File$sd2),]))
  # Heteroscedasticity, balanced designs ("Het_bal")
  id_sdSD_bal=as.numeric(rownames(File[(File$n1==File$n2)&(File$sd1 < File$sd2),]))
  id_SDsd_bal=as.numeric(rownames(File[(File$n1==File$n2)&(File$sd1 > File$sd2),]))
  
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
    names<-expand.grid("CP_",c("Hedges_g","Glass_g1","Glass_g2","Shieh_g","Hedges_gprime"))
    colnames(res) <- c("n1","n2",paste0(names[,1],names[,2]))
    res[,1:2] <- cbind(combi[,1],combi[,2])
    res[,3] <- tapply(Sel$C.P_hedges.g,list(Sel$n1,Sel$n2),mean)[1:9]
    res[,4] <- tapply(Sel$C.P_glass1.g,list(Sel$n1,Sel$n2),mean)[1:9]
    res[,5] <- tapply(Sel$C.P_glass2.g,list(Sel$n1,Sel$n2),mean)[1:9]
    res[,6] <- tapply(Sel$C.P_shieh.g,list(Sel$n1,Sel$n2),mean)[1:9]
    res[,7] <- tapply(Sel$C.P_hedges.gprime,list(Sel$n1,Sel$n2),mean)[1:9]
    # Select only rows with no "NA"  
    res <- subset(res,res[,3] != "NA") 
    # Select only bias columns
    res.CP <- t(res[,3:7])
    colnames(res.CP) <- paste0(res[,1],":",res[,2])
    
    param <- str_extract_all(list.files(Path)[j], "[[:digit:]]+\\.*[[:digit:]]*")
    if(param[[1]][2]==2.08){
      G1 <- -as.numeric(param[[1]][2])
    } else {G1 <- as.numeric(param[[1]][2])}
    G2 <- as.numeric(param[[1]][4])
    
    setwd("D:/Documents/Github_projects/Effect-sizes/Scripts outputs/Confidence intervals/Graphs/Unbiased estimators/")
    
    # plot for CP
    
    if (j==1){ylabelbias="Coverage probability"
    } else {ylabelbias=""}

    if(G1==-2.08){
      g1=2.08
    } else {g1=G1}
    
    png(file=paste0("G1=",g1, " & G2=",G2,";",names(Conditions_id)[i],".png"),width=1400,height=1700, units = "px", res = 300)  
    ylim=c(0,1)  
    mar=c(2,2,1,1)
    barplot(res.CP, 
            ylim=ylim,
            col = c("black","grey40","grey60","grey80","white"),
            beside = TRUE,
            main=paste0("G1=",G1,"; G2=",G2),
            #xaxt="n",
            cex.lab=1.5,
            cex.main=1.5,
            ylab=ylabelbias,
            cex.names=.8
    )

    abline(h=.95,lty=2,lwd=2)    
    par(xpd=TRUE)
    text(-1,-.08,labels=expression(bold(paste("n"["1"],":","n"["2"],":"))),cex=1)

    dev.off()
    
  }
  
}

#     Condition c, as a function of the SD-ratio
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Because the comparison pattern is very similar whatever n = 20, 50 or 100
# We will only plot the results as a function of the SD-ratio when n = 20

Path <-  "D:/Documents/Github_projects/Effect-sizes/Scripts outputs/Confidence intervals/Data summary/.txt files/"

for (j in seq_len(length(list.files(Path)))){

  File=read.table(paste0(Path,list.files(Path)[j]),header=T,sep=";",dec=",")  
  
  # Extract "Het_bal" subcategory 

  id_sdSD_bal=as.numeric(rownames(File[(File$n1==File$n2)&(File$sd1 < File$sd2),]))
  id_SDsd_bal=as.numeric(rownames(File[(File$n1==File$n2)&(File$sd1 > File$sd2),]))
  id_Het_bal <- c(id_sdSD_bal,id_SDsd_bal)
  
  Sel <- File[id_Het_bal,]
  
  # Possible SD-ratios
  sd2val <- as.numeric(levels(factor(Sel$sd2)))
  sd1val<- as.numeric(levels(factor(Sel$sd1)))
  combi <- expand.grid(sd1val, sd2val)
  
  # Matrix containing biases information
  res <- matrix(0,6,7)  
  names<-expand.grid("CP_",c("Hedges_g","Glass_g1","Glass_g2","Shieh_g","Hedges_gprime"))
  colnames(res) <- c("sd1","sd2",paste0(names[,1],names[,2]))
  res[,1:2] <- cbind(combi[,1],combi[,2])
  res[,3] <- tapply(Sel$C.P_hedges.g,list(Sel$sd1,Sel$sd2),mean)[1:6]
  res[,4] <- tapply(Sel$C.P_glass1.g,list(Sel$sd1,Sel$sd2),mean)[1:6]
  res[,5] <- tapply(Sel$C.P_glass2.g,list(Sel$sd1,Sel$sd2),mean)[1:6]
  res[,6] <- tapply(Sel$C.P_shieh.g,list(Sel$sd1,Sel$sd2),mean)[1:6]
  res[,7] <- tapply(Sel$C.P_hedges.gprime,list(Sel$sd1,Sel$sd2),mean)[1:6]
  # Select only rows with no "NA"  
  res <- subset(res,res[,3] != "NA") 
  # Select only bias columns
  res.CP <- t(res[,3:7])
  colnames(res.CP) <- paste0(res[,1],":",res[,2])

  param <- str_extract_all(list.files(Path)[j], "[[:digit:]]+\\.*[[:digit:]]*")
  if(param[[1]][2]==2.08){
    G1 <- -as.numeric(param[[1]][2])
  } else {G1 <- as.numeric(param[[1]][2])}
  G2 <- as.numeric(param[[1]][4])
  
  setwd("D:/Documents/Github_projects/Effect-sizes/Scripts outputs/Confidence intervals/Graphs/Unbiased estimators/")
  if(G1==-2.08){
    g1=2.08
  } else {g1=G1}

  png(file=paste0("SDR_G1=",g1, " & G2=",G2,";",names(Conditions_id)[i],".png"),width=1400,height=1700, units = "px", res = 300)  
  ylim=c(0,1)    
  mar=c(2,2,1,1)
  if (j==1){ylabelbias="Coverage probability"
  } else {ylabelbias=""}
  
  barplot(res.CP, 
          ylim=ylim,
          col = c("black","grey40","grey60","grey80","white"),
          beside = TRUE,
          main=paste0("G1=",G1,"; G2=",G2),
          #xaxt="n",
          cex.lab=1.5,
          cex.main=1.5,
          cex.names=.8,
          ylab=ylabelbias
  )
  
  abline(h=.95,lty=2,lwd=2)    
  par(xpd=TRUE)
  text(-2,-.08,labels=expression(bold(paste("\u03c3"["1"],":","\u03c3"["2"],":"))),cex=1)
  
  dev.off()
  

}

#     Condition d and e
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

# If sdSD & nN or SDsd & Nn: positive correlation between n and sd (rpos; condition d)
# If sdSD & Nn or SDsd & nN: negative correlation between n and sd (rneg; condition e)

# Note: There are too many conditions to represent them all in a single Figure. 
# We will therefore generate 3 figures per condition, as a function of the Total sample size (n1+n2)

### For a constant SD, effect of N
Path <-  "D:/Documents/Github_projects/Effect-sizes/Scripts outputs/Confidence intervals/Data summary/.txt files/"

#Ratiolarger/smaller
plot_hetr <- function(ratio){for (j in seq_len(length(list.files(Path)))){
  
  File=read.table(paste0(Path,list.files(Path)[j]),header=T,sep=";",dec=",")  
  
  # Extract "Het_unbal" subcategory (for a specific total sample size)
  
  if (ratio == 10){
    
    # Heteroscedasticity
    id_sdSD_nN=as.numeric(rownames(File[(File$n1<File$n2)&(File$sd2 == 10),]))
    id_sdSD_Nn=as.numeric(rownames(File[(File$n1>File$n2)&(File$sd2 == 10),]))
    id_SDsd_Nn=as.numeric(rownames(File[(File$n1>File$n2)&(File$sd2 == .1),]))
    id_SDsd_nN=as.numeric(rownames(File[(File$n1<File$n2)&(File$sd2 == .1),]))
    
    Conditions_id <- list(id_Het_firstsmaller=c(id_sdSD_nN,id_sdSD_Nn),id_Het_firstlarger=c(id_SDsd_nN,id_SDsd_Nn))
    
  } else if (ratio == 4){

    id_sdSD_nN=as.numeric(rownames(File[(File$n1<File$n2)&(File$sd2 == 4),]))
    id_sdSD_Nn=as.numeric(rownames(File[(File$n1>File$n2)&(File$sd2 == 4),]))
    id_SDsd_Nn=as.numeric(rownames(File[(File$n1>File$n2)&(File$sd2 == .25),]))
    id_SDsd_nN=as.numeric(rownames(File[(File$n1<File$n2)&(File$sd2 == .25),]))
    
    Conditions_id <- list(id_Het_firstsmaller=c(id_sdSD_nN,id_sdSD_Nn),id_Het_firstlarger=c(id_SDsd_nN,id_SDsd_Nn))
    
  } else if (ratio == 2){

    id_sdSD_nN=as.numeric(rownames(File[(File$n1<File$n2)&(File$sd2 == 2),]))
    id_sdSD_Nn=as.numeric(rownames(File[(File$n1>File$n2)&(File$sd2 == 2),]))
    id_SDsd_Nn=as.numeric(rownames(File[(File$n1>File$n2)&(File$sd2 == .5),]))
    id_SDsd_nN=as.numeric(rownames(File[(File$n1<File$n2)&(File$sd2 == .5),]))
    
    Conditions_id <- list(id_Het_firstsmaller=c(id_sdSD_nN,id_sdSD_Nn),id_Het_firstlarger=c(id_SDsd_nN,id_SDsd_Nn))
    
  }  
  
  
  for (i in seq_len(length(Conditions_id))){ 
    
    Sel <- File[Conditions_id[[i]],]
    
    # Possible SD-ratios and nratios
    sd2val <- as.numeric(levels(factor(Sel$sd2)))
    sd1val<- as.numeric(levels(factor(Sel$sd1)))
    n2val <- as.numeric(levels(factor(Sel$n2)))
    n1val<- as.numeric(levels(factor(Sel$n1)))
    combi <- expand.grid(n1val, n2val,sd1val, sd2val)
    
    # Matrix containing biases information
    res <- matrix(0,9,9)  
    names<-expand.grid("CP_",c("Hedges","Glass1","Glass2","Shieh","Hedges_deltaprime"))
    colnames(res) <- c("sd1","sd2","n1","n2",paste0(names[,1],names[,2]))
    res[,1:4] <- cbind(combi[,3],combi[,4],combi[,1],combi[,2])
    res[,5] <- tapply(Sel$C.P_hedges.g,list(Sel$sd1,Sel$sd2,Sel$n1,Sel$n2),mean)[1:9]
    res[,6] <- tapply(Sel$C.P_glass1.g,list(Sel$sd1,Sel$sd2,Sel$n1,Sel$n2),mean)[1:9]
    res[,7] <- tapply(Sel$C.P_glass2.g,list(Sel$sd1,Sel$sd2,Sel$n1,Sel$n2),mean)[1:9]
    res[,8] <- tapply(Sel$C.P_shieh.g,list(Sel$sd1,Sel$sd2,Sel$n1,Sel$n2),mean)[1:9]
    res[,9] <- tapply(Sel$C.P_hedges.gprime,list(Sel$sd1,Sel$sd2,Sel$n1,Sel$n2),mean)[1:9]
    # Select only rows with no "NA"  
    res <- subset(res,res[,5] != "NA") 
    # Select only bias columns
    res.CP <- t(res[,5:9])
    colnames(res.CP) <- paste0(res[,3],":",res[,4],"\n",res[,1],":",res[,2])
    
    setwd("D:/Documents/Github_projects/Effect-sizes/Scripts outputs/Confidence intervals/Graphs/Unbiased estimators/")

    param <- str_extract_all(list.files(Path)[j], "[[:digit:]]+\\.*[[:digit:]]*")
    if(param[[1]][2]==2.08){
      G1 <- -as.numeric(param[[1]][2])
    } else (G1 <- as.numeric(param[[1]][2]))
    G2 <- as.numeric(param[[1]][4])
    
    setwd(paste0("D:/Documents/Github_projects/Effect-sizes/Scripts outputs/Confidence intervals/Graphs/Unbiased estimators/",names(Conditions_id)[i]))
    if(G1==-2.08){
      g1=2.08
    } else {g1=G1}
    png(file=paste0("G1=",g1, " & G2=",G2,";",names(Conditions_id)[i],"; SDR=",ratio,".png"),width=1400,height=1700, units = "px", res = 300)  
    

    # plot for the CP
    
    if (j==1){ylabelbias="Coverage probability"
    } else {ylabelbias=""}
    
    ylim=c(0,1)    
    mar=c(2,2,1,1)
    barplot(res.CP, 
            col = c("black","grey40","grey60","grey80","white"),
            beside = TRUE,
            main=paste0("G1=",G1,"; G2=",G2),
            #xaxt="n",
            ylim=ylim,
            cex.lab=1.5,
            cex.name=.6,
            cex.main=1.5,
            ylab=ylabelbias
    )
    abline(h=.95,lty=2,lwd=2)    
    par(xpd=TRUE)
    text(-2,-.2,labels=expression(bold(paste("\u03c3"["1"],":","\u03c3"["2"],":"))),cex=.8)
    text(-2,-.1,labels=expression(bold(paste(n["1"],":",n["2"],":"))),cex=.8)
    
    dev.off()
    
  }
}
}

plot_hetr(ratio=10)
plot_hetr(ratio=4)
plot_hetr(ratio=2)