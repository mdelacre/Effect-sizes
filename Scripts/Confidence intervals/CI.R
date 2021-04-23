#Required packages

#devtools::install_github("mdelacre/deffectsize",force=T)
library(deffectsize)
library(stringr)

Mainfolder="D:/Documents/ES MEASURES/"
subfolder=list.files(Mainfolder)
Folder=paste0(Mainfolder,subfolder)

for (i in seq_len(length(Folder))){ 

  for (j in seq_len(length(list.files(Folder[i])))){ 

  filepath = paste0(Folder[i],"/",list.files(Folder[i])[j])
  file=readRDS(filepath)
  lev_biased <- c("tstat","tdf","tncp") 
  estimator_biased <- c("Bcohen","Bglass1","Bglass2","Bcohen d'","Bshieh")
  lev_unbiased <- c("est","LL","UL") 
  estimator_unbiased <- c("Bcohen","Bglass1","Bglass2","Bcohen d'","Bshieh","Uhedge","Uglass1","Uglass2","Ucohen d'","Ushieh")

  # B = biased estimators; U = unbiased estimators
  col.res <- c(do.call(paste, c(expand.grid(lev_biased,estimator_biased), sep = "_")),
               do.call(paste, c(expand.grid(lev_unbiased,estimator_unbiased), sep = "_")))
  res<-matrix(0,length(file[,1]),length(col.res)) 
  colnames(res) <- col.res
  
  # Extracting population parameters values and sample sizes from file names
  param <- str_extract_all(filepath, "[[:digit:]]+\\.*[[:digit:]]*")
  n1 <- as.numeric(param[[1]][9])
  n2 <- as.numeric(param[[1]][10])
  N <- n1+n2
  m1 <- as.numeric(param[[1]][11])
  m2 <- as.numeric(param[[1]][12])
  sd1 <- as.numeric(param[[1]][13])
  sd2 <- as.numeric(param[[1]][14])
  nratio <- n1/n2

  # Extracting sample parameters (except for n1 and n2) from file
  mean1 <- file[,1] 
  mean2 <- file[,2]
  s1 <- file[,3]
  s2 <- file[,4]

  ## Information about t-statistics and t-distributions, for each estimator
  
  ### Cohen's ds
  #### t stat
  s_pooled <- sqrt(((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-2))
  res[,1]<- (mean1-mean2)/(s_pooled*sqrt(1/n1+1/n2))
  #### df
  res[,2]<-n1+n2-2
  #### ncp
  sd_pooled <- sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
  res[,3]<- (m1-m2)/(sd_pooled*sqrt(1/n1+1/n2))
  
  ### Glass's ds using sd1 as standardizer
  
  #### t stat
  res[,4] <- ((mean1-mean2)/(s1*sqrt(1/n1+s2^2/(n2*s1^2))))
  #### df
  res[,5]<-n1-1
  #### ncp
  res[,6] <- ((m1-m2)/(sd1*sqrt(1/n1+sd2^2/(n2*sd1^2))))
  
  ### Glass's ds using sd2 as standardizer
  
  #### t stat
  res[,7] <- ((mean2-mean1)/(s2*sqrt(1/n2+s1^2/(n1*s2^2))))
  #### df
  res[,8]<-n2-1
  #### ncp
  res[,9] <- ((m2-m1)/(sd2*sqrt(1/n2+sd1^2/(n1*sd2^2))))
  
  ### Cohen's d's
  
  #### t stat
  res[,10] <- (mean1-mean2)/sqrt(s1^2/n1+s2^2/n2) #(sqrt(n1*n2)*(mean1-mean2))/sqrt(n2*s1^2+n1*s2^2)
  #### df
  res[,11]<-((n1-1)*(n2-1)*(s1^2+s2^2)^2)/((n2-1)*s1^4+(n1-1)*s2^4)
  #### ncp
  res[,12] <- (sqrt(n1*n2)*(m1-m2))/sqrt(n2*sd1^2+n1*sd2^2)
  
  ### Shieh's ds
  
  #### t stat
  res[,13] <- (mean1-mean2)/sqrt(s1^2/n1+s2^2/n2) 
  #### df
  res[,14]<-((s1^2/n1+s2^2/n2)^2)/((s1^2/n1)^2/(n1-1)+(s2^2/n2)^2/(n2-1))
  #### ncp
  res[,15] <- (m1-m2)/sqrt(sd1^2/n1+sd2^2/n2) 
  
  ## Estimates for each estimator

  ### Cohen's ds
  res[,16]<- file[,9]
  ### Glass's ds using sd1 as standardizer
  res[,19]<- file[,11]
  ### Glass's ds using sd2 as standardizer
  res[,22]<- -file[,12] # minus because I computed (m1-m2) in ES MEASURES
                        # while I use (m2-m1) when computing confidence intervals around Glass's ds using sd2 as standardizer.
  ### Cohen's d's
  res[,25]<- file[,17]
  ### Shieh's d's
  res[,28]<- file[,15]
  ### Hedges' gs
  res[,31] <- file[,10]
  ### Glass's gs using sd1 as standardizer
  res[,34] <- file[,13]
  ### Glass's gs using sd2 as standardizer
  res[,37] <- -file[,14] # same explanation as for Glass's ds using sd2 as standardizer
  ### Hedges g's
  res[,40] <- file[,18]
  ### Shieh's gs
  res[,43] <- file[,16]

  ## Lower and upper confidence limits, using deffectsize package
  
    for(k in seq_len(nrow(file))){
    CI_Bcohen <- cohen_CI(m1=mean1[k],m2=mean2[k],sd1=s1[k],sd2=s2[k],n1=n1,n2=n2,conf.level=.95,var.equal=T,unbiased=F, alternative="two.sided")
    res[k,17:18] <- CI_Bcohen$CI

    CI_Bglass1<-glass_CI(m1=mean1[k],m2=mean2[k],sd1=s1[k],sd2=s2[k],n1=n1,n2=n2,conf.level=.95,unbiased=F, alternative="two.sided")
    res[k,20:21] <- CI_Bglass1$CI

    CI_Bglass2<-glass_CI(m1=mean2[k],m2=mean1[k],sd1=s2[k],sd2=s1[k],n1=n2,n2=n1,conf.level=.95,unbiased=F, alternative="two.sided")
    res[k,23:24] <- CI_Bglass2$CI
    
    CI_Bcohenprime <- cohen_CI(m1=mean1[k],m2=mean2[k],sd1=s1[k],sd2=s2[k],n1=n1,n2=n2,conf.level=.95,var.equal=F,unbiased=F, alternative="two.sided")
    res[k,26:27] <- CI_Bcohenprime$CI
    
    CI_Bshieh<-shieh_CI(m1=mean1[k],m2=mean2[k],sd1=s1[k],sd2=s2[k],n1=n1,n2=n2,conf.level=.95,unbiased=F, alternative="two.sided")
    res[k,29:30] <- CI_Bshieh$CI
    
    CI_Ucohen <- cohen_CI(m1=mean1[k],m2=mean2[k],sd1=s1[k],sd2=s2[k],n1=n1,n2=n2,conf.level=.95,var.equal=T,unbiased=T, alternative="two.sided")
    res[k,32:33] <- CI_Ucohen$CI
    
    CI_Uglass1<-glass_CI(m1=mean1[k],m2=mean2[k],sd1=s1[k],sd2=s2[k],n1=n1,n2=n2,conf.level=.95,unbiased=T, alternative="two.sided")
    res[k,35:36] <- CI_Uglass1$CI
    
    CI_Uglass2<-glass_CI(m1=mean2[k],m2=mean1[k],sd1=s2[k],sd2=s1[k],n1=n2,n2=n1,conf.level=.95,unbiased=T, alternative="two.sided")
    res[k,38:39] <- CI_Uglass2$CI
    
    CI_Ucohenprime <- cohen_CI(m1=mean1[k],m2=mean2[k],sd1=s1[k],sd2=s2[k],n1=n1,n2=n2,conf.level=.95,var.equal=F,unbiased=T, alternative="two.sided")
    res[k,41:42] <- CI_Ucohenprime$CI
    
    CI_Ushieh<-shieh_CI(m1=mean1[k],m2=mean2[k],sd1=s1[k],sd2=s2[k],n1=n1,n2=n2,conf.level=.95,unbiased=T, alternative="two.sided")
    res[k,44:45] <- CI_Ushieh$CI
    

     }

chem <- "D:/Documents/CI.MEASURES/"
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

