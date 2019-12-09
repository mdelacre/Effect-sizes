library(stringr)

Mainfolder="C:/Users/Marie/Documents/ES MEASURES/"
subfolder=list.files(Mainfolder)
Folder=paste0(Mainfolder,subfolder)


for (i in seq_len(length(Folder))){

  # set up empty container for all estimated parameters 
  matr <-matrix(0,length(list.files(Folder[i])),24)
  indic_param <- c("m","sd", "skew","kurt")
  sample <- 1:2
  estimator <- c("Average_","Stdev_")
  columns_names=expand.grid(estimator,sample,indic_param)
  colnames(matr) <- c("n1","n2","m1","m2","sd1","sd2","skew","kurt",paste0(columns_names[,1],columns_names[,3],columns_names[,2]))

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
    skewness <- as.numeric(param[[1]][2])
    kurtosis <- as.numeric(param[[1]][4])

    # Compute bias
    
    colnames(matr) <- c("n1","n2","m1","m2","sd1","sd2","skew","kurt",paste0(columns_names[,1],columns_names[,3],columns_names[,2]))
    
    matr[j,1]<- n1
    matr[j,2]<- n2
    matr[j,3]<- m1
    matr[j,4]<- m2
    matr[j,5]<- sd1
    matr[j,6]<- sd2
    matr[j,7]<- skewness
    matr[j,8]<- kurtosis
    
    matr[j,9]<- mean(file[,1])
    matr[j,10]<- sd(file[,1])
    matr[j,11]<- mean(file[,2])
    matr[j,12]<- sd(file[,2])
    matr[j,13]<- mean(file[,3])
    matr[j,14]<- sd(file[,3])
    matr[j,15]<- mean(file[,4])
    matr[j,16]<- sd(file[,4])
    matr[j,17]<- mean(file[,5])
    matr[j,18]<- sd(file[,5])
    matr[j,19]<- mean(file[,6])
    matr[j,20]<- sd(file[,6])
    matr[j,21]<- mean(file[,7])
    matr[j,22]<- sd(file[,7])
    matr[j,23]<- mean(file[,8])
    matr[j,24]<- sd(file[,8])
    
  }
  
  setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs")
  shapeparam<-str_extract_all(Folder[i], "[[:digit:]]+\\.*[[:digit:]]*")
  sub<-paste0("G1=",shapeparam[[1]][2],",G2=",shapeparam[[1]][4]) # Sign is not extracted but it's not a big deal
  write.table(matr,paste0(sub,",SampDistr.txt"),sep=";",dec=",")
  
}

