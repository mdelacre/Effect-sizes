library(stringr)

  # set up empty container for all estimated parameters
  good_mes <-matrix(0,length(list.files(Folder)),5+3*6)
  goodness_indic <- c("bias_","eff_","MSE_")
  estimator <- c("Cohen","Hedge","Glass1","Glass2","Shieh", "Shieh_corr")
  columns_names=expand.grid(estimator,goodness_indic)
  colnames(good_mes) <- c("n1","n2","n1/n2","m1-m2","sd1/sd2",paste0(columns_names[,2],columns_names[,1]))
  
  Mainfolder="C:/Users/Marie/Documents/ES MEASURES/"
  subfolder=list.files(Mainfolder)
  Folder=paste0(Mainfolder,subfolder)
  
  for (j in seq_len(length(list.files(Folder)))){


  filepath = paste0(Folder,"/",list.files(Folder)[j])[2]
  file=readRDS(filepath)
  param <- str_extract_all(list.files(Folder)[j], "[[:digit:]]+\\.*[[:digit:]]*")
  n1 <- as.numeric(param[[1]][1])
  n2 <- as.numeric(param[[1]][2])
  m1 <- as.numeric(param[[1]][3])
  m2 <- as.numeric(param[[1]][4])
  sd1 <- as.numeric(param[[1]][5])
  sd2 <- as.numeric(param[[1]][6])
    
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
  
  setwd("C:/Users/Marie/Documents/ES MEASURES")
  write.table(good_mes,"good_mes.txt",sep=";",dec=",")
  
}
