library(stringr)

Folder="D:/ES MEASURES/G1=0,G2=0"

# 1) Considering the normality assumption, we know the theoretical distribution of different estimators

for (i in 2){   # i = 2 = le dossier contenant les simulations de distributions normales

  for (j in seq_len(length(list.files(Folder)))){ 

    filepath = paste0(Folder,"/",list.files(Folder)[j])
    file=readRDS(filepath)
    
  # Le titre des fichiers contenant les valeurs de tous les paramètres de la population, 
  # j'extrais ces différentes valeurs via str_extract_all
    param <- str_extract_all(list.files(Folder)[j], "[[:digit:]]+\\.*[[:digit:]]*")
    n1 <- as.numeric(param[[1]][5])
    n2 <- as.numeric(param[[1]][6])
    m1 <- as.numeric(param[[1]][7])
    m2 <- as.numeric(param[[1]][8])
    sd1 <- as.numeric(param[[1]][9])
    sd2 <- as.numeric(param[[1]][10])
    N <- n1+n2
    nratio <- n1/n2
    
    if(sd1==sd2){
    par(mfrow=c(2,2))

# distribution théorique attendue pour le d de Cohen
      #    Cohen_ds <- file[,9]
      #    delta=(m1-m2)/sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
      #    df=n1+n2-2 
      #    ncp=delta*sqrt((n1*n2)/(n1+n2))
      #    theo_cohen <- sqrt((n1+n2)/(n1*n2))*rt(1000000,df=df,ncp=ncp)
      #    plot(density(Cohen_ds),main=paste("mu1-mu2=",m1-m2,"\n sd-ratio=",sd1/sd2,"\n nratio=",n1/n2),xlab=paste("delta=",delta))
      #    lines(density(theo_cohen),col="blue")    
    
#distribution théorique attendue pour le d de glass, utilisant sd1 comme standardiseur

      #    delta1=(m1-m2)/sd1    
      #    glass1<-file[,11] #file[,11] = la colonne contenant les valeurs de d de glass utilisant sd1 comme standardiseur, pour chaque itération de la simulation
      #    df1=n1-1
      #    ncp1=delta1/sqrt((1/n1)+(sd2^2/(n2*sd1^2)))
      #    theo1 <- sqrt((1/n1)+(sd2^2/(n2*sd1^2)))*rt(1000000,df=df1,ncp=ncp1) 
      #    plot(density(glass1)) # distributions de glass1t = distribution du d de glass utilisant sd1 comme standardiseur, multiplié par sqrt((n1*n2)/(n1+n2))
      #    lines(density(theo1),col="blue")    

#distribution théorique attendue pour le d de glass, utilisant sd2 comme standardiseur
    
    #    delta2=(m1-m2)/sd2    
    #    glass2<-file[,12] #file[,11] = la colonne contenant les valeurs de d de glass utilisant sd1 comme standardiseur, pour chaque itération de la simulation
    #    df2=n2-1
    #    ncp2=delta2/sqrt((1/n2)+(sd1^2/(n1*sd2^2)))
    #    theo2 <- sqrt((1/n2)+(sd1^2/(n1*sd2^2)))*rt(1000000,df=df2,ncp=ncp2) 
    #    plot(density(glass2))
    #    lines(density(theo2),col="blue")    
    
# distribution théorique attendue pour le d de Shieh
    #welch_t<-file[,15]*sqrt(N) #file[,15] = la colonne contenant les valeurs de d de Shieh, pour chaque itération de la simulation
    #delta_W=(m1-m2)/sqrt(sd1^2/(n1/N)+sd2^2/(n2/N))
    #df_W=(sd1^2/n1+sd2^2/n2)^2/((sd1^2/n1)^2/(n1-1)+(sd2^2/n2)^2/(n2-1)) 
    #ncp_W=delta_W*sqrt(N)
    #theo_W <- rt(1000000,df=df_W,ncp=ncp_W)     
    #plot(density(welch_t),main=paste("mu1-mu2=",m1-m2,"\n sd-ratio=",sd1/sd2,"\n nratio=",n1/n2),xlab=paste("delta=",delta_W))
    #lines(density(theo_W),col="red")    
    
# distribution théorique attendue pour le d cohen avec variances non poolées: SOLUTION PROPOSEE DANS L ANNEXE
    cohendprime_t <- file[,17]  
    delta_t=(m1-m2)/sqrt((sd1^2+sd2^2)/2)
    df_t=((n1-1)*(n2-1)*(sd1^2+sd2^2)^2)/((n2-1)*sd1^4+(n1-1)*sd2^4)
    ncp_t=delta_t*sqrt((n1*n2*(sd1^2+sd2^2))/(2*(n2*sd1^2+n1*sd2^2)))
    theo_cohendprime <- sqrt((2*(n2*sd1^2+n1*sd2^2))/(n1*n2*(sd1^2+sd2^2)))*rt(1000000,df=df_t,ncp=ncp_t)    
    plot(density(cohendprime_t),main=paste("mu1-mu2=",m1-m2,"\n sd-ratio=",sd1/sd2,"\n nratio=",n1/n2),xlab=paste("delta=",delta_t))
    lines(density(theo_cohendprime),col="red")    

# distribution théorique attendue pour le d cohen avec variances non poolées: SOLUTION PROPOSEE DANS CORROLARY 5
    cohendprime_t <- file[,17]  
    delta_t=(m1-m2)/sqrt((sd1^2+sd2^2)/2)
    df_t=((n1-1)*(n2-1)*(sd1^2+sd2^2)^2)/((n2-1)*sd1^4+(n1-1)*sd2^4)
    ncp_t=delta_t*((n2*sd1^2+n1*sd2^2)/(n1*n2))/sqrt((sd1^2+sd2^2)/2)
    theo_cohendprime <- ((n2*sd1^2+n1*sd2^2)/(n1*n2))/sqrt((sd1^2+sd2^2)/2)*rt(1000000,df=df_t,ncp=ncp_t)
      #(sqrt(2)*(n2*sd1^2+n1*sd2^2))/(sqrt(sd1^2+sd2^2)*(n1*n2))*rt(1000000,df=df_t,ncp=ncp_t)    
    plot(density(cohendprime_t),main=paste("mu1-mu2=",m1-m2,"\n sd-ratio=",sd1/sd2,"\n nratio=",n1/n2),xlab=paste("delta=",delta_t))
    lines(density(theo_cohendprime),col="red")    
    
# TENTATIVE PERSO équivalente pour le d de Cohen avec variances non poolées
#    welch_t <- file[,17]*sqrt((n1*n2*(file[,3]^2+file[,4]^2))/(2*(n2*file[,3]^2+n1*file[,4]^2))) 
#    delta_w=(m1-m2)/sqrt((sd1^2+sd2^2)/2)
#    df_w=(sd1^2/n1+sd2^2/n2)^2/((sd1^2/n1)^2/(n1-1)+(sd2^2/n2)^2/(n2-1))#quand n1=n2, c'est égal à df_w=((n1-1)*(n2-1)*(sd1^2+sd2^2)^2)/((n1-1)*sd1^4+(n2-1)*sd2^4)
#    ncp_w=delta_w*sqrt((n1*n2*(sd1^2+sd2^2))/(2*(n2*sd1^2+n1*sd2^2)))#delta_w*(sigma_eq*sqrt(N*nratio))/((nratio+1)*sigma_uneq)
#    theo_w <- rt(1000000,df=df_w,ncp=ncp_w)    
#    plot(density(welch_t),main=paste("mu1-mu2=",m1-m2,"\n sd-ratio=",sd1/sd2,"\n nratio=",n1/n2),xlab=paste("delta=",delta_w))
#    lines(density(theo_w),col="red")    

# distribution théorique attendue pour le d cohen avec variances non poolées VERSION SANS BIAIS
#    cohendprime_unbiased <- file[,18] 
#    delta=(m1-m2)/sqrt((sd1^2+sd2^2)/2)
#    df_unbiased=((n1-1)*(n2-1)*(sd1^2+sd2^2)^2)/((n2-1)*sd1^4+(n1-1)*sd2^4)
#    ncp_unbiased=delta*sqrt((sd1^2+sd2^2)/2)/sqrt(sd1^2/n1+sd2^2/n2)
#    k <- sqrt(sd1^2/n1+sd2^2/n2)/sqrt((sd1^2+sd2^2)/2)
#    cf <- gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2))
#    theo_cohendprime_unbiased <- k*cf*rt(1000000,df=df_unbiased,ncp=ncp_unbiased)    
#    plot(density(cohendprime_unbiased),main=paste("mu1-mu2=",m1-m2,"\n sd-ratio=",sd1/sd2,"\n nratio=",n1/n2),xlab=paste("delta=",delta))
#    lines(density(theo_cohendprime_unbiased),col="red")    

} else {
      par(mfrow=c(1,2))

  
# distribution théorique attendue pour le d de Cohen
  #    Cohen_ds <- file[,9]
  #    delta=(m1-m2)/sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
  #    df=n1+n2-2 
  #    ncp=delta*sqrt((n1*n2)/(n1+n2))
  #    theo_cohen <- sqrt((n1+n2)/(n1*n2))*rt(1000000,df=df,ncp=ncp)
  #    plot(density(Cohen_ds),main=paste("mu1-mu2=",m1-m2,"\n sd-ratio=",sd1/sd2,"\n nratio=",n1/n2),xlab=paste("delta=",delta))
  #    lines(density(theo_cohen),col="blue")    

# distribution théorique attendue pour le d de glass, utilisant sd1 comme standardiseur
  
  #    delta1=(m1-m2)/sd1    
  #    glass1<-file[,11] #file[,11] = la colonne contenant les valeurs de d de glass utilisant sd1 comme standardiseur, pour chaque itération de la simulation
  #    df1=n1-1
  #    ncp1=delta1/sqrt((1/n1)+(sd2^2/(n2*sd1^2)))
  #    theo1 <- sqrt((1/n1)+(sd2^2/(n2*sd1^2)))*rt(1000000,df=df1,ncp=ncp1) 
  #    plot(density(glass1)) # distributions de glass1t = distribution du d de glass utilisant sd1 comme standardiseur, multiplié par sqrt((n1*n2)/(n1+n2))
  #    lines(density(theo1),col="blue")    

# distribution théorique attendue pour le d de glass, utilisant sd2 comme standardiseur
  
  #    delta2=(m1-m2)/sd2    
  #    glass2<-file[,12] #file[,11] = la colonne contenant les valeurs de d de glass utilisant sd1 comme standardiseur, pour chaque itération de la simulation
  #    df2=n2-1
  #    ncp2=delta2/sqrt((1/n2)+(sd1^2/(n1*sd2^2)))
  #    theo2 <- sqrt((1/n2)+(sd1^2/(n1*sd2^2)))*rt(1000000,df=df2,ncp=ncp2) 
  #    plot(density(glass2))
  #    lines(density(theo2),col="blue")    
  
# distribution théorique attendue pour le d de Shieh
#welch_t<-file[,15]*sqrt(N) #file[,15] = la colonne contenant les valeurs de d de Shieh, pour chaque itération de la simulation
#delta_W=(m1-m2)/sqrt(sd1^2/(n1/N)+sd2^2/(n2/N))
#df_W=(sd1^2/n1+sd2^2/n2)^2/((sd1^2/n1)^2/(n1-1)+(sd2^2/n2)^2/(n2-1)) 
#ncp_W=delta_W*sqrt(N)
#theo_W <- rt(1000000,df=df_W,ncp=ncp_W)     
#plot(density(welch_t),main=paste("mu1-mu2=",m1-m2,"\n sd-ratio=",sd1/sd2,"\n nratio=",n1/n2),xlab=paste("delta=",delta_W))
#lines(density(theo_W),col="blue")    

    # distribution théorique attendue pour le d cohen avec variances non poolées: SOLUTION PROPOSEE DANS L ANNEXE
    cohendprime_t <- file[,17] 
    delta_t=(m1-m2)/sqrt((sd1^2+sd2^2)/2)
    df_t=((n1-1)*(n2-1)*(sd1^2+sd2^2)^2)/((n2-1)*sd1^4+(n1-1)*sd2^4)
    ncp_t=delta_t*sqrt((n1*n2*(sd1^2+sd2^2))/(2*(n2*sd1^2+n1*sd2^2)))
    theo_cohendprime <- sqrt((2*(n2*sd1^2+n1*sd2^2))/(n1*n2*(sd1^2+sd2^2)))*rt(1000000,df=df_t,ncp=ncp_t)    
    plot(density(cohendprime_t),main=paste("mu1-mu2=",m1-m2,"\n sd-ratio=",sd1/sd2,"\n nratio=",n1/n2),xlab=paste("delta=",delta_t))
    lines(density(theo_cohendprime),col="red")    

    # distribution théorique attendue pour le d cohen avec variances non poolées: SOLUTION PROPOSEE DANS CORROLARY 5
    cohendprime_t <- file[,17]  
    delta_t=(m1-m2)/sqrt((sd1^2+sd2^2)/2)
    df_t=((n1-1)*(n2-1)*(sd1^2+sd2^2)^2)/((n2-1)*sd1^4+(n1-1)*sd2^4)
    ncp_t=delta_t*((n2*sd1^2+n1*sd2^2)/(n1*n2))/sqrt((sd1^2+sd2^2)/2)
    theo_cohendprime <- ((n2*sd1^2+n1*sd2^2)/(n1*n2))/sqrt((sd1^2+sd2^2)/2)*rt(1000000,df=df_t,ncp=ncp_t)
    plot(density(cohendprime_t),main=paste("mu1-mu2=",m1-m2,"\n sd-ratio=",sd1/sd2,"\n nratio=",n1/n2),xlab=paste("delta=",delta_t))
    lines(density(theo_cohendprime),col="red")    
    
# TENTATIVE PERSO équivalente pour le d de Cohen avec variances non poolées

#      welch_t <- file[,17]*sqrt((n1*n2*(file[,3]^2+file[,4]^2))/(2*(n2*file[,3]^2+n1*file[,4]^2))) 
#      delta_w=(m1-m2)/sqrt((sd1^2+sd2^2)/2)
#      df_w=(sd1^2/n1+sd2^2/n2)^2/((sd1^2/n1)^2/(n1-1)+(sd2^2/n2)^2/(n2-1))#quand n1=n2, c'est égal à df_w=((n1-1)*(n2-1)*(sd1^2+sd2^2)^2)/((n1-1)*sd1^4+(n2-1)*sd2^4)
#      ncp_w=delta_w*sqrt((n1*n2*(sd1^2+sd2^2))/(2*(n2*sd1^2+n1*sd2^2)))#delta_w*(sigma_eq*sqrt(N*nratio))/((nratio+1)*sigma_uneq)
#      theo_w <- rt(1000000,df=df_w,ncp=ncp_w)    
#      plot(density(welch_t),main=paste("mu1-mu2=",m1-m2,"\n sd-ratio=",sd1/sd2,"\n nratio=",n1/n2),xlab=paste("delta=",delta_w))
#      lines(density(theo_w),col="blue")    

  # distribution théorique attendue pour le d cohen avec variances non poolées VERSION SANS BIAIS
#  cohendprime_unbiased <- file[,18] 
#  delta=(m1-m2)/sqrt((sd1^2+sd2^2)/2)
#  df_unbiased=((n1-1)*(n2-1)*(sd1^2+sd2^2)^2)/((n2-1)*sd1^4+(n1-1)*sd2^4)
#  ncp_unbiased=delta*sqrt((sd1^2+sd2^2)/2)/sqrt(sd1^2/n1+sd2^2/n2)
#  k <- sqrt(sd1^2/n1+sd2^2/n2)/sqrt((sd1^2+sd2^2)/2)
#  cf <- gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2))
#  theo_cohendprime_unbiased <- k*cf*rt(1000000,df=df_unbiased,ncp=ncp_unbiased)    
#  plot(density(cohendprime_unbiased),main=paste("mu1-mu2=",m1-m2,"\n sd-ratio=",sd1/sd2,"\n nratio=",n1/n2),xlab=paste("delta=",delta))
  #  lines(density(theo_cohendprime_unbiased),col="blue")    
  
       }
  
  
   }
}

###############################################
####                                       ####
####      When distributions are normal    ####
####                                       ####
###############################################

# 2) if we know the theoretical distribution, one can compute the mean and variance

  param <- c("emp","theo") # empirical (based on simulations) vs. theoretical (based on theoretical distribution)
  lev <- c("mean","bias","relbias","var") 
  estimator <- c("cohen","hedge","glass1","glass2","unbiasedglass1","unbiasedglass2","cohen d'","unbiased cohen d'","shieh","unbiasedshieh")

  col.res <- do.call(paste, c(expand.grid(param,lev, estimator), sep = "_"))
  res<-matrix(0,length(list.files(Folder)),length(col.res)+3) # +3 for mean diff, nratio and sample ratio 
  colnames(res) <- c("m.diff","sd.ratio","n.ratio",col.res)
  
  for (j in seq_len(length(list.files(Folder)))){ 
    
    filepath = paste0(Folder,"/",list.files(Folder)[j])
    file=readRDS(filepath)
    # Le titre des fichiers contenant les valeurs de tous les paramètres de la population, 
    # j'extrais ces différentes valeurs via str_extract_all
    param <- str_extract_all(list.files(Folder)[j], "[[:digit:]]+\\.*[[:digit:]]*")
    n1 <- as.numeric(param[[1]][5])
    n2 <- as.numeric(param[[1]][6])
    N <- n1+n2
    m1 <- as.numeric(param[[1]][7])
    m2 <- as.numeric(param[[1]][8])
    sd1 <- as.numeric(param[[1]][9])
    sd2 <- as.numeric(param[[1]][10])
    nratio <- n1/n2

    ### Mean difference, sd-ratio and sample sizes ratio
    res[j,1] <- m1-m2
    res[j,2] <- sd1/sd2
    res[j,3] <- n1/n2
      
    ### 1A) Mean, bias, relative bias and variance for Cohen's ds
    
      # Mean
      res[j,4] <- mean(file[,9]) # empirical
      df=n1+n2-2 
      cohen_delta=(m1-m2)/sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/df)
      res[j,5] <- (cohen_delta*sqrt(df/2)*gamma((df-1)/2))/gamma(df/2) # theo

      # bias 
      res[j,6] <-  mean(file[,9]) - cohen_delta
      res[j,7] <-  cohen_delta*(-1+(sqrt((n1+n2-2)/2)*gamma((n1+n2-3)/2))/gamma((n1+n2-2)/2))

      #relative bias
      res[j,8] <-  (mean(file[,9]) - cohen_delta)/cohen_delta
      res[j,9] <-  (-1+(sqrt((n1+n2-2)/2)*gamma((n1+n2-3)/2))/gamma((n1+n2-2)/2))

      # Variance
      res[j,10] <-  var(file[,9])
      res[j,11] <- (N*df)/(n1*n2*(df-2))+cohen_delta^2*(df/(df-2)-(sqrt(df/2)*gamma((df-1)/2)/gamma(df/2))^2)

      ### 1B) Mean, bias and variance for Hedge's gs
      
      # Mean
      res[j,12] <- mean(file[,10]) # empirical
      df=n1+n2-2 
      res[j,13] <- cohen_delta
      
      # bias 
      res[j,14] <-  mean(file[,10]) - cohen_delta
      res[j,15] <-  0
      
      #relative bias
      res[j,16] <-  (mean(file[,10]) - cohen_delta)/cohen_delta
      res[j,17] <-  0
      
      # Variance
      res[j,18] <-  var(file[,10])
      res[j,19] <- ((N*df)/(n1*n2*(df-2))+cohen_delta^2*(df/(df-2)-(sqrt(df/2)*gamma((df-1)/2)/gamma(df/2))^2))*(gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2)))^2
        
    ### 2A) Mean, bias and variance for Glass's ds using sd1 as standardizer

    # Mean
    res[j,20] <- mean(file[,11]) # empirical
    glass_delta1= (m1-m2)/sd1
    df=n1-1 
    res[j,21] <- (glass_delta1*sqrt(df/2)*gamma((df-1)/2))/gamma(df/2) # theo
      
    # bias
    res[j,22] <-  mean(file[,11]) - glass_delta1 
    res[j,23] <-  glass_delta1*(-1+(sqrt((n1-1)/2)*gamma((n1-2)/2))/gamma((n1-1)/2))
    
    # relative bias
    res[j,24] <-  (mean(file[,11]) - glass_delta1)/glass_delta1 
    res[j,25] <-  (-1+(sqrt((n1-1)/2)*gamma((n1-2)/2))/gamma((n1-1)/2))
    
    # Variance
    res[j,26] <-  var(file[,11])
    res[j,27] <-  ((n1-1)/(n1-3))*(1/n1+sd2^2/(n2*sd1^2)+glass_delta1^2)-glass_delta1^2*((sqrt((n1-1)/2)*gamma((n1-2)/2))/gamma((n1-1)/2))^2       
      
    ### 3A) Mean, bias and variance for Glass's ds using sd2 as standardizer        
      
    # Mean
    res[j,28] <- mean(file[,12]) # empirical
    glass_delta2= (m1-m2)/sd2
    df=n2-1 
    res[j,29] <- (glass_delta2*sqrt(df/2)*gamma((df-1)/2))/gamma(df/2) # theo
    
    # bias
    res[j,30] <-  mean(file[,12]) - glass_delta2
    res[j,31] <-  glass_delta2*(-1+(sqrt((n2-1)/2)*gamma((n2-2)/2))/gamma((n2-1)/2))

    # relative bias
    res[j,32] <-  (mean(file[,12]) - glass_delta2)/glass_delta2
    res[j,33] <-  (-1+(sqrt((n2-1)/2)*gamma((n2-2)/2))/gamma((n2-1)/2))
    
    # Variance
    res[j,34] <-  var(file[,12])
    res[j,35] <-  ((n2-1)/(n2-3))*(1/n2+sd1^2/(n1*sd2^2)+glass_delta2^2)-glass_delta2^2*((sqrt((n2-1)/2)*gamma((n2-2)/2))/gamma((n2-1)/2))^2   

    ### 2B) Mean, bias and variance for unbiased Glass's ds using sd1 as standardizer
    
    # Mean
    res[j,36] <- mean(file[,13]) # empirical
    glass_delta1= (m1-m2)/sd1
    df=n1-1 
    res[j,37] <- glass_delta1
    
    # bias
    res[j,38] <-  mean(file[,13]) - glass_delta1 
    res[j,39] <-  0
    
    # relative bias
    res[j,40] <-  (mean(file[,13]) - glass_delta1)/glass_delta1 
    res[j,41] <-  0
    
    # Variance
    res[j,42] <-  var(file[,13])
    res[j,43] <-  (((n1-1)/(n1-3))*(1/n1+sd2^2/(n2*sd1^2)+glass_delta1^2)-glass_delta1^2*((sqrt((n1-1)/2)*gamma((n1-2)/2))/gamma((n1-1)/2))^2)*(gamma((n1-1)/2)/(sqrt((n1-1)/2)*gamma((n1-2)/2)))^2
    
    ### 3B) Mean, bias and variance for unbiased Glass's ds using sd2 as standardizer        
    
    # Mean
    res[j,44] <- mean(file[,14]) # empirical
    glass_delta2= (m1-m2)/sd2
    df=n2-1 
    res[j,45] <- glass_delta2 
    
    # bias
    res[j,46] <-  mean(file[,14]) - glass_delta2
    res[j,47] <-  0
    
    # relative bias
    res[j,48] <-  (mean(file[,14]) - glass_delta2)/glass_delta2
    res[j,49] <-  0
    
    # Variance
    res[j,50] <-  var(file[,14])
    res[j,51] <-  (((n2-1)/(n2-3))*(1/n2+sd1^2/(n1*sd2^2)+glass_delta2^2)-glass_delta2^2*((sqrt((n2-1)/2)*gamma((n2-2)/2))/gamma((n2-1)/2))^2)*(gamma((n2-1)/2)/(sqrt((n2-1)/2)*gamma((n2-2)/2)))^2
    
    ### 4A) Mean, bias, relative bias and variance for Cohen's d's
    # Mean 
    res[j,52] <- mean(file[,17])
    cohen_deltaprime <- (m1-m2)/sqrt((sd1^2+sd2^2)/2)    
    df <- ((n1-1)*(n2-1)*(sd1^2+sd2^2)^2)/((n2-1)*sd1^4+(n1-1)*sd2^4)
    
    res[j,53]   <- cohen_deltaprime*(sqrt(df/2)*gamma((df-1)/2))/gamma(df/2)
    #cohen_deltaprime*(4*df-1)/(4*(df-1))
    
    # Rem.: l'espérance de (file[,17])*sqrt((file[,3]^2+file[,4]^2)/(n2*file[,3]^2+n1*file[,4]^2))
    # est égale à cohen_deltaprime*sqrt((sd1^2+sd2^2)/(n2*sd1^2+n1*sd2^2))*((sqrt(df/2)*gamma((df-1)/2))/gamma(df/2))
    # avec df = (sd1^2/n1+sd2^2/n2)^2/((sd1^2/n1)^2/(n1-1)+(sd2^2/n2)^2/(n2-1))
    # Déduit sur base de l'équation de l'espérance de Shieh, et du lien entre Shieh's ds et Cohen's d's...
    
    # bias 
    res[j,54] <-  mean(file[,17]) - cohen_deltaprime
    res[j,55] <-  cohen_deltaprime*((sqrt(df/2)*gamma((df-1)/2))/gamma(df/2)-1)
    #cohen_deltaprime*((4*df-1)/(4*(df-1))-1)
    
    #relative bias
    res[j,56] <-  (mean(file[,17]) - cohen_deltaprime)/cohen_deltaprime
    res[j,57] <- (cohen_deltaprime*((sqrt(df/2)*gamma((df-1)/2))/gamma(df/2)-1))/cohen_deltaprime 
    # (cohen_deltaprime*((sqrt(df/2)*gamma((df-1)/2))/gamma(df/2)-1))/cohen_deltaprime
    
    # Variance
    res[j,58] <-  var(file[,17])
    res[j,59] <- (df/(df-2))*(2*(sd1^2/n1+sd2^2/n2)/(sd1^2+sd2^2)+cohen_deltaprime^2)-cohen_deltaprime^2*(sqrt(df/2)*gamma((df-1)/2)/gamma(df/2))^2
    

    ### 4B) Mean, bias, relative bias and variance for unbiased Cohen's d's
    # Mean 
    res[j,60] <- mean(file[,18])
    cohen_deltaprime <- (m1-m2)/sqrt((sd1^2+sd2^2)/2)    
    df <- ((n1-1)*(n2-1)*(sd1^2+sd2^2)^2)/((n2-1)*sd1^4+(n1-1)*sd2^4)
    
    res[j,61]   <- cohen_deltaprime
    
    # bias 
    res[j,62] <-  mean(file[,18]) - cohen_deltaprime
    res[j,63] <-  0
    
    #relative bias
    res[j,64] <-  (mean(file[,18]) - cohen_deltaprime)/cohen_deltaprime
    res[j,65] <- 0
    
    # Variance
    res[j,66] <-  var(file[,18])
    
    k <- sqrt(sd1^2/n1+sd2^2/n2)/sqrt((sd1^2+sd2^2)/2)
    cf <- gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2))
    
    
    res[j,67] <- ((df/(df-2))*(2*(sd1^2/n1+sd2^2/n2)/(sd1^2+sd2^2)+cohen_deltaprime^2)-cohen_deltaprime^2*(sqrt(df/2)*gamma((df-1)/2)/gamma(df/2))^2)*(gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2)))^2
      #(df/(df-2))*cf^2*(2*(sd1^2/n1+sd2^2/n2)/(sd1^2+sd2^2))+cohen_deltaprime^2*(df/(df-2)*cf^2-1)

    ### 5A) Mean, bias and variance for Shieh's ds
    # Mean
    
    res[j,68] <- mean(file[,15]) # empirical
    q1 <- n1/(n1+n2)
    q2 <- n2/(n1+n2)
    N <- n1+n2
    shieh_delta <- (m1-m2)/sqrt(sd1^2/q1+sd2^2/q2)    
    df <- (sd1^2/n1+sd2^2/n2)^2/((sd1^2/n1)^2/(n1-1)+(sd2^2/n2)^2/(n2-1))
    res[j,69]   <- (shieh_delta*sqrt(df/2)*gamma((df-1)/2))/gamma(df/2) # theo
    # = (cohen_deltaprime*sqrt((n1*n2*(sd1^2+sd2^2))/(2*N*(n2*sd1^2+n1*sd2^2)))*sqrt(df/2)*gamma((df-1)/2))/gamma(df/2) # theo
    
    # bias 
    res[j,70] <-  mean(file[,15]) - shieh_delta
    res[j,71] <-  shieh_delta*((sqrt(df/2)*gamma((df-1)/2))/gamma(df/2)-1) # theo

    # relative bias 
    res[j,72] <-  (mean(file[,15]) - shieh_delta)/shieh_delta
    res[j,73] <-  ((sqrt(df/2)*gamma((df-1)/2))/gamma(df/2)-1) # theo
    
    # Variance
    res[j,74] <-  var(file[,15])
    res[j,75] <- df/((df-2)*N)*(1+shieh_delta^2*N)-shieh_delta^2*((sqrt(df/2)*gamma((df-1)/2))/gamma(df/2))^2
      
    #  Note: when n1=n2 and sd1=sd2 simultaneously , var(shieh's d)= (n1+n2-2)/(4*(n1+n2-4)*n1*n2/(n1+n2))*(1+(n1*n2/(n1+n2))*4*shieh_delta^2)-shieh_delta^2*((sqrt((n1+n2-2)/2)*gamma(((n1+n2-2)-1)/2))/gamma((n1+n2-2)/2))^2

    ### 5B) Mean, bias and variance for unbiased Shieh's ds
    # Mean
    
    res[j,76] <- mean(file[,16]) # empirical
    q1 <- n1/(n1+n2)
    q2 <- n2/(n1+n2)
    N <- n1+n2
    shieh_delta <- (m1-m2)/sqrt(sd1^2/q1+sd2^2/q2)    
    df <- (sd1^2/n1+sd2^2/n2)^2/((sd1^2/n1)^2/(n1-1)+(sd2^2/n2)^2/(n2-1))
    res[j,77]   <- shieh_delta # theo

    # bias 
    res[j,78] <-  mean(file[,16]) - shieh_delta
    res[j,79] <-  0
    
    # relative bias 
    res[j,80] <-  (mean(file[,16]) - shieh_delta)/shieh_delta
    res[j,81] <-  0
    
    # Variance
    res[j,82] <-  var(file[,16])
    res[j,83] <- (df/((df-2)*N)*(1+shieh_delta^2*N)-shieh_delta^2*((sqrt(df/2)*gamma((df-1)/2))/gamma(df/2))^2)*(gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2)))^2
      #df/((df-2)*N)*(gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2)))^2+shieh_delta^2*(df/(df-2)*(gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2)))^2-1)
    

      }
  

homosc_eqn <- subset(res,res[,2] == 1 & res[,3]==1) # designs balancés, homoscédasticité
homosc_uneqn <- subset(res,res[,2] == 1 & res[,3]!=1) # designs non balancés, homoscédasticité
heterosc_eqn <- subset(res,res[,2] != 1 & res[,3]==1) # designs balancés, hétéroscédasticité   
heterosc_uneqn <- subset(res,res[,2] != 1 & res[,3]!=1) # designs non balancés, hétéroscédasticité   

nonnull_homosc_eqn <- subset(res,res[,2] == 1 & res[,3]==1 & res[,1]!=0) # designs balancés, homoscédasticité
nonnull_homosc_uneqn <- subset(res,res[,2] == 1 & res[,3]!=1 & res[,1]!=0) # designs non balancés, homoscédasticité
nonnull_heterosc_eqn <- subset(res,res[,2] != 1 & res[,3]==1 & res[,1]!=0) # designs balancés, hétéroscédasticité   
nonnull_heterosc_uneqn <- subset(res,res[,2] != 1 & res[,3]!=1 & res[,1]!=0) # designs non balancés, hétéroscédasticité   


##################### Cohen ##################### 

### Mean
# equal variances and sample sizes (15 scenarios)
round(max(abs(homosc_eqn[,4]-homosc_eqn[,5])),3) # .012 (je trouve .05 en utilisant l'estimation de la moyenne pour tous les cas où n1=n2)
round(sum(abs(homosc_eqn[,4]-homosc_eqn[,5])),3) # .036 (je trouve .181 en utilisant l'estimation de la moyenne pour tous les cas où n1=n2) 
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(homosc_uneqn[,4]-homosc_uneqn[,5])),3) # .005
round(sum(abs(homosc_uneqn[,4]-homosc_uneqn[,5])),3) # .038
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(heterosc_eqn[,4]-heterosc_eqn[,5])),3) # .08 (je trouve .03 en utilisant l'estimation de la moyenne pour tous les cas où n1=n2)
round(sum(abs(heterosc_eqn[,4]-heterosc_eqn[,5])),3) # .885 (je trouve .404 en utilisant l'estimation de la moyenne pour tous les cas où n1=n2)
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(heterosc_uneqn[,4]-heterosc_uneqn[,5])),3) # .252 (rem.: j'obtiens un meilleur résultat quand j'utilise la formule théorique pour les cas où sd1!=sd2 mais n1=n2... enfin moins pourri que l'équation de base :) 
round(sum(abs(heterosc_uneqn[,4]-heterosc_uneqn[,5])),3) # 2.654

### Bias
# equal variances and sample sizes (15 scenarios)
round(max(abs(homosc_eqn[,6]-homosc_eqn[,7])),3) # .012
round(sum(abs(homosc_eqn[,6]-homosc_eqn[,7])),3) # .036
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(homosc_uneqn[,6]-homosc_uneqn[,7])),3) # .005
round(sum(abs(homosc_uneqn[,6]-homosc_uneqn[,7])),3) # .038
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(heterosc_eqn[,6]-heterosc_eqn[,7])),3) # .08
round(sum(abs(heterosc_eqn[,6]-heterosc_eqn[,7])),3) # .885
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(heterosc_uneqn[,6]-heterosc_uneqn[,7])),3) # .252
round(sum(abs(heterosc_uneqn[,6]-heterosc_uneqn[,7])),3) # 2.654

### Relative Bias
# equal variances and sample sizes (15 scenarios)
round(max(abs(nonnull_homosc_eqn[,8]-nonnull_homosc_eqn[,9])),3) # .004
round(sum(abs(nonnull_homosc_eqn[,8]-nonnull_homosc_eqn[,9])),3) # .015
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(nonnull_homosc_uneqn[,8]-nonnull_homosc_uneqn[,9])),3) # .002
round(sum(abs(nonnull_homosc_uneqn[,8]-nonnull_homosc_uneqn[,9])),3) # .017
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(nonnull_heterosc_eqn[,8]-nonnull_heterosc_eqn[,9])),3) # .018
round(sum(abs(nonnull_heterosc_eqn[,8]-nonnull_heterosc_eqn[,9])),3) # .42
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(nonnull_heterosc_uneqn[,8]-nonnull_heterosc_uneqn[,9])),3) # .028
round(sum(abs(nonnull_heterosc_uneqn[,8]-nonnull_heterosc_uneqn[,9])),3) # .963

### Variance
# equal variances and sample sizes (15 scenarios)
round(max(homosc_eqn[,10]/homosc_eqn[,11]),3) # 1.006
round(min(homosc_eqn[,10]/homosc_eqn[,11]),3) # .91
# equal variances, unequal sample sizes (30 scenarios) 
round(max(homosc_uneqn[,10]/homosc_uneqn[,11]),3) # 1.017
round(min(homosc_uneqn[,10]/homosc_uneqn[,11]),3) # .951
# unequal variances, equal sample sizes (90 scenarios) 
round(max(heterosc_eqn[,10]/heterosc_eqn[,11]),3) # 1.753
round(min(heterosc_eqn[,10]/heterosc_eqn[,11]),3) # 1.005
# unequal variances and sample sizes (90 scenarios) 
round(max(heterosc_uneqn[,10]/heterosc_uneqn[,11]),3) # 5.624
round(min(heterosc_uneqn[,10]/heterosc_uneqn[,11]),3) # 0.208

##################### Hedges ##################### 

### Mean
# equal variances and sample sizes (15 scenarios)
round(max(abs(homosc_eqn[,12]-homosc_eqn[,13])),3) # .011
round(sum(abs(homosc_eqn[,12]-homosc_eqn[,13])),3) # .035 
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(homosc_uneqn[,12]-homosc_uneqn[,13])),3) # .005
round(sum(abs(homosc_uneqn[,12]-homosc_uneqn[,13])),3) # .037
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(heterosc_eqn[,12]-heterosc_eqn[,13])),3) # .079
round(sum(abs(heterosc_eqn[,12]-heterosc_eqn[,13])),3) # .873
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(heterosc_uneqn[,12]-heterosc_uneqn[,13])),3) # .25
round(sum(abs(heterosc_uneqn[,12]-heterosc_uneqn[,13])),3) # 2.634

### Bias
# equal variances and sample sizes (15 scenarios)
round(max(abs(homosc_eqn[,14]-homosc_eqn[,15])),3) # .011
round(sum(abs(homosc_eqn[,14]-homosc_eqn[,15])),3) # .035
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(homosc_uneqn[,14]-homosc_uneqn[,15])),3) # .005
round(sum(abs(homosc_uneqn[,14]-homosc_uneqn[,15])),3) # .037
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(heterosc_eqn[,14]-heterosc_eqn[,15])),3) # .079
round(sum(abs(heterosc_eqn[,14]-heterosc_eqn[,15])),3) # .873
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(heterosc_uneqn[,14]-heterosc_uneqn[,15])),3) # .25
round(sum(abs(heterosc_uneqn[,14]-heterosc_uneqn[,15])),3) # 2.634

### relative bias
# equal variances and sample sizes (15 scenarios)
round(max(abs(nonnull_homosc_eqn[,16]-nonnull_homosc_eqn[,17])),3) # .004
round(sum(abs(nonnull_homosc_eqn[,16]-nonnull_homosc_eqn[,17])),3) # .014
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(nonnull_homosc_uneqn[,16]-nonnull_homosc_uneqn[,17])),3) # .002
round(sum(abs(nonnull_homosc_uneqn[,16]-nonnull_homosc_uneqn[,17])),3) # .016
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(nonnull_heterosc_eqn[,16]-nonnull_heterosc_eqn[,17])),3) # .017
round(sum(abs(nonnull_heterosc_eqn[,16]-nonnull_heterosc_eqn[,17])),3) # .414
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(nonnull_heterosc_uneqn[,16]-nonnull_heterosc_uneqn[,17])),3) # .027
round(sum(abs(nonnull_heterosc_uneqn[,16]-nonnull_heterosc_uneqn[,17])),3) # .956

### variance
# equal variances and sample sizes (15 scenarios)
round(max(homosc_eqn[,18]/homosc_eqn[,19]),3) # 1.006
round(min(homosc_eqn[,18]/homosc_eqn[,19]),3) # .91 
# equal variances, unequal sample sizes (30 scenarios) 
round(max(homosc_uneqn[,18]/homosc_uneqn[,19]),3) # 1.017
round(min(homosc_uneqn[,18]/homosc_uneqn[,19]),3) # .951
# unequal variances, equal sample sizes (90 scenarios) 
round(max(heterosc_eqn[,18]/heterosc_eqn[,19]),3) # 1.753
round(min(heterosc_eqn[,18]/heterosc_eqn[,19]),3) # 1.005
# unequal variances and sample sizes (90 scenarios) 
round(max(heterosc_uneqn[,18]/heterosc_uneqn[,19]),3) # 5.624
round(min(heterosc_uneqn[,18]/heterosc_uneqn[,19]),3) # .208

##################### Glass's ds (with sd1) ##################### 

### Mean

# equal variances and sample sizes (15 scenarios)
round(max(abs(homosc_eqn[,20]-homosc_eqn[,21])),3) # .022
round(sum(abs(homosc_eqn[,20]-homosc_eqn[,21])),3) # .064
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(homosc_uneqn[,20]-homosc_uneqn[,21])),3) # .019
round(sum(abs(homosc_uneqn[,20]-homosc_uneqn[,21])),3) # .127
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(heterosc_eqn[,20]-heterosc_eqn[,21])),3) # .037
round(sum(abs(heterosc_eqn[,20]-heterosc_eqn[,21])),3) # .455
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(heterosc_uneqn[,20]-heterosc_uneqn[,21])),3) # .026
round(sum(abs(heterosc_uneqn[,20]-heterosc_uneqn[,21])),3) # .835

### Bias
# equal variances and sample sizes (15 scenarios)
round(max(abs(homosc_eqn[,22]-homosc_eqn[,23])),3) # .022
round(sum(abs(homosc_eqn[,22]-homosc_eqn[,23])),3) # .064
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(homosc_uneqn[,22]-homosc_uneqn[,23])),3) # .019
round(sum(abs(homosc_uneqn[,22]-homosc_uneqn[,23])),3) # .127
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(heterosc_eqn[,22]-heterosc_eqn[,23])),3) # .037
round(sum(abs(heterosc_eqn[,22]-heterosc_eqn[,23])),3) # .455
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(heterosc_uneqn[,22]-heterosc_uneqn[,23])),3) # .026
round(sum(abs(heterosc_uneqn[,22]-heterosc_uneqn[,23])),3) # .835

### relative bias
# equal variances and sample sizes (15 scenarios)
round(max(abs(nonnull_homosc_eqn[,24]-nonnull_homosc_eqn[,25])),3) # .006
round(sum(abs(nonnull_homosc_eqn[,24]-nonnull_homosc_eqn[,25])),3) # .025
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(nonnull_homosc_uneqn[,24]-nonnull_homosc_uneqn[,25])),3) # .007
round(sum(abs(nonnull_homosc_uneqn[,24]-nonnull_homosc_uneqn[,25])),3) # .053
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(nonnull_heterosc_eqn[,24]-nonnull_heterosc_eqn[,25])),3) # .012
round(sum(abs(nonnull_heterosc_eqn[,24]-nonnull_heterosc_eqn[,25])),3) # .168
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(nonnull_heterosc_uneqn[,24]-nonnull_heterosc_uneqn[,25])),3) # .009
round(sum(abs(nonnull_heterosc_uneqn[,24]-nonnull_heterosc_uneqn[,25])),3) # .316

### Variance
# equal variances and sample sizes (15 scenarios)
round(max(homosc_eqn[,26]/homosc_eqn[,27]),3) # 1.006
round(min(homosc_eqn[,26]/homosc_eqn[,27]),3) # .897
# equal variances, unequal sample sizes (30 scenarios) 
round(max(homosc_uneqn[,26]/homosc_uneqn[,27]),3) # 1.006
round(min(homosc_uneqn[,26]/homosc_uneqn[,27]),3) # .891
# unequal variances, equal sample sizes (90 scenarios) 
round(max(heterosc_eqn[,26]/heterosc_eqn[,27]),3) # 1.004
round(min(heterosc_eqn[,26]/heterosc_eqn[,27]),3) # .888
# unequal variances and sample sizes (90 scenarios) 
round(max(heterosc_uneqn[,26]/heterosc_uneqn[,27]),3) # 1.009
round(min(heterosc_uneqn[,26]/heterosc_uneqn[,27]),3) # .881

##################### Glass's ds (with sd2) ##################### 

### Mean
# equal variances and sample sizes (15 scenarios)
round(max(abs(homosc_eqn[,28]-homosc_eqn[,29])),3) # .023
round(sum(abs(homosc_eqn[,28]-homosc_eqn[,29])),3) # .069
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(homosc_uneqn[,28]-homosc_uneqn[,29])),3) # .027
round(sum(abs(homosc_uneqn[,28]-homosc_uneqn[,29])),3) # .137
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(heterosc_eqn[,28]-heterosc_eqn[,29])),3) # .23
round(sum(abs(heterosc_eqn[,28]-heterosc_eqn[,29])),3) # 1.083
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(heterosc_uneqn[,28]-heterosc_uneqn[,29])),3) # .219
round(sum(abs(heterosc_uneqn[,28]-heterosc_uneqn[,29])),3) # 2.191

### bias
# equal variances and sample sizes (15 scenarios)
round(max(abs(homosc_eqn[,30]-homosc_eqn[,31])),3) # .023
round(sum(abs(homosc_eqn[,30]-homosc_eqn[,31])),3) # .069 
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(homosc_uneqn[,30]-homosc_uneqn[,31])),3) # .027
round(sum(abs(homosc_uneqn[,30]-homosc_uneqn[,31])),3) # .137
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(heterosc_eqn[,30]-heterosc_eqn[,31])),3) # .23
round(sum(abs(heterosc_eqn[,30]-heterosc_eqn[,31])),3) # 1.083
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(heterosc_uneqn[,30]-heterosc_uneqn[,31])),3) # .219
round(sum(abs(heterosc_uneqn[,30]-heterosc_uneqn[,31])),3) # 2.191

### relative bias
# equal variances and sample sizes (15 scenarios)
round(max(abs(nonnull_homosc_eqn[,32]-nonnull_homosc_eqn[,33])),3) # .007
round(sum(abs(nonnull_homosc_eqn[,32]-nonnull_homosc_eqn[,33])),3) # .028 
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(nonnull_homosc_uneqn[,32]-nonnull_homosc_uneqn[,33])),3) # .007
round(sum(abs(nonnull_homosc_uneqn[,32]-nonnull_homosc_uneqn[,33])),3) # .056
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(nonnull_heterosc_eqn[,32]-nonnull_heterosc_eqn[,33])),3) # .012
round(sum(abs(nonnull_heterosc_eqn[,32]-nonnull_heterosc_eqn[,33])),3) # .168
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(nonnull_heterosc_uneqn[,32]-nonnull_heterosc_uneqn[,33])),3) # .009
round(sum(abs(nonnull_heterosc_uneqn[,32]-nonnull_heterosc_uneqn[,33])),3) # .313

### Variance
# equal variances and sample sizes (15 scenarios)
round(max(homosc_eqn[,34]/homosc_eqn[,35]),3) # 1.005
round(min(homosc_eqn[,34]/homosc_eqn[,35]),3) # .889
# equal variances, unequal sample sizes (30 scenarios) 
round(max(homosc_uneqn[,34]/homosc_uneqn[,35]),3) # 1.015
round(min(homosc_uneqn[,34]/homosc_uneqn[,35]),3) # .881
# unequal variances, equal sample sizes (90 scenarios) 
round(max(heterosc_eqn[,34]/heterosc_eqn[,35]),3) # 1.008
round(min(heterosc_eqn[,34]/heterosc_eqn[,35]),3) # .883
# unequal variances and sample sizes (90 scenarios) 
round(max(heterosc_uneqn[,34]/heterosc_uneqn[,35]),3) # 1.011
round(min(heterosc_uneqn[,34]/heterosc_uneqn[,35]),3) # .872

##################### Unbiased Glass's ds (with sd1) ##################### 

### Mean

# equal variances and sample sizes (15 scenarios)
round(max(abs(homosc_eqn[,36]-homosc_eqn[,37])),3) # .021
round(sum(abs(homosc_eqn[,36]-homosc_eqn[,37])),3) # .062
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(homosc_uneqn[,36]-homosc_uneqn[,37])),3) # .018
round(sum(abs(homosc_uneqn[,36]-homosc_uneqn[,37])),3) # .123
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(heterosc_eqn[,36]-heterosc_eqn[,37])),3) # .036
round(sum(abs(heterosc_eqn[,36]-heterosc_eqn[,37])),3) # .44
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(heterosc_uneqn[,36]-heterosc_uneqn[,37])),3) # .025
round(sum(abs(heterosc_uneqn[,36]-heterosc_uneqn[,37])),3) # .807

### Bias
# equal variances and sample sizes (15 scenarios)
round(max(abs(homosc_eqn[,38]-homosc_eqn[,39])),3) # .021
round(sum(abs(homosc_eqn[,38]-homosc_eqn[,39])),3) # .062
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(homosc_uneqn[,38]-homosc_uneqn[,39])),3) # .018
round(sum(abs(homosc_uneqn[,38]-homosc_uneqn[,39])),3) # .123
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(heterosc_eqn[,38]-heterosc_eqn[,39])),3) # .036
round(sum(abs(heterosc_eqn[,38]-heterosc_eqn[,39])),3) # .44
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(heterosc_uneqn[,38]-heterosc_uneqn[,39])),3) # .025
round(sum(abs(heterosc_uneqn[,38]-heterosc_uneqn[,39])),3) # .807

### relative bias
# equal variances and sample sizes (15 scenarios)
round(max(abs(nonnull_homosc_eqn[,40]-nonnull_homosc_eqn[,41])),3) # .006
round(sum(abs(nonnull_homosc_eqn[,40]-nonnull_homosc_eqn[,41])),3) # .024
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(nonnull_homosc_uneqn[,40]-nonnull_homosc_uneqn[,41])),3) # .006
round(sum(abs(nonnull_homosc_uneqn[,40]-nonnull_homosc_uneqn[,41])),3) # .052
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(nonnull_heterosc_eqn[,40]-nonnull_heterosc_eqn[,41])),3) # .012
round(sum(abs(nonnull_heterosc_eqn[,40]-nonnull_heterosc_eqn[,41])),3) # .163
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(nonnull_heterosc_uneqn[,40]-nonnull_heterosc_uneqn[,41])),3) # .009
round(sum(abs(nonnull_heterosc_uneqn[,40]-nonnull_heterosc_uneqn[,41])),3) # .305

### Variance
# equal variances and sample sizes (15 scenarios)
round(max(homosc_eqn[,42]/homosc_eqn[,43]),3) # 1.006
round(min(homosc_eqn[,42]/homosc_eqn[,43]),3) # .897
# equal variances, unequal sample sizes (30 scenarios) 
round(max(homosc_uneqn[,42]/homosc_uneqn[,43]),3) # 1.006
round(min(homosc_uneqn[,42]/homosc_uneqn[,43]),3) # .891
# unequal variances, equal sample sizes (90 scenarios) 
round(max(heterosc_eqn[,42]/heterosc_eqn[,43]),3) # 1.004
round(min(heterosc_eqn[,42]/heterosc_eqn[,43]),3) # .888
# unequal variances and sample sizes (90 scenarios) 
round(max(heterosc_uneqn[,42]/heterosc_uneqn[,43]),3) # 1.009
round(min(heterosc_uneqn[,42]/heterosc_uneqn[,43]),3) # .881

##################### Unbiased glass's ds (with sd2) ##################### 

### Mean
# equal variances and sample sizes (15 scenarios)
round(max(abs(homosc_eqn[,44]-homosc_eqn[,45])),3) # .022
round(sum(abs(homosc_eqn[,44]-homosc_eqn[,45])),3) # .066
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(homosc_uneqn[,44]-homosc_uneqn[,45])),3) # .026
round(sum(abs(homosc_uneqn[,44]-homosc_uneqn[,45])),3) # .133
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(heterosc_eqn[,44]-heterosc_eqn[,45])),3) # .221
round(sum(abs(heterosc_eqn[,44]-heterosc_eqn[,45])),3) # 1.044
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(heterosc_uneqn[,44]-heterosc_uneqn[,45])),3) # .21
round(sum(abs(heterosc_uneqn[,44]-heterosc_uneqn[,45])),3) # 2.115

### bias
# equal variances and sample sizes (15 scenarios)
round(max(abs(homosc_eqn[,46]-homosc_eqn[,47])),3) # .022
round(sum(abs(homosc_eqn[,46]-homosc_eqn[,47])),3) # .066 
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(homosc_uneqn[,46]-homosc_uneqn[,47])),3) # .026
round(sum(abs(homosc_uneqn[,46]-homosc_uneqn[,47])),3) # .133
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(heterosc_eqn[,46]-heterosc_eqn[,47])),3) # .221
round(sum(abs(heterosc_eqn[,46]-heterosc_eqn[,47])),3) # 1.044
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(heterosc_uneqn[,46]-heterosc_uneqn[,47])),3) # .21
round(sum(abs(heterosc_uneqn[,46]-heterosc_uneqn[,47])),3) # 2.115

### relative bias
# equal variances and sample sizes (15 scenarios)
round(max(abs(nonnull_homosc_eqn[,48]-nonnull_homosc_eqn[,49])),3) # .006
round(sum(abs(nonnull_homosc_eqn[,48]-nonnull_homosc_eqn[,49])),3) # .027 
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(nonnull_homosc_uneqn[,48]-nonnull_homosc_uneqn[,49])),3) # .007
round(sum(abs(nonnull_homosc_uneqn[,48]-nonnull_homosc_uneqn[,49])),3) # .054
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(nonnull_heterosc_eqn[,48]-nonnull_heterosc_eqn[,49])),3) # .011
round(sum(abs(nonnull_heterosc_eqn[,48]-nonnull_heterosc_eqn[,49])),3) # .163
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(nonnull_heterosc_uneqn[,48]-nonnull_heterosc_uneqn[,49])),3) # .009
round(sum(abs(nonnull_heterosc_uneqn[,48]-nonnull_heterosc_uneqn[,49])),3) # .302

### Variance
# equal variances and sample sizes (15 scenarios)
round(max(homosc_eqn[,50]/homosc_eqn[,51]),3) # 1.005
round(min(homosc_eqn[,50]/homosc_eqn[,51]),3) # .889
# equal variances, unequal sample sizes (30 scenarios) 
round(max(homosc_uneqn[,50]/homosc_uneqn[,51]),3) # 1.015
round(min(homosc_uneqn[,50]/homosc_uneqn[,51]),3) # .881
# unequal variances, equal sample sizes (90 scenarios) 
round(max(heterosc_eqn[,50]/heterosc_eqn[,51]),3) # 1.008
round(min(heterosc_eqn[,50]/heterosc_eqn[,51]),3) # .883
# unequal variances and sample sizes (90 scenarios) 
round(max(heterosc_uneqn[,50]/heterosc_uneqn[,51]),3) # 1.011
round(min(heterosc_uneqn[,50]/heterosc_uneqn[,51]),3) # .872

##################### Cohen's d's ##################### 

### Mean
# equal variances and sample sizes (15 scenarios)
round(max(abs(homosc_eqn[,52]-homosc_eqn[,53])),3) # .012
round(sum(abs(homosc_eqn[,52]-homosc_eqn[,53])),3) # .036 
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(homosc_uneqn[,52]-homosc_uneqn[,53])),3) # 0.01
round(sum(abs(homosc_uneqn[,52]-homosc_uneqn[,53])),3) # 0.076
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(heterosc_eqn[,52]-heterosc_eqn[,53])),3) # 0.036
round(sum(abs(heterosc_eqn[,52]-heterosc_eqn[,53])),3) # .303
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(heterosc_uneqn[,52]-heterosc_uneqn[,53])),3) # .03
round(sum(abs(heterosc_uneqn[,52]-heterosc_uneqn[,53])),3) # .586

### Bias
# equal variances and sample sizes (15 scenarios)
round(max(abs(homosc_eqn[,54]-homosc_eqn[,55])),3) # .012
round(sum(abs(homosc_eqn[,54]-homosc_eqn[,55])),3) # .036
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(homosc_uneqn[,54]-homosc_uneqn[,55])),3) # .01
round(sum(abs(homosc_uneqn[,54]-homosc_uneqn[,55])),3) # .076
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(heterosc_eqn[,54]-heterosc_eqn[,55])),3) # .036
round(sum(abs(heterosc_eqn[,54]-heterosc_eqn[,55])),3) # .303
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(heterosc_uneqn[,54]-heterosc_uneqn[,55])),3) # .03
round(sum(abs(heterosc_uneqn[,54]-heterosc_uneqn[,55])),3) # .586

### relative bias
# equal variances and sample sizes (15 scenarios)
round(max(abs(nonnull_homosc_eqn[,56]-nonnull_homosc_eqn[,57])),3) # .004
round(sum(abs(nonnull_homosc_eqn[,56]-nonnull_homosc_eqn[,57])),3) # .015
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(nonnull_homosc_uneqn[,56]-nonnull_homosc_uneqn[,57])),3) # .003
round(sum(abs(nonnull_homosc_uneqn[,56]-nonnull_homosc_uneqn[,57])),3) # .032
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(nonnull_heterosc_eqn[,56]-nonnull_heterosc_eqn[,57])),3) # .012
round(sum(abs(nonnull_heterosc_eqn[,56]-nonnull_heterosc_eqn[,57])),3) # .153
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(nonnull_heterosc_uneqn[,56]-nonnull_heterosc_uneqn[,57])),3) # .009
round(sum(abs(nonnull_heterosc_uneqn[,56]-nonnull_heterosc_uneqn[,57])),3) # .278

### variance
# equal variances and sample sizes (15 scenarios)
round(max(homosc_eqn[,58]/homosc_eqn[,59]),3) # 1.006
round(min(homosc_eqn[,58]/homosc_eqn[,59]),3) # .91 
# equal variances, unequal sample sizes (30 scenarios) 
round(max(homosc_uneqn[,58]/homosc_uneqn[,59]),3) # 1.007
round(min(homosc_uneqn[,58]/homosc_uneqn[,59]),3) # .902
# unequal variances, equal sample sizes (90 scenarios) 
round(max(heterosc_eqn[,58]/heterosc_eqn[,59]),3) # 1.007
round(min(heterosc_eqn[,58]/heterosc_eqn[,59]),3) # .874
# unequal variances and sample sizes (90 scenarios) 
round(max(heterosc_uneqn[,58]/heterosc_uneqn[,59]),3) # 1.011
round(min(heterosc_uneqn[,58]/heterosc_uneqn[,59]),3) # .86


##################### unbiased cohen's d's ##################### 

### Mean
# equal variances and sample sizes (15 scenarios)
round(max(abs(homosc_eqn[,60]-homosc_eqn[,61])),3) # .015
round(sum(abs(homosc_eqn[,60]-homosc_eqn[,61])),3) # .046 
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(homosc_uneqn[,60]-homosc_uneqn[,61])),3) # 0.01
round(sum(abs(homosc_uneqn[,60]-homosc_uneqn[,61])),3) # 0.079
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(heterosc_eqn[,60]-heterosc_eqn[,61])),3) # 0.034
round(sum(abs(heterosc_eqn[,60]-heterosc_eqn[,61])),3) # .277
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(heterosc_uneqn[,60]-heterosc_uneqn[,61])),3) # .029
round(sum(abs(heterosc_uneqn[,60]-heterosc_uneqn[,61])),3) # .514

### Bias
# equal variances and sample sizes (15 scenarios)
round(max(abs(homosc_eqn[,62]-homosc_eqn[,63])),3) # .015
round(sum(abs(homosc_eqn[,62]-homosc_eqn[,63])),3) # .046
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(homosc_uneqn[,62]-homosc_uneqn[,63])),3) # .01
round(sum(abs(homosc_uneqn[,62]-homosc_uneqn[,63])),3) # .079
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(heterosc_eqn[,62]-heterosc_eqn[,63])),3) # .034
round(sum(abs(heterosc_eqn[,62]-heterosc_eqn[,63])),3) # .277
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(heterosc_uneqn[,62]-heterosc_uneqn[,63])),3) # .029
round(sum(abs(heterosc_uneqn[,62]-heterosc_uneqn[,63])),3) # .514

### relative bias
# equal variances and sample sizes (15 scenarios)
round(max(abs(nonnull_homosc_eqn[,64]-nonnull_homosc_eqn[,65])),3) # .004
round(sum(abs(nonnull_homosc_eqn[,64]-nonnull_homosc_eqn[,65])),3) # .018
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(nonnull_homosc_uneqn[,64]-nonnull_homosc_uneqn[,65])),3) # .003
round(sum(abs(nonnull_homosc_uneqn[,64]-nonnull_homosc_uneqn[,65])),3) # .033
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(nonnull_heterosc_eqn[,64]-nonnull_heterosc_eqn[,65])),3) # .011
round(sum(abs(nonnull_heterosc_eqn[,64]-nonnull_heterosc_eqn[,65])),3) # .141
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(nonnull_heterosc_uneqn[,64]-nonnull_heterosc_uneqn[,65])),3) # .008
round(sum(abs(nonnull_heterosc_uneqn[,64]-nonnull_heterosc_uneqn[,65])),3) # .246

### variance
# equal variances and sample sizes (15 scenarios)
round(max(homosc_eqn[,66]/homosc_eqn[,67]),3) # 1.006
round(min(homosc_eqn[,66]/homosc_eqn[,67]),3) # .908 
# equal variances, unequal sample sizes (30 scenarios) 
round(max(homosc_uneqn[,66]/homosc_uneqn[,67]),3) # 1.007
round(min(homosc_uneqn[,66]/homosc_uneqn[,67]),3) # .925
# unequal variances, equal sample sizes (90 scenarios) 
round(max(heterosc_eqn[,66]/heterosc_eqn[,67]),3) # 1.008
round(min(heterosc_eqn[,66]/heterosc_eqn[,67]),3) # .89
# unequal variances and sample sizes (90 scenarios) 
round(max(heterosc_uneqn[,66]/heterosc_uneqn[,67]),3) # 1.011
round(min(heterosc_uneqn[,66]/heterosc_uneqn[,67]),3) # .882

##################### Shieh ##################### 

### Mean
# equal variances and sample sizes (15 scenarios)
round(max(abs(homosc_eqn[,68]-homosc_eqn[,69])),3) # .006
round(sum(abs(homosc_eqn[,68]-homosc_eqn[,69])),3) # .018 
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(homosc_uneqn[,68]-homosc_uneqn[,69])),3) # .008
round(sum(abs(homosc_uneqn[,68]-homosc_uneqn[,69])),3) # .065
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(heterosc_eqn[,68]-heterosc_eqn[,69])),3) # .018
round(sum(abs(heterosc_eqn[,68]-heterosc_eqn[,69])),3) # .152
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(heterosc_uneqn[,68]-heterosc_uneqn[,69])),3) # .009
round(sum(abs(heterosc_uneqn[,68]-heterosc_uneqn[,69])),3) # .254

### Bias
# equal variances and sample sizes (15 scenarios)
round(max(abs(homosc_eqn[,70]-homosc_eqn[,71])),3) # .006
round(sum(abs(homosc_eqn[,70]-homosc_eqn[,71])),3) # .018
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(homosc_uneqn[,70]-homosc_uneqn[,71])),3) # .008
round(sum(abs(homosc_uneqn[,70]-homosc_uneqn[,71])),3) # .065
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(heterosc_eqn[,70]-heterosc_eqn[,71])),3) # .018
round(sum(abs(heterosc_eqn[,70]-heterosc_eqn[,71])),3) # .152
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(heterosc_uneqn[,70]-heterosc_uneqn[,71])),3) # .009
round(sum(abs(heterosc_uneqn[,70]-heterosc_uneqn[,71])),3) # .254

### relative bias
round(max(abs(nonnull_homosc_eqn[,72]-nonnull_homosc_eqn[,73])),3) # .004
round(sum(abs(nonnull_homosc_eqn[,72]-nonnull_homosc_eqn[,73])),3) # .015
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(nonnull_homosc_uneqn[,72]-nonnull_homosc_uneqn[,73])),3) # .005
round(sum(abs(nonnull_homosc_uneqn[,72]-nonnull_homosc_uneqn[,73])),3) # .064
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(nonnull_heterosc_eqn[,72]-nonnull_heterosc_eqn[,73])),3) # .012
round(sum(abs(nonnull_heterosc_eqn[,72]-nonnull_heterosc_eqn[,73])),3) # .153
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(nonnull_heterosc_uneqn[,72]-nonnull_heterosc_uneqn[,73])),3) # .009
round(sum(abs(nonnull_heterosc_uneqn[,72]-nonnull_heterosc_uneqn[,73])),3) # .31

### Variance
# equal variances and sample sizes (15 scenarios)
round(max(homosc_eqn[,74]/homosc_eqn[,75]),3) # 1.006
round(min(homosc_eqn[,74]/homosc_eqn[,75]),3) # .91 
# equal variances, unequal sample sizes (30 scenarios) 
round(max(homosc_uneqn[,74]/homosc_uneqn[,75]),3) # 1.005
round(min(homosc_uneqn[,74]/homosc_uneqn[,75]),3) # .865
# unequal variances, equal sample sizes (90 scenarios) 
round(max(heterosc_eqn[,74]/heterosc_eqn[,75]),3) # 1.007
round(min(heterosc_eqn[,74]/heterosc_eqn[,75]),3) # .874
# unequal variances and sample sizes (90 scenarios) 
round(max(heterosc_uneqn[,74]/heterosc_uneqn[,75]),3) # 1.011
round(min(heterosc_uneqn[,74]/heterosc_uneqn[,75]),3) # .867

##################### Unbiased Shieh ##################### 

### Mean
# equal variances and sample sizes (15 scenarios)
round(max(abs(homosc_eqn[,76]-homosc_eqn[,77])),3) # .008
round(sum(abs(homosc_eqn[,76]-homosc_eqn[,77])),3) # .023 
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(homosc_uneqn[,76]-homosc_uneqn[,77])),3) # .007
round(sum(abs(homosc_uneqn[,76]-homosc_uneqn[,77])),3) # .052
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(heterosc_eqn[,76]-heterosc_eqn[,77])),3) # .017
round(sum(abs(heterosc_eqn[,76]-heterosc_eqn[,77])),3) # .139
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(heterosc_uneqn[,76]-heterosc_uneqn[,77])),3) # .008
round(sum(abs(heterosc_uneqn[,76]-heterosc_uneqn[,77])),3) # .237

### Bias
# equal variances and sample sizes (15 scenarios)
round(max(abs(homosc_eqn[,78]-homosc_eqn[,79])),3) # .008
round(sum(abs(homosc_eqn[,78]-homosc_eqn[,79])),3) # .023
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(homosc_uneqn[,78]-homosc_uneqn[,79])),3) # .007
round(sum(abs(homosc_uneqn[,78]-homosc_uneqn[,79])),3) # .052
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(heterosc_eqn[,78]-heterosc_eqn[,79])),3) # .017
round(sum(abs(heterosc_eqn[,78]-heterosc_eqn[,79])),3) # .139
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(heterosc_uneqn[,78]-heterosc_uneqn[,79])),3) # .008
round(sum(abs(heterosc_uneqn[,78]-heterosc_uneqn[,79])),3) # .237

### relative bias
round(max(abs(nonnull_homosc_eqn[,80]-nonnull_homosc_eqn[,81])),3) # .004
round(sum(abs(nonnull_homosc_eqn[,80]-nonnull_homosc_eqn[,81])),3) # .018
# equal variances, unequal sample sizes (30 scenarios) 
round(max(abs(nonnull_homosc_uneqn[,80]-nonnull_homosc_uneqn[,81])),3) # .004
round(sum(abs(nonnull_homosc_uneqn[,80]-nonnull_homosc_uneqn[,81])),3) # .051
# unequal variances, equal sample sizes (90 scenarios) 
round(max(abs(nonnull_heterosc_eqn[,80]-nonnull_heterosc_eqn[,81])),3) # .011
round(sum(abs(nonnull_heterosc_eqn[,80]-nonnull_heterosc_eqn[,81])),3) # .141
# unequal variances and sample sizes (90 scenarios) 
round(max(abs(nonnull_heterosc_uneqn[,80]-nonnull_heterosc_uneqn[,81])),3) # .009
round(sum(abs(nonnull_heterosc_uneqn[,80]-nonnull_heterosc_uneqn[,81])),3) # .288

### Variance
# equal variances and sample sizes (15 scenarios)
round(max(homosc_eqn[,82]/homosc_eqn[,83]),3) # 1.006
round(min(homosc_eqn[,82]/homosc_eqn[,83]),3) # .908
# equal variances, unequal sample sizes (30 scenarios) 
round(max(homosc_uneqn[,82]/homosc_uneqn[,83]),3) # 1.007
round(min(homosc_uneqn[,82]/homosc_uneqn[,83]),3) # .9
# unequal variances, equal sample sizes (90 scenarios) 
round(max(heterosc_eqn[,82]/heterosc_eqn[,83]),3) # 1.008
round(min(heterosc_eqn[,82]/heterosc_eqn[,83]),3) # .89
# unequal variances and sample sizes (90 scenarios) 
round(max(heterosc_uneqn[,82]/heterosc_uneqn[,83]),3) # 1.011
round(min(heterosc_uneqn[,82]/heterosc_uneqn[,83]),3) # .881

# 3) On peut alors étudier les paramètres dont dépendent le biais et la variance de chaque estimateur, et l'impact de leur évolution


