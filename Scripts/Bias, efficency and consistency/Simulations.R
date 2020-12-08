for (package in c("PearsonDS","gsl")) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}


##### Function to load packages

# Note: when kurtosis = 3 and skewness = 0, one has a normal distribution

get_write <- function(object,
                      kurt=0,
                      skew=0,
                      n1=40,n2=40,
                      sd1=2,sd2=2,
                      m1=1,m2=0){
  
  # compute mean and standard deviation
  Skewness <-  paste0(" G1=",skew)
  Kurtosis <-  paste0(" G2=",kurt)
  nobs <- paste0(" n=[",n1,",",n2,"]")
  means <-  paste0(" means=[",m1,",",m2,"]")
  stdevs <- paste0(" sds=[",sd1,",",sd2,"]")
  fname <-  paste(Skewness, Kurtosis,nobs, means, stdevs, sep=",")
  fname <-  paste0(fname, ".rds")
  saveRDS(object, file = fname)
  
}

##### Function to generate dataset,realize statistical test and compute (and extract) p-value

get_simu     <- function(nSims=1000000,n1=50,n2=50,  
                         kurt=3,skew=0,
                         sd1=2, sd2=2,  
                         m1=1,m2=0){

  # set up empty container for all estimated parameters
  ES <-matrix(0,nSims,18+4*2) # Six colums to store the Cohen's d, Shieh's d and modified shieh's d measures
  # 8 columns to store de mean, sd, skewness and kurtosis of each measures
  estimator <- c("Cohen's d","Hedge's g","Glass's sd1","Glass's sd2","Unbiased Glass1","Unbiased Glass2","Shieh's d", "Unbiased Shieh's d","cohen's d'","unbiased cohen's d'")
  descr=expand.grid(paste(1:2),c("mean","sd", "skewness","kurtosis"))
  colnames(ES) <- c(paste0(descr[,2],descr[,1]),estimator)
  
  for (i in 1:nSims){
    
    # Dataset as a function of the number of groups (k)
    y1 <- rpearson(n1,moments=c(m1,sd1^2,skewness=skew*(n1-2)/sqrt(n1*(n1-1)),kurtosis=(kurt*(n1-2)*(n1-3)-6*(n1-1))/(n1^2-1)+3)) # Note: the function takes the variance                                                                              # as an argument (not the sd)
    y2 <- rpearson(n2,moments=c(m2,sd2^2,skewness=skew*(n2-2)/sqrt(n2*(n2-1)),kurtosis=(kurt*(n2-2)*(n2-3)-6*(n2-1))/(n2^2-1)+3))

    # For the explanation about the transformation of skew1,skew2,kurt1 and kurt2,
    ### Descriptives
    mean1 <- mean(y1)
    mean2 <- mean(y2)
    sdev1 <- sd(y1)
    sdev2 <- sd(y2)

    M2_y1 <- sum((y1-mean1)^2)/n1 # variance non corrigée
    M3_y1 <- sum((y1-mean1)^3)/n1
    M4_y1 <- sum((y1-mean1)^4)/n1
    kurtosis1 <- ((n1-1)/((n1-2)*(n1-3)))*((n1+1)*(M4_y1/(M2_y1^2)-3)+6)
    skewness1 <- sqrt(n1*(n1-1))/(n1-2)*(M3_y1/M2_y1^(3/2))

    M2_y2 <- sum((y2-mean2)^2)/n2 # variance non corrigée
    M3_y2 <- sum((y2-mean2)^3)/n2
    M4_y2 <- sum((y2-mean2)^4)/n2
    kurtosis2 <- ((n2-1)/((n2-2)*(n2-3)))*((n2+1)*(M4_y2/(M2_y2^2)-3)+6)
    skewness2 <- sqrt(n2*(n2-1))/(n2-2)*(M3_y2/M2_y2^(3/2))

    ES[i,1:8] <- c(mean1,mean2,sdev1,sdev2,kurtosis1,kurtosis2,skewness1,skewness2)
    
    ### Effect sizes measures
    
    # Cohen's d 
    pooled_sd <-sqrt(((n1-1)*sdev1^2+(n2-1)*sdev2^2)/(n1+n2-2))
    cohen_d <- (mean1-mean2)/pooled_sd
    
    # Hedges's g
    N <- n1+n2
    hedge_g <- cohen_d*gamma((N-2)/2)/(sqrt((N-2)/2)*gamma((N-3)/2)) #(1-3/(4*(n1+n2)-9))   
    
    # Glass's delta
    glass_sd1 <- (mean1-mean2)/sdev1
    glass_sd2 <- (mean1-mean2)/sdev2

    # Unbiased Glass's delta
    unbiased_glass_sd1 <- glass_sd1*gamma((n1-1)/2)/(sqrt((n1-1)/2)*gamma((n1-2)/2))
    unbiased_glass_sd2 <- glass_sd2*gamma((n2-1)/2)/(sqrt((n2-1)/2)*gamma((n2-2)/2))
    
    # Shieh's d
    q1 <- n1/N
    q2 <- n2/N
    shieh_d <- (mean1-mean2)/sqrt(sdev1^2/q1+sdev2^2/q2)       

    # Unbiased Shieh's d
    df <- (sdev1^2/n1+sdev2^2/n2)^2/((sdev1^2/n1)^2/(n1-1)+(sdev2^2/n2)^2/(n2-1))
    correction <- gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2))
    unbiased_shieh <- shieh_d*correction

    # Cohen's d's

    cohen_d_prime <- (mean1-mean2)/sqrt((sdev1^2+sdev2^2)/2)
    df <- (sdev1^2/n1+sdev2^2/n2)^2/((sdev1^2/n1)^2/(n1-1)+(sdev2^2/n2)^2/(n2-1))
    unbiased_cohen_d_prime <- shieh_d_corr*gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2))

    ES[i,9:18] <- c(cohen_d,hedge_g,glass_sd1,glass_sd2,unbiased_glass_sd1,unbiased_glass_sd2,shieh_d,unbiased_shieh,cohen_d_prime,unbiased_cohen_d_prime)
    
  }
  
  # Extraction of the ES matrix 
  chem <- "C:/Users/Marie/Documents/ES MEASURES/"
  sschem<- paste0("G1=",skew,",G2=",kurt) 
  setwd(dir=paste0(chem,sschem)) # destination file  
  get_write(ES,kurt,skew,n1,n2,sd1,sd2,m1,m2)

}

#####################################################################################
#############                       APPLICATIONS                        #############
#####################################################################################

#### Normal case

n1 <- c(20,50,100)
n2 <- c(20,50,100)
m1 <- seq(0,4,1)
m2 <- 0
Kurt <- 0 # when asymetry, kurtosis > 0, always!
Skew <- 0
sd1 <- 1
sd2 <- c(.1,.25,.5,1,2,4,10)

Simu=expand.grid(Skew,Kurt,n1,n2,sd1,sd2,m1,m2)
length(Simu[,1])
colnames(Simu)<-c("skewness","kurtosis","n1","n2","sd1","sd2","m1","m2")

#length(Simu[,1])

# create subfolders
#fold<-expand.grid(Skew,Kurt)
#for (j in seq_len(length(fold[,1]))){
#  setwd("C:/Users/Marie/Documents/ES MEASURES")
#  dir.create(paste0("G1=",fold[j,1],",G2=",fold[j,2]))
#}

# performing simulations  
for (i in seq_len(length(Simu[,1]))){
  get_simu(nSims=100000,n1=Simu[i,3],n2=Simu[i,4],sd1=Simu[i,5],sd2=Simu[i,6],m1=Simu[i,7],m2=Simu[i,8],skew=Simu[i,1],kurt=Simu[i,2])  
}

#### Abnormal cases

n1 <- c(20,50,100)
n2 <- c(20,50,100)
m1 <- seq(0,4,1)
m2 <- 0
Kurt <- 95.75
Skew <- c(-2.08,0,6.32) # skewness = 0 is the symmetric case
sd1 <- 1
sd2 <- c(.1,.25,.5,1,2,4,10)

Simu=expand.grid(Skew,Kurt,n1,n2,sd1,sd2,m1,m2)
length(Simu[,1])
colnames(Simu)<-c("skewness","kurtosis","n1","n2","sd1","sd2","m1","m2")

#length(Simu[,1])

length(Simu$skewness[Simu$kurtosis==95.75])

# create subfolders
#fold<-expand.grid(Skew,Kurt)
#for (j in seq_len(length(fold[,1]))){  
#  setwd("C:/Users/Marie/Documents/ES MEASURES")
#  dir.create(paste0("G1=",fold[j,1],",G2=",fold[j,2]))
#}

#View(Simu)

setwd("C:/Users/Marie/Documents/ES MEASURES/G1=0,G2=95.75")

# performing simulations  
  for (i in seq_len(length(Simu[,1]))){ 
    get_simu(nSims=100000,n1=Simu[i,3],n2=Simu[i,4],sd1=Simu[i,5],sd2=Simu[i,6],m1=Simu[i,7],m2=Simu[i,8],skew=Simu[i,1],kurt=Simu[i,2])  
  }

