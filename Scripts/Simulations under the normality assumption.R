for (package in c("moments")) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}

##### Function to load packages

get_write <- function(object,n1,n2,
                      sd1,sd2,
                      m1,m2){

# compute mean and standard deviation
  distr <-  "Distr=[normal,normal]"
  nobs <- paste0(" n=[",n1,",",n2,"]")
  means <-  paste0(" means=[",m1,",",m2,"]")
  stdevs <- paste0(" sds=[",sd1,",",sd2,"]")
  fname <-  paste(distr, nobs, means, stdevs, sep=",")
  fname <-  paste0(fname, ".rds")
  saveRDS(object, file = fname)
  
}

##### Function to generate dataset,realize statistical test and compute (and extract) p-value

get_simu     <- function(nSims=1000000,n1=50,n2=50,  
                           sd1=2, sd2=2,  
                           m1=1,m2=0){

# set up empty container for all estimated parameters
     ES <-matrix(0,nSims,6+4*2) # Six colums to store the Cohen's d, Shieh's d and modified shieh's d measures
                                # 8 columns to store de mean, sd, skewness and kurtosis of each measures
     estimator <- c("Cohen's d","Hedge's d","Glass's sd1","Glass's sd2","Shieh's d", "Shieh's d corr")
     descr=expand.grid(paste(1:2),c("mean","sd", "skewness","kurtosis"))
     colnames(ES) <- c(paste0(descr[,2],descr[,1]),estimator)

    for (i in 1:nSims){

      # Dataset as a function of the number of groups (k)
          y1 <- rnorm(n1,m1,sd1)
          y2 <- rnorm(n2,m2,sd2)
          
      ### Descriptives

       mean1 <- mean(y1)
       mean2 <- mean(y2)
       sdev1 <- sd(y1)
       sdev2 <- sd(y2)
       kurt1 <- kurtosis(y1)
       kurt2 <- kurtosis(y2)
       skew1 <- skewness(y1)
       skew2 <- skewness(y2)
       
       ES[i,1:8] <- c(mean1,mean2,sdev1,sdev2,kurt1,kurt2,skew1,skew2)
       
      ### Effect sizes measures
       
       # Cohen's d
       pooled_sd <-sqrt(((n1-1)*sdev1^2+(n2-1)*sdev2^2)/(n1+n2-2))
       cohen_d <- (mean1-mean2)/pooled_sd

       # Hedges's g
       hedge_g <- cohen_d*(1-3/(4*(n1+n2-9)))   
       
       # Glass's delta
       glass_sd1 <- (mean1-mean2)/sdev1
       glass_sd2 <- (mean1-mean2)/sdev2
       
       # Shieh's d
       N <- n1+n2
       q1 <- n1/N
       q2 <- n2/N
       shieh_d <- (mean1-mean2)/sqrt(sdev1^2/q1+sdev2^2/q2)       

      # Shieh's d after correction
       nratio <- n1/n2
       shieh_d_corr <- (mean1-mean2)/sqrt(sdev1^2/q1+sdev2^2/q2)*((nratio+1)/sqrt(nratio))
       
       ES[i,9:14] <- c(cohen_d,hedge_g,glass_sd1,glass_sd2,shieh_d,shieh_d_corr)
       
       }

   # Extraction of the ES matrix 
    setwd(dir="C:/Users/Marie/Documents/ES MEASURES/Normal") # destination file  
    get_write(ES,n1,n2,sd1,sd2,m1,m2)
    
}

#####################################################################################
#############                       APPLICATIONS                        #############
#####################################################################################

n1=c(20,50,100)
n2=c(20,50,100)
m1=seq(0,4,1)
m2=0
sd1=1
sd2=c(.1,.25,.5,1,2,4,10)

Simu=expand.grid(n1,n2,sd1,sd2,m1,m2)
colnames(Simu)<-c("n1","n2","sd1","sd2","m1","m2")
# View(Simu)

for (i in seq_len(length(A[,1]))){
  get_simu(nSims=100000,Simu[i,1],Simu[i,2],Simu[i,3],Simu[i,4],Simu[i,5],Simu[i,6])  
}


