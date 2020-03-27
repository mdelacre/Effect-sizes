for (package in c("PearsonDS","gsl")) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}


##### Function to load packages

# Note: when kurtosis = 3 and skewness = 0, one has a normal distribution

get_write <- function(object,
                      conf.level=.95,
                      alternative="two.sided",
                      kurt=0,
                      skew=0,
                      n1=40,n2=40,
                      sd1=2,sd2=2,
                      m1=1,m2=0){
   
  # compute mean and standard deviation
  Conf.level <-  paste0(" conf.level=",conf.level) 
  Alt <- alternative
  Skewness <-  paste0(" G1=",skew)
  Kurtosis <-  paste0(" G2=",kurt)
  nobs <- paste0(" n=[",n1,",",n2,"]")
  means <-  paste0(" means=[",m1,",",m2,"]")
  stdevs <- paste0(" sds=[",sd1,",",sd2,"]")
  fname <-  paste(Conf.level,Alt,Skewness, Kurtosis,nobs, means, stdevs, sep=",")
  fname <-  paste0(fname, ".rds")
  saveRDS(object, file = fname)
  
}

##### Function to generate dataset,realize statistical test and compute (and extract) p-value

get_simu     <- function(nSims=1000000,conf.level=.95,alternative="two.sided",n1=50,n2=50,  
                         kurt=3,skew=0,
                         sd1=2, sd2=2,  
                         m1=1,m2=0){

  # set up empty container for all estimated parameters
  CL <-matrix(0,nSims,12) 
  # 5 columnms dedicated to Student's t-test: t stat, cohen's d, hedge's g and confidence limits of Cohen's d
  # 7 columnms dedicated to Welch's t-test: t stat, Shieh's d and confidence limits of Shieh's d
  #                                         Transformed Shieh's d and confidence limits of Transformed Shieh's d

# Student's t-test
  ST <- c("student_p","cohen_ds","hedges_gs","low_lim_cohen","up_lim_cohen")  
  # Welch's t-test
  WT <- c("welch_p","shieh_ds","low_lim_shieh","up_lim_shieh","transfoshieh_ds","low_lim__transfoshieh","up_lim__transfoshieh")  
  colnames(CL) <- c(ST,WT)
  
  for (i in 1:nSims){
    
    # Dataset as a function of the number of groups (k)
    y1 <- rpearson(n1,moments=c(m1,sd1^2,skewness=skew*(n1-2)/sqrt(n1*(n1-1)),kurtosis=(kurt*(n1-2)*(n1-3)-6*(n1-1))/(n1^2-1)+3))                                                                             # as an argument (not the sd)
    y2 <- rpearson(n2,moments=c(m2,sd2^2,skewness=skew*(n2-2)/sqrt(n2*(n2-1)),kurtosis=(kurt*(n2-2)*(n2-3)-6*(n2-1))/(n2^2-1)+3))

    ### Descriptives
    mean1 <- mean(y1)
    mean2 <- mean(y2)
    sdev1 <- sd(y1)
    sdev2 <- sd(y2)
    
    # Student's t-test
    
    if(alternative=="two.sided"){
      student_p <-t.test(y1,y2,var.equal=TRUE,alternative="two.sided")$p.value       
    } if(alternative=="less"){
      student_p <-t.test(y1,y2,var.equal=TRUE,alternative="less")$p.value 
    } if(alternative=="greater"){
      student_p <-t.test(y1,y2,var.equal=TRUE,alternative="greater")$p.value 
    }

        pooled_sd <-sqrt(((n1-1)*sdev1^2+(n2-1)*sdev2^2)/(n1+n2-2))
    cohen_ds <- (mean1-mean2)/pooled_sd
    hedges_gs<- cohen_ds*(1-3/(4*(n1+n2-9)))
    
    cohen_lim <- Cohen.CI(Group.1=y1, Group.2=y2,conf.level=conf.level,alternative=alternative)
    
    CL[i,1:5] <- c(student_p,cohen_ds,hedges_gs,cohen_lim)
    
    ### Welch's t-test
    if(alternative=="two.sided"){
      welch_p  <-t.test(y1,y2,var.equal=FALSE,alternative="two.sided")$p.value
    } if(alternative=="less"){
      welch_p  <-t.test(y1,y2,var.equal=FALSE,alternative="less")$p.value
    } if(alternative=="greater"){
      welch_p  <-t.test(y1,y2,var.equal=FALSE,alternative="greater")$p.value
    }
    
      N <- n1+n2
      q1 <- n1/N
      q2 <- n2/N
    shieh_ds <- (mean1-mean2)/sqrt(sdev1^2/q1+sdev2^2/q2)
    
    shieh_lim <- Shieh.CI(Group.1=y1, Group.2=y2,conf.level=conf.level,alternative=alternative)
    
    CL[i,6:9] <- c(welch_p,shieh_ds,shieh_lim)

      nratio <- n1/n2
      sigma_bal <- sqrt((sdev1^2+sdev2^2)/2)
      sigma_unbal <- sqrt((1-q1)*sdev1^2+(1-q2)*sdev2^2) # we give more weight to the variance of the smallest group
    corr_shieh_ds <- (mean1-mean2)/sqrt(sdev1^2/q1+sdev2^2/q2)*(((nratio+1)*sigma_unbal)/(2*sigma_bal*sqrt(nratio))) # what value of Shieh's delta would be obtain if n1=n2?
    corr_shieh_lim <- shieh_lim*(((nratio+1)*sigma_unbal)/(2*sigma_bal*sqrt(nratio)))

    
    CL[i,10:12] <- c(corr_shieh_ds,corr_shieh_lim)

    
  }

  # Extraction of the ES matrix 
  chem <- "C:/Users/Marie/Documents/CI MEASURES/"
  sschem<- paste0(alternative,"/") 
  sschem2<- paste0("G1=",skew,",G2=",kurt) 
  setwd(dir=paste0(chem,sschem,sschem2)) # destination file  
  get_write(CL,conf.level,alternative,kurt,skew,n1,n2,sd1,sd2,m1,m2)

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
  get_simu(nSims=10,conf.level=.95,alternative="two.sided",n1=Simu[i,3],n2=Simu[i,4],sd1=Simu[i,5],sd2=Simu[i,6],m1=Simu[i,7],m2=Simu[i,8],skew=Simu[i,1],kurt=Simu[i,2])  
}

for (i in seq_len(length(Simu[,1]))){
  get_simu(nSims=100000,conf.level=.95,alternative="greater",n1=Simu[i,3],n2=Simu[i,4],sd1=Simu[i,5],sd2=Simu[i,6],m1=Simu[i,7],m2=Simu[i,8],skew=Simu[i,1],kurt=Simu[i,2])  
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
  for (i in c(379,381)){ #il manque 382,383,384
    get_simu(nSims=100000,conf.level=.95,n1=Simu[i,3],n2=Simu[i,4],sd1=Simu[i,5],sd2=Simu[i,6],m1=Simu[i,7],m2=Simu[i,8],skew=Simu[i,1],kurt=Simu[i,2])  
  }


# performing simulations  
for (i in seq_len(length(Simu[,1]))){ #
  get_simu(nSims=100000,conf.level=.95,alternative="greater",n1=Simu[i,3],n2=Simu[i,4],sd1=Simu[i,5],sd2=Simu[i,6],m1=Simu[i,7],m2=Simu[i,8],skew=Simu[i,1],kurt=Simu[i,2])  
}


