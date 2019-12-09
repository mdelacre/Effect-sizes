

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

get_write(1,n1=10,n2=10,sd1=2,sd2=2,m1=1,m2=0)  

##### Function to generate dataset,realize statistical test and compute (and extract) p-value

get_simu     <- function(nSims=100000,n1=50,n2=50,  
                           sd1=2, sd2=2,  
                           m1=1,m2=0){

# set up empty container for all estimated parameters
     ES <-matrix(0,nSims,3) # Three colums to store the Cohen's d, Shieh's and modified shieh's d measures
     colnames(ES) <- c("Cohen's d", "Shieh's", "Modified Shieh's")

    for (i in 1:nSims){

      # Dataset as a function of the number of groups (k)
          y1 <- rnorm(n1,m1,sd1)
          y2 <- rnorm(n2,m2,sd2)

       # Effect sizes measures

       mean1 <- mean(y1)
       mean2 <- mean(y2)
       sdev1 <- sd(y1)
       sdev2 <- sd(y2)
       pooled_sd<-sqrt(((n1-1)*sdev1^2+(n2-1)*sdev2^2)/(n1+n2-2))
       ES[i,1] <- (mean1-mean2)/pooled_sd
       
       N <- n1+n2
       q1 <- n1/N
       q2 <- n2/N
       ES[i,2] <- (mean1-mean2)/sqrt(sdev1^2/q1+sdev2^2/q2)
       
       nratio <- n1/n2 
       ES[i,3] <- (mean1-mean2)/sqrt(sdev1^2/q1+sdev2^2/q2)*((nratio+1)/sqrt(nratio))
       }

   # p-values extraction in a specified file 
    setwd(dir="C:/Users/Marie/Documents/transfo_shieh test") # destination file  
    get_write(ES,n1,n2,sd1,sd2,m1,m2)
    
}

#####################################################################################
#############                       APPLICATIONS                        #############
#####################################################################################

#############                    Type 1 error rate                      #############

### normal distributions

get_simu(n1=20,n2=10,m1=0,m2=1,sd1=2,sd2=1)
get_simu(n1=20,n2=10,m1=0,m2=1,sd1=2,sd2=2)
get_simu(n1=20,n2=10,m1=0,m2=1,sd1=2,sd2=4)
get_simu(n1=20,n2=10,m1=0,m2=1,sd1=2,sd2=8)
get_simu(n1=20,n2=20,m1=0,m2=1,sd1=2,sd2=1)
get_simu(n1=20,n2=20,m1=0,m2=1,sd1=2,sd2=2)
get_simu(n1=20,n2=20,m1=0,m2=1,sd1=2,sd2=4)
get_simu(n1=20,n2=20,m1=0,m2=1,sd1=2,sd2=8)
get_simu(n1=20,n2=30,m1=0,m2=1,sd1=2,sd2=1)
get_simu(n1=20,n2=30,m1=0,m2=1,sd1=2,sd2=2)
get_simu(n1=20,n2=30,m1=0,m2=1,sd1=2,sd2=4)
get_simu(n1=20,n2=30,m1=0,m2=1,sd1=2,sd2=8)
get_simu(n1=20,n2=40,m1=0,m2=1,sd1=2,sd2=1)
get_simu(n1=20,n2=40,m1=0,m2=1,sd1=2,sd2=2)
get_simu(n1=20,n2=40,m1=0,m2=1,sd1=2,sd2=4)
get_simu(n1=20,n2=40,m1=0,m2=1,sd1=2,sd2=8)
get_simu(n1=30,n2=15,m1=0,m2=1,sd1=2,sd2=1)
get_simu(n1=30,n2=15,m1=0,m2=1,sd1=2,sd2=2)
get_simu(n1=30,n2=15,m1=0,m2=1,sd1=2,sd2=4)
get_simu(n1=30,n2=15,m1=0,m2=1,sd1=2,sd2=8)
get_simu(n1=30,n2=30,m1=0,m2=1,sd1=2,sd2=1)
get_simu(n1=30,n2=30,m1=0,m2=1,sd1=2,sd2=2)
get_simu(n1=30,n2=30,m1=0,m2=1,sd1=2,sd2=4)
get_simu(n1=30,n2=30,m1=0,m2=1,sd1=2,sd2=8)
get_simu(n1=30,n2=45,m1=0,m2=1,sd1=2,sd2=1)
get_simu(n1=30,n2=45,m1=0,m2=1,sd1=2,sd2=2)
get_simu(n1=30,n2=45,m1=0,m2=1,sd1=2,sd2=4)
get_simu(n1=30,n2=45,m1=0,m2=1,sd1=2,sd2=8)
get_simu(n1=30,n2=60,m1=0,m2=1,sd1=2,sd2=1)
get_simu(n1=30,n2=60,m1=0,m2=1,sd1=2,sd2=2)
get_simu(n1=30,n2=60,m1=0,m2=1,sd1=2,sd2=4)
get_simu(n1=30,n2=60,m1=0,m2=1,sd1=2,sd2=8)
get_simu(n1=40,n2=20,m1=0,m2=1,sd1=2,sd2=1)
get_simu(n1=40,n2=20,m1=0,m2=1,sd1=2,sd2=2)
get_simu(n1=40,n2=20,m1=0,m2=1,sd1=2,sd2=4)
get_simu(n1=40,n2=20,m1=0,m2=1,sd1=2,sd2=8)
get_simu(n1=40,n2=40,m1=0,m2=1,sd1=2,sd2=1)
get_simu(n1=40,n2=40,m1=0,m2=1,sd1=2,sd2=2)
get_simu(n1=40,n2=40,m1=0,m2=1,sd1=2,sd2=4)
get_simu(n1=40,n2=40,m1=0,m2=1,sd1=2,sd2=8)
get_simu(n1=40,n2=60,m1=0,m2=1,sd1=2,sd2=1)
get_simu(n1=40,n2=60,m1=0,m2=1,sd1=2,sd2=2)
get_simu(n1=40,n2=60,m1=0,m2=1,sd1=2,sd2=4)
get_simu(n1=40,n2=60,m1=0,m2=1,sd1=2,sd2=8)
get_simu(n1=40,n2=80,m1=0,m2=1,sd1=2,sd2=1)
get_simu(n1=40,n2=80,m1=0,m2=1,sd1=2,sd2=2)
get_simu(n1=40,n2=80,m1=0,m2=1,sd1=2,sd2=4)
get_simu(n1=40,n2=80,m1=0,m2=1,sd1=2,sd2=8)
get_simu(n1=50,n2=25,m1=0,m2=1,sd1=2,sd2=1)
get_simu(n1=50,n2=25,m1=0,m2=1,sd1=2,sd2=2)
get_simu(n1=50,n2=25,m1=0,m2=1,sd1=2,sd2=4)
get_simu(n1=50,n2=25,m1=0,m2=1,sd1=2,sd2=8)
get_simu(n1=50,n2=50,m1=0,m2=1,sd1=2,sd2=1)
get_simu(n1=50,n2=50,m1=0,m2=1,sd1=2,sd2=2)
get_simu(n1=50,n2=50,m1=0,m2=1,sd1=2,sd2=4)
get_simu(n1=50,n2=50,m1=0,m2=1,sd1=2,sd2=8)
get_simu(n1=50,n2=75,m1=0,m2=1,sd1=2,sd2=1)
get_simu(n1=50,n2=75,m1=0,m2=1,sd1=2,sd2=2)
get_simu(n1=50,n2=75,m1=0,m2=1,sd1=2,sd2=4)
get_simu(n1=50,n2=75,m1=0,m2=1,sd1=2,sd2=8)
get_simu(n1=50,n2=100,m1=0,m2=1,sd1=2,sd2=1)
get_simu(n1=50,n2=100,m1=0,m2=1,sd1=2,sd2=2)
get_simu(n1=50,n2=100,m1=0,m2=1,sd1=2,sd2=4)
get_simu(n1=50,n2=100,m1=0,m2=1,sd1=2,sd2=8)
get_simu(n1=100,n2=50,m1=0,m2=1,sd1=2,sd2=1)
get_simu(n1=100,n2=50,m1=0,m2=1,sd1=2,sd2=2)
get_simu(n1=100,n2=50,m1=0,m2=1,sd1=2,sd2=4)
get_simu(n1=100,n2=50,m1=0,m2=1,sd1=2,sd2=8)
get_simu(n1=100,n2=100,m1=0,m2=1,sd1=2,sd2=1)
get_simu(n1=100,n2=100,m1=0,m2=1,sd1=2,sd2=2)
get_simu(n1=100,n2=100,m1=0,m2=1,sd1=2,sd2=4)
get_simu(n1=100,n2=100,m1=0,m2=1,sd1=2,sd2=8)
get_simu(n1=100,n2=150,m1=0,m2=1,sd1=2,sd2=1)
get_simu(n1=100,n2=150,m1=0,m2=1,sd1=2,sd2=2)
get_simu(n1=100,n2=150,m1=0,m2=1,sd1=2,sd2=4)
get_simu(n1=100,n2=150,m1=0,m2=1,sd1=2,sd2=8)
get_simu(n1=100,n2=200,m1=0,m2=1,sd1=2,sd2=1)
get_simu(n1=100,n2=200,m1=0,m2=1,sd1=2,sd2=2)
get_simu(n1=100,n2=200,m1=0,m2=1,sd1=2,sd2=4)
get_simu(n1=100,n2=200,m1=0,m2=1,sd1=2,sd2=8)



