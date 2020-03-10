for (package in "rootSolve") {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}


#########################################
##  Simulation - Confidence intervals  ##
#########################################




#  Obtain confidence limits for mu1-mu2 (see Method 2 in CI Reminder)
#--------------------------------------------------------------

meandiff.CI <- function(Group.1, Group.2,conf.level)
{ 
  n1 <- length(Group.1)
  n2 <- length(Group.2)
  mean1 <- mean(Group.1)
  mean2 <- mean(Group.2)
  sd1 <- sd(Group.1)
  sd2 <- sd(Group.2)
  
  pooled_sd <- sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
  SE <- pooled_sd*sqrt(1/n1+1/n2)
  

  # Lower limit = limit of mu1-mu2 such as 1-pt(q=t_obs, df=n1+n2-2) = (1-conf.level)/2 = alpha/2
  # with t_obs = ((mean1-mean2)-(theo_mudiff))/SE

  f=function(theo_mudiff,rep) 1-pt(q=((mean1-mean2)-(theo_mudiff))/SE, df=n1+n2-2)-rep
  out=uniroot(f,c(-100,100),rep=(1-conf.level)/2)
  theo_mudiff.1 <- out$root
  
  # Upper limit = limit of mu1-mu2 such as pt(q=t_obs, df=n1+n2-2) = (1-conf.level)/2 = alpha/2
  # with t_obs = ((mean1-mean2)-(theo_mudiff))/SE
  
  f=function(theo_mudiff,rep) pt(q=((mean1-mean2)-(theo_mudiff))/SE, df=n1+n2-2)-rep
  out=uniroot(f,c(-100,100),rep=(1-conf.level)/2)
  theo_mudiff.2 <- out$root

  result <- c(theo_mudiff.1, theo_mudiff.2) 
  return(result)
}

# Application
Group.1 <- round(rnorm(15,5,2))
Group.2 <- round(rnorm(12,4,1))
meandiff.CI(Group.1,Group.2,.99)

#alpha <- 1 - conf.level
#n1 <- length(Group.1)
#n2 <- length(Group.2)
#mean1 <- mean(Group.1)
#mean2 <- mean(Group.2)
#sd1 <- sd(Group.1)
#sd2 <- sd(Group.2)
#pooled_sd <- sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
#SE <- pooled_sd*sqrt(1/n1+1/n2)

#(mean1-mean2)-qt(1-alpha/2,df=n1+n2-2)*SE
#(mean1-mean2)+qt(1-alpha/2,df=n1+n2-2)*SE

#--------------------------------------------------------------
#  Obtain Student confidence limits following Cohen's procedure
#--------------------------------------------------------------
Group.1 <- round(rnorm(15,5,2))
Group.2 <- round(rnorm(12,4,1))

Cohen.CI <- function(Group.1, Group.2,conf.level)
{ 
  n1 <- length(Group.1)
  n2 <- length(Group.2)
  
  #  perform two-sample t-test assuming equal variance (same assumptions of cohen's d)
  t_obs <- t.test(Group.1, Group.2, alternative = "two.sided", var.equal = TRUE)$statistic
  
  #  find limits of lambda at pt = alpha/2 or 1-alpha/2

  # lower limit = limit of lambda such as 1-pt(q=t_obs, df=n1+n2-2, ncp = lambda) = (1-conf.level)/2 = alpha/2
  f=function(lambda,rep) 1-pt(q=t_obs, df=n1+n2-2, ncp = lambda)-rep
  out=uniroot(f,c(-100,100),rep=(1-conf.level)/2)
  lambda.1 <- out$root
  delta.1 <- lambda.1*sqrt(1/n1+1/n2)   # See explanation for delta.2

    
  # upper limit = limit of lambda such as pt(q=t_obs, df=n1+n2-2, ncp = lambda) = (1-conf.level)/2 = alpha/2
  f=function(lambda,rep) pt(q=t_obs, df=n1+n2-2, ncp = lambda)-rep
  out=uniroot(f,c(-100,100),rep=(1-conf.level)/2)
  lambda.2 <- out$root
  delta.2 <- lambda.2*sqrt(1/n1+1/n2)   # Pr[t_obs(alpha/2,df=df,ncp=lambda) <= lambda <= t_obs(1-alpha/2,df=df,ncp=lambda)=.025
                                        # because lambda = delta * sqrt[n1n2/(n1+n2)] :
                                        # Pr[t_obs(alpha/2,df=df,ncp=lambda) <= delta * sqrt[n1n2/(n1+n2)] <= t_obs(1-alpha/2,df=df,ncp=lambda)=.025
                                        # Pr[t_obs(alpha/2,df=df,ncp=lambda)*sqrt[(n1+n2)/n1n2] <= delta <= t_obs(1-alpha/2,df=df,ncp=lambda)*sqrt[(n1+n2)/n1n2]=.025
                                        # and sqrt[(n1+n2)/n1n2]=sqrt(1/n1+1/n2)
  result <- c(delta.1, delta.2) 
  return(result)
}

#--------------------------------------------------------------
#  Obtain Welch confidence limits following Shieh's procedure
#--------------------------------------------------------------
 
Shieh.CL <- function(Group.1, Group.2)
{ 
  n1 <- length(Group.1)
  n2 <- length(Group.2)
  s1 <- sd(Group.1)
  s2 <- sd(Group.2)
  
  #  perform two-sample Welch t-test (same assumptions of Shieh's d)
  V0 <- t.test(Group.1, Group.2, alternative = "two.sided", var.equal = FALSE)$statistic
  
  #  sample estimates for degrees of freedom nu of noncentral t distribution
  nu <- (s1^2/n1 + s2^2/n2)^2 / ((s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1))
  
  #  obtain the point estimator of the standardized mean difference
  G.nu <- gamma(nu/2)/(sqrt((n1+n2)*nu/2)*gamma((nu-1)/2))
  delta.nu <- G.nu*V0
  
  #  find limits of lambda at pt = 0.025 or 0.975
  lambda <- 0.01
  
  if (pt(q=V0, df=nu, ncp = lambda) > 0.025) {
    while (pt(q = V0, df = nu, ncp = lambda) - 0.025 > 0.0001) {
      lambda <- lambda + 0.01
    }} else if (pt(q=V0, df= nu, ncp = lambda) < 0.025) {
      while (0.025 - pt(q = V0, df = nu, ncp = lambda) > 0.0001){
        lambda <- lambda - 0.01
      }} else {lambda <- lambda}
  
  lambda.1 <- lambda
  delta.1 <- lambda.1/sqrt(n1+n2)
  
  lambda <- 0.01
  
  if (pt(q=V0, df=nu, ncp = lambda) > 0.975) {
    while (pt(q = V0, df = nu, ncp = lambda) - 0.975 > 0.0001) {
      lambda <- lambda + 0.01
    }} else if (pt(q = V0, df = nu, ncp = lambda) < 0.975) {
      while (0.975 - pt(q = V0, df = nu, ncp = lambda) > 0.0001){
        lambda <- lambda - 0.01
      }} else {lambda <- lambda}  
  
  lambda.2 <- lambda
  delta.2 <- lambda.2/sqrt(n1+n2)
  
  delta.low <- min(delta.1, delta.2)
  delta.upp <- max(delta.1, delta.2)
  
  result <- c(delta.low, delta.upp, delta.nu) 
  return(result)
}


###################
#  Normal groups  #
###################

CI.norm <- function(Nsim=100000,n1,n2,mu1,mu2,sigma1,sigma2) 
{
  #  Non-Central CI of ds
  low.NC.ds <- c()
  upp.NC.ds <- c()
  w.NC.ds <- c()
  count.NC.ds <- 0
  count.NC.ep.ds <- 0    # to calculate empirical power
  
  #  Non-Central CI of Shieh's d
  low.shieh <- c()
  upp.shieh <- c()
  w.shieh <- c()
  count.shieh <- 0
  count.ep.shieh <- 0
  
  #  BCs CIs for ds
  low.ds <- c()
  upp.ds <- c()
  w.ds <- c()
  count.ds <- 0
  count.ep.ds <- 0
  
  #  BCs CIs for du
  low.du <- c()
  upp.du <- c()
  w.du <- c()
  count.du <- 0
  count.ep.du <- 0

  # ----- Population parameters -----
  ds.theo <- (mu1-mu2) / sqrt( ((n1-1)*sigma1^2 + (n2-1)*sigma2^2) / (n1+n2-2) )
  
  du.theo <- ds.theo * (1 - 3/(4*(n1+n2-2)-1))
  
  shieh.theo <- Unbiased.shieh(mu1,mu2,sigma1,sigma2,n1,n2)
    
  for  (j in 1:Nsim) {
    
    sample1 <- rnorm(n=n1, mean=mu1, sd=sigma1)
    sample2 <- rnorm(n=n2, mean=mu2, sd=sigma2)
    
    
    # ----- Confidence intervals -----
    
    #  NC of ds
    low.NC.ds[j] <- Cohen.CI(sample1, sample2)[1]
    upp.NC.ds[j] <- Cohen.CI(sample1, sample2)[2]
    w.NC.ds[j] <- upp.NC.ds[j] - low.NC.ds[j] 
    if ((low.NC.ds[j] < ds.theo) & (ds.theo < upp.NC.ds[j])){
      count.NC.ds <- count.NC.ds + 1} 
    if ((low.NC.ds[j] > 0 ) | (upp.NC.ds[j] < 0)){
      count.NC.ep.ds <- count.NC.ep.ds + 1} 
    
    #  NC of shieh
    low.shieh[j] <- Shieh.CL(Group.1 = sample1, Group.2 = sample2)[1]
    upp.shieh[j] <- Shieh.CL(Group.1 = sample1, Group.2 = sample2)[2]
    w.shieh[j] <- upp.shieh[j] - low.shieh[j] 
    if ((low.shieh[j] < shieh.theo) & (shieh.theo < upp.shieh[j])) {
      count.shieh <- count.shieh + 1} 
    if ((low.shieh[j] > 0 ) | (upp.shieh[j] < 0)){
      count.ep.shieh <- count.ep.shieh + 1}
    
    #  BCa of ds
    low.ds[j] <- BCa.CL.ds(Group.1 = sample1, Group.2 = sample2)[1]
    upp.ds[j] <- BCa.CL.ds(Group.1 = sample1, Group.2 = sample2)[2]
    w.ds[j] <- upp.ds[j] - low.ds[j] 
    if ((low.ds[j] < ds.theo) & (ds.theo < upp.ds[j])) {
      count.ds <- count.ds + 1} 
    if ((low.ds[j] > 0 ) | (upp.ds[j] < 0)){
      count.ep.ds <- count.ep.ds + 1}
    
    #  BCa of du/Hedge's gs
    low.du[j] <- BCa.CL.du(Group.1 = sample1, Group.2 = sample2)[1]
    upp.du[j] <- BCa.CL.du(Group.1 = sample1, Group.2 = sample2)[2]
    w.du[j] <- upp.du[j] - low.du[j] 
    if ((low.du[j] < du.theo) & (du.theo < upp.du[j])) {
      count.du <- count.du + 1} 
    if ((low.du[j] > 0 ) | (upp.du[j] < 0)){
      count.ep.du <- count.ep.du + 1}
  }
  

  # coverage probability
  cp.NC.ds <- count.NC.ds/Nsim
  cp.shieh <- count.shieh/Nsim
  cp.ds <- count.ds/Nsim
  cp.du <- count.du/Nsim

  
  # empirical power
  ep.NC.ds <- ifelse(ds.theo==0,1-count.NC.ep.ds/Nsim,count.NC.ep.ds/Nsim)
  ep.shieh <- ifelse(ds.theo==0,1-count.ep.shieh/Nsim,count.ep.shieh/Nsim)
  ep.ds <- ifelse(ds.theo==0,1-count.ep.ds/Nsim,count.ep.ds/Nsim)
  ep.du <- ifelse(ds.theo==0,1-count.ep.du/Nsim,count.ep.du/Nsim)

  
  result <- c(ds.theo,cp.NC.ds,mean(low.NC.ds),median(low.NC.ds),sd(low.NC.ds),
              mean(upp.NC.ds),median(upp.NC.ds),sd(upp.NC.ds),
              mean(w.NC.ds),median(w.NC.ds),sd(w.NC.ds),ep.NC.ds,
              
              shieh.theo,cp.shieh,mean(low.shieh),median(low.shieh),sd(low.shieh),
              mean(upp.shieh),median(upp.shieh),sd(upp.shieh),
              mean(w.shieh),median(w.shieh),sd(w.shieh),ep.shieh,
              
              ds.theo,cp.ds,mean(low.ds),median(low.ds),sd(low.ds),
              mean(upp.ds),median(upp.ds),sd(upp.ds),
              mean(w.ds),median(w.ds),sd(w.ds),ep.ds,
              
              du.theo,cp.du,mean(low.du),median(low.du),sd(low.du),
              mean(upp.du),median(upp.du),sd(upp.du),
              mean(w.du),median(w.du),sd(w.du),ep.du)
  
  out <- matrix(result,nrow=12,ncol=4) 
  colnames(out) <- c("NC_ds","NC_shieh","BCa_ds","BCa_du")
  rownames(out) <- c("Theo_value","CP",
                     "mean.low","median.low","sd.low",
                     "mean.upp","median.upp","sd.upp",
                     "mean.w","median.w","sd.w",
                     "EP")
  return(out)
}


##  Get simulations cross settings
set.seed(625)

#  N=c(10,10), Sigma = c(1,1)
cal.mu1(d=0,sigma1=1,n1=10,n2=10) # mu1=0 
CI.1 <- CI.norm(n1=10,n2=10,mu1=0,mu2=0,sigma1=1,sigma2=1)
write.csv(CI.1,file="CI.1")

cal.mu1(d=0.2,sigma1=1,n1=10,n2=10) # mu1=0.2
CI.2 <- CI.norm(n1=10,n2=10,mu1=0.2,mu2=0,sigma1=1,sigma2=1)
write.csv(CI.2,file="CI.2")

cal.mu1(d=0.5,sigma1=1,n1=10,n2=10) # mu1=0.5
CI.3 <- CI.norm(n1=10,n2=10,mu1=0.5,mu2=0,sigma1=1,sigma2=1)
write.csv(CI.3,file="CI.3")

cal.mu1(d=1,sigma1=1,n1=10,n2=10) # mu1=1
CI.4 <- CI.norm(n1=10,n2=10,mu1=1,mu2=0,sigma1=1,sigma2=1)
write.csv(CI.4,file="CI.4")

cal.mu1(d=2,sigma1=1,n1=10,n2=10) # mu1=2
CI.5 <- CI.norm(n1=10,n2=10,mu1=2,mu2=0,sigma1=1,sigma2=1)
write.csv(CI.5,file="CI.5")

#  N=c(10,10), Sigma = c(2,1)
CI.6 <- CI.norm(n1=10,n2=10,mu1=0,mu2=0,sigma1=2,sigma2=1)
write.csv(CI.6,file="CI.6")

mu1 <- cal.mu1(d=0.2,sigma1=2,n1=10,n2=10) 
CI.7 <- CI.norm(n1=10,n2=10,mu1=mu1,mu2=0,sigma1=2,sigma2=1)
write.csv(CI.7,file="CI.7")

mu1 <- cal.mu1(d=0.5,sigma1=2,n1=10,n2=10) 
CI.8 <- CI.norm(n1=10,n2=10,mu1=mu1,mu2=0,sigma1=2,sigma2=1)
write.csv(CI.8,file="CI.8")

mu1 <- cal.mu1(d=1,sigma1=2,n1=10,n2=10) 
CI.9 <- CI.norm(n1=10,n2=10,mu1=mu1,mu2=0,sigma1=2,sigma2=1)
write.csv(CI.9,file="CI.9")

mu1 <- cal.mu1(d=2,sigma1=2,n1=10,n2=10) 
CI.10 <- CI.norm(n1=10,n2=10,mu1=mu1,mu2=0,sigma1=2,sigma2=1)
write.csv(CI.10,file="CI.10")

#  N=c(20,10), Sigma = c(1,1)
mu1 <- cal.mu1(d=0,sigma1=1,n1=20,n2=10)
CI.11 <- CI.norm(n1=20,n2=10,mu1=mu1,mu2=0,sigma1=1,sigma2=1)
write.csv(CI.11,file="CI.11")

mu1 <- cal.mu1(d=0.2,sigma1=1,n1=20,n2=10)
CI.12 <- CI.norm(n1=20,n2=10,mu1=mu1,mu2=0,sigma1=1,sigma2=1)
write.csv(CI.12,file="CI.12")

mu1 <- cal.mu1(d=0.5,sigma1=1,n1=20,n2=10)
CI.13 <- CI.norm(n1=20,n2=10,mu1=mu1,mu2=0,sigma1=1,sigma2=1)
write.csv(CI.13,file="CI.13")

mu1 <- cal.mu1(d=1,sigma1=1,n1=20,n2=10)
CI.14 <- CI.norm(n1=20,n2=10,mu1=mu1,mu2=0,sigma1=1,sigma2=1)
write.csv(CI.14,file="CI.14")

mu1 <- cal.mu1(d=2,sigma1=1,n1=20,n2=10)
CI.15 <- CI.norm(n1=20,n2=10,mu1=mu1,mu2=0,sigma1=1,sigma2=1)
write.csv(CI.15,file="CI.15")

#  N=c(20,10), Sigma = c(2,1)
mu1 <- cal.mu1(d=0,sigma1=2,n1=20,n2=10)
CI.16 <- CI.norm(n1=20,n2=10,mu1=mu1,mu2=0,sigma1=2,sigma2=1)
write.csv(CI.16,file="CI.16")

mu1 <- cal.mu1(d=0.2,sigma1=2,n1=20,n2=10) 
CI.17 <- CI.norm(n1=20,n2=10,mu1=mu1,mu2=0,sigma1=2,sigma2=1)
write.csv(CI.17,file="CI.17")

mu1 <- cal.mu1(d=0.5,sigma1=2,n1=20,n2=10) 
CI.18 <- CI.norm(n1=20,n2=10,mu1=mu1,mu2=0,sigma1=2,sigma2=1)
write.csv(CI.18,file="CI.18")

mu1 <- cal.mu1(d=1,sigma1=2,n1=20,n2=10) 
CI.19 <- CI.norm(n1=20,n2=10,mu1=mu1,mu2=0,sigma1=2,sigma2=1)
write.csv(CI.19,file="CI.19")

mu1 <- cal.mu1(d=2,sigma1=2,n1=20,n2=10) 
CI.20 <- CI.norm(n1=20,n2=10,mu1=mu1,mu2=0,sigma1=2,sigma2=1)
write.csv(CI.20,file="CI.20")

#  N=c(75,50), sigma = c(1,1)
mu1 <- cal.mu1(d=0,sigma1=1,n1=75,n2=50)
CI.21 <- CI.norm(n1=75,n2=50,mu1=mu1,mu2=0,sigma1=1,sigma2=1)
write.csv(CI.21,file="CI.21")

mu1 <- cal.mu1(d=0.2,sigma1=1,n1=75,n2=50)
CI.22 <- CI.norm(n1=75,n2=50,mu1=mu1,mu2=0,sigma1=1,sigma2=1)
write.csv(CI.22,file="CI.22")

mu1 <- cal.mu1(d=0.5,sigma1=1,n1=75,n2=50)
CI.23 <- CI.norm(n1=75,n2=50,mu1=mu1,mu2=0,sigma1=1,sigma2=1)
write.csv(CI.23,file="CI.23")

mu1 <- cal.mu1(d=1,sigma1=1,n1=75,n2=50)
CI.24 <- CI.norm(n1=75,n2=50,mu1=mu1,mu2=0,sigma1=1,sigma2=1)
write.csv(CI.24,file="CI.24")

mu1 <- cal.mu1(d=2,sigma1=1,n1=75,n2=50)
CI.25 <- CI.norm(n1=75,n2=50,mu1=mu1,mu2=0,sigma1=1,sigma2=1)
write.csv(CI.25,file="CI.25")

#  N=c(75,50), sigma = c(2,1)
mu1 <- cal.mu1(d=0,sigma1=2,n1=75,n2=50)
CI.26 <- CI.norm(n1=75,n2=50,mu1=mu1,mu2=0,sigma1=2,sigma2=1)
write.csv(CI.26,file="CI.26")

mu1 <- cal.mu1(d=0.2,sigma1=2,n1=75,n2=50)
CI.27 <- CI.norm(n1=75,n2=50,mu1=mu1,mu2=0,sigma1=2,sigma2=1)
write.csv(CI.27,file="CI.27")

mu1 <- cal.mu1(d=0.5,sigma1=2,n1=75,n2=50)
CI.28 <- CI.norm(n1=75,n2=50,mu1=mu1,mu2=0,sigma1=2,sigma2=1)
write.csv(CI.28,file="CI.28")

mu1 <- cal.mu1(d=1,sigma1=2,n1=75,n2=50)
CI.29 <- CI.norm(n1=75,n2=50,mu1=mu1,mu2=0,sigma1=2,sigma2=1)
write.csv(CI.29,file="CI.29")

set.seed(625)
mu1 <- cal.mu1(d=2,sigma1=2,n1=75,n2=50)
CI.30 <- CI.norm(n1=75,n2=50,mu1=mu1,mu2=0,sigma1=2,sigma2=1)
write.csv(CI.30,file="CI.30")

#  read files 
#  table 4.5 in appendix
CI.1 <- read.csv("CI.1")
CI.2 <- read.csv("CI.2")
CI.3 <- read.csv("CI.3")
CI.4 <- read.csv("CI.4")
CI.5 <- read.csv("CI.5")

#  table 4.6 in appendix
CI.6 <- read.csv("CI.6")
CI.7 <- read.csv("CI.7")
CI.8 <- read.csv("CI.8")
CI.9 <- read.csv("CI.9")
CI.10 <- read.csv("CI.10")

#  table 4.7 in appendix
CI.11 <- read.csv("CI.11")
CI.12 <- read.csv("CI.12")
CI.13 <- read.csv("CI.13")
CI.14 <- read.csv("CI.14")
CI.15 <- read.csv("CI.15")

#  table 4.8 in appendix
CI.16 <- read.csv("CI.16")
CI.17 <- read.csv("CI.17")
CI.18 <- read.csv("CI.18")
CI.19 <- read.csv("CI.19")
CI.20 <- read.csv("CI.20")

#  table 4.9 in appendix
CI.21 <- read.csv("CI.21")
CI.22 <- read.csv("CI.22")
CI.23 <- read.csv("CI.23")
CI.24 <- read.csv("CI.24")
CI.25 <- read.csv("CI.25")

#  table 4.10 in appendix
CI.26 <- read.csv("CI.26")
CI.27 <- read.csv("CI.27")
CI.28 <- read.csv("CI.28")
CI.29 <- read.csv("CI.29")
CI.30 <- read.csv("CI.30")



#######################
#  Non-Normal groups  #
#######################

CI.snorm <- function(Nsim=100000,n1,n2,mu1,mu2,sigma1,sigma2,xi1,xi2) 
{
  #  Non-Central CI of ds
  low.NC.ds <- c()
  upp.NC.ds <- c()
  w.NC.ds <- c()
  count.NC.ds <- 0
  count.NC.ep.ds <- 0    # to calculate empirical power
  
  #  Non-Central CI of Shieh's d
  low.shieh <- c()
  upp.shieh <- c()
  w.shieh <- c()
  count.shieh <- 0
  count.ep.shieh <- 0
  
  #  BCs CIs for ds
  low.ds <- c()
  upp.ds <- c()
  w.ds <- c()
  count.ds <- 0
  count.ep.ds <- 0
  
  #  BCs CIs for du
  low.du <- c()
  upp.du <- c()
  w.du <- c()
  count.du <- 0
  count.ep.du <- 0
  
  # ----- Population parameters -----
  ds.theo <- (mu1-mu2) / sqrt( ((n1-1)*sigma1^2 + (n2-1)*sigma2^2) / (n1+n2-2) )
  
  du.theo <- ds.theo * (1 - 3/(4*(n1+n2-2)-1))
  
  shieh.theo <- Unbiased.shieh(mu1,mu2,sigma1,sigma2,n1,n2)
  
  for  (j in 1:Nsim) {
    
    sample1 <- rsnorm(n=n1, mean=mu1, sd=sigma1, xi=xi1)
    sample2 <- rsnorm(n=n2, mean=mu2, sd=sigma2, xi=xi2)
    
    
    # ----- Confidence intervals -----
    
    #  NC of ds
    low.NC.ds[j] <- Cohen.CI(sample1, sample2)[1]
    upp.NC.ds[j] <- Cohen.CI(sample1, sample2)[2]
    w.NC.ds[j] <- upp.NC.ds[j] - low.NC.ds[j] 
    if ((low.NC.ds[j] < ds.theo) & (ds.theo < upp.NC.ds[j])){
      count.NC.ds <- count.NC.ds + 1} 
    if ((low.NC.ds[j] > 0 ) | (upp.NC.ds[j] < 0)){
      count.NC.ep.ds <- count.NC.ep.ds + 1} 
    
    #  NC of shieh
    low.shieh[j] <- Shieh.CL(Group.1 = sample1, Group.2 = sample2)[1]
    upp.shieh[j] <- Shieh.CL(Group.1 = sample1, Group.2 = sample2)[2]
    w.shieh[j] <- upp.shieh[j] - low.shieh[j] 
    if ((low.shieh[j] < shieh.theo) & (shieh.theo < upp.shieh[j])) {
      count.shieh <- count.shieh + 1} 
    if ((low.shieh[j] > 0 ) | (upp.shieh[j] < 0)){
      count.ep.shieh <- count.ep.shieh + 1}
    
    #  BCa of ds
    low.ds[j] <- BCa.CL.ds(Group.1 = sample1, Group.2 = sample2)[1]
    upp.ds[j] <- BCa.CL.ds(Group.1 = sample1, Group.2 = sample2)[2]
    w.ds[j] <- upp.ds[j] - low.ds[j] 
    if ((low.ds[j] < ds.theo) & (ds.theo < upp.ds[j])) {
      count.ds <- count.ds + 1} 
    if ((low.ds[j] > 0 ) | (upp.ds[j] < 0)){
      count.ep.ds <- count.ep.ds + 1}
    
    #  BCa of du
    low.du[j] <- BCa.CL.du(Group.1 = sample1, Group.2 = sample2)[1]
    upp.du[j] <- BCa.CL.du(Group.1 = sample1, Group.2 = sample2)[2]
    w.du[j] <- upp.du[j] - low.du[j] 
    if ((low.du[j] < du.theo) & (du.theo < upp.du[j])) {
      count.du <- count.du + 1} 
    if ((low.du[j] > 0 ) | (upp.du[j] < 0)){
      count.ep.du <- count.ep.du + 1}
  }
  
  
  # coverage probability
  cp.NC.ds <- count.NC.ds/Nsim
  cp.shieh <- count.shieh/Nsim
  cp.ds <- count.ds/Nsim
  cp.du <- count.du/Nsim
  
  
  # empirical power
  ep.NC.ds <- ifelse(ds.theo==0,1-count.NC.ep.ds/Nsim,count.NC.ep.ds/Nsim)
  ep.shieh <- ifelse(ds.theo==0,1-count.ep.shieh/Nsim,count.ep.shieh/Nsim)
  ep.ds <- ifelse(ds.theo==0,1-count.ep.ds/Nsim,count.ep.ds/Nsim)
  ep.du <- ifelse(ds.theo==0,1-count.ep.du/Nsim,count.ep.du/Nsim)
  
  
  result <- c(ds.theo,cp.NC.ds,mean(low.NC.ds),median(low.NC.ds),sd(low.NC.ds),
              mean(upp.NC.ds),median(upp.NC.ds),sd(upp.NC.ds),
              mean(w.NC.ds),median(w.NC.ds),sd(w.NC.ds),ep.NC.ds,
              
              shieh.theo,cp.shieh,mean(low.shieh),median(low.shieh),sd(low.shieh),
              mean(upp.shieh),median(upp.shieh),sd(upp.shieh),
              mean(w.shieh),median(w.shieh),sd(w.shieh),ep.shieh,
              
              ds.theo,cp.ds,mean(low.ds),median(low.ds),sd(low.ds),
              mean(upp.ds),median(upp.ds),sd(upp.ds),
              mean(w.ds),median(w.ds),sd(w.ds),ep.ds,
              
              du.theo,cp.du,mean(low.du),median(low.du),sd(low.du),
              mean(upp.du),median(upp.du),sd(upp.du),
              mean(w.du),median(w.du),sd(w.du),ep.du)
  
  out <- matrix(result,nrow=12,ncol=4) 
  colnames(out) <- c("NC_ds","NC_shieh","BCa_ds","BCa_du")
  rownames(out) <- c("Theo_value","CP",
                     "mean.low","median.low","sd.low",
                     "mean.upp","median.upp","sd.upp",
                     "mean.w","median.w","sd.w",
                     "EP")
  return(out)
}


##  Get simulations cross settings
set.seed(625)

#  N=c(10,10), Sigma = c(1,1)
CI.snorm.1 <- CI.snorm(n1=10,n2=10,mu1=0,mu2=0,sigma1=1,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.1,file="CI.snorm.1")

CI.snorm.2 <- CI.snorm(n1=10,n2=10,mu1=0.2,mu2=0,sigma1=1,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.2,file="CI.snorm.2")

CI.snorm.3 <- CI.snorm(n1=10,n2=10,mu1=0.5,mu2=0,sigma1=1,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.3,file="CI.snorm.3")

CI.snorm.4 <- CI.snorm(n1=10,n2=10,mu1=1,mu2=0,sigma1=1,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.4,file="CI.snorm.4")

CI.snorm.5 <- CI.snorm(n1=10,n2=10,mu1=2,mu2=0,sigma1=1,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.5,file="CI.snorm.5")

#  N=c(10,10), Sigma = c(2,1)
CI.snorm.6 <- CI.snorm(n1=10,n2=10,mu1=0,mu2=0,sigma1=2,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.6,file="CI.snorm.6")

mu1 <- cal.mu1(d=0.2,sigma1=2,n1=10,n2=10) 
CI.snorm.7 <- CI.snorm(n1=10,n2=10,mu1=mu1,mu2=0,sigma1=2,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.7,file="CI.snorm.7")

mu1 <- cal.mu1(d=0.5,sigma1=2,n1=10,n2=10) 
CI.snorm.8 <- CI.snorm(n1=10,n2=10,mu1=mu1,mu2=0,sigma1=2,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.8,file="CI.snorm.8")

mu1 <- cal.mu1(d=1,sigma1=2,n1=10,n2=10) 
CI.snorm.9 <- CI.snorm(n1=10,n2=10,mu1=mu1,mu2=0,sigma1=2,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.9,file="CI.snorm.9")

mu1 <- cal.mu1(d=2,sigma1=2,n1=10,n2=10) 
CI.snorm.10 <- CI.snorm(n1=10,n2=10,mu1=mu1,mu2=0,sigma1=2,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.10,file="CI.snorm.10")

#  N=c(20,10), Sigma = c(1,1)
mu1 <- cal.mu1(d=0,sigma1=1,n1=20,n2=10)
CI.snorm.11 <- CI.snorm(n1=20,n2=10,mu1=mu1,mu2=0,sigma1=1,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.11,file="CI.snorm.11")

mu1 <- cal.mu1(d=0.2,sigma1=1,n1=20,n2=10)
CI.snorm.12 <- CI.snorm(n1=20,n2=10,mu1=mu1,mu2=0,sigma1=1,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.12,file="CI.snorm.12")

mu1 <- cal.mu1(d=0.5,sigma1=1,n1=20,n2=10)
CI.snorm.13 <- CI.snorm(n1=20,n2=10,mu1=mu1,mu2=0,sigma1=1,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.13,file="CI.snorm.13")

mu1 <- cal.mu1(d=1,sigma1=1,n1=20,n2=10)
CI.snorm.14 <- CI.snorm(n1=20,n2=10,mu1=mu1,mu2=0,sigma1=1,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.14,file="CI.snorm.14")

mu1 <- cal.mu1(d=2,sigma1=1,n1=20,n2=10)
CI.snorm.15 <- CI.snorm(n1=20,n2=10,mu1=mu1,mu2=0,sigma1=1,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.15,file="CI.snorm.15")

#  N=c(20,10), Sigma = c(2,1)
mu1 <- cal.mu1(d=0,sigma1=2,n1=20,n2=10)
CI.snorm.16 <- CI.snorm(n1=20,n2=10,mu1=mu1,mu2=0,sigma1=2,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.16,file="CI.snorm.16")

mu1 <- cal.mu1(d=0.2,sigma1=2,n1=20,n2=10) 
CI.snorm.17 <- CI.snorm(n1=20,n2=10,mu1=mu1,mu2=0,sigma1=2,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.17,file="CI.snorm.17")

mu1 <- cal.mu1(d=0.5,sigma1=2,n1=20,n2=10) 
CI.snorm.18 <- CI.snorm(n1=20,n2=10,mu1=mu1,mu2=0,sigma1=2,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.18,file="CI.snorm.18")

mu1 <- cal.mu1(d=1,sigma1=2,n1=20,n2=10) 
CI.snorm.19 <- CI.snorm(n1=20,n2=10,mu1=mu1,mu2=0,sigma1=2,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.19,file="CI.snorm.19")

mu1 <- cal.mu1(d=2,sigma1=2,n1=20,n2=10) 
CI.snorm.20 <- CI.snorm(n1=20,n2=10,mu1=mu1,mu2=0,sigma1=2,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.20,file="CI.snorm.20")

#  N=c(75,50), sigma = c(1,1)
mu1 <- cal.mu1(d=0,sigma1=1,n1=75,n2=50)
CI.snorm.21 <- CI.snorm(n1=75,n2=50,mu1=mu1,mu2=0,sigma1=1,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.21,file="CI.snorm.21")

mu1 <- cal.mu1(d=0.2,sigma1=1,n1=75,n2=50)
CI.snorm.22 <- CI.snorm(n1=75,n2=50,mu1=mu1,mu2=0,sigma1=1,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.22,file="CI.snorm.22")

set.seed(625)
mu1 <- cal.mu1(d=0.5,sigma1=1,n1=75,n2=50)
CI.snorm.23 <- CI.snorm(n1=75,n2=50,mu1=mu1,mu2=0,sigma1=1,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.23,file="CI.snorm.23")

set.seed(625)
mu1 <- cal.mu1(d=1,sigma1=1,n1=75,n2=50)
CI.snorm.24 <- CI.snorm(n1=75,n2=50,mu1=mu1,mu2=0,sigma1=1,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.24,file="CI.snorm.24")

set.seed(625)
mu1 <- cal.mu1(d=2,sigma1=1,n1=75,n2=50)
CI.snorm.25 <- CI.snorm(n1=75,n2=50,mu1=mu1,mu2=0,sigma1=1,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.25,file="CI.snorm.25")

#  N=c(75,50), sigma = c(2,1)
set.seed(625)
mu1 <- cal.mu1(d=0,sigma1=2,n1=75,n2=50)
CI.snorm.26 <- CI.snorm(n1=75,n2=50,mu1=mu1,mu2=0,sigma1=2,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.26,file="CI.snorm.26")

set.seed(625)
mu1 <- cal.mu1(d=0.2,sigma1=2,n1=75,n2=50)
CI.snorm.27 <- CI.snorm(n1=75,n2=50,mu1=mu1,mu2=0,sigma1=2,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.27,file="CI.snorm.27")

set.seed(625)
mu1 <- cal.mu1(d=0.5,sigma1=2,n1=75,n2=50)
CI.snorm.28 <- CI.snorm(n1=75,n2=50,mu1=mu1,mu2=0,sigma1=2,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.28,file="CI.snorm.28")

set.seed(625)
mu1 <- cal.mu1(d=1,sigma1=2,n1=75,n2=50)
CI.snorm.29 <- CI.snorm(n1=75,n2=50,mu1=mu1,mu2=0,sigma1=2,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.29,file="CI.snorm.29")

set.seed(625)
mu1 <- cal.mu1(d=2,sigma1=2,n1=75,n2=50)
CI.snorm.30 <- CI.snorm(n1=75,n2=50,mu1=mu1,mu2=0,sigma1=2,sigma2=1,xi1=-1,xi2=2)
write.csv(CI.snorm.30,file="CI.snorm.30")


#  read files
#  table 4.11 in appendix
CI.snorm.1 <- read.csv("CI.snorm.1")
CI.snorm.2 <- read.csv("CI.snorm.2")
CI.snorm.3 <- read.csv("CI.snorm.3")
CI.snorm.4 <- read.csv("CI.snorm.4")
CI.snorm.5 <- read.csv("CI.snorm.5")

#  table 4.12 in appendix
CI.snorm.6 <- read.csv("CI.snorm.6")
CI.snorm.7 <- read.csv("CI.snorm.7")
CI.snorm.8 <- read.csv("CI.snorm.8")
CI.snorm.9 <- read.csv("CI.snorm.9")
CI.snorm.10 <- read.csv("CI.snorm.10")

#  table 4.13 in appendix
CI.snorm.11 <- read.csv("CI.snorm.11")
CI.snorm.12 <- read.csv("CI.snorm.12")
CI.snorm.13 <- read.csv("CI.snorm.13")
CI.snorm.14 <- read.csv("CI.snorm.14")
CI.snorm.15 <- read.csv("CI.snorm.15")

#  table 4.14 in appendix
CI.snorm.16 <- read.csv("CI.snorm.16")
CI.snorm.17 <- read.csv("CI.snorm.17")
CI.snorm.18 <- read.csv("CI.snorm.18")
CI.snorm.19 <- read.csv("CI.snorm.19")
CI.snorm.20 <- read.csv("CI.snorm.20")

#  table 4.15 in appendix
CI.snorm.21 <- read.csv("CI.snorm.21")
CI.snorm.22 <- read.csv("CI.snorm.22")
CI.snorm.23 <- read.csv("CI.snorm.23")
CI.snorm.24 <- read.csv("CI.snorm.24")
CI.snorm.25 <- read.csv("CI.snorm.25")

#  table 4.16 in appendix
CI.snorm.26 <- read.csv("CI.snorm.26")
CI.snorm.27 <- read.csv("CI.snorm.27")
CI.snorm.28 <- read.csv("CI.snorm.28")
CI.snorm.29 <- read.csv("CI.snorm.29")
CI.snorm.30 <- read.csv("CI.snorm.30")



###################################
#  Plot the coverage - figure 4.5 #
###################################

coverage <- function(data){
  out <- data[[1]][2,]
  for (i in 1:(length(data)-1)){
    out <- rbind(out,data[[i+1]][2,])
    i <- i+1
  }
  return(out)
}

#  plot for norm
# noquote(paste0("CI.",c(1:30),sep=","))
data <- list(CI.1,  CI.2,  CI.3,  CI.4,  CI.5,  
             CI.11, CI.12, CI.13, CI.14, CI.15, 
             CI.21, CI.22, CI.23, CI.24, CI.25, 
             CI.6,  CI.7,  CI.8,  CI.9,  CI.10,
             CI.16, CI.17, CI.18, CI.19, CI.20, 
             CI.26, CI.27, CI.28, CI.29, CI.30)
ci_norm <- coverage(data)[,-1]

pdf("ci.norm.pdf",width = 6, height = 4)
barplot(t(ci_norm-0.85), beside=TRUE,
        ylim=c(0,0.15))
dev.off()


#  plot for snorm
# noquote(paste0("CI.snorm.",c(1:30),sep=","))
data <- list(CI.snorm.1,  CI.snorm.2,  CI.snorm.3,  CI.snorm.4,  CI.snorm.5,  
             CI.snorm.11, CI.snorm.12, CI.snorm.13, CI.snorm.14, CI.snorm.15, 
             CI.snorm.21, CI.snorm.22, CI.snorm.23, CI.snorm.24, CI.snorm.25, 
             CI.snorm.6,  CI.snorm.7,  CI.snorm.8,  CI.snorm.9,  CI.snorm.10, 
             CI.snorm.16, CI.snorm.17, CI.snorm.18, CI.snorm.19, CI.snorm.20,
             CI.snorm.26, CI.snorm.27, CI.snorm.28, CI.snorm.29, CI.snorm.30)
ci_snorm <-  coverage(data)[,-1]

pdf("ci.snorm.pdf",width = 6, height = 4)
barplot(t(ci_snorm-0.85), beside=TRUE,
        ylim=c(0,0.15))
dev.off()


#####################################
#  Plot the mean width - figure 4.6 #
#####################################

width <- function(data){
  out <- data[[1]][9,]
  for (i in 1:(length(data)-1)){
    out <- rbind(out,data[[i+1]][9,])
    i <- i+1
  }
  return(out)
}

#  plot for norm
data <- list(CI.1,  CI.2,  CI.3,  CI.4,  CI.5,  
             CI.11, CI.12, CI.13, CI.14, CI.15, 
             CI.21, CI.22, CI.23, CI.24, CI.25, 
             CI.6,  CI.7,  CI.8,  CI.9,  CI.10,
             CI.16, CI.17, CI.18, CI.19, CI.20, 
             CI.26, CI.27, CI.28, CI.29, CI.30)
width_norm <- width(data)[,-1]

pdf("width.norm.pdf",width = 6, height = 4)
barplot(t(width_norm), beside=TRUE,
        ylim=c(0,3))
dev.off()


#  plot for snorm
# noquote(paste0("CI.snorm.",c(1:30),sep=","))
data <- list(CI.snorm.1,  CI.snorm.2,  CI.snorm.3,  CI.snorm.4,  CI.snorm.5,  
             CI.snorm.11, CI.snorm.12, CI.snorm.13, CI.snorm.14, CI.snorm.15, 
             CI.snorm.21, CI.snorm.22, CI.snorm.23, CI.snorm.24, CI.snorm.25, 
             CI.snorm.6,  CI.snorm.7,  CI.snorm.8,  CI.snorm.9,  CI.snorm.10, 
             CI.snorm.16, CI.snorm.17, CI.snorm.18, CI.snorm.19, CI.snorm.20,
             CI.snorm.26, CI.snorm.27, CI.snorm.28, CI.snorm.29, CI.snorm.30)
width_snorm <-  width(data)[,-1]

pdf("width.snorm.pdf",width = 6, height = 4)
barplot(t(width_snorm), beside=TRUE,
        ylim=c(0,3))
dev.off()





















