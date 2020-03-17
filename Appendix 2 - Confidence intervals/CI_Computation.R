for (package in "rootSolve") {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}


#########################################
##    Compute Confidences intervals    ##
#########################################

#  Obtain confidence limits for mu1-mu2 under the assumption of homoscedasticity
#-------------------------------------------------------------------------------

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

# Check: the method returns approximately same CI as the classical method based on pivotal quantity
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

Group.1 <- rnorm(20,5,2)
Group.2 <- rnorm(24,4,1)
Cohen.CI(Group.1, Group.2,conf.level=.95)
  
#--------------------------------------------------------------
#  Obtain Welch confidence limits following Shieh's procedure
#--------------------------------------------------------------
 
Shieh.CL <- function(Group.1, Group.2,conf.level)
{ 
  n1 <- length(Group.1)
  n2 <- length(Group.2)
  s1 <- sd(Group.1)
  s2 <- sd(Group.2)
  
  #  perform two-sample Welch t-test (same assumptions of Shieh's d)
  test <- t.test(Group.1, Group.2, alternative = "two.sided", var.equal = FALSE)
  w_obs <- test$statistic
  #  sample estimates for degrees of freedom DF of noncentral t distribution
  DF <- test$parameter # = (sd1^2/n1 + sd2^2/n2)^2 / ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))

  #  find limits of lambda at pt = alpha/2 or 1-alpha/2
  
  # lower limit = limit of lambda such as 1-pt(q=t_obs, df=n1+n2-2, ncp = lambda) = (1-conf.level)/2 = alpha/2
  f=function(lambda,rep) 1-pt(q=w_obs, df=DF, ncp = lambda)-rep
  out=uniroot(f,c(-100,100),rep=(1-conf.level)/2)
  lambda.1 <- out$root
  delta.1 <- lambda.1/sqrt(n1+n2)

  # upper limit = limit of lambda such as pt(q=t_obs, df=n1+n2-2, ncp = lambda) = (1-conf.level)/2 = alpha/2
  f=function(lambda,rep) pt(q=w_obs, df=DF, ncp = lambda)-rep
  out=uniroot(f,c(-100,100),rep=(1-conf.level)/2)
  lambda.2 <- out$root
  delta.2 <- lambda.2/sqrt(n1+n2)

  result <- c(delta.1, delta.2) 
  return(result)
}


















