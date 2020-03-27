for (package in "rootSolve") {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}
 

#########################################
##    Compute Confidences intervals    ##
#########################################

#  Obtain confidence limits for mu1-mu2 

meandiff.CI <- function(Group.1, Group.2,conf.level=.95,var.equal=FALSE,alternative){ # alternative="two.sided", less or "greater"
  n1 <- length(Group.1)
  n2 <- length(Group.2)
  mean1 <- mean(Group.1)
  mean2 <- mean(Group.2)
  sd1 <- sd(Group.1)
  sd2 <- sd(Group.2)
  
  if(alternative=="two.sided"){
    
    if (var.equal==TRUE){
      pooled_sd <- sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
      SE <- pooled_sd*sqrt(1/n1+1/n2) # standard error
      
      # Lower limit = limit of mu1-mu2 such as 1-pt(q=t_obs, df=n1+n2-2) = (1-conf.level)/2 = alpha/2
      # with t_obs = ((mean1-mean2)-(theo_mudiff))/SE
      
      f=function(theo_mudiff,rep) 1-pt(q=((mean1-mean2)-(theo_mudiff))/SE, df=n1+n2-2)-rep
      out=uniroot(f,c(0,2),rep=(1-conf.level)/2,extendInt = "yes")
      theo_mudiff.1 <- out$root
      
      # upper limit = limit of mu1-mu2 such as pt(q=t_obs, df=n1+n2-2) = (1-conf.level)/2 = alpha/2
      # with t_obs = ((mean1-mean2)-(theo_mudiff))/SE
      
      f=function(theo_mudiff,rep) pt(q=((mean1-mean2)-(theo_mudiff))/SE, df=n1+n2-2)-rep
      out=uniroot(f,c(0,2),rep=(1-conf.level)/2,extendInt = "yes")
      theo_mudiff.2 <- out$root
      
      result <- c(theo_mudiff.1, theo_mudiff.2) 
      
    } else if (var.equal==FALSE){
      
      SE <- sqrt(sd1^2/n1+sd2^2/n2) # standard error
      
      test <- t.test(Group.1, Group.2, alternative = "two.sided", var.equal = FALSE)
      DF <- test$parameter 
      
      # Lower limit = limit of mu1-mu2 such as 1-pt(q=t_obs, df=DF) = (1-conf.level)/2 = alpha/2
      # with t_obs = ((mean1-mean2)-(theo_mudiff))/SE
      # and DF = (sd1^2/n1 + sd2^2/n2)^2 / ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))
      
      f=function(theo_mudiff,rep) 1-pt(q=((mean1-mean2)-(theo_mudiff))/SE, df=DF)-rep
      out=uniroot(f,c(0,2),rep=(1-conf.level)/2,extendInt = "yes")
      
      theo_mudiff.1 <- out$root
      
      # upper limit = limit of mu1-mu2 such as pt(q=t_obs, df=DF) = (1-conf.level)/2 = alpha/2
      # with t_obs = ((mean1-mean2)-(theo_mudiff))/SE
      # and DF = (sd1^2/n1 + sd2^2/n2)^2 / ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))
      
      f=function(theo_mudiff,rep) pt(q=((mean1-mean2)-(theo_mudiff))/SE, df=DF)-rep
      out=uniroot(f,c(0,2),rep=(1-conf.level)/2,extendInt = "yes")
      theo_mudiff.2 <- out$root
      
      result <- c(theo_mudiff.1, theo_mudiff.2) 
    }
    
  } else if (alternative == "greater"){
    
    if (var.equal==TRUE){
      pooled_sd <- sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
      SE <- pooled_sd*sqrt(1/n1+1/n2)
      
      
      # Lower limit = limit of mu1-mu2 such as 1-pt(q=t_obs, df=n1+n2-2) = (1-conf.level) = alpha
      # with t_obs = ((mean1-mean2)-(theo_mudiff))/SE
      
      f=function(theo_mudiff,rep) 1-pt(q=((mean1-mean2)-(theo_mudiff))/SE, df=n1+n2-2)-rep
      out=uniroot(f,c(0,2),rep=1-conf.level,extendInt = "yes")
      
      theo_mudiff.1 <- out$root
      
      # upper limit = +Inf
      
      theo_mudiff.2 <- Inf
      
      result <- c(theo_mudiff.1, theo_mudiff.2) 
      
    } else if (var.equal==FALSE) {
      SE <- sqrt(sd1^2/n1+sd2^2/n2)
      
      test <- t.test(Group.1, Group.2, alternative = "greater", var.equal = FALSE)
      DF <- test$parameter
      
      # Lower limit = limit of mu1-mu2 such as 1-pt(q=t_obs, df=DF) = (1-conf.level) = alpha
      # with t_obs = ((mean1-mean2)-(theo_mudiff))/SE
      # and DF = (sd1^2/n1 + sd2^2/n2)^2 / ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))
      
      f=function(theo_mudiff,rep) 1-pt(q=((mean1-mean2)-(theo_mudiff))/SE, df=DF)-rep
      out=uniroot(f,c(0,2),rep=1-conf.level,extendInt = "yes")
      theo_mudiff.1 <- out$root
      
      # Upper limit = +Inf
      
      theo_mudiff.2 <-  Inf
      
      result <- c(theo_mudiff.1, theo_mudiff.2) 
    }
    
    
  } else if (alternative == "less"){
    
    if (var.equal==TRUE){
      pooled_sd <- sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
      SE <- pooled_sd*sqrt(1/n1+1/n2)
      
      
      # Lower limit = - inf
      
      theo_mudiff.1 <- -Inf
      
      # Upper limit = limit of mu1-mu2 such as pt(q=t_obs, df=n1+n2-2) = (1-conf.level) = alpha
      # with t_obs = ((mean1-mean2)-(theo_mudiff))/SE
      
      f=function(theo_mudiff,rep) pt(q=((mean1-mean2)-(theo_mudiff))/SE, df=n1+n2-2)-rep
      out=uniroot(f,c(0,2),rep=1-conf.level,extendInt = "yes")
      theo_mudiff.2 <- out$root
      
      result <- c(theo_mudiff.1, theo_mudiff.2) 
      
    } else if (var.equal==FALSE) {
      SE <- sqrt(sd1^2/n1+sd2^2/n2)
      
      test <- t.test(Group.1, Group.2, alternative = "less", var.equal = FALSE)
      DF <- test$parameter    
      
      # Lower limit = - inf
      
      theo_mudiff.1 <- -Inf
      
      # Upper limit = limit of mu1-mu2 such as pt(q=t_obs, df=DF) = (1-conf.level) = alpha
      # with t_obs = ((mean1-mean2)-(theo_mudiff))/SE
      # and DF = (sd1^2/n1 + sd2^2/n2)^2 / ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))
      
      f=function(theo_mudiff,rep) pt(q=((mean1-mean2)-(theo_mudiff))/SE, df=DF)-rep
      out=uniroot(f,c(0,2),rep=1-conf.level,extendInt = "yes")
      theo_mudiff.2 <- out$root
      
      result <- c(theo_mudiff.1, theo_mudiff.2) 
    }
    
  }
  
  return(result)
}

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

#SE <- sqrt(sd1^2/n1+sd2^2/n2)
#(mean1-mean2)-qt(1-alpha/2,df=DF)*SE
#(mean1-mean2)+qt(1-alpha/2,df=DF)*SE


#--------------------------------------------------------------
#  Obtain Student confidence limits following Cohen's procedure
#--------------------------------------------------------------

Cohen.CI <- function(Group.1, Group.2,conf.level,alternative="two.sided") #alternative = "two.sided", "greater" or "less"
{ 
  n1 <- length(Group.1)
  n2 <- length(Group.2)
  
  if(alternative=="two.sided"){

    #  perform two-sample t-test assuming equal variance (same assumptions of cohen's d)
    t_obs <- t.test(Group.1, Group.2, alternative = "two.sided", var.equal = TRUE)$statistic
    
    # lower limit = limit of lambda such as 1-pt(q=t_obs, df=n1+n2-2, ncp = lambda) = (1-conf.level)/2 = alpha/2
    f=function(lambda,rep) 1-pt(q=t_obs, df=n1+n2-2, ncp = lambda)-rep
    out=uniroot(f,c(0,2),rep=(1-conf.level)/2,extendInt = "yes")
    lambda.1 <- out$root
    delta.1 <- lambda.1*sqrt(1/n1+1/n2)   # See explanation for delta.2
    
    
    # upper limit = limit of lambda such as pt(q=t_obs, df=n1+n2-2, ncp = lambda) = (1-conf.level)/2 = alpha/2
    f=function(lambda,rep) pt(q=t_obs, df=n1+n2-2, ncp = lambda)-rep
    out=uniroot(f,c(0,2),rep=(1-conf.level)/2,extendInt = "yes")
    lambda.2 <- out$root
    delta.2 <- lambda.2*sqrt(1/n1+1/n2)   # Pr[t_obs(alpha/2,df=df,ncp=lambda) <= lambda <= t_obs(1-alpha/2,df=df,ncp=lambda)=.025
    # because lambda = delta * sqrt[n1n2/(n1+n2)] :
    # Pr[t_obs(alpha/2,df=df,ncp=lambda) <= delta * sqrt[n1n2/(n1+n2)] <= t_obs(1-alpha/2,df=df,ncp=lambda)=.025
    # Pr[t_obs(alpha/2,df=df,ncp=lambda)*sqrt[(n1+n2)/n1n2] <= delta <= t_obs(1-alpha/2,df=df,ncp=lambda)*sqrt[(n1+n2)/n1n2]=.025
    # and sqrt[(n1+n2)/n1n2]=sqrt(1/n1+1/n2)
    result <- c(delta.1, delta.2) 
    
  } else if (alternative == "greater"){
    
    t_obs <- t.test(Group.1, Group.2, alternative = "greater", var.equal = TRUE)$statistic
    
    # lower limit = limit of lambda such as 1-pt(q=t_obs, df=n1+n2-2, ncp = lambda) = (1-conf.level) = alpha
    f=function(lambda,rep) 1-pt(q=t_obs, df=n1+n2-2, ncp = lambda)-rep
    out=uniroot(f,c(0,2),rep=1-conf.level,extendInt = "yes")

    lambda.1 <- out$root
    delta.1 <- lambda.1*sqrt(1/n1+1/n2)   # See explanation for delta.2
    
    # upper limit = +Inf
    delta.2 <- +Inf
    result <- c(delta.1, delta.2) 
    
  } else if (alternative == "less"){
    
    #  perform two-sample t-test assuming equal variance (same assumptions of cohen's d)
    t_obs <- t.test(Group.1, Group.2, alternative = "less", var.equal = TRUE)$statistic
    
    # lower limit = -Inf
    delta.1 <- -Inf
    
    # upper limit = limit of lambda such as pt(q=t_obs, df=n1+n2-2, ncp = lambda) = (1-conf.level) = alpha
    f=function(lambda,rep) pt(q=t_obs, df=n1+n2-2, ncp = lambda)-rep
    out=uniroot(f,c(0,2),rep=1-conf.level,extendInt = "yes")
    lambda.2 <- out$root
    delta.2 <- lambda.2*sqrt(1/n1+1/n2)
    
    result <- c(delta.1, delta.2) 
    
  }
  
  return(result)
}

#--------------------------------------------------------------
#  Obtain Welch confidence limits following Shieh's procedure
#--------------------------------------------------------------
 
Shieh.CI <- function(Group.1, Group.2,conf.level,alternative="two.sided")
{ 
  n1 <- length(Group.1)
  n2 <- length(Group.2)
  s1 <- sd(Group.1)
  s2 <- sd(Group.2)
  
  if(alternative=="two.sided"){

    #  perform two-sample Welch t-test (same assumptions of Shieh's d)
    test <- t.test(Group.1, Group.2, alternative = "two.sided", var.equal = FALSE)
    w_obs <- test$statistic
    #  sample estimates for degrees of freedom DF of noncentral t distribution
    DF <- test$parameter 
    
    # lower limit = limit of lambda such as 1-pt(q=t_obs, df=DF, ncp = lambda) = (1-conf.level)/2 = alpha/2
    # with DF = (sd1^2/n1 + sd2^2/n2)^2 / ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))    
    
    f=function(lambda,rep) 1-pt(q=w_obs, df=DF, ncp = lambda)-rep
    out=uniroot(f,c(0,2),rep=(1-conf.level)/2,extendInt = "yes")
    
    lambda.1 <- out$root
    delta.1 <- lambda.1/sqrt(n1+n2)
    
    # upper limit = limit of lambda such as pt(q=t_obs, df=DF, ncp = lambda) = (1-conf.level)/2 = alpha/2
    # with DF = (sd1^2/n1 + sd2^2/n2)^2 / ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))    
    
    f=function(lambda,rep) pt(q=w_obs, df=DF, ncp = lambda)-rep
    out=uniroot(f,c(0,2),rep=(1-conf.level)/2,extendInt = "yes")
    lambda.2 <- out$root
    delta.2 <- lambda.2/sqrt(n1+n2)
    
    result <- c(delta.1, delta.2) 
    
  } else if (alternative == "greater"){

    #  perform two-sample Welch t-test (same assumptions of Shieh's d)
    test <- t.test(Group.1, Group.2, alternative = "greater", var.equal = FALSE)
    w_obs <- test$statistic
    #  sample estimates for degrees of freedom DF of noncentral t distribution
    DF <- test$parameter 
    
    # lower limit = limit of lambda such as 1-pt(q=t_obs, df=DF, ncp = lambda) = (1-conf.level) = alpha
    # with DF = (sd1^2/n1 + sd2^2/n2)^2 / ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))    
    
    f=function(lambda,rep) 1-pt(q=w_obs, df=DF, ncp = lambda)-rep
    out=uniroot(f,c(0,2),rep=1-conf.level,extendInt = "yes")
    
    lambda.1 <- out$root
    delta.1 <- lambda.1/sqrt(n1+n2)
    
    # upper limit = limit of lambda such as pt(q=t_obs, df=DF, ncp = lambda) = (1-conf.level) = alpha
    # with DF = (sd1^2/n1 + sd2^2/n2)^2 / ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))    
    
    delta.2 <- +Inf
    
    result <- c(delta.1, delta.2)    
    
  } else if (alternative == "less"){
    
    #  perform two-sample Welch t-test (same assumptions of Shieh's d)
    test <- t.test(Group.1, Group.2, alternative = "less", var.equal = FALSE)
    w_obs <- test$statistic
    #  sample estimates for degrees of freedom DF of noncentral t distribution
    DF <- test$parameter 
    
    # lower limit = limit of lambda such as 1-pt(q=t_obs, df=DF, ncp = lambda) = (1-conf.level) = alpha
    # with DF = (sd1^2/n1 + sd2^2/n2)^2 / ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))    
    
    delta.1 <- -Inf
    
    # upper limit = limit of lambda such as pt(q=t_obs, df=DF, ncp = lambda) = (1-conf.level) = alpha
    # with DF = (sd1^2/n1 + sd2^2/n2)^2 / ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))    
    
    f=function(lambda,rep) pt(q=w_obs, df=DF, ncp = lambda)-rep
    out=uniroot(f,c(0,2),rep=1-conf.level,extendInt = "yes")
    lambda.2 <- out$root
    delta.2 <- lambda.2/sqrt(n1+n2)
    
    result <- c(delta.1, delta.2) 
    
  }
  
  return(result)
}















