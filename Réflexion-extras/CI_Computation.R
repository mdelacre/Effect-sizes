for (package in "rootSolve") {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
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















