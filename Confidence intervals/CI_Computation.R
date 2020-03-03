for (package in "rootSolve") {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}


#########################################
##  Simulation - Confidence intervals  ##
#########################################

#--------------------------------------------------------------
#  Obtain Student confidence limits following Cohen's procedure
#--------------------------------------------------------------

Cohen.CI <- function(Group.1, Group.2,alpha)
{ 
  n1 <- length(Group.1)
  n2 <- length(Group.2)
  
  #  perform two-sample t-test assuming equal variance (same assumptions of cohen's d)
  t_obs <- t.test(Group.1, Group.2, alternative = "two.sided", var.equal = TRUE)$statistic
  
  #  find limits of lambda at pt = alpha/2 or 1-alpha/2
  
  # Limit of lambda such as pt(q=t_obs, df=n1+n2-2, ncp = lambda) = alpha/2
  f=function(lambda,rep) pt(q=t_obs, df=n1+n2-2, ncp = lambda)-rep
  out=uniroot(f,c(-100,100),rep=alpha/2)
  lambda.1 <- out$root
  delta.1 <- lambda.1*sqrt(1/n1+1/n2)   # Pr[t_obs(alpha/2,df=df,ncp=lambda) <= lambda <= t_obs(1-alpha/2,df=df,ncp=lambda)=.025
                                        # because lambda = delta * sqrt[n1n2/(n1+n2)] :
                                        # Pr[t_obs(alpha/2,df=df,ncp=lambda) <= delta * sqrt[n1n2/(n1+n2)] <= t_obs(1-alpha/2,df=df,ncp=lambda)=.025
                                        # Pr[t_obs(alpha/2,df=df,ncp=lambda)*sqrt[(n1+n2)/n1n2] <= delta <= t_obs(1-alpha/2,df=df,ncp=lambda)*sqrt[(n1+n2)/n1n2]=.025
                                        # and sqrt[(n1+n2)/n1n2]=sqrt(1/n1+1/n2)

  # Limit of lambda such as 1-pt(q=t_obs, df=n1+n2-2, ncp = lambda) = alpha/2
  f=function(lambda,rep) 1-pt(q=t_obs, df=n1+n2-2, ncp = lambda)-rep
  out=uniroot(f,c(-100,100),rep=alpha/2)
  lambda.2 <- out$root
  delta.2 <- lambda.2*sqrt(1/n1+1/n2)   # See explanation for delta.1

  delta.low <- min(delta.1, delta.2)
  delta.upp <- max(delta.1, delta.2)
  
  result <- c(delta.low, delta.upp) 
  return(result)
}

















